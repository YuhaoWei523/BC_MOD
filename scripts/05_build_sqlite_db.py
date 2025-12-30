import scanpy as sc
import sqlite3
import os
import warnings

warnings.filterwarnings("ignore")

# --- Configuration ---
INPUT_FILE = "../results/processed_breast_cancer.h5ad"
DB_NAME = "../dbs/scrna.db"  # Conforming to your app's path


def build_database():
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Input file not found: {INPUT_FILE}")
        return

    print(f"📂 Loading: {INPUT_FILE} ...")
    adata = sc.read_h5ad(INPUT_FILE)

    # Check annotations
    if 'subtype' not in adata.obs or 'cell_type_major' not in adata.obs:
        raise ValueError("Missing 'subtype' or 'cell_type_major' in adata.obs")

    # Create Group Keys
    print("🔄 Creating group keys...")
    adata.obs['group_key'] = adata.obs['subtype'].astype(str) + "::" + adata.obs['cell_type_major'].astype(str)
    unique_groups = adata.obs['group_key'].unique()
    print(f"   Identified {len(unique_groups)} biological groups.")

    # Calculate Pseudobulk
    print("📊 Calculating Pseudobulk Means...")

    # Use RAW data if available (contains all genes)
    use_adata = adata.raw.to_adata() if adata.raw else adata

    records = []

    for i, group in enumerate(unique_groups):
        if "Unknown" in group: continue
        subtype, cell_type = group.split("::")
        print(f"   Processing [{i + 1}/{len(unique_groups)}]: {subtype} - {cell_type}")

        # Mask & Extract
        cells_mask = adata.obs['group_key'] == group
        chunk_X = use_adata[cells_mask].X

        # Mean Expression
        if hasattr(chunk_X, "toarray"):
            mean_expression = chunk_X.mean(axis=0).A1
        else:
            mean_expression = chunk_X.mean(axis=0)

        gene_names = use_adata.var_names

        # Filter & Append (Expression > 0.01)
        for gene_idx, exp_val in enumerate(mean_expression):
            if exp_val > 0.01:
                records.append((
                    gene_names[gene_idx],
                    subtype,
                    cell_type,
                    float(f"{exp_val:.4f}")
                ))

    # Write to SQLite
    print(f"💾 Writing {len(records)} records to DB...")
    os.makedirs(os.path.dirname(DB_NAME), exist_ok=True)
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()

    cursor.execute('DROP TABLE IF EXISTS Table_Expression')
    cursor.execute('''
        CREATE TABLE Table_Expression (
            Gene TEXT,
            Subtype TEXT,
            CellType TEXT,
            Avg_Expression REAL
        )
    ''')

    cursor.executemany('INSERT INTO Table_Expression VALUES (?,?,?,?)', records)
    conn.commit()

    # Indices
    print("⚡ Creating Indices...")
    cursor.execute('CREATE INDEX idx_gene ON Table_Expression (Gene)')
    cursor.execute('CREATE INDEX idx_subtype_cell ON Table_Expression (Subtype, CellType)')

    conn.close()
    print(f"\n🎉 Database Built: {os.path.abspath(DB_NAME)}")


if __name__ == "__main__":
    build_database()