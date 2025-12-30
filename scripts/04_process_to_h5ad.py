import scanpy as sc
import os
import warnings
import gc
import matplotlib
import pandas as pd

matplotlib.use('Agg')  # Backend for non-GUI plotting
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- Configuration ---
DATA_DIR = "../data"
OUTPUT_FILE = "../results/processed_breast_cancer.h5ad"
FIGURES_DIR = "../results/figures"
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ⚡ Memory Protection Limit ⚡
MAX_CELLS_PER_SAMPLE = 3000
JOIN_METHOD = 'inner'

# Subtype Mappings (Metadata)
GSE161529_GSM_MAP = {
    # ... (Please paste your full GSM mapping dictionary here) ...
    # For brevity in this response, assume it is populated
}

DATASET_LEVEL_MAP = {
    "GSE240112": "Luminal_A",
    "GSE262288": "Luminal_B",
    "GSE274139": "Luminal_B",
    "GSE306201": "Luminal_A",
    "GSE289825": "TNBC"
}


def get_sample_subtype(gse_id, gsm_id):
    """Retrieve subtype based on GSE/GSM ID."""
    if gse_id == "GSE161529":
        return GSE161529_GSM_MAP.get(gsm_id, "Unknown")
    return DATASET_LEVEL_MAP.get(gse_id, "Unknown")


def safe_load_10x(directory):
    """Safely load 10x Genomics matrix with error handling for headers."""
    try:
        return sc.read_10x_mtx(directory, cache=False)
    except ValueError as e:
        if "Length of passed value" in str(e) and "var_names" in str(e):
            print(f"    ⚠️ Fixing features header mismatch...")
            mtx_path = os.path.join(directory, "matrix.mtx.gz")
            feats_path = os.path.join(directory, "features.tsv.gz")
            if not os.path.exists(feats_path):
                feats_path = os.path.join(directory, "genes.tsv.gz")

            adata = sc.read_mtx(mtx_path).T
            features = pd.read_csv(feats_path, header=None, sep='\t')

            # Handle header row mismatch
            if len(features) == adata.n_vars + 1:
                features = features.iloc[1:]

            if features.shape[1] > 1:
                adata.var_names = features.iloc[:, 1].values
                adata.var['gene_ids'] = features.iloc[:, 0].values
            else:
                adata.var_names = features.iloc[:, 0].values
            return adata
        else:
            raise e


def load_and_merge_incrementally():
    """Incrementally load and merge samples to save memory."""
    merged_adata = None
    count = 0

    print(f"🚀 Starting Incremental Merge (Limit: {MAX_CELLS_PER_SAMPLE}, Join={JOIN_METHOD})...")

    for gse_id in os.listdir(DATA_DIR):
        gse_path = os.path.join(DATA_DIR, gse_id)
        if not os.path.isdir(gse_path) or gse_id == "GSE274139": continue

        for root, dirs, files in os.walk(gse_path):
            if "matrix.mtx.gz" in files and "barcodes.tsv.gz" in files:
                try:
                    sample_id = os.path.basename(root)
                    subtype = get_sample_subtype(gse_id, sample_id)

                    # 1. Load single sample
                    current_adata = safe_load_10x(root)
                    current_adata.var_names_make_unique()

                    # 2. Downsample
                    if MAX_CELLS_PER_SAMPLE is not None and current_adata.n_obs > MAX_CELLS_PER_SAMPLE:
                        sc.pp.subsample(current_adata, n_obs=MAX_CELLS_PER_SAMPLE, random_state=42)
                        print(f"  [{count + 1}] Loaded: {gse_id}/{sample_id} -> Subsampled to {MAX_CELLS_PER_SAMPLE}")
                    else:
                        print(f"  [{count + 1}] Loaded: {gse_id}/{sample_id} ({current_adata.n_obs} cells)")

                    # 3. Add Metadata
                    current_adata.obs['batch'] = gse_id
                    current_adata.obs['sample_id'] = sample_id
                    current_adata.obs['subtype'] = subtype

                    # 4. Merge
                    if merged_adata is None:
                        merged_adata = current_adata
                    else:
                        merged_adata = sc.concat(
                            [merged_adata, current_adata],
                            join=JOIN_METHOD,
                            index_unique=None
                        )

                    # 5. Garbage Collection
                    del current_adata
                    gc.collect()
                    count += 1

                except Exception as e:
                    print(f"❌ Error loading {sample_id}: {e}")

    if merged_adata is None: return None
    merged_adata.obs_names_make_unique()
    print(f"✅ Merge Complete! Final Shape: {merged_adata.shape}")
    return merged_adata


def run_preprocess(adata):
    print("\n⚡ Starting Preprocessing...")

    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 25, :]
    sc.pp.filter_cells(adata, min_genes=200)
    print(f"  Cells after QC: {adata.n_obs}")

    # Normalize & Log
    print("  Normalizing & Log1p...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG
    print("  Selecting HVGs (Top 3000)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False)

    # Backup Raw
    adata.raw = adata

    # Slice (Keep only HVG for PCA to save memory)
    print("  Slicing to HVGs...")
    adata = adata[:, adata.var.highly_variable]
    gc.collect()

    # PCA & Scale
    print("  Scaling & PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')

    # Integration (Harmony or PCA)
    try:
        import harmony
        print("  Running Harmony Integration...")
        sc.external.pp.harmony_integrate(adata, 'batch')
        use_rep = 'X_pca_harmony'
    except:
        print("  Harmony not found, using PCA...")
        use_rep = 'X_pca'

    # Clustering
    print("  Clustering (Leiden)...")
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # Annotation
    print("  Annotating Cell Types...")
    marker_genes = {
        "T_cells": ["CD3D", "CD3E", "CD2"],
        "B_cells": ["CD79A", "MS4A1"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18"],
        "Macrophage": ["CD68", "CD163", "AIF1"],
        "Fibroblast": ["COL1A1", "DCN"],
        "Endothelial": ["PECAM1", "VWF"]
    }

    scored_types = []
    for cell_type, genes in marker_genes.items():
        # Look up genes in RAW (Full genome)
        valid_genes = [g for g in genes if g in adata.raw.var_names]
        if valid_genes:
            sc.tl.score_genes(adata, valid_genes, score_name=cell_type, use_raw=True)
            scored_types.append(cell_type)

    if scored_types:
        scores = adata.obs[scored_types]
        adata.obs['cell_type_major'] = scores.idxmax(axis=1)
    else:
        adata.obs['cell_type_major'] = "Unknown"

    # Plotting
    print("  Saving UMAP...")
    sc.pl.umap(
        adata,
        color=['cell_type_major', 'batch', 'subtype'],
        save="_all_cells_integrated.png",
        show=False
    )
    # Move plot to figures dir
    if os.path.exists("figures/umap_all_cells_integrated.png"):
        shutil.move("figures/umap_all_cells_integrated.png", os.path.join(FIGURES_DIR, "umap.png"))

    return adata


if __name__ == "__main__":
    adata_merged = load_and_merge_incrementally()
    if adata_merged:
        adata_processed = run_preprocess(adata_merged)
        print(f"\n💾 Saving to: {OUTPUT_FILE}")
        adata_processed.write(OUTPUT_FILE, compression="gzip")
        print("\n🎉 Preprocessing Complete!")