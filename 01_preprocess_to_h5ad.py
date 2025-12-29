import scanpy as sc
import os
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- 配置 ---
DATA_DIR = "./data"
OUTPUT_FILE = "./results/processed_breast_cancer.h5ad"  # 中间结果保存路径
os.makedirs("./results", exist_ok=True)

# ==========================================
# 1. 定义精确的 GSM 映射字典 (GSE161529 专用)
# ==========================================
GSE161529_GSM_MAP = {
    # --- Normal / Control ---
    'GSM4909253': 'Normal', 'GSM4909254': 'Normal', 'GSM4909255': 'Normal',
    'GSM4909256': 'Normal', 'GSM4909257': 'Normal', 'GSM4909258': 'Normal',
    'GSM4909259': 'Normal', 'GSM4909260': 'Normal', 'GSM4909261': 'Normal',
    'GSM4909262': 'Normal', 'GSM4909263': 'Normal', 'GSM4909264': 'Normal',
    'GSM4909265': 'Normal', 'GSM4909266': 'Normal', 'GSM4909267': 'Normal',
    'GSM4909268': 'Normal', 'GSM4909269': 'Normal', 'GSM4909270': 'Normal',
    'GSM4909271': 'Normal', 'GSM4909272': 'Normal', 'GSM4909273': 'Normal',
    'GSM4909274': 'Normal', 'GSM4909275': 'Normal', 'GSM4909276': 'Normal',
    # 将 BRCA1 pre-neoplastic (癌前病变) 归为 Normal/Control 组，或可标记为 Normal
    'GSM4909277': 'Normal', 'GSM4909278': 'Normal', 'GSM4909279': 'Normal', 'GSM4909280': 'Normal',

    # --- TNBC (三阴性) ---
    'GSM4909281': 'TNBC', 'GSM4909282': 'TNBC', 'GSM4909283': 'TNBC', 'GSM4909284': 'TNBC',
    'GSM4909285': 'TNBC', 'GSM4909286': 'TNBC', 'GSM4909287': 'TNBC', 'GSM4909288': 'TNBC',

    # --- HER2 Positive ---
    'GSM4909289': 'HER2_Positive', 'GSM4909290': 'HER2_Positive', 'GSM4909291': 'HER2_Positive',
    'GSM4909292': 'HER2_Positive', 'GSM4909293': 'HER2_Positive', 'GSM4909294': 'HER2_Positive',

    # --- Luminal (ER+/PR+) -> 映射为 Luminal_A ---
    'GSM4909295': 'Luminal_A',  # PR+
    'GSM4909296': 'Luminal_A', 'GSM4909297': 'Luminal_A', 'GSM4909298': 'Luminal_A',
    'GSM4909299': 'Luminal_A', 'GSM4909300': 'Luminal_A', 'GSM4909301': 'Luminal_A',
    'GSM4909302': 'Luminal_A', 'GSM4909303': 'Luminal_A', 'GSM4909304': 'Luminal_A',
    'GSM4909305': 'Luminal_A', 'GSM4909306': 'Luminal_A', 'GSM4909307': 'Luminal_A',
    'GSM4909308': 'Luminal_A',  # Lymph-node 也算作该病人的亚型
    'GSM4909309': 'Luminal_A', 'GSM4909310': 'Luminal_A', 'GSM4909311': 'Luminal_A',
    'GSM4909312': 'Luminal_A', 'GSM4909313': 'Luminal_A', 'GSM4909314': 'Luminal_A',
    'GSM4909315': 'Luminal_A', 'GSM4909316': 'Luminal_A', 'GSM4909317': 'Luminal_A',
    'GSM4909318': 'Luminal_A', 'GSM4909319': 'Luminal_A', 'GSM4909320': 'Luminal_A',
    'GSM4909321': 'Luminal_A'
}

# 其他单一亚型数据集的映射 (Dataset Level Mapping)
DATASET_LEVEL_MAP = {
    "GSE240112": "Luminal_A",
    "GSE262288": "Luminal_B",
    "GSE274139": "Luminal_B",  # 注意：此数据集如果不兼容建议跳过
    "GSE306201": "Luminal_A",
    "GSE289825": "TNBC"
}


# ==========================================
# 2. 修改后的加载函数
# ==========================================
def get_sample_subtype(gse_id, gsm_id):
    """根据 GSE 和 GSM 返回准确的亚型"""
    # 优先检查是否存在于 GSM 级精细映射中 (针对 GSE161529)
    if gse_id == "GSE161529":
        return GSE161529_GSM_MAP.get(gsm_id, "Unknown")

    # 否则使用数据集级别的统一映射
    return DATASET_LEVEL_MAP.get(gse_id, "Unknown")


def load_data_recursive():
    adatas = []
    print(f"正在加载数据...")

    for gse_id in os.listdir(DATA_DIR):
        gse_path = os.path.join(DATA_DIR, gse_id)
        if not os.path.isdir(gse_path) or gse_id == "GSE274139": continue

        # 递归查找
        for root, dirs, files in os.walk(gse_path):
            if "matrix.mtx.gz" in files and "barcodes.tsv.gz" in files:
                try:
                    sample_id = os.path.basename(root)  # 获取 GSM 编号

                    # 🔍 核心修改：获取精准亚型
                    subtype = get_sample_subtype(gse_id, sample_id)

                    print(f"  - {gse_id} / {sample_id} -> {subtype}")

                    adata = sc.read_10x_mtx(root, cache=True)
                    adata.var_names_make_unique()

                    # 注入元数据
                    adata.obs['batch'] = gse_id
                    adata.obs['sample_id'] = sample_id
                    adata.obs['subtype'] = subtype  # 写入精准亚型

                    adatas.append(adata)
                except Exception as e:
                    print(f"❌ Error {sample_id}: {e}")

    if not adatas: return None
    print("正在合并所有样本...")
    return sc.concat(adatas, join='outer')

def run_preprocess(adata):
    print("\n开始预处理流水线...")
    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 25, :]
    sc.pp.filter_cells(adata, min_genes=200)

    # 归一化与降维
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata  # 备份全基因数据
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')

    # Harmony 去批次
    try:
        import harmony
        print("  Running Harmony...")
        sc.external.pp.harmony_integrate(adata, 'batch')
        use_rep = 'X_pca_harmony'
    except:
        print("  Harmony not found, using PCA")
        use_rep = 'X_pca'

    # 聚类
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # 自动注释 (Marker)
    print("  正在注释细胞类型...")
    marker_genes = {
        "T_cells": ["CD3D", "CD3E", "CD2"],
        "B_cells": ["CD79A", "MS4A1"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18"],
        "Macrophage": ["CD68", "CD163", "AIF1"],
        "Fibroblast": ["COL1A1", "DCN"],
        "Endothelial": ["PECAM1", "VWF"]
    }

    # 计算得分
    for cell_type, genes in marker_genes.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if valid_genes:
            sc.tl.score_genes(adata, valid_genes, score_name=cell_type)

    # 赋予最大分数的类型
    scores = adata.obs[list(marker_genes.keys())]
    adata.obs['cell_type_major'] = scores.idxmax(axis=1)

    return adata


if __name__ == "__main__":
    adata_merged = load_data_recursive()
    if adata_merged:
        adata_processed = run_preprocess(adata_merged)

        # ⚠️ 保存结果 ⚠️
        print(f"\n正在保存中间文件到: {OUTPUT_FILE}")
        adata_processed.write(OUTPUT_FILE, compression="gzip")
        print("✅ 预处理完成并保存！")