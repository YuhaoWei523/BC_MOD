import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

# 1. 设置参数
DATA_DIR = "./data"
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')


def load_and_preprocess():
    adatas = []

    # --- A. 加载数据 ---
    for gse_id in os.listdir(DATA_DIR):
        gse_path = os.path.join(DATA_DIR, gse_id)
        if os.path.isdir(gse_path) and "matrix.mtx.gz" in os.listdir(gse_path):
            try:
                print(f"Loading {gse_id}...")
                adata = sc.read_10x_mtx(gse_path, cache=True)
                adata.obs['batch'] = gse_id  # 标记来源，用于后续去批次效应

                # 下面这一步需要你手动准备一个 CSV，或者简单地全部标记为 "Unknown" 后续填补
                # adata.obs['cancer_subtype'] = 'TNBC'

                adata.var_names_make_unique()
                adatas.append(adata)
            except Exception as e:
                print(f"跳过 {gse_id}: {e}")

    if not adatas:
        print("未加载到数据，请检查文件名！")
        return None

    # 合并所有数据集
    adata_full = sc.concat(adatas, join='outer')
    print(f"合并后规模: {adata_full.shape}")

    # --- B. 质量控制 (QC) ---
    # 按照中期报告标准：过滤低质量细胞 [cite: 10]
    adata_full.var['mt'] = adata_full.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_full, qc_vars=['mt'], inplace=True)

    # 过滤：基因数 > 200, 线粒体 < 20%
    adata_full = adata_full[adata_full.obs.n_genes_by_counts > 200, :]
    adata_full = adata_full[adata_full.obs.pct_counts_mt < 20, :]

    # --- C. 归一化与降维 ---
    # 使用 LogNormalize (SCTransform 较慢，对于课程设计 LogNorm 足够且稳健)
    sc.pp.normalize_total(adata_full, target_sum=1e4)
    sc.pp.log1p(adata_full)
    sc.pp.highly_variable_genes(adata_full, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # 保存原始数据到 raw (用于后续差异分析)
    adata_full.raw = adata_full

    # 仅保留高变基因进行计算
    adata_full = adata_full[:, adata_full.var.highly_variable]
    sc.pp.scale(adata_full, max_value=10)
    sc.tl.pca(adata_full, svd_solver='arpack')

    # --- D. 批次校正 (Harmony) ---
    # 既然中期报告提到了 Harmony ，这里必须使用
    # pip install harmony-pytorch
    try:
        import harmony
        print("Running Harmony integration...")
        sc.external.pp.harmony_integrate(adata_full, 'batch')
        use_rep = 'X_pca_harmony'
    except:
        print("Harmony未安装，使用标准PCA")
        use_rep = 'X_pca'

    # --- E. 聚类与注释 ---
    sc.pp.neighbors(adata_full, n_neighbors=10, n_pcs=40, use_rep=use_rep)
    sc.tl.umap(adata_full)
    sc.tl.leiden(adata_full, resolution=0.5)  # 分群

    # *关键步骤*：你需要为这些数字 cluster (0, 1, 2...) 赋予生物学名称
    # 在课程中，可以使用简化的 Marker 自动注释
    marker_genes = {
        "T_cells": ["CD3D", "CD3E", "CD2"],
        "B_cells": ["CD79A", "MS4A1"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18"],
        "Macrophage": ["CD68", "CD163"],
        "Fibroblast": ["COL1A1", "DCN"]
    }

    # 简单的打分注释 (实际科研需人工核对)
    sc.tl.score_genes_cell_cycle(adata_full, s_genes=[], g2m_genes=[])  # 初始化
    for cell_type, genes in marker_genes.items():
        valid_genes = [g for g in genes if g in adata_full.var_names]
        if valid_genes:
            sc.tl.score_genes(adata_full, valid_genes, score_name=cell_type)

    # 将得分为最大值的类型赋给细胞
    scores = adata_full.obs[marker_genes.keys()]
    adata_full.obs['cell_type_major'] = scores.idxmax(axis=1)

    return adata_full

# 运行预处理
# processed_adata = load_and_preprocess()
# processed_adata.write("./data/integrated_breast_cancer.h5ad")