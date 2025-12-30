import os
import shutil
import re

DATA_DIR = "./data"


def organize_gse_folder(gse_id):
    gse_path = os.path.join(DATA_DIR, gse_id)
    if not os.path.isdir(gse_path):
        return

    print(f"📂 正在整理: {gse_id} ...")

    # 1. 扫描所有文件
    files = [f for f in os.listdir(gse_path) if os.path.isfile(os.path.join(gse_path, f))]

    for filename in files:
        # 跳过压缩包本身
        if filename.endswith('.tar'):
            continue

        # 2. 提取 GSM 编号 (例如从 GSM4909319_... 中提取 GSM4909319)
        match = re.search(r'(GSM\d+)', filename)
        if not match:
            print(f"  ⚠️ 跳过无 GSM 编号文件: {filename}")
            continue

        gsm_id = match.group(1)

        # 3. 创建 GSM 子文件夹
        gsm_dir = os.path.join(gse_path, gsm_id)
        os.makedirs(gsm_dir, exist_ok=True)

        # 4. 移动并重命名
        src_path = os.path.join(gse_path, filename)

        # 判断目标文件名
        lower_name = filename.lower()
        target_name = None

        if "matrix" in lower_name and "mtx" in lower_name:
            target_name = "matrix.mtx.gz"
        elif "barcodes" in lower_name and "tsv" in lower_name:
            target_name = "barcodes.tsv.gz"
        elif ("features" in lower_name or "genes" in lower_name) and "tsv" in lower_name:
            target_name = "features.tsv.gz"

        # 5. 执行移动
        if target_name:
            dst_path = os.path.join(gsm_dir, target_name)
            # 防止自己覆盖自己
            if src_path != dst_path:
                shutil.move(src_path, dst_path)
                print(f"  ✅ 归档: {gsm_id} -> {target_name}")
        else:
            # 如果是 quant.sf 或其他非 10x 文件，直接原样移入文件夹，不改名
            dst_path = os.path.join(gsm_dir, filename)
            shutil.move(src_path, dst_path)
            print(f"  ℹ️ 移动非标准文件: {filename} -> {gsm_dir}")

    # 清理空文件夹逻辑可按需添加


# 执行所有 GSE
for gse_id in os.listdir(DATA_DIR):
    if gse_id.startswith("GSE"):
        organize_gse_folder(gse_id)