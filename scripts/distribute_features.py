import os
import shutil
import glob

DATA_DIR = "./data"

# 你手动下载的文件可能叫不同名字，这里定义好
# 格式: "GSE_ID": "你现在的features文件名"
MISSING_MAP = {
    "GSE161529": "GSE161529_features.tsv.gz",
    "GSE306201": "GSE306201_Features.tsv.gz"  # 注意大小写，根据你实际下载的文件名修改
}


def distribute_shared_features():
    for gse_id, feature_filename in MISSING_MAP.items():
        gse_path = os.path.join(DATA_DIR, gse_id)
        source_file = os.path.join(gse_path, feature_filename)

        # 1. 检查源文件是否存在
        if not os.path.exists(source_file):
            # 尝试模糊搜索，防止文件名大小写或后缀微小差异
            candidates = glob.glob(os.path.join(gse_path, "*eatures.tsv.gz"))
            if candidates:
                source_file = candidates[0]
                print(f"⚠️ 在 {gse_id} 中找到了 {os.path.basename(source_file)}，将使用它。")
            else:
                print(f"❌ 错误: 在 {gse_path} 下找不到 {feature_filename}，请确认你已放入该文件。")
                continue

        print(f"📂 正在为 {gse_id} 分发 features 文件...")

        # 2. 遍历所有 GSM 子文件夹
        count = 0
        for gsm_id in os.listdir(gse_path):
            gsm_path = os.path.join(gse_path, gsm_id)

            # 确保是文件夹且是 GSM 开头
            if os.path.isdir(gsm_path) and gsm_id.startswith("GSM"):
                target_path = os.path.join(gsm_path, "features.tsv.gz")

                # 如果目标文件夹里没有 features，就复制进去
                if not os.path.exists(target_path):
                    try:
                        shutil.copy(source_file, target_path)
                        count += 1
                    except Exception as e:
                        print(f"  ❌ 复制失败 {gsm_id}: {e}")
                else:
                    # 如果已经有了，可能是之前遗留的或者别的名字，可以选覆盖或跳过
                    pass

        print(f"✅ {gse_id}: 已将 features 文件分发到 {count} 个样本文件夹中。")


if __name__ == "__main__":
    distribute_shared_features()