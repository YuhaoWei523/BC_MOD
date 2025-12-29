import os
import pandas as pd
import requests
from bs4 import BeautifulSoup
import time
import tarfile
from tqdm import tqdm


class GEOMatrixDownloaderV2:
    def __init__(self, output_dir="./data"):
        self.output_dir = output_dir
        self.session = requests.Session()
        # 模拟浏览器头，防止被拒
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        os.makedirs(output_dir, exist_ok=True)

    def download_file(self, url, dest_path):
        """流式下载文件"""
        try:
            with self.session.get(url, stream=True, timeout=60) as r:
                if r.status_code == 404:
                    return False
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))

                with open(dest_path, 'wb') as f, tqdm(
                        desc=os.path.basename(dest_path),
                        total=total_size,
                        unit='iB',
                        unit_scale=True,
                        unit_divisor=1024,
                ) as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        bar.update(len(chunk))
            return True
        except Exception as e:
            if os.path.exists(dest_path):
                os.remove(dest_path)
            return False

    def construct_ftp_url(self, gse_id):
        """
        策略A：直接构造 FTP/HTTP 路径（最稳健，不需要解析网页）
        规律：https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE12345/suppl/GSE12345_RAW.tar
        """
        # 提取数字部分
        num_part = gse_id.replace("GSE", "")
        # GEO 的存储桶逻辑：GSE161529 -> GSE161nnn
        if len(num_part) < 3:
            stub = "GSE" + num_part
        else:
            stub = "GSE" + num_part[:-3] + "nnn"

        base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse_id}/suppl"

        # 尝试几种常见的文件名组合
        candidates = [
            f"{base_url}/{gse_id}_RAW.tar",
            f"{base_url}/{gse_id}_RAW.tar.gz"
        ]
        return candidates

    def extract_tar(self, file_path, extract_path):
        try:
            print(f"  📦 解压中: {file_path}")
            if file_path.endswith('.tar'):
                with tarfile.open(file_path, 'r:') as tar:
                    tar.extractall(path=extract_path)
            return True
        except Exception as e:
            print(f"  ❌ 解压错误: {e}")
            return False

    def process_gse(self, gse_id):
        # 修正输入可能带有的重复 GSE 前缀
        if gse_id.startswith("GSEGSE"):
            gse_id = gse_id.replace("GSEGSE", "GSE")

        print(f"\n{'=' * 60}")
        print(f"处理数据集: {gse_id}")

        gse_dir = os.path.join(self.output_dir, gse_id)
        os.makedirs(gse_dir, exist_ok=True)

        # --- 🔍 更新后的检查逻辑: 递归查找 .mtx.gz ---
        # 原因：解压后或整理后的文件通常在 GSMxxx 子文件夹中，os.listdir 看不到
        has_data = False
        for root, dirs, files in os.walk(gse_dir):
            if any(f.endswith('matrix.mtx.gz') or f.endswith('.mtx.gz') for f in files):
                has_data = True
                break

        if has_data:
            print("  ✅ 检测到子目录中已存在 .mtx.gz 文件，跳过下载。")
            return 'skipped'
        # ------------------------------------------------

        # --- 策略 A: 暴力构造链接 (优先尝试) ---
        print("  🔄 尝试直接构造下载链接...")
        candidates = self.construct_ftp_url(gse_id)
        success = False

        for url in candidates:
            filename = url.split('/')[-1]
            save_path = os.path.join(gse_dir, filename)

            # 额外检查：如果 tar 包已经存在但没解压（可能是上次中断了）
            if os.path.exists(save_path):
                print(f"  📦 发现已存在的压缩包: {filename}，尝试直接解压...")
                if self.extract_tar(save_path, gse_dir):
                    success = True
                    break
                else:
                    print("  ⚠️ 现有压缩包解压失败（可能已损坏），重新下载...")
                    os.remove(save_path)  # 删除坏文件

            print(f"  ➡️ 尝试下载: {filename}")

            if self.download_file(url, save_path):
                print("  ✅ 下载成功！")
                self.extract_tar(save_path, gse_dir)
                success = True
                break

        if success:
            return 'success'

        # --- 策略 B: 爬虫解析 (备用) ---
        print("  ⚠️ 直接构造失败，尝试解析网页...")
        print(f"  ❌ {gse_id} 自动下载失败，请手动检查: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}")
        return 'failed'

    def run_from_csv(self, csv_path):
        df = pd.read_csv(csv_path)
        # 确保去重
        gse_ids = df['gse_id'].unique()
        print(f"待处理数据集: {len(gse_ids)} 个")

        for gse_id in gse_ids:
            self.process_gse(gse_id)


if __name__ == "__main__":
    # 使用之前生成的 CSV
    csv_file = "./data/breast_cancer_scRNA_seq_GSE.csv"
    downloader = GEOMatrixDownloaderV2(output_dir="./data")
    if os.path.exists(csv_file):
        downloader.run_from_csv(csv_file)
    else:
        print("未找到 CSV 文件")