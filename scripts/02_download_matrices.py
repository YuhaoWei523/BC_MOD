import os
import pandas as pd
import requests
import tarfile
from tqdm import tqdm

# --- Configuration ---
INPUT_CSV = "../data/breast_cancer_scRNA_seq_GSE.csv"
DATA_DIR = "../data"


class GEOMatrixDownloader:
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        os.makedirs(output_dir, exist_ok=True)

    def download_file(self, url, dest_path):
        """Stream download with progress bar."""
        try:
            with self.session.get(url, stream=True, timeout=60) as r:
                if r.status_code == 404: return False
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))

                with open(dest_path, 'wb') as f, tqdm(
                        desc=os.path.basename(dest_path),
                        total=total_size,
                        unit='iB',
                        unit_scale=True
                ) as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        bar.update(size)
            return True
        except Exception as e:
            print(f"Error downloading {url}: {e}")
            return False

    def extract_tar(self, tar_path, extract_path):
        """Extract tar/tar.gz files."""
        try:
            if tarfile.is_tarfile(tar_path):
                with tarfile.open(tar_path) as tar:
                    tar.extractall(path=extract_path)
                print(f"Extracted: {tar_path}")
                return True
        except Exception as e:
            print(f"Extraction failed: {e}")
        return False

    def process_dataset(self, gse_id):
        """Construct URL and download."""
        base_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file"
        # Standard filenames for 10x data usually stored in RAW.tar
        filename = f"{gse_id}_RAW.tar"
        save_path = os.path.join(self.output_dir, gse_id, filename)

        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        if os.path.exists(save_path):
            print(f"File exists: {save_path}, skipping download.")
            self.extract_tar(save_path, os.path.dirname(save_path))
            return

        print(f"Downloading {gse_id}...")
        if self.download_file(base_url, save_path):
            self.extract_tar(save_path, os.path.dirname(save_path))


if __name__ == "__main__":
    if os.path.exists(INPUT_CSV):
        df = pd.read_csv(INPUT_CSV)
        downloader = GEOMatrixDownloader(DATA_DIR)

        # Example loop
        for gse_id in df['gse_id'].unique():
            downloader.process_dataset(gse_id)
    else:
        print(f"CSV not found at {INPUT_CSV}")