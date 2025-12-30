import os
import shutil
import re
import glob

# --- Configuration ---
DATA_DIR = "../data"

# Mapping for datasets that share a single features file across all samples
SHARED_FEATURES_MAP = {
    "GSE161529": "GSE161529_features.tsv.gz",
    "GSE306201": "GSE306201_Features.tsv.gz"
}


class FileManager:
    def __init__(self, data_dir):
        self.data_dir = data_dir

    def standardize_filenames(self, gse_id):
        """
        Renames raw files (e.g., GSM123_barcodes.tsv.gz) to
        standard 10x names (barcodes.tsv.gz).
        """
        gse_path = os.path.join(self.data_dir, gse_id)
        if not os.path.isdir(gse_path): return

        print(f"📂 Organizing folder: {gse_id} ...")
        files = [f for f in os.listdir(gse_path) if os.path.isfile(os.path.join(gse_path, f))]

        for filename in files:
            if filename.endswith('.tar'): continue

            # Extract GSM ID
            match = re.search(r'(GSM\d+)', filename)
            if not match: continue
            gsm_id = match.group(1)

            # Create Sample Directory
            gsm_dir = os.path.join(gse_path, gsm_id)
            os.makedirs(gsm_dir, exist_ok=True)

            # Determine Target Name
            lower_name = filename.lower()
            target_name = None
            if "matrix" in lower_name and "mtx" in lower_name:
                target_name = "matrix.mtx.gz"
            elif "barcodes" in lower_name and "tsv" in lower_name:
                target_name = "barcodes.tsv.gz"
            elif ("features" in lower_name or "genes" in lower_name) and "tsv" in lower_name:
                target_name = "features.tsv.gz"

            # Move File
            if target_name:
                src = os.path.join(gse_path, filename)
                dst = os.path.join(gsm_dir, target_name)
                if src != dst:
                    shutil.move(src, dst)
                    # print(f"  Moved: {filename} -> {gsm_id}/{target_name}")

    def distribute_features(self):
        """
        Copies the shared features file to individual GSM folders
        for specific datasets.
        """
        print("\n🔄 Distributing shared feature files...")
        for gse_id, feat_name in SHARED_FEATURES_MAP.items():
            gse_path = os.path.join(self.data_dir, gse_id)
            if not os.path.exists(gse_path): continue

            # Locate source feature file (allow fuzzy matching)
            source_file = os.path.join(gse_path, feat_name)
            if not os.path.exists(source_file):
                candidates = glob.glob(os.path.join(gse_path, "*eatures.tsv.gz"))
                if candidates:
                    source_file = candidates[0]
                else:
                    print(f"  ⚠️ Warning: Shared feature file not found for {gse_id}")
                    continue

            # Distribute to GSM folders
            for gsm_id in os.listdir(gse_path):
                gsm_path = os.path.join(gse_path, gsm_id)
                if os.path.isdir(gsm_path) and gsm_id.startswith("GSM"):
                    target = os.path.join(gsm_path, "features.tsv.gz")
                    if not os.path.exists(target):
                        shutil.copy(source_file, target)
                        print(f"  ✅ Copied features to {gse_id}/{gsm_id}")


if __name__ == "__main__":
    manager = FileManager(DATA_DIR)

    # 1. First, organize filenames into folders
    for gse in os.listdir(DATA_DIR):
        if gse.startswith("GSE"):
            manager.standardize_filenames(gse)

    # 2. Then, distribute missing feature files
    manager.distribute_features()
    print("\n🎉 File organization complete!")