import os
import pandas as pd
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup
import time
from pathlib import Path
import subprocess
import sys
import re
import io  # Import io for text streaming
import xml.etree.ElementTree as ET  # Import the XML parser


class GEODataDownloader:
    """Download raw FASTQ data from GEO datasets using SRA Toolkit"""

    def __init__(self, output_dir="./data"):
        self.output_dir = output_dir
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        os.makedirs(output_dir, exist_ok=True)
        self.check_sra_toolkit()

    def check_sra_toolkit(self):
        """Check if SRA Toolkit is installed"""
        try:
            result = subprocess.run(['prefetch', '--version'],
                                    capture_output=True, text=True, timeout=5)
            print(f"✓ SRA Toolkit found: {result.stdout.strip()}")
            self.sra_available = True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            print("✗ SRA Toolkit not found!")
            print("\nTo download FASTQ files, please install SRA Toolkit:")
            print("  Ubuntu/Debian: sudo apt-get install sra-toolkit")
            print("  macOS: brew install sratoolkit")
            print("  Or download from: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit")
            print("\nAlternatively, the script will provide download commands you can run manually.")
            self.sra_available = False

    def get_sra_project_ids_from_gse(self, gse_id):
        """
        Gets the SRA Project (SRP) or BioProject (PRJNA) IDs from the
        main GSE series page.
        """
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"
        project_ids = []
        try:
            response = self.session.get(url, timeout=15)
            response.raise_for_status()
            soup = BeautifulSoup(response.content, 'html.parser')

            for link in soup.find_all('a', href=True):
                href = link.get('href', '')

                if 'bioproject/' in href:
                    matches = re.findall(r'(PRJNA\d+)', href)
                    project_ids.extend(matches)

                if 'Traces/study' in href or '/sra?term=SRP' in href:
                    matches = re.findall(r'(SRP\d+)', href)
                    project_ids.extend(matches)

            if not project_ids:
                print(f"    No BioProject (PRJNA) or SRA Study (SRP) links found on {gse_id} page.")

            return list(set(project_ids))

        except Exception as e:
            print(f"  Error getting project IDs for {gse_id}: {e}")
            return []

    def get_srr_from_eutils(self, project_id):
        """
        Gets all SRR run accessions using the NCBI E-utilities API.
        This version correctly parses the XML response.
        """
        srr_ids = []
        try:
            # Step 1: ESearch - Find SRA database IDs for this project
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            esearch_params = {
                'db': 'sra',
                'term': project_id,
                'retmode': 'json',
                'retmax': 1000
            }
            response = self.session.get(esearch_url, params=esearch_params, timeout=20)
            response.raise_for_status()
            esearch_data = response.json()

            idlist = esearch_data.get('esearchresult', {}).get('idlist', [])
            if not idlist:
                print(f"    ESearch found no SRA IDs for {project_id}")
                return []

            print(f"    ESearch found {len(idlist)} SRA records for {project_id}.")

            # Step 2: EFetch - Get the RunInfo (XML) for these SRA IDs
            sra_ids_str = ",".join(idlist)
            efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            efetch_params = {
                'db': 'sra',
                'id': sra_ids_str,
                'rettype': 'runinfo'
            }

            response = self.session.get(efetch_url, params=efetch_params, timeout=45)
            response.raise_for_status()

            # --- XML PARSING LOGIC ---
            try:
                tree = ET.fromstring(response.text)
                srr_ids = [run.text for run in tree.findall('.//Run')]

                if not srr_ids:
                    print(f"    XML received, but no <Run> tags found for {project_id}")

                return srr_ids

            except ET.ParseError as e:
                print(f"    Failed to parse XML response: {e}")
                print(f"    Raw text (first 1000 chars): {response.text[:1000]}")
                return []

        except Exception as e:
            print(f"    Error during E-utilities query for {project_id}: {e}")
            return []

    def get_all_srr_ids(self, gse_id):
        """
        Get all SRR run IDs for a GSE dataset.
        """
        print(f"  Finding SRA Project (SRP) or BioProject (PRJNA) IDs for {gse_id}...")

        project_ids = self.get_sra_project_ids_from_gse(gse_id)

        if not project_ids:
            print(f"  ✗ No SRA projects found for {gse_id}")
            return []

        print(f"  Found project(s): {', '.join(project_ids)}")
        print(f"  Retrieving all SRR runs via E-utilities API...")

        all_srr_ids = []
        for project_id in project_ids:
            srr_ids = self.get_srr_from_eutils(project_id)
            if srr_ids:
                print(f"    {project_id}: Found {len(srr_ids)} SRR run(s).")
                all_srr_ids.extend(srr_ids)
            else:
                print(f"    {project_id}: No SRR runs found via API.")
            time.sleep(0.5)

        final_ids = list(set(all_srr_ids))
        if not final_ids:
            print(f"  ✗ No SRA data found for {gse_id}")

        return final_ids

    # --- THIS FUNCTION IS MODIFIED ---
    def download_fastq_with_prefetch(self, srr_id, output_dir):
        """
        Download FASTQ using prefetch and fasterq-dump.
        This version prints live progress from the SRA tools.
        """
        if not self.sra_available:
            print(f"    ⊙ SRA Toolkit not available. Run manually:")
            print(f"       prefetch {srr_id} && fasterq-dump {srr_id} -O {output_dir}")
            return False

        try:
            print(f"    Prefetching {srr_id}...")
            prefetch_cmd = ['prefetch', srr_id, '-O', output_dir]

            # Run prefetch, allowing its output (progress) to print to the console
            # We remove capture_output=True and text=True
            result = subprocess.run(prefetch_cmd, timeout=600)

            if result.returncode != 0:
                print(f"    ✗ Prefetch failed with return code {result.returncode}")
                # Don't return False here, fasterq-dump can often handle downloading
                # if prefetch fails.

            print(f"    Converting to FASTQ...")
            fasterq_cmd = ['fasterq-dump', srr_id, '-O', output_dir, '-e', '4', '--split-files']

            # Run fasterq-dump, allowing its output (progress) to print to the console
            result = subprocess.run(fasterq_cmd, timeout=1800)

            if result.returncode != 0:
                print(f"    ✗ fasterq-dump failed with return code {result.returncode}")
                return False

            print(f"    ✓ Successfully downloaded {srr_id}")

            # Clean up SRA file to save space
            local_sra_file = os.path.join(output_dir, srr_id + '.sra')
            if os.path.exists(local_sra_file):
                try:
                    os.remove(local_sra_file)
                except:
                    pass

            sra_dir = os.path.join(output_dir, srr_id)
            if os.path.isdir(sra_dir) and not os.listdir(sra_dir):
                try:
                    os.rmdir(sra_dir)
                except:
                    pass

            return True

        except subprocess.TimeoutExpired:
            print(f"    ✗ Download timeout for {srr_id}")
            return False
        except Exception as e:
            print(f"    ✗ Error downloading {srr_id}: {e}")
            return False

    def download_gse_fastq(self, gse_id, max_samples=None):
        """
        Download FASTQ files for a GSE dataset
        """
        print(f"\n{'=' * 60}")
        print(f"Processing: {gse_id}")
        print(f"{'=' * 60}")

        gse_dir = os.path.join(self.output_dir, gse_id)
        os.makedirs(gse_dir, exist_ok=True)

        result = {
            'gse_id': gse_id, 'status': 'pending', 'total_runs': 0,
            'downloaded_runs': 0, 'failed_runs': 0, 'srr_ids': []
        }

        srr_ids = self.get_all_srr_ids(gse_id)

        if not srr_ids:
            result['status'] = 'no_sra_data'
            return result

        result['total_runs'] = len(srr_ids)
        result['srr_ids'] = srr_ids

        srr_to_download = srr_ids
        if max_samples:
            srr_to_download = srr_ids[:max_samples]
            print(f"  Limiting to first {max_samples} SRA runs")

        print(f"\n  Found {len(srr_ids)} total SRA runs.")
        print(f"  Attempting to download {len(srr_to_download)} runs.")
        print(f"  Output directory: {gse_dir}")
        print()

        for idx, srr_id in enumerate(srr_to_download, 1):
            print(f"  [{idx}/{len(srr_to_download)}] Downloading {srr_id}...")

            fastq_files = [
                os.path.join(gse_dir, f"{srr_id}.fastq"),
                os.path.join(gse_dir, f"{srr_id}_1.fastq"),
                os.path.join(gse_dir, f"{srr_id}_2.fastq")
            ]

            if any(os.path.exists(f) for f in fastq_files):
                print(f"    ⊙ Already exists, skipping")
                result['downloaded_runs'] += 1
                continue

            success = self.download_fastq_with_prefetch(srr_id, gse_dir)

            if success:
                result['downloaded_runs'] += 1
            else:
                result['failed_runs'] += 1
            time.sleep(1)

        if result['downloaded_runs'] == len(srr_to_download):
            result['status'] = 'success'
        elif result['downloaded_runs'] > 0:
            result['status'] = 'partial'
        else:
            result['status'] = 'failed'

        print(f"\n  Summary for {gse_id}:")
        print(f"    Total runs found: {result['total_runs']}")
        print(f"    Attempted to download: {len(srr_to_download)}")
        print(f"    ✓ Successfully downloaded: {result['downloaded_runs']}")
        print(f"    ✗ Failed: {result['failed_runs']}")
        return result

    def download_from_csv(self, csv_path, max_samples_per_gse=None):
        """
        Download FASTQ data for all GSE IDs in a CSV file
        """
        try:
            df = pd.read_csv(csv_path)
            print(f"Found {len(df)} datasets in {csv_path}")
        except Exception as e:
            print(f"Error reading CSV: {e}")
            return None

        if 'gse_id' not in df.columns:
            print("Error: CSV must have 'gse_id' column")
            return None

        results = []
        for idx, row in df.iterrows():
            gse_id = row['gse_id']
            result = self.download_gse_fastq(gse_id, max_samples=max_samples_per_gse)
            results.append(result)

            progress_df = pd.DataFrame(results)
            progress_file = os.path.join(self.output_dir, 'download_progress.csv')
            progress_df.to_csv(progress_file, index=False)
            time.sleep(2)

        print(f"\n{'=' * 60}")
        print("DOWNLOAD SUMMARY")
        print(f"{'=' * 60}")

        total_runs = sum(r['total_runs'] for r in results)
        downloaded_runs = sum(r['downloaded_runs'] for r in results)
        failed_runs = sum(r['failed_runs'] for r in results)

        print(f"Total datasets processed: {len(results)}")
        print(f"Total SRA runs found: {total_runs}")
        print(f"  ✓ Successfully downloaded: {downloaded_runs}")
        print(f"  ✗ Failed: {failed_runs}")
        print(f"\nData saved in: {self.output_dir}")
        return results


# Main execution
if __name__ == "__main__":
    csv_path = "./data/breast_cancer_scRNA_seq_GSE.csv"
    downloader = GEODataDownloader(output_dir="./data")
    results = downloader.download_from_csv(csv_path, max_samples_per_gse=None)

    print("\n" + "=" * 60)
    print("Download complete! Check ./data/<GSE_ID>/ for FASTQ files")
    print("=" * 60)