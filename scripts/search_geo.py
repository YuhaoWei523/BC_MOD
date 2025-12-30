import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
from urllib.parse import urlencode


class GEOCrawler:
    """Crawler for retrieving scRNA-seq dataset IDs from GEO related to breast cancer"""

    def __init__(self):
        self.base_url = "https://www.ncbi.nlm.nih.gov/geo/browse/"
        self.search_url = "https://www.ncbi.nlm.nih.gov/gds/"
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })

    def search_geo(self, query, max_results=100):
        """
        Search GEO database for datasets matching the query

        Args:
            query: Search terms (e.g., "breast cancer single cell RNA-seq")
            max_results: Maximum number of results to retrieve

        Returns:
            List of dataset IDs (GSE numbers)
        """
        datasets = []
        retstart = 0
        retmax = 20  # GEO returns 20 results per page

        print(f"Searching GEO for: {query}")

        while len(datasets) < max_results:
            params = {
                'term': query,
                'retstart': retstart,
                'retmax': retmax
            }

            url = f"{self.search_url}?{urlencode(params)}"

            try:
                response = self.session.get(url, timeout=10)
                response.raise_for_status()

                soup = BeautifulSoup(response.content, 'html.parser')

                # Find all dataset links (GSE accessions)
                links = soup.find_all('a', href=True)

                page_datasets = []
                for link in links:
                    href = link.get('href', '')
                    if 'acc=GSE' in href:
                        gse_id = href.split('acc=')[1].split('&')[0]
                        if gse_id.startswith('GSE') and gse_id not in datasets:
                            page_datasets.append(gse_id)

                if not page_datasets:
                    print(f"No more results found at offset {retstart}")
                    break

                datasets.extend(page_datasets)
                print(f"Found {len(page_datasets)} datasets (Total: {len(datasets)})")

                retstart += retmax
                time.sleep(0.5)  # Be polite to the server

            except Exception as e:
                print(f"Error during search: {e}")
                break

        return datasets[:max_results]

    def get_dataset_info(self, gse_id):
        """
        Retrieve detailed information for a specific GSE dataset

        Args:
            gse_id: GEO Series accession (e.g., "GSE123456")

        Returns:
            Dictionary with dataset information
        """
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"

        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()

            soup = BeautifulSoup(response.content, 'html.parser')

            info = {
                'gse_id': gse_id,
                'title': '',
                'organism': '',
                'summary': '',
                'platform': '',
                'samples': 0
            }

            # Extract title
            title_tag = soup.find('td', text='Title')
            if title_tag:
                info['title'] = title_tag.find_next_sibling('td').get_text(strip=True)

            # Extract organism
            organism_tag = soup.find('td', text='Organism')
            if organism_tag:
                info['organism'] = organism_tag.find_next_sibling('td').get_text(strip=True)

            # Extract summary
            summary_tag = soup.find('td', text='Summary')
            if summary_tag:
                info['summary'] = summary_tag.find_next_sibling('td').get_text(strip=True)

            # Extract platform
            platform_links = soup.find_all('a', href=lambda x: x and 'acc=GPL' in x)
            if platform_links:
                info['platform'] = platform_links[0].get_text(strip=True)

            # Count samples
            sample_links = soup.find_all('a', href=lambda x: x and 'acc=GSM' in x)
            info['samples'] = len(sample_links)

            return info

        except Exception as e:
            print(f"Error fetching info for {gse_id}: {e}")
            return None

    # --- THIS FUNCTION HAS BEEN UPDATED ---
    def filter_scrna_datasets(self, datasets):
        """
        Filter datasets to find those likely to be scRNA-seq AND about human breast cancer

        Args:
            datasets: List of GSE IDs

        Returns:
            List of filtered GSE IDs with metadata
        """

        # Keywords for scRNA-seq (can be in title or summary)
        scrna_keywords = [
            'single cell', 'single-cell', 'scrna', 'sc-rna',
            '10x', 'drop-seq', 'smart-seq', 'single-nucleus', 'snrna',
            'single cell transcriptome', 'single-cell transcriptome'
        ]

        # Keywords for topic (MUST be in title for high relevancy)
        topic_keywords = [
            'breast cancer', 'breast tumor', 'mammary carcinoma', 'dcis'
        ]

        filtered_data = []

        for i, gse_id in enumerate(datasets):
            print(f"Processing {i + 1}/{len(datasets)}: {gse_id}")

            info = self.get_dataset_info(gse_id)

            if not info:
                print(f"  ✗ Skipping {gse_id}: Could not fetch info.")
                continue

            # --- FILTER 1: Check Organism ---
            # This is the most important filter to remove mouse models
            if 'homo sapiens' not in info['organism'].lower():
                print(f"  ✗ Skipping {gse_id}: Organism is {info['organism']}, not Homo sapiens.")
                continue

            # Prepare text
            title_text = info['title'].lower()
            full_text = (title_text + ' ' + info['summary']).lower()

            # --- FILTER 2: Check for scRNA-seq keywords (in title OR summary) ---
            is_scrna = any(keyword.lower() in full_text for keyword in scrna_keywords)

            if not is_scrna:
                print(f"  ✗ Skipping {gse_id}: No scRNA-seq keywords found.")
                continue

            # --- FILTER 3: Check for topic keywords (MUST be in title) ---
            # This ensures the dataset is *about* breast cancer, not just
            # mentioning it in the summary.
            is_topic_in_title = any(keyword.lower() in title_text for keyword in topic_keywords)

            if is_topic_in_title:
                filtered_data.append(info)
                print(f"  ✓ scRNA-seq dataset found: {info['title'][:60]}...")
            else:
                print(f"  ✗ Skipping {gse_id}: Title does not contain breast cancer keywords.")

            time.sleep(0.3)  # Be polite to the server

        return filtered_data

    def search_multiple_queries(self, queries, max_results_per_query=50):
        """
        Search GEO with multiple queries and combine results

        Args:
            queries: List of search queries
            max_results_per_query: Maximum results per query

        Returns:
            List of unique GSE IDs
        """
        all_datasets = set()  # Use set to avoid duplicates

        for query in queries:
            print(f"\n{'=' * 60}")
            print(f"Searching with query: {query}")
            print(f"{'=' * 60}")

            datasets = self.search_geo(query, max_results=max_results_per_query)
            all_datasets.update(datasets)
            print(f"Unique datasets so far: {len(all_datasets)}")
            time.sleep(1)  # Delay between queries

        return list(all_datasets)


# Example usage
if __name__ == "__main__":
    import os

    crawler = GEOCrawler()

    # --- QUERIES HAVE BEEN UPDATED ---
    # Use NCBI's advanced search syntax to be more specific
    # We now filter for "Homo sapiens" at the search level
    queries = [
        "(breast cancer[Title/Abstract] OR breast tumor[Title/Abstract]) AND (scRNA-seq[Title/Abstract] OR single-cell[Title/Abstract]) AND (Homo sapiens[Organism])",
        "(mammary[Title/Abstract] OR DCIS[Title/Abstract]) AND (scRNA-seq[Title/Abstract] OR single-cell[Title/Abstract]) AND (Homo sapiens[Organism])",
        "breast cancer single cell transcriptome (Homo sapiens[Organism])",
        "breast tumor single nucleus (Homo sapiens[Organism])"
    ]
    csv_path = "./data/breast_cancer_scRNA_seq_GSE.csv"

    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)

    # Search with all queries and get unique datasets
    print("Starting multi-query search...")
    all_datasets = crawler.search_multiple_queries(queries, max_results_per_query=50)

    print(f"\n{'=' * 60}")
    print(f"Total unique datasets found: {len(all_datasets)}")
    print(f"{'=' * 60}\n")

    # Filter and get detailed info
    scrna_datasets = crawler.filter_scrna_datasets(all_datasets)

    # Save to CSV and remove duplicates based on GSE ID
    if scrna_datasets:
        df = pd.DataFrame(scrna_datasets)
        # Remove duplicates based on gse_id, keep first occurrence
        df = df.drop_duplicates(subset=['gse_id'], keep='first')
        df = df.sort_values('gse_id')  # Sort by GSE ID

        df.to_csv(csv_path, index=False)
        print(f"\n{'=' * 60}")
        print(f"Saved {len(df)} unique scRNA-seq datasets to {csv_path}")
        print(f"{'=' * 60}\n")

        # Display summary
        print("\nSummary of found datasets:")
        for i, (_, ds) in enumerate(df.iterrows(), 1):
            print(f"\n{i}. {ds['gse_id']}")
            print(f"   Title: {ds['title'][:70]}...")
            print(f"   Samples: {ds['samples']}")
    else:
        print("\nNo scRNA-seq datasets found matching the criteria.")