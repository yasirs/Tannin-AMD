import os
import requests
import shutil

# Configuration
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
OUTPUT_DIR = os.path.join(BASE_DIR, "data/external/geo/GSE99287")
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99287/suppl/"

FILES_TO_DOWNLOAD = [
    "GSE99287_RPE_ATACSeq_peak_counts.txt.gz",
    "GSE99287_peak_counts_annotation.txt.gz"
]

def download_file(filename):
    url = BASE_URL + filename
    local_path = os.path.join(OUTPUT_DIR, filename)
    
    if os.path.exists(local_path):
        print(f"File already exists: {local_path}")
        return

    print(f"Downloading {url}...")
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Successfully downloaded {filename}")
    except Exception as e:
        print(f"Error downloading {filename}: {e}")

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    for f in FILES_TO_DOWNLOAD:
        download_file(f)

if __name__ == "__main__":
    main()
