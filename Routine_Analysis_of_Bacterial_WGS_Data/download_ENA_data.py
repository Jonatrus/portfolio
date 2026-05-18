import urllib.request
from pathlib import Path
import requests

OUTPUT_DIR = Path("data/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PROJECT = "PRJNA891285"

url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={PROJECT}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
response = requests.get(url)
lines = response.text.strip().split("\n")[1:]

for line in lines:
    run_accession, ftp_paths = line.split("\t")
    
    for ftp_path in ftp_paths.split(";"):
        url = f"ftp://{ftp_path}"
        filename = OUTPUT_DIR / ftp_path.split("/")[-1]
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, filename)

print("Done!")