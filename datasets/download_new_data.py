import os
import urllib.request
import sys

def download_file(url, filepath):
    print(f"Downloading {url}...")
    print(f"Saving to {filepath}")
    
    if os.path.exists(filepath):
        print(f"File {filepath} already exists. Skipping download.")
        return

    def reporthook(blocknum, blocksize, totalsize):
        readsofar = blocknum * blocksize
        if totalsize > 0:
            percent = readsofar * 1e2 / totalsize
            sys.stdout.write(f"\r{percent:.2f}%  -  {readsofar / (1024*1024):.2f} MB / {totalsize / (1024*1024):.2f} MB")
            sys.stdout.flush()
        else:
            sys.stdout.write(f"\rRead {readsofar / (1024*1024):.2f} MB")
            sys.stdout.flush()

    try:
        urllib.request.urlretrieve(url, filepath, reporthook)
        print("\nDownload Complete!")
    except Exception as e:
        print(f"\nError downloading {url}: {e}")

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    ncrna_dir = os.path.join(base_dir, 'human_ncrna')
    
    if not os.path.exists(ncrna_dir):
        os.makedirs(ncrna_dir)

    print("="*60)
    print("      ASMSG PHASE 1: DATA ACQUISITION (HUMAN ONLY)")
    print("="*60)

    # 1. GENCODE Human lncRNA FASTA (Very fast, ~25MB)
    gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.lncRNA_transcripts.fa.gz"
    gencode_path = os.path.join(ncrna_dir, "gencode.v44.lncRNA_transcripts.fa.gz")
    
    # 2. Ensembl Human general ncRNA FASTA (Very fast, ~15MB)
    ensembl_url = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
    ensembl_path = os.path.join(ncrna_dir, "Homo_sapiens.GRCh38.ncrna.fa.gz")

    print("\n[1/3] Downloading GENCODE Human lncRNA sequences (~25 MB)...")
    download_file(gencode_url, gencode_path)
    
    print("\n[2/3] Downloading Ensembl Human general ncRNA sequences (~15 MB)...")
    download_file(ensembl_url, ensembl_path)

    print("\n[3/3] RNADisease v4.0 Download Instructions")
    print("-" * 60)
    print("Please download RNADisease v4.0 manually:")
    print("1. Go to: http://www.rnadisease.org/download")
    print("2. Scroll to 'Experimental Data' section.")
    print("3. Download the 'All RNA-disease' file.")
    print(f"4. Save it to: {os.path.join(base_dir, 'rnadisease', 'alldata.csv')}")
    print("="*60)
