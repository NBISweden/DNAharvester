import os
import argparse
import pandas as pd
from tqdm import tqdm
from Bio import Entrez
import requests
import zipfile
import io

# Configure your NCBI email
Entrez.email = "your_email@example.com"

def clean_fasta_sequence(fasta_file):
    """Modify sequence headers: remove ', complete genome' and replace spaces with underscores."""
    try:
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
        with open(fasta_file, 'w') as f:
            for line in lines:
                if line.startswith('>'):
                    header = line.strip()
                    header = header.replace(", complete genome", "")
                    header = header.replace(" ", "_")
                    f.write(header + '\n')
                else:
                    f.write(line)
        return True
    except Exception as e:
        print(f"[ERROR] Could not clean {fasta_file}: {e}")
    return False


def download_nuccore_fasta(acc, outdir):
    """Download FASTA from NCBI nuccore using Entrez."""
    try:
        handle = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text")
        seq_data = handle.read()
        handle.close()
        if seq_data.strip():
            outfile = os.path.join(outdir, f"{acc}.fasta")
            with open(outfile, "w") as f:
                f.write(seq_data)
            clean_fasta_sequence(outfile)
            return outfile
    except Exception as e:
        print(f"[ERROR] Could not fetch {acc} via Entrez: {e}")
    return None


def download_genome_gca_gcf(acc, outdir):
    """Download genome FASTA using NCBI Datasets API for GCA/GCF accessions, unzip and extract fasta."""
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
    try:
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with zipfile.ZipFile(io.BytesIO(response.content)) as z:
                for file in z.namelist():
                    if file.endswith(".fna") or file.endswith(".fasta"):
                        fasta_data = z.read(file)
                        outfile = os.path.join(outdir, f"{acc}.fasta")
                        with open(outfile, "wb") as f:
                            f.write(fasta_data)
                        clean_fasta_sequence(outfile)
                        return outfile
            print(f"[ERROR] No FASTA found in archive for {acc}")
        else:
            print(f"[ERROR] Failed to fetch {acc}, status {response.status_code}")
    except Exception as e:
        print(f"[ERROR] Exception for {acc}: {e}")
    return None


def append_to_previous_db(previous_db, fasta_file):
    """Append a new FASTA sequence into an existing database file."""
    try:
        with open(fasta_file, "r") as f_in, open(previous_db, "a") as f_out:
            f_out.write("\n")
            f_out.write(f_in.read())
        return True
    except Exception as e:
        print(f"[ERROR] Could not append {fasta_file} to {previous_db}: {e}")
    return False


def main():
    parser = argparse.ArgumentParser(description="Download FASTA sequences from NCBI.")
    parser.add_argument("--table", required=True, help="Input table file (CSV/TSV) with column 'Accessions'")
    parser.add_argument("--outpath", required=True, help="Output directory for FASTA files")
    parser.add_argument("--previous_db", required=False, help="Path to an existing FASTA database (.fna) to append new sequences")
    args = parser.parse_args()

    infile = args.table
    outdir = args.outpath
    previous_db = args.previous_db

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Read table (auto-detect CSV/TSV)
    if infile.endswith(".csv"):
        df = pd.read_csv(infile)
    else:
        df = pd.read_csv(infile, sep="\t")

    if "Accessions" not in df.columns:
        raise ValueError("Input table must have a column named 'Accessions'")

    # Iterate with tqdm
    for acc in tqdm(df["Accessions"].dropna().unique(), desc="Downloading FASTA"):
        if acc.startswith("GCA_") or acc.startswith("GCF_"):
            fasta_file = download_genome_gca_gcf(acc, outdir)
        else:
            fasta_file = download_nuccore_fasta(acc, outdir)

        if fasta_file:
            print(f"[OK] {acc} saved at {fasta_file}.")
            if previous_db:
                appended = append_to_previous_db(previous_db, fasta_file)
                if appended:
                    print(f"[OK] {acc} appended to {previous_db}.")
                    # Remove temporary fasta file after appending
                    try:
                        os.remove(fasta_file)
                        print(f"[CLEANUP] Removed temporary file {fasta_file}.")
                    except Exception as e:
                        print(f"[ERROR] Could not remove {fasta_file}: {e}")
        else:
            print(f"[FAILED] {acc} not retrieved.")


if __name__ == "__main__":
    main()
