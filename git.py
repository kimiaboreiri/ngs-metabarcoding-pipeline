import os
import subprocess
from glob import glob


# PATHS

raw_data = " "
trimmed = "/mnt/c/Users/Kimia/Documents/fall/trimmed/trimmed"
fastqc_out = "/mnt/c/Users/Kimia/Documents/fall/trimmed/trimmed/fastqc"

os.makedirs(trimmed, exist_ok=True)
os.makedirs(fastqc_out, exist_ok=True)
'''

# PRIMERS

# step 1 – reverse-complement primers (3')
RC_FORWARD = "CTAGAGGAGCCTGTTCTA"
RC_REVERSE = "GGGGTATCTAATCCCAGT"

# step 2 – main primers (5')
FORWARD_PRIMER = "ACTGGGATTAGATACCCC"
REVERSE_PRIMER = "TAGAACAGGCTCCTCTAG"


# FIND FASTQ FILES
=
r1_files = sorted(
    glob(os.path.join(raw_data, "**", "*_R1_001.fastq.gz"), recursive=True)
)

print(f"Found {len(r1_files)} R1 files")

if len(r1_files) == 0:
    raise RuntimeError("No R1 FASTQ files found. Check paths.")


# FASTQC 

print("Running FastQC...")
for fq in glob(os.path.join(raw_data, "**", "*.fastq.gz"), recursive=True):
    subprocess.run(
        ["fastqc", fq, "-o", fastqc_out],
        check=True
    )


# CUTADAPT

print("Running cutadapt...")

for r1 in r1_files:
    base = r1.replace("_R1_001.fastq.gz", "")
    r2 = base + "_R2_001.fastq.gz"

    if not os.path.exists(r2):
        print("WARNING: Missing R2 for", r1)
        continue

    sample = os.path.basename(base)

    step1_r1 = os.path.join(trimmed, sample + "_R1_step1.fastq.gz")
    step1_r2 = os.path.join(trimmed, sample + "_R2_step1.fastq.gz")

    final_r1 = os.path.join(trimmed, sample + "_R1_trimmed.fastq.gz")
    final_r2 = os.path.join(trimmed, sample + "_R2_trimmed.fastq.gz")

    print("Processing:", sample)

    # ---- STEP 1: remove RC primers (3')
    subprocess.run([
        "cutadapt",
        "-a", RC_FORWARD,
        "-A", RC_REVERSE,
        "-o", step1_r1,
        "-p", step1_r2,
        r1, r2
    ], check=True)

    # ---- STEP 2: remove main primers (5')
    subprocess.run([
        "cutadapt",
        "-g", "^" + FORWARD_PRIMER,
        "-G", "^" + REVERSE_PRIMER,
        "--discard-untrimmed",
        "-o", final_r1,
        "-p", final_r2,
        step1_r1, step1_r2
    ], check=True)

    os.remove(step1_r1)
    os.remove(step1_r2)
'''

# RUN DADA2

print("Running DADA2...")
subprocess.run(["Rscript", "dada2.R"], check=True)

print("Pipeline finished successfully ")

