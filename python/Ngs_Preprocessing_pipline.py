import os
import subprocess
import sys

def run_command(command):
    """Run a shell command and handle errors."""
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"Error: {command} failed.")
        sys.exit(1)

def main(sample_id, genome_id):
    print("Pipe Start")

    # Setup directories
    SRA_DATA = "/mnt/c/Users/Hitesh/Desktop/Unix/SRA_data"
    file_path = os.path.join(SRA_DATA, f"{sample_id}.fastq")
    QC_output = "/mnt/c/Users/Hitesh/Desktop/Unix/QC-out"
    Trim_out = "/mnt/c/Users/Hitesh/Desktop/Unix/Trim-Out"
    Trim_out_html = os.path.join(Trim_out, "html")
    Trim_out_json = os.path.join(Trim_out, "json")
    Genom = os.path.join("/mnt/c/Users/Hitesh/Desktop/Unix/Genom", genome_id, "*.fna")
    Genom_out = os.path.join("/mnt/c/Users/Hitesh/Desktop/Unix/Genom", genome_id)
    Sam_out = "/mnt/c/Users/Hitesh/Desktop/Unix/SAM"
    Bam_out = os.path.join(Sam_out, "Bam")

    # Create directories if they don't exist
    os.makedirs(SRA_DATA, exist_ok=True)
    os.makedirs(QC_output, exist_ok=True)
    os.makedirs(Trim_out, exist_ok=True)
    os.makedirs(Trim_out_html, exist_ok=True)
    os.makedirs(Trim_out_json, exist_ok=True)
    os.makedirs(Sam_out, exist_ok=True)
    os.makedirs(Bam_out, exist_ok=True)

    print("Directories setup complete.")

    # Download data
    run_command(f"fastq-dump -O {SRA_DATA} {sample_id}")
    print("fastq-dump complete.")

    # Quality control
    run_command(f"fastqc -t 2 -o {QC_output} {file_path}")
    print("FastQC complete.")

    # Trimming
    run_command(f"fastp --detect_adapter_for_pe -i {file_path} -o {os.path.join(Trim_out, f'{sample_id}_trimmed.fastq')} -j {os.path.join(Trim_out_json, f'{sample_id}_json.json')} -h {os.path.join(Trim_out_html, f'{sample_id}_html.html')}")
    print("Trimming complete.")

    # Build genome index with bwa
    run_command(f"bwa index {Genom}")
    print("Genome indexing complete.")

    # Alignment with bwa
    run_command(f"bwa mem {Genom} {os.path.join(Trim_out, f'{sample_id}_trimmed.fastq')} > {os.path.join(Sam_out, f'{sample_id}.sam')}")
    print("Alignment complete.")

    # Convert SAM to BAM
    run_command(f"samtools view -S -b {os.path.join(Sam_out, f'{sample_id}.sam')} > {os.path.join(Bam_out, f'{sample_id}.bam')}")
    print("Bam conversion complete.")

    # Sort BAM
    run_command(f"samtools sort -o {os.path.join(Bam_out, f'{sample_id}_sort.bam')} {os.path.join(Bam_out, f'{sample_id}.bam')}")
    print("Bam sort complete.")

    # Index BAM
    run_command(f"samtools index {os.path.join(Bam_out, f'{sample_id}_sort.bam')}")
    print("Bam index complete.")

    # HTSeq count
    run_command(f"htseq-count --format bam {os.path.join(Bam_out, f'{sample_id}_sort.bam')} {os.path.join(Genom_out, 'genomic.gtf')} > {os.path.join(Bam_out, f'{sample_id}_hts.ham')}")
    print("htseq ham complete.")

    # Convert ham to txt
    with open(os.path.join(Bam_out, f"{sample_id}_hts.ham"), "r") as ham_file:
        with open(os.path.join(Bam_out, f"{sample_id}_hts.txt"), "w") as txt_file:
            txt_file.write(ham_file.read())

    print("Txt done")
    print("Pipe End")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <sample_id> <genome_id>")
        sys.exit(1)

    sample_id = sys.argv[1]
    genome_id = sys.argv[2]
    main(sample_id, genome_id)
