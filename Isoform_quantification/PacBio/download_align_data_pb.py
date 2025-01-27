import os

### Downloading, Alignment
def create_prg_file(url, down_path, fasta_path, sam_path, bam_path, ref, prg_file_path):
   
    header = f"#!/bin/bash\n" + \
    "module load minimap2/2.17\n" + \
    "module load samtools/1.9\n" + \
    f"minimap2 -ax splice -uf -C5 {ref} {fasta_path} > {sam_path}\n" + \
    f"samtools view -bS {sam_path} > {bam_path}\n"
    #f"wget -O {down_path} {url}\n"
    #f"samtools fastq {down_path} > {fasta_path}\n"


    with open(prg_file_path, "w") as f:
        f.write(header)
    return prg_file_path



def create_slurm_file(prg_file_path, job_name, slurm_file_path):

    # header = f"#!/bin/bash\n" + \
    # f"#SBATCH --job-name={job_name}                     # Set the job name\n" + \
    # "#SBATCH --mem=20G                                  # Request 20G\n" + \
    # f"#SBATCH --output={log_dir}/out_{job_name}.%j      # Send stdout/err to\n" + \
    # "#SBATCH --mail-user=sraghavendra@nygenome.org      # Where to send mail\n" + \
    # f"{prg_file_path}"
    header = f"#!/bin/bash\n" + \
    f"#SBATCH --job-name={job_name}                     # Set the job name to \"JobExample1\"\n" + \
    "#SBATCH --partition=pe2                            # Partition Name\n" + \
    "#SBATCH --mail-type=END,FAIL                       # Mail events (NONE, BEGIN, END, FAIL, ALL)\n" + \
    "#SBATCH --mem=50G                                  # Request 20G\n" + \
    f"#SBATCH --output={log_dir}/out_{job_name}.%j      # Send stdout/err to\n" + \
    "#SBATCH --mail-user=sraghavendra@nygenome.org      # Where to send mail\n" + \
    f"{prg_file_path}"

    with open (slurm_file_path, "w") as f:
        f.write(header)
    return slurm_file_path


def get_aln_file_name(day, rep, ext):
    file_name = f"aln_{day}{rep}.{ext}"
    return file_name

def get_reads_file_name(day, rep, ext):
    file_name = f"flnc_{day}{rep}.{ext}"
    return file_name

def get_prg_file_name(day, rep):
    file_name = f"prg_file_{day}{rep}.sh"
    return file_name

def get_slurm_file_name(day, rep):
    file_name = f"slurm_{day}{rep}.sh"
    return file_name

reference = '/gpfs/commons/home/sraghavendra/PacBio/reference/transcriptome.fna'
down_dir = '/gpfs/commons/home/sraghavendra/PacBio/reads/long'
sam_dir = '/gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome'
bam_dir = '/gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome'
log_dir = '/gpfs/commons/home/sraghavendra/PacBio/logs/long'
scripts_path = '/gpfs/commons/home/sraghavendra/PacBio/scripts/long'

day_rep = [
#    ('0', '1'),
#    ('0', '2'),
    ('0', '3'),
    ('1', '1'),
    ('2', '1'),
    ('3', '1'),
    ('3', '2'),
    ('4', '1'),
    ('5', '1'), 
    ('5', '2'),
    ('5', '3')
]


def gen_combination():


    for day, rep in day_rep:                                 

        prg_file_path = os.path.join(scripts_path, get_prg_file_name(day, rep))
        slurm_file_path = os.path.join(scripts_path, get_slurm_file_name(day, rep))

        job_name = f"aln_pb_{day}{rep}"

        url = f"https://downloads.pacbcloud.com/public/dataset/MAS-Seq-bulk-2023-WTC11/day{day}-rep{rep}/2-FLNC/flnc.bam"
        down_path = os.path.join(down_dir, get_reads_file_name(day, rep, 'bam'))
        fasta_path = os.path.join(down_dir, get_reads_file_name(day, rep, 'fastq'))
        sam_path = os.path.join(sam_dir, get_aln_file_name(day, rep, 'sam'))
        bam_path = os.path.join(bam_dir, get_aln_file_name(day, rep, 'bam'))
        ref = reference

        create_prg_file(url, down_path, fasta_path, sam_path, bam_path, ref, prg_file_path=prg_file_path)
        create_slurm_file(prg_file_path=prg_file_path, job_name=job_name, slurm_file_path=slurm_file_path)

        os.system(f"chmod u+x {prg_file_path}")
        os.system(f"chmod u+x {slurm_file_path}")
        os.system(f"sbatch {slurm_file_path}")


def main():
    gen_combination()

if __name__ == "__main__":
    main()
