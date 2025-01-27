import os

### Downloading, Alignment
def create_prg_file(url1, url2, down_path_r1, down_path_r2, bam_path, ref_dir, prg_file_path):
   
    # header = f"#!/bin/bash\n" + \
    # "module load minimap2/2.17\n" + \
    # "module load samtools/1.9\n" + \
    # f"wget {url1} > {down_path_r1}\n" + \
    # f"wget {url2} > {down_path_r2}\n"
    # #f"minimap2 -ax sr -t4 {ref} {down_path_r1} {down_path_r2} > {sam_path}\n"
    # #f"samtools view -bS {sam_path} > {bam_path}\n"

    header = f"#!/bin/bash\n" + \
    "module load star/2.5.2a\n" + \
    "module load samtools/1.9\n" + \
    f"STAR --genomeDir {ref_dir} " + \
    "--runThreadN 6 " + \
    f"--readFilesIn {down_path_r1} {down_path_r2} " + \
    f"--outFileNamePrefix {bam_path} " + \
    "--outSAMtype BAM SortedByCoordinate\n"
    #f"samtools flagstat {bam_path} > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/spliced/report_01_short.txt"

    with open(prg_file_path, "w") as f:
        f.write(header)
    return prg_file_path



def create_slurm_file(prg_file_path, job_name, slurm_file_path):

    header = f"#!/bin/bash\n" + \
    f"#SBATCH --job-name={job_name}                     # Set the job name to \"JobExample1\"\n" + \
    "#SBATCH --partition=pe2                            # Partition Name\n" + \
    "#SBATCH --mail-type=END,FAIL                       # Mail events (NONE, BEGIN, END, FAIL, ALL)\n" + \
    "#SBATCH --mem=100G                                  # Request 50G\n" + \
    f"#SBATCH --output={log_dir}/out_{job_name}.%j      # Send stdout/err to\n" + \
    "#SBATCH --mail-user=sraghavendra@nygenome.org      # Where to send mail\n" + \
    f"{prg_file_path}"

    with open (slurm_file_path, "w") as f:
        f.write(header)
    return slurm_file_path


def get_aln_file_name(day, rep):
    file_name = f"aln_{day}{rep}_short"#.{ext}"
    return file_name

def get_reads_file_name(day, rep, pair, ext):
    file_name = f"ill_{day}{rep}_{pair}.{ext}"
    return file_name

def get_prg_file_name(day, rep):
    file_name = f"prg_file_ill_{day}{rep}.sh"
    return file_name

def get_slurm_file_name(day, rep):
    file_name = f"slurm_ill_{day}{rep}.sh"
    return file_name

ref_dir = '/gpfs/commons/home/sraghavendra/PacBio/reference/transc/'
down_dir = '/gpfs/commons/home/sraghavendra/PacBio/reads/short/'
#sam_dir = '/gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome'
#bam_dir = 'gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/Illumina_Liz'
bam_dir = '/gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome'
log_dir = '/gpfs/commons/home/sraghavendra/PacBio/logs/short/'
scripts_path = '/gpfs/commons/home/sraghavendra/PacBio/scripts/short'


file_pairs = [
    ['20230519_Day0-1_LS11085_S130_R1_001.fastq.gz',
    '20230519_Day0-1_LS11085_S130_R2_001.fastq.gz'],
    ['20230519_Day1-1_LS11085_S131_R1_001.fastq.gz',
    '20230519_Day1-1_LS11085_S131_R2_001.fastq.gz'],
    ['20230519_Day2-1_LS11085_S132_R1_001.fastq.gz',
    '20230519_Day2-1_LS11085_S132_R2_001.fastq.gz'],
    ['20230519_Day4-1_LS11085_S134_R1_001.fastq.gz',
    '20230519_Day4-1_LS11085_S134_R2_001.fastq.gz'],
    ['20230519_Day5-1_LS11085_S135_R1_001.fastq.gz',
    '20230519_Day5-1_LS11085_S135_R2_001.fastq.gz'],
    ['20230610_Day0-2_LS11170_S32_R1_001.fastq.gz',
    '20230610_Day0-2_LS11170_S32_R2_001.fastq.gz'],
    ['20230610_Day0-3_LS11170_S33_R1_001.fastq.gz',
    '20230610_Day0-3_LS11170_S33_R2_001.fastq.gz'],
    ['20230610_Day3-2_LS11170_S34_R1_001.fastq.gz',
    '20230610_Day3-2_LS11170_S34_R2_001.fastq.gz'],
    ['20230610_Day5-2_LS11170_S35_R1_001.fastq.gz',
    '20230610_Day5-2_LS11170_S35_R2_001.fastq.gz'],
    ['20230610_Day5-3_LS11170_S36_R1_001.fastq.gz',
    '20230610_Day5-3_LS11170_S36_R2_001.fastq.gz'],
    ['Day3-1_LS11085_R1_001.fastq.gz',
    'Day3-1_LS11085_R2_001.fastq.gz']
]


def gen_combination():


    for file in file_pairs:  

        parts = file[0].split('_')
        if file==file_pairs[-1]:
            dr = parts[0]
            day, rep = dr.split('-')
            day = day[-1]
        else:
            dr = parts[1]
            day, rep = dr.split('-')
            day = day[-1]

        prg_file_path = os.path.join(scripts_path, get_prg_file_name(day, rep))
        slurm_file_path = os.path.join(scripts_path, get_slurm_file_name(day, rep))

        job_name = f"down_aln_ill_{day}{rep}"

        url1 = f"https://downloads.pacbcloud.com/public/dataset/MAS-Seq-bulk-2023-WTC11/illumina/{file[0]}"
        url2 = f"https://downloads.pacbcloud.com/public/dataset/MAS-Seq-bulk-2023-WTC11/illumina/{file[1]}"
        down_path1 = os.path.join(down_dir, get_reads_file_name(day, rep, 'R1', 'fastq'))
        down_path2 = os.path.join(down_dir, get_reads_file_name(day, rep, 'R2', 'fastq'))
        #sam_path = os.path.join(sam_dir, get_aln_file_name(day, rep, 'sam'))
        bam_path = os.path.join(bam_dir, get_aln_file_name(day, rep))

        create_prg_file(url1, url2, down_path1, down_path2, bam_path, ref_dir, prg_file_path=prg_file_path)
        create_slurm_file(prg_file_path=prg_file_path, job_name=job_name, slurm_file_path=slurm_file_path)

        os.system(f"chmod u+x {prg_file_path}")
        os.system(f"chmod u+x {slurm_file_path}")
        os.system(f"sbatch {slurm_file_path}")


def main():
    gen_combination()

if __name__ == "__main__":
    main()
