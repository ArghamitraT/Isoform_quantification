import os

### Downloading, Alignment
def create_prg_file(orig_path, ds_path1, ds_path2, ds_path3, per, prg_file_path):
   
    header = f"#!/bin/bash\n" + \
    "module load minimap2/2.17\n" + \
    "module load samtools/1.9\n" + \
    f"samtools view -b -s 0.{per} {orig_path} > {ds_path1}\n" + \
    f"samtools view -b -s 1.{per} {orig_path} > {ds_path2}\n" + \
    f"samtools view -b -s 2.{per} {orig_path} > {ds_path3}\n"

    with open(prg_file_path, "w") as f:
        f.write(header)
    return prg_file_path



def create_slurm_file(prg_file_path, job_name, slurm_file_path):

    header = f"#!/bin/bash\n" + \
    f"#SBATCH --job-name={job_name}                     # Set the job name to \"JobExample1\"\n" + \
    "#SBATCH --partition=pe2                            # Partition Name\n" + \
    "#SBATCH --mail-type=END,FAIL                       # Mail events (NONE, BEGIN, END, FAIL, ALL)\n" + \
    "#SBATCH --mem=100G                                 # Request 20G\n" + \
    f"#SBATCH --output={log_dir}/out_{job_name}.%j      # Send stdout/err to\n" + \
    "#SBATCH --mail-user=sraghavendra@nygenome.org      # Where to send mail\n" + \
    f"{prg_file_path}"

    with open (slurm_file_path, "w") as f:
        f.write(header)
    return slurm_file_path


def get_aln_file_name(day, rep, ext):
    file_name = f"long_aln_{day}{rep}.{ext}"
    return file_name

def get_ds_aln_file_name(day, rep, num, per):
    file_name = f"ds_{per}_num{num}_aln_{day}{rep}_long.bam"
    return file_name

def get_prg_file_name(day, rep, per):
    file_name = f"prg_file_{per}_{day}{rep}.sh"
    return file_name

def get_slurm_file_name(day, rep, per):
    file_name = f"slurm_{per}_{day}{rep}.sh"
    return file_name

orig_dir = '/gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome'
ds_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/other_days'
log_dir = '/gpfs/commons/home/sraghavendra/PacBio/logs/long/downsample'
scripts_path = '/gpfs/commons/home/sraghavendra/PacBio/scripts/long/downsample'

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

        for per in [5, 10]:                                 

            prg_file_path = os.path.join(scripts_path, get_prg_file_name(day, rep, per))
            slurm_file_path = os.path.join(scripts_path, get_slurm_file_name(day, rep, per))

            job_name = f"ds_{per}_pb_{day}{rep}"

            orig_path = os.path.join(orig_dir, get_aln_file_name(day, rep, 'bam'))
            ds_path1 = os.path.join(ds_dir, get_ds_aln_file_name(day, rep, 1, per))
            ds_path2 = os.path.join(ds_dir, get_ds_aln_file_name(day, rep, 2, per))
            ds_path3 = os.path.join(ds_dir, get_ds_aln_file_name(day, rep, 3, per))

            create_prg_file(orig_path, ds_path1, ds_path2, ds_path3, per, prg_file_path=prg_file_path)
            create_slurm_file(prg_file_path=prg_file_path, job_name=job_name, slurm_file_path=slurm_file_path)

            os.system(f"chmod u+x {prg_file_path}")
            os.system(f"chmod u+x {slurm_file_path}")
            os.system(f"sbatch {slurm_file_path}")


def main():
    gen_combination()

if __name__ == "__main__":
    main()
