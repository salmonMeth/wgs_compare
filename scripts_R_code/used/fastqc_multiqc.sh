#runs FastQC and MultiQC on the fasta files and saves the resulting reports in new directories


# in the `cd /scratch/project_2019524/PinkSalmon_lowinput` commands
#change the paths to the paths of folders containing the fq.gz files

  
module load bio-apps/v202603
module load fastqc/0.12.1
module load multiqc/1.35

#for lowinput
cd /scratch/project_2019524/PinkSalmon_lowinput
mkdir -p qc/fastqc
mkdir -p qc/multiqc
fastqc *.fq.gz --outdir qc/fastqc --threads 8
multiqc qc/fastqc --outdir qc/multiqc

#for lowinput P100
cd /scratch/project_2019524/PinkSalmon_lowinput_P100
mkdir -p qc/fastqc
mkdir -p qc/multiqc
fastqc *.fq.gz --outdir qc/fastqc --threads 8
multiqc qc/fastqc --outdir qc/multiqc

#for normal input
cd /scratch/project_2019524/PinkSalmon_normalinput
mkdir -p qc/fastqc
mkdir -p qc/multiqc
fastqc *.fq.gz --outdir qc/fastqc --threads 8
multiqc qc/fastqc --outdir qc/multiqc
