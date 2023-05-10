#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --job-name=cellranger
#SBATCH --mail-type=END
#SBATCH --mail-user=ernouha@gmail.com
#SBATCH --output=/scratch/slurm_cellranger_count_virilis.out




# ########################################################
# # Making a reference genome compatible with cellranger #
# ########################################################

# #Needs to be done only once, then the same ref genome can be used each time
# 

cd /scratch/

cellranger mkref \
--genome=Drosophila_Virilis_DVirRS2_extended1 \
--fasta=/scratch/GCF_003285735.1_DvirRS2_genomic.fna \
--genes=/scratch/GCF_003285735.1_DvirRS2_genomic_extended_stranded_reads_210323.gtf

# #Specify the folder where the analysis should be done (and the output folder stored)
analysis_folder=/scratch/
echo $analysis_folder
#Specify the ID of your sample (will be the name of the folder where the output of cell ranger is stored)
sample_id='Virilis_larva1_output'

#Specify where your reference transcriptome is stored
ref_transcriptome=/scratch/Drosophila_Virilis_DVirRS2_extended1
echo $ref_transcriptome 
#Specify where your reads are stored
reads_location=/shared/10Xdata/sequencing_data/Dvirilis/fastq
echo $reads_location
#Name of your sample
sample_name='Dvirilis_larva1'

#Number of cells you expect to recover
cell_number='15000'

# ###########
# # Mapping #
# ###########

# #module purge
# #module load cellranger/2.1.0


# R
# #Moves to the folder where the analysis should be done
cd $analysis_folder


#CHANGED TO CELLRANGER V7 ----- include introns default true ==> parameter removed

# #Does the mapping
# #--nosecondary: optional. Remove it if you want cell ranger to do the secondary analysis (dimensionality reduction, clustering and visualization)
cellranger count --id=$sample_id \
                   --transcriptome=$ref_transcriptome \
                   --fastqs=$reads_location \
                   --sample=$sample_name \
                   --expect-cells=$cell_number \
                   --nosecondary
