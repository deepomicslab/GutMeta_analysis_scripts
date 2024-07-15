#!/bin/bash
#SBATCH --job-name=gutmeta
#SBATCH --cpus-per-task=30
#SBATCH --time=24:00:00
#SBATCH --output=/data2/platform/gutmeta_v2_platform/task_logs/output/output_%j.output
#SBATCH --error=/data2/platform/gutmeta_v2_platform/task_logs/error/error_%j.error
#SBATCH --mem=50G

# module load  GCCcore/11.2.0 Python/3.9.6 
#module load GCC/11.2.0  OpenMPI/4.1.1
#module load MMseqs2/13-45111
# module load scikit-bio/0.5.6-Python-3.9.6
echo 'module load complete'

#export PATH=/home/platform/phage_db/tools/prodigal-gv:$PATH
#export EGGNOG_DATA_DIR=/home/platform/phage_db/phage_data/data/tools_data
#export PATH=/home/platform/phage_db/tools/eggnog-mapper:/home/platform/phage_db/tools/eggnog-mapper/eggnog-mapper/bin:"$PATH"


#mkdir -p $2
#cd /home/platform/phage_db/phage_api/workspace/analysis_script/annotation_v2
#prodigal -i $1 -f gff -o $2/gene.gff3 -d $2/gene.fna -a $2/gene.faa -p meta
#emapper.py -i $2/gene.faa --output $2/emapper_out -m diamond --cpu 30
#mkdir $2/iterate_annot
#python iterate_annot.py $2/gene.faa $2/emapper_out.emapper.annotations $2/iterate_annot 100 30
#python complete_gff.py $2/gene.gff3 $2/emapper_out.emapper.annotations $2/iterate_annot $2/sequence.gff3 $2/acc_list.txt
#eval "$(/apps/software/Miniconda3/4.9.2/bin/conda shell.bash hook)"
#/apps/software/Miniconda3/4.9.2/bin/conda activate /data2/platform/.conda/envs/gutmeta
cd /data2/platform/gutmeta_v2_platform/analysis_scripts/test_script
echo $1, $2, $3
python main.py --abdf=$1 --ann=$2 --output=$3
echo 'task complete'
