#SBATCH --mail-type=fail
#SBATCH --job-name="QC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=25G
#SBATCH --partition=pcourse80

source /data/users/tagaoglu/BC7107_20/scripts/module.sh

#create and go to the TP directory
mkdir /data/users/${USER}/QCres
cd /data/users/${USER}/QCres

ln -s /data/users/tagaoglu/BC7107_20/QC/*R1*.fastq.gz .
ln -s /data/users/tagaoglu/BC7107_20/QC/*R2*.fastq.gz .
ln -s /data/users/tagaoglu/BC7107_20/QC/TruSeq3-SE .

for k in `ls -1 *.fastq.gz`; 
do fastqc -t 2 ${k};
done

trimmomatic SE -phred33 -threads 4 JK2*_R1*.fastq.gz JK2_cleaned.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51

trimmomatic PE -phred33 -threads 4 D-4*_R1*.fastq.gz D-4*_R2*.fastq.gz D-4_1trim.fastq.gz D-4_1unpaired.fastq.gz D-4_2trim.fastq.gz D-4_2unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

fastqc -t 2 JK2_cleaned.fastq.gz
fastqc -t 2 D-4_*trim.fastq.gz

multiqc .
