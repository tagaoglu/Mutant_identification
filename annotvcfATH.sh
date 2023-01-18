#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name="vcfannot"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH --mem=25G
#SBATCH --partition=pcourse80
#SBATCH -o vcfannot-output.txt
#SBATCH -e vcfannot-error.txt
#SBATCH -p pcourse80

#load modules
source /data/users/tagaoglu/BC7107_20/scripts/module.sh

#create and go to the TP directory
cd /data/users/$USER/3_genotyping

#the file vitE_allchr.vcf.gz must be in the curren dir
if [ ! -f "./vitaminE_allchr.vcf.gz" ]
then
        echo "missing vcf file: vitaminE_allchr.vcf.gz"
	exit 0
fi

#annotate vcf file using snpEff
#first if needed download and install snpEff locally
if [ ! -d "./snpEff" ]
then
        echo "downloading snpEff" 
        wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
        unzip snpEff_latest_core.zip
        rm snpEff_latest_core.zip
fi

# get the database for the Athaliana genome.
#java -Xmx4g -jar ./snpEff/snpEff.jar download Arabidopsis_thaliana

#annotate the VCF file
java -Xmx4g -jar ./snpEff/snpEff.jar -no-upstream -no-downstream Arabidopsis_thaliana vitaminE_allchr.vcf.gz > vitE_annot.vcf

#remove synonymous and intergenic variants 
cat vitE_annot.vcf | java -jar ./snpEff/SnpSift.jar filter "(( ANN[*].EFFECT != 'synonymous_variant') & ( ANN[*].EFFECT != 'intergenic_region'))"   > vitE_coding.vcf

# keep only variants that are found in less than 5 strains and where last genome is reference
cat vitE_coding.vcf | java -jar ./snpEff/SnpSift.jar filter  "((countHom()>2) & (countHet()=0) & (isRef(GEN[6]))) " > vitE_coding_notRef.vcf

# keep only mutations due to EMS
cat vitE_coding_notRef.vcf | java -jar ./snpEff/SnpSift.jar filter  "(REF='G' & ALT='A') | (REF='C' & ALT='T')" > vitE_coding_notRef_EMS.vcf

Part 2:
#Check presence of mutants in a list:

#get reference annotation and genome
ln -s /data/users/tagaoglu/BC7107_20/reference/ATH.gff3 Athal_annot.gff3
ln -s /data/users/tagaoglu/BC7107_20/reference/ATH.fa Athal_genome.fa
ln -s /data/users/tagaoglu/BC7107_20/reference/geneList.txt .

# extract bed file from gene list & annotations
grep -f geneList.txt Athal_annot.gff3 | grep mRNA | perl -ane '$F[8] =~ m/(AT\w+\.\d)\;/; $ID=$1; print "$F[0]\t$F[3]\t$F[4]\t\.\t$ID\n"' > geneList.bed

#get mutation falling into the gene list...
bedtools intersect -a vitE_coding_notRef_EMS.vcf -b geneList.bed > known_variants.txt

# Tip: modify before importing into excel
perl -pe 's/;ANN=/;\tANN=/g' vitE_coding_notRef_EMS.vcf > vitE_coding_notRef_EMS_forexcel.vcf  
perl -i -pe 's/FORMAT\t/ANN\tFORMAT\t/' vitE_coding_notRef_EMS_forexcel.vcf

# Visualize results with IGV
