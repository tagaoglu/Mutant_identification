# Pipeline to identify potential causative mutations for a certain phenotype

Below you can find descriptions of the scripts in this repository:

* **all required modules**              
module.sh

* **Performing quality control**      
doQC.sh

* **Sequence genome, remap (and identify variants)**       
remap_script.sh        
BQSR-script.sh

* **Finalize SNV calling**        
finalize_SNV_calling.sh      

* **Identify variants (and Annotation)**       
annotvcfATH.sh


## Summary of the goals of this project:
- Remapping of reads of 6 mutants and one parental strain on a reference genome (TAIR10 Arabidopsis thaliana)
- QC of the reads with FastQC and MultiQC
- Trimming and cleaning the reads if necessary
- Remapping the reads with bwa
- Base recalibration and SNV calling with GATK recommended pipeline to get a joint VCF file
- SNV effect annotation of the VCF file with SNPeffector
- Selection of the interesting SNVs by applying several filters with SnpSift:
  - keep SNVs in coding regions,
  - non-synonymous,
  - homozygotes for the mutants and absent in the parental strain,
  - originating from EMS mutations (G->A and C->T)
  
A list of known genes was checked for the presence of variants    

The final goal is to identify new SNVs explaining the Low Vitamine E (LOVE) phenotype of the 6 mutants.
