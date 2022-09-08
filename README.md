

# **Exploring the Genetic and Genomic Connection Underlying Neurodegeneration with Brain Iron Accumulation and the Risk for Parkinson’s Disease**
 **Authors:** Pilar Alvarez Jerez, Sara Bandres Ciga, Anni Moore, Mary Makarious
 
 **Background:**
 - Neurodegeneration with brain iron accumulation (NBIA) represents a group of inherited heterogeneous neurodegenerative disorders characterized by iron accumulation and the presence of axonal spheroids in the basal ganglia and other brain areas. 
- In PD, iron accumulation is a cardinal feature of degenerating regions in the brain and seems to be a key player in mechanisms that precipitate cell death.

**Aims:** 
- To explore the relationship between NBIA associated genes and PD etiology. 
- Examine whether a genetic burden of variants in NBIA genes could contribute to the risk of developing PD by performing gene-based and single-variant association analyses.


## Structure of ReadMe
### [1) Locating and extracting data for our genes of interest](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#1-locating-and-extracting-data)

### [2) Annotation and finding variants of interest](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#3-annotation-and-variants-of-interest)

### [3) Association analyses](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#2-association-analyses-1)

### [4) Burden analyses](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#4-burden-analyses-1)


### [5) Inspection of compound heterozygotes](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#5-compound-heterozygotes)

### [6) Other data](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#6-other-data-1)



### Lets get started....


## 1) Locating and extracting data
Genes of interest (positions in hg38)

**_ATP13A2_**  1: 16,985,958-17,011,928  
**_DCAF17_** 2: 171,434,217-171,485,052
**_CP_**  3: 149,162,410-149,221,829  
**_GTPBP2_**  6: 43,605,316-43,629,264  
**_VAC14_**   16: 70,687,439-70,801,160 
**_FA2H_**  16: 74,712,955-74,774,831  
**_COASY_**. 17: 42,561,467-42,566,277  
**_C19orf12_**  19: 29,698,886-29,715,789  
**_FTL_**  19: 48,965,309-48,966,879  
**_PANK2_**  20: 3,888,839-3,929,882  
**_PLA2G6_**  22: 38,111,495-38,214,778  

For the purpose of this ReadMe, we will be working with WGS AMP_PD v 2.5 data. The whole analysis was repeated with WES data from UKBiobank.
Additionally, when appropriate, we will only be showing example code for one gene to keep it clean. The code can be run individually per gene of interest or in loops.
 
#### Working directory for analysis
	cd /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/
 Modules needed
```
module load plink
module load annovar
module load samtools
module load rvtests
```

 Extract genes into plink files, add pheno and update sex
	 

    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO \
    --chr 1 --from-bp 16985958 --to-bp 17011928 \ 
    --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --make-bed --out 

     # Pheno and sex files made from /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt



 
    
## 2) Annotation and variants of interest
 

	# Creating VCF
    plink --bfile AMP_ATP13A2 --recode 'vcf-fid' --out AMP_ATP13A2
    bgzip AMP_ATP13A2.vcf
	tabix -f -p vcf AMP_ATP13A2.vcf.gz
    
    # Annotating variants
    table_annovar.pl AMP_ATP13A2.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_ATP13A2.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    # Trim and format annotation
    head -1 AMP_ATP13A2.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_ATP13A2.annovar.hg38_multianno.txt > AMP_ATP13A2.trimmed.annotation.txt


For further analysis, we filtered annotated variants per gene by exonic/splicing and removed synonymous variants.  
Variants were then merged with known variants from HGMD via R. 

Example R script

    module load R/3.6.0 
    R
    require(data.table)
    library(dplyr)

	# Note the "forR" files are those filtered for coding variants
    A <- fread("AMP_ATP13A2.trimmed.annotation.forR.txt", header = T)
    B <- fread("/Users/alvarezjerezp2/Desktop/HGMD_info/HGMD_AT13A2.txt", header=T)
    colnames(B) <- c("CHR", "Start","ID","REF","ALT","HGMD_Description","HGVS_nucleotide", "HGVS_protein","Variant_class","Variant_type", "Reported_phenotype")
    total <- merge(A,B, by="Start")
    write.table(total, file="AMP_ATP13A2_HGMD.txt",  quote=F, row.names = F, sep = "\t") 


We then created a list with all the combined variants of interest for all genes. 
List format was as follows:

    # head variants_of_interest_AMP.txt
    1	17000494	17000494	ATP13A2
    1	17005754	17005754	ATP13A2
    2	171449955	171449955	DCAF17
    2	171481070	171481070	DCAF17
    3	149178609	149178609	CP
    3	149183513	149183513	CP

This list was used for association analyses and compound heterozygous insight. 

Separately creating a one column list of those genes names that have variants in the variants_of_interest_AMP.txt file to extract variants from plink files.

	# head gene_list.txt
	ATP13A2
    DCAF17
    CP
    FA2H
    C19orf12
    FTL
    PANK2
    PLA2G6

Loop to extract variants of interest from plink files.

    cat genes_list.txt | while read LINE
    do
    plink --bfile AMP_$LINE --extract range variants_of_interest_AMP.txt \
    --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --make-bed --out TEMP_$LINE
    done

Create list of temps to merge together minus first temp file (in our case  no TEMP_ATP13A2).

    # vim temp_outputs
    # head temp_outputs.txt
    TEMP_DCAF17
    TEMP_CP
    TEMP_FA2H
    TEMP_C19orf12
    TEMP_FTL
    TEMP_PANK2
    TEMP_PLA2G6
    
Merge plink files

    plink --bfile TEMP_ATP13A2 --merge-list temp_outputs.txt \
    --make-bed --out ALL_AMP

Recoding to .raw file for downstream analysis

    plink -bfile ALL_AMP --recode A --out ALL_AMP
    
   ## 3) Association analyses

Plink association with variants of interest

    plink --bfile ALL_AMP --assoc --out  ALL_AMP_freqs
Fisher association and logistic regression per gene

    # Fisher
    plink --bfile AMP_ATP13A2 --fisher \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
     --out AMP_ATP13A2
    # Logistic regression 
    plink --bfile AMP_ATP13A2 --logistic \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
    --out AMP_ATP13A2 --ci 0.95
# 4) Burden analyses
Here we will be using the rvtests package to run analyses on all variants and coding variants.

Running individual gene burdens for both all and coding variants, as well as a gene-set burden analysis.

First we will create our coding files from the annotation. These are different than the variants_of_interest_AMP.txt file as these have not been joined with HGMD data.

    # Pull out exonic and splicing variants from annotation files
    awk '$6=="exonic" {print}' AMP_ATP13A2.trimmed.annotation.txt > AMP_ATP13A2.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_ATP13A2.trimmed.annotation.txt > AMP_ATP13A2.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_ATP13A2.trimmed.annotation.txt > AMP_ATP13A2.trimmed.annotation.exonicsplicing.variants.txt
    # Joing all the coding variants
    cat AMP_ATP13A2.trimmed.annotation.exonic.variants.txt AMP_ATP13A2.trimmed.annotation.splicing.variants.txt AMP_ATP13A2.trimmed.annotation.exonicsplicing.variants.txt > AMP_ATP13A2.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_ATP13A2.trimmed.annotation.coding.variants.txt > AMP_ATP13A2.trimmed.annotation.coding.variants.SNPs.txt
    
    # Make coding plink bfiles
    plink --bfile AMP_ATP13A2 --extract range AMP_ATP13A2.trimmed.annotation.coding.variants.txt \
    --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --make-bed --out AMP_ATP13A2_CODING
    
    # Make coding VCF
    plink --bfile AMP_ATP13A2_CODING --recode 'vcf-fid' --out AMP_ATP13A2_CODING
    
    bgzip AMP_ATP13A2_CODING.vcf
    tabix -f -p vcf AMP_ATP13A2_CODING.vcf.gz
    
    
Need to make pheno file with rvtests requirements, which include MAT and PAT columns.

    cp /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    pheno_rvtests_AMP.txt
    
    # Added PAT and MAT col with 0 values
    head pheno_rvtests_AMP.txt
    # FID IID PAT MAT SEX AGE_ANALYSIS PC1 PC2 PC3 PC4 PC5 PD_PHENO
    # 001_10 001_10 0 0 1 62 -0.04517551 -0.006875955 0.01393995 0.04467011 -0.02931715 1

Example rvtest analyses for ATP13A2 at MAF 0.01 and 0.03
   
   

	module load rvtests
	
	# All variants, MAF 0.03
    rvtest --noweb --inVcf AMP_ATP13A2.vcf.gz \
    --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO \
    --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
	--kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
	--geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 \
	--out AMP_ATP13A2_PD_BURDEN_maf003
	
	# All variants, MAF 0.01
    rvtest --noweb --inVcf AMP_ATP13A2.vcf.gz \
    --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO \
    --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 \
    --out AMP_ATP13A2_PD_BURDEN_maf001
	
	# Coding variants, MAF 0.03
    rvtest --noweb --inVcf AMP_ATP13A2_CODING.vcf.gz \
    --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO \
    --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 \
    --out AMP_ATP13A2_PD_BURDEN_maf003_CODING
    
	# Coding variants, MAF 0.01
    rvtest --noweb --inVcf AMP_ATP13A2_CODING.vcf.gz \
    --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO \
    --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 \
    --out AMP_ATP13A2_PD_BURDEN_maf001_CODING
    
    # This example only has for ATP13A2, repeat commands for all other genes


Next running gene-set burden analysis for all variants and coding variants.

Like earlier, need to combine gene files. For all variants, combine the files from the very beginning (AMP_ATP13A2).

    vim gene_all_outputs.txt
    # Has AMP_$GENE list for all genes except ATP13A2
    plink --bfile AMP_ATP13A2 --merge-list gene_all_outputs.txt \
    --make-bed --out GENE_SET_ALL
    
    # Update sex and pheno for that joined list
    plink --bfile GENE_SET_ALL --update-sex  /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out GENE_SET_ALL
    
	# Make VCF
    plink --bfile GENE_SET_ALL --recode 'vcf-fid' --out GENE_SET_ALL
    bgzip GENE_SET_ALL.vcf
    tabix -f -p vcf GENE_SET_ALL.vcf.gz

Do same thing for coding files created above.

    vim coding_all_outputs_AMP.txt
    # Has AMP_$GENE_CODING list for all genes except ATP13A2.
    plink --bfile AMP_ATP13A2_CODING --merge-list coding_all_outputs_AMP.txt \
    --make-bed --out GENE_SET_CODING_AMP
    
    plink --bfile GENE_SET_CODING_AMP --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out GENE_SET_CODING_AMP
	
	plink --bfile GENE_SET_CODING_AMP --recode 'vcf-fid' --out GENE_SET_CODING_AMP
    bgzip GENE_SET_CODING_AMP.vcf
    tabix -f -p vcf GENE_SET_CODING_AMP.vcf.gz

   Run gene-set burden analysis
  

     sh run_burden_pathways_NBIA_AMP.sh GENE_SET_ALL.vcf.gz 0.05
     sh run_burden_pathways_NBIA_AMP.sh GENE_SET_CODING_AMP.vcf.gz 0.05
   Look inside run_burden_pathways_NBIA_AMP.sh
	     
	    # !/bin/sh
	    #sh run_burden_pathways_NBIA_AMP.sh GENE_SET_ALL.vcf.gz 0.05
	    module load rvtests
	    FILENAME=$1
	    MAF=$2
	    OUTNAME=${FILENAME/".vcf.gz"/""}
	    COV_NAME=pheno_rvtests_AMP.txt
	    rvtest --noweb --hide-covar --out ${OUTNAME}_freqUpper${MAF}_PATHWAY --kernel skato \
	    --inVcf /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/${FILENAME} \
	    --pheno ${COV_NAME} \
	    --pheno-name PD_PHENO \
	    --covar ${COV_NAME} \
	    --freqUpper $MAF \
	    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
	    --setFile /data/CARD/projects/GBA_PILAR/NBIA/set_file.txt
       
# 5) Compound heterozygotes

    python
    import pandas as pd
    
    # Load in genotype file and get genotypes as 2 (homoz for alternate allele), 1 (heterozygous) and 0 homoz for reference allele
    df = pd.read_csv('ALL_AMP.raw', sep = ' ')
    df.head()
    # Remove no phenotype samples
    df = df[df.PHENOTYPE != -9]
    df.head()
    df.shape
    # Making subset table of sample info
    info = df[['FID','IID','PAT','MAT','SEX','PHENOTYPE']]
    df.columns
    # Making subset table without info
    df_sub = df[['FID', 'IID', 'chr1:16988161:T:A_A',
           'chr1:16991787:G:A_A', 'chr1:17000494:G:A_A', 'chr1:17005754:G:A_A',
           'chr2:171449955:C:T_T', 'chr2:171481070:T:C_C', 'chr3:149178609:C:G_G',
           'chr3:149183513:C:T_T', 'chr3:149185366:G:A_A', 'chr3:149199783:G:A_A',
           'chr3:149212616:C:G_G', 'chr16:74740048:C:T_T', 'chr16:74740049:G:A_A',
           'chr16:74774524:C:T_T', 'chr16:74774662:G:C_C', 'chr19:29702747:T:C_C',
           'chr19:29702888:C:T_T', 'chr19:29702977:C:A_A', 'chr19:29702977:C:G_G',
           'chr19:48966341:G:T_T', 'chr20:3889237:A:T_T', 'chr20:3907997:A:G_G',
           'chr20:3910812:A:G_G', 'chr20:3918695:G:A_A', 'chr20:3918719:A:G_G',
           'chr22:38112212:A:C_C', 'chr22:38112541:G:A_A', 'chr22:38126390:T:C_C',
           'chr22:38132890:C:T_T', 'chr22:38132917:C:T_T', 'chr22:38132952:G:A_A',
           'chr22:38133007:G:A_A', 'chr22:38143275:C:T_T', 'chr22:38169240:T:C_C',
           'chr22:38169326:G:A_A', 'chr22:38169336:C:T_T', 'chr22:38169411:G:A_A',
           'chr22:38169426:T:C_C']]
    df_sub.head()
    df_sub.shape
    # Identifying patients homozygous for alt allele
    homo = df_sub[df_sub.values == 2]
    homo = homo.drop_duplicates()
    print(homo.shape)
    homo.head()
    homo.columns
    homo_info = pd.merge(homo,info, how='left', left_on=['FID','IID'], right_on = ['FID','IID'])
    homo_info = homo_info[['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE', 'chr1:16988161:T:A_A',
           'chr1:16991787:G:A_A', 'chr1:17000494:G:A_A', 'chr1:17005754:G:A_A',
           'chr2:171449955:C:T_T', 'chr2:171481070:T:C_C', 'chr3:149178609:C:G_G',
           'chr3:149183513:C:T_T', 'chr3:149185366:G:A_A', 'chr3:149199783:G:A_A',
           'chr3:149212616:C:G_G', 'chr16:74740048:C:T_T', 'chr16:74740049:G:A_A',
           'chr16:74774524:C:T_T', 'chr16:74774662:G:C_C', 'chr19:29702747:T:C_C',
           'chr19:29702888:C:T_T', 'chr19:29702977:C:A_A', 'chr19:29702977:C:G_G',
           'chr19:48966341:G:T_T', 'chr20:3889237:A:T_T', 'chr20:3907997:A:G_G',
           'chr20:3910812:A:G_G', 'chr20:3918695:G:A_A', 'chr20:3918719:A:G_G',
           'chr22:38112212:A:C_C', 'chr22:38112541:G:A_A', 'chr22:38126390:T:C_C',
           'chr22:38132890:C:T_T', 'chr22:38132917:C:T_T', 'chr22:38132952:G:A_A',
           'chr22:38133007:G:A_A', 'chr22:38143275:C:T_T', 'chr22:38169240:T:C_C',
           'chr22:38169326:G:A_A', 'chr22:38169336:C:T_T', 'chr22:38169411:G:A_A',
           'chr22:38169426:T:C_C']]
    homo_info.head()
    # Option to drop NA samples here 
    homo_info_control = homo_info.loc[homo_info['PHENOTYPE'] == 1]
    print(homo_info_control.shape)
    homo_info_case = homo_info.loc[homo_info['PHENOTYPE'] == 2]
    print(homo_info_case.shape)
    homo_info.to_csv('AMP_coding.homo_MAF001.txt', sep = '\t', index=False)
    homo_info_control.to_csv('AMP_coding_control.homo_MAF001.txt', sep = '\t', index=False)
    homo_info_case.to_csv('AMP_coding_case.homo_MAF001.txt', sep = '\t', index=False)
    # Identifying patients heterozygous for alt allele
    het = df_sub[df_sub.values == 1]
    het = het.drop_duplicates()
    print(het.shape)
    het.head()
    het.columns
    het['count'] = het [['chr1:16988161:T:A_A', 'chr1:16991787:G:A_A',
           'chr1:17000494:G:A_A', 'chr1:17005754:G:A_A', 'chr2:171449955:C:T_T',
           'chr2:171481070:T:C_C', 'chr3:149178609:C:G_G', 'chr3:149183513:C:T_T',
           'chr3:149185366:G:A_A', 'chr3:149199783:G:A_A', 'chr3:149212616:C:G_G',
           'chr16:74740048:C:T_T', 'chr16:74740049:G:A_A', 'chr16:74774524:C:T_T',
           'chr16:74774662:G:C_C', 'chr19:29702747:T:C_C', 'chr19:29702888:C:T_T',
           'chr19:29702977:C:A_A', 'chr19:29702977:C:G_G', 'chr19:48966341:G:T_T',
           'chr20:3889237:A:T_T', 'chr20:3907997:A:G_G', 'chr20:3910812:A:G_G',
           'chr20:3918695:G:A_A', 'chr20:3918719:A:G_G', 'chr22:38112212:A:C_C',
           'chr22:38112541:G:A_A', 'chr22:38126390:T:C_C', 'chr22:38132890:C:T_T',
           'chr22:38132917:C:T_T', 'chr22:38132952:G:A_A', 'chr22:38133007:G:A_A',
           'chr22:38143275:C:T_T', 'chr22:38169240:T:C_C', 'chr22:38169326:G:A_A',
           'chr22:38169336:C:T_T', 'chr22:38169411:G:A_A', 'chr22:38169426:T:C_C']].sum(axis=1)
    het.head()
    # Finding samples with more than one alt allele (compount het)
    het_comp = het.loc[het['count'] > 1]
    het_comp = het_comp.drop_duplicates()
    print(het_comp.shape)
    het_comp.head()
    new_df = pd.merge(het_comp, info,  how='left', left_on=['FID','IID'], right_on = ['FID','IID'])
    new_df = new_df[['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE', 'chr1:16988161:T:A_A',
           'chr1:16991787:G:A_A', 'chr1:17000494:G:A_A', 'chr1:17005754:G:A_A',
           'chr2:171449955:C:T_T', 'chr2:171481070:T:C_C', 'chr3:149178609:C:G_G',
           'chr3:149183513:C:T_T', 'chr3:149185366:G:A_A', 'chr3:149199783:G:A_A',
           'chr3:149212616:C:G_G', 'chr16:74740048:C:T_T', 'chr16:74740049:G:A_A',
           'chr16:74774524:C:T_T', 'chr16:74774662:G:C_C', 'chr19:29702747:T:C_C',
           'chr19:29702888:C:T_T', 'chr19:29702977:C:A_A', 'chr19:29702977:C:G_G',
           'chr19:48966341:G:T_T', 'chr20:3889237:A:T_T', 'chr20:3907997:A:G_G',
           'chr20:3910812:A:G_G', 'chr20:3918695:G:A_A', 'chr20:3918719:A:G_G',
           'chr22:38112212:A:C_C', 'chr22:38112541:G:A_A', 'chr22:38126390:T:C_C',
           'chr22:38132890:C:T_T', 'chr22:38132917:C:T_T', 'chr22:38132952:G:A_A',
           'chr22:38133007:G:A_A', 'chr22:38143275:C:T_T', 'chr22:38169240:T:C_C',
           'chr22:38169326:G:A_A', 'chr22:38169336:C:T_T', 'chr22:38169411:G:A_A',
           'chr22:38169426:T:C_C', 'count']]
    new_df.head()
    # Option here to drop NA samples
    print(new_df.shape)
    controls = new_df.loc[new_df['PHENOTYPE'] == 1]
    cases = new_df.loc[new_df['PHENOTYPE'] == 2]
    
    print(controls.shape)
    print(cases.shape)
    new_df.to_csv('AMP_coding_all.compound.hets_MAF001.txt', sep = '\t', index=True)
    new_df.head()
    cases.to_csv('AMP_coding_case.compound.hets_MAF001.txt', sep = '\t', index=True)
    controls.to_csv('AMP_coding_control.compound.hets_MAF001.txt', sep = '\t', index=True)


Compound heterozygotes were looked at for all variants/gene. Need to later go into results and only count those alt alleles that are within same gene. 



# 6) Other data

The analysis was then repeated with whole exome sequencing UKBiobank data. Follows same pattern as above.

Where to find data:
	 

    # My working directory 
    cd /data/CARD/projects/GBA_PILAR/NBIA/UKB
    # Split up variants by gene (hg38)
    # ATP13A2 1: 16,985,958-17,011,928	
    awk '$4 > 16985958' /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/UKBexomeOQFE_chr1.bim | awk '$4 < 17011928' > ATP13A2_variants_UKB.txt
        
    plink --bed /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/ukb23155_c1_b0_v1.bed \
    --bim /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/UKBexomeOQFE_chr1.bim \
    --fam /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/ukb23155_c1_b0_v1_s200632.fam \
    --extract ATP13A2_variants_UKB.txt --out ATP13A2_UKB --make-bed
    
    
    	
    # Covariate file found here: /data/CARD/projects/GBA_PILAR/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt 
	# To format pheno file
    awk '{print $1, $2, $9}' /data/CARD/projects/GBA_PILAR/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt > pheno_UKB_NBIA.txt

For transcriptomics script please contact us.



Done!! 
