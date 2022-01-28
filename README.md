
# **Neurodegeneration with brain iron accumulation (NBIA) and the risk for Parkinson’s disease (PD)**
 **Authors:** Sara Bandres Ciga, Pilar Alvarez Jerez, Raquel Duran etc
 
**Aim:** Neurodegeneration with brain iron accumulation (NBIA) represents a group of inherited heterogeneous neurodegenerative disorders characterized by iron accumulation and the presence of axonal spheroids in the basal ganglia and other brain areas. Given the rarity of the disease, it can often go unrecognized or misdiagnosed. In PD, iron accumulation is a cardinal feature of degenerating regions in the brain and seems to be a key player in mechanisms that precipitate cell death. The aim of this study was to further explore the relationship between NBIA associated genes and PD etiology. We examined whether a genetic burden of variants in NBIA genes could contribute to the risk of developing PD by performing gene-based and single variant association analyses.

## Structure of ReadMe
[1) Locating and extracting data for our genes of interest](https://github.com/neurogenetics/NBIA_PD/blob/main/README.md#1-locating-and-extracting-data)
For the purpose of this ReadMe, we will be working with WGS AMP_PD v 2.5 data. The whole analysis is repeated with WES data from UKBiobank and will be briefly shown at the end.
### 2) Association analyses
### 3) Annotation and finding variants of interest
Annotated variants will be filtered by coding variants and those reported in HGMD. 
### 4) Burden analyses
Running individual gene burdens for both all and coding variants, as well as a gene-set burden analysis.
### 5) Inspection of compound heterozygotes


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
*WDR45 was initially included but limited data for X chromosome

	# Working directory for analysis
	cd /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/
	# Modules needed
	module load plink
	module load annovar
	module load samtools
	module load rvtests

 Extract genes into plink files, add pheno and update sex
	 

    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 1 --from-bp 16985958 --to-bp 17011928 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_ATP13A2
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 2 --from-bp 171434217 --to-bp 171485052 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_DCAF17
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 3 --from-bp 149162410 --to-bp 149221829 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_CP
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 6 --from-bp 43605316 --to-bp 43629264 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_GTPBP2
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 16 --from-bp 70687439 --to-bp 70801160 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_VAC14
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 16 --from-bp 74712955 --to-bp 74774831 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed  --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_FA2H
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 17 --from-bp 42561467 --to-bp 42566277 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_COASY    
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 19 --from-bp 29698886 --to-bp 29715789 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_C19orf12
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 19 --from-bp 48965309 --to-bp 48966879 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_FTL
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 20 --from-bp 3888839 --to-bp 3929882 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_PANK2
    plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/euro_king_pca_v2.5_July2021/AMPv2.5_sampleQC_EURO --chr 22 --from-bp 38111495 --to-bp 38214778 --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/AMP_PLA2G6   

     # Pheno and sex files made from /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt

## 2) Association analyses

Fisher association and logistic regression

    plink --bfile AMP_ATP13A2 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_ATP13A2
    plink --bfile AMP_ATP13A2 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_ATP13A2 --ci 0.95
    
    plink --bfile AMP_C19orf12 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_C19orf12
    plink --bfile AMP_C19orf12 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_C19orf12 --ci 0.95
    
    plink --bfile AMP_COASY --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_COASY
    plink --bfile AMP_COASY --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_COASY --ci 0.95
    
    plink --bfile AMP_CP --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_CP
    plink --bfile AMP_CP --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_CP --ci 0.95
    
    plink --bfile AMP_DCAF17 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_DCAF17
    plink --bfile AMP_DCAF17 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_DCAF17 --ci 0.95
    
    plink --bfile AMP_FA2H --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_FA2H
    plink --bfile AMP_FA2H --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_FA2H --ci 0.95
    
    plink --bfile AMP_FTL --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_FTL
    plink --bfile AMP_FTL --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_FTL --ci 0.95
    
    plink --bfile AMP_GTPBP2 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_GTPBP2
    plink --bfile AMP_GTPBP2 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_GTPBP2 --ci 0.95
    
    plink --bfile AMP_PANK2 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_PANK2
    plink --bfile AMP_PANK2 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_PANK2 --ci 0.95
    
    plink --bfile AMP_PLA2G6 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_PLA2G6
    plink --bfile AMP_PLA2G6 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_PLA2G6 --ci 0.95
    
    plink --bfile AMP_VAC14 --fisher --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --out AMP_VAC14
    plink --bfile AMP_VAC14 --logistic --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt \
    --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
    --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --out AMP_VAC14 --ci 0.95

## 3) Annotation and variants of interest

        plink --bfile AMP_ATP13A2 --recode 'vcf-fid' --out AMP_ATP13A2
    bgzip AMP_ATP13A2.vcf
    tabix -f -p vcf AMP_ATP13A2.vcf.gz
    
    table_annovar.pl AMP_ATP13A2.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_ATP13A2.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_ATP13A2.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_ATP13A2.annovar.hg38_multianno.txt > AMP_ATP13A2.trimmed.annotation.txt
    
    
    plink --bfile AMP_C19orf12 --recode 'vcf-fid' --out AMP_C19orf12
    bgzip AMP_C19orf12.vcf
    tabix -f -p vcf AMP_C19orf12.vcf.gz
    
    table_annovar.pl AMP_C19orf12.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_C19orf12.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_C19orf12.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_C19orf12.annovar.hg38_multianno.txt > AMP_C19orf12.trimmed.annotation.txt
    
    
    plink --bfile AMP_COASY --recode 'vcf-fid' --out AMP_COASY
    bgzip AMP_COASY.vcf
    tabix -f -p vcf AMP_COASY.vcf.gz
    
    table_annovar.pl AMP_COASY.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_COASY.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_COASY.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_COASY.annovar.hg38_multianno.txt > AMP_COASY.trimmed.annotation.txt
    
    
    plink --bfile AMP_CP --recode 'vcf-fid' --out AMP_CP
    bgzip AMP_CP.vcf
    tabix -f -p vcf AMP_CP.vcf.gz
    
    table_annovar.pl AMP_CP.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_CP.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_CP.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_CP.annovar.hg38_multianno.txt > AMP_CP.trimmed.annotation.txt
    
    
    
    
    plink --bfile AMP_DCAF17 --recode 'vcf-fid' --out AMP_DCAF17
    bgzip AMP_DCAF17.vcf
    tabix -f -p vcf AMP_DCAF17.vcf.gz
    
    table_annovar.pl AMP_DCAF17.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_DCAF17.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_DCAF17.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_DCAF17.annovar.hg38_multianno.txt > AMP_DCAF17.trimmed.annotation.txt
    
    
    
    plink --bfile AMP_FA2H --recode 'vcf-fid' --out AMP_FA2H
    bgzip AMP_FA2H.vcf
    tabix -f -p vcf AMP_FA2H.vcf.gz
    
    table_annovar.pl AMP_FA2H.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_FA2H.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_FA2H.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_FA2H.annovar.hg38_multianno.txt > AMP_FA2H.trimmed.annotation.txt
    
    
    
    
    plink --bfile AMP_FTL --recode 'vcf-fid' --out AMP_FTL
    bgzip AMP_FTL.vcf
    tabix -f -p vcf AMP_FTL.vcf.gz
    
    table_annovar.pl AMP_FTL.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_FTL.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_FTL.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_FTL.annovar.hg38_multianno.txt > AMP_FTL.trimmed.annotation.txt
    
    
    
    plink --bfile AMP_GTPBP2 --recode 'vcf-fid' --out AMP_GTPBP2
    bgzip AMP_GTPBP2.vcf
    tabix -f -p vcf AMP_GTPBP2.vcf.gz
    
    table_annovar.pl AMP_GTPBP2.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_GTPBP2.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_GTPBP2.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_GTPBP2.annovar.hg38_multianno.txt > AMP_GTPBP2.trimmed.annotation.txt
    
    
    plink --bfile AMP_PANK2 --recode 'vcf-fid' --out AMP_PANK2
    bgzip AMP_PANK2.vcf
    tabix -f -p vcf AMP_PANK2.vcf.gz
    
    table_annovar.pl AMP_PANK2.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_PANK2.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_PANK2.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_PANK2.annovar.hg38_multianno.txt > AMP_PANK2.trimmed.annotation.txt
    
    
    
    plink --bfile AMP_PLA2G6 --recode 'vcf-fid' --out AMP_PLA2G6
    bgzip AMP_PLA2G6.vcf
    tabix -f -p vcf AMP_PLA2G6.vcf.gz
    
    table_annovar.pl AMP_PLA2G6.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_PLA2G6.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_PLA2G6.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_PLA2G6.annovar.hg38_multianno.txt > AMP_PLA2G6.trimmed.annotation.txt
    
    
    plink --bfile AMP_VAC14 --recode 'vcf-fid' --out AMP_VAC14
    bgzip AMP_VAC14.vcf
    tabix -f -p vcf AMP_VAC14.vcf.gz
    
    table_annovar.pl AMP_VAC14.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
    --thread 16 \
    -out AMP_VAC14.annovar \
    -remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
    
    head -1 AMP_VAC14.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct AMP_VAC14.annovar.hg38_multianno.txt > AMP_VAC14.trimmed.annotation.txt



For further analysis we filtered annotated variants per gene by exonic or splicing and removed synonymous.  After we merged those remaining variants with known variants from HGMD via R. 

Example R script

    module load R/3.6.0 
    R
    require(data.table)
    library(dplyr)


    A <- fread("AMP_FA2H.trimmed.annotation.forR.txt", header = T)
    B <- fread("/Users/alvarezjerezp2/Desktop/HGMD_info/HGMD_FA2H.txt", header=T)
    colnames(B) <- c("CHR", "Start","ID","REF","ALT","HGMD_Description","HGVS_nucleotide", "HGVS_protein","Variant_class","Variant_type", "Reported_phenotype")
    total <- merge(A,B, by="Start")
    write.table(total, file="AMP_FA2H_HGMD.txt",  quote=F, row.names = F, sep = "\t") 

Created a combined variants_of_interest_AMP.txt list with: chr start stop gene_name, with no header. Using this list for association analysis and compound heterozygous insight. Separately creating a one column list of those genes names that have variants in the variants_of_interest_AMP.txt file.

Example

    # head variants_of_interest_AMP.txt
    1	17000494	17000494	ATP13A2
    1	17005754	17005754	ATP13A2
    2	171449955	171449955	DCAF17
    2	171481070	171481070	DCAF17
    3	149178609	149178609	CP
    3	149183513	149183513	CP
	
	# head gene_list.txt
	ATP13A2
    DCAF17
    CP
    FA2H
    C19orf12
    FTL
    PANK2
    PLA2G6

Loop to extract variants of interest from plink files

    cat genes_list.txt | while read LINE
    do
    plink --bfile AMP_$LINE --extract range variants_of_interest_AMP.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt \
    --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out TEMP_$LINE
    done

Create list of temps to merge together minus first temp (in our case  no TEMP_ATP13A2)

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

    plink --bfile TEMP_ATP13A2 --merge-list temp_outputs.txt --make-bed --out ALL_AMP

Association with variants of interest

    plink --bfile ALL_AMP --assoc --out  ALL_AMP_freqs

Recoding to .raw file for downstream analysis

    plink -bfile ALL_AMP --recode A --out ALL_AMP
    
# 4) Burden analyses
Here we will be using the rvtests package to run analyses on all variants and coding variants.

First we will create our coding files from the annotation. These are different than the variants_of_interest_AMP.txt file as these have not been joined with HGMD data.

    awk '$6=="exonic" {print}' AMP_ATP13A2.trimmed.annotation.txt > AMP_ATP13A2.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_ATP13A2.trimmed.annotation.txt > AMP_ATP13A2.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_ATP13A2.trimmed.annotation.txt > AMP_ATP13A2.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_ATP13A2.trimmed.annotation.exonic.variants.txt AMP_ATP13A2.trimmed.annotation.splicing.variants.txt AMP_ATP13A2.trimmed.annotation.exonicsplicing.variants.txt > AMP_ATP13A2.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_ATP13A2.trimmed.annotation.coding.variants.txt > AMP_ATP13A2.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_ATP13A2 --extract range AMP_ATP13A2.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_ATP13A2_CODING
    plink --bfile AMP_ATP13A2_CODING --recode 'vcf-fid' --out AMP_ATP13A2_CODING
    
    bgzip AMP_ATP13A2_CODING.vcf
    tabix -f -p vcf AMP_ATP13A2_CODING.vcf.gz
    
    
    awk '$6=="exonic" {print}' AMP_C19orf12.trimmed.annotation.txt > AMP_C19orf12.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_C19orf12.trimmed.annotation.txt > AMP_C19orf12.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_C19orf12.trimmed.annotation.txt > AMP_C19orf12.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_C19orf12.trimmed.annotation.exonic.variants.txt AMP_C19orf12.trimmed.annotation.splicing.variants.txt AMP_C19orf12.trimmed.annotation.exonicsplicing.variants.txt > AMP_C19orf12.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_C19orf12.trimmed.annotation.coding.variants.txt > AMP_C19orf12.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_C19orf12 --extract range AMP_C19orf12.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_C19orf12_CODING
    plink --bfile AMP_C19orf12_CODING --recode 'vcf-fid' --out AMP_C19orf12_CODING
    
    bgzip AMP_C19orf12_CODING.vcf
    tabix -f -p vcf AMP_C19orf12_CODING.vcf.gz
    
    
    
    awk '$6=="exonic" {print}' AMP_COASY.trimmed.annotation.txt > AMP_COASY.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_COASY.trimmed.annotation.txt > AMP_COASY.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_COASY.trimmed.annotation.txt > AMP_COASY.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_COASY.trimmed.annotation.exonic.variants.txt AMP_COASY.trimmed.annotation.splicing.variants.txt AMP_COASY.trimmed.annotation.exonicsplicing.variants.txt > AMP_COASY.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_COASY.trimmed.annotation.coding.variants.txt > AMP_COASY.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_COASY --extract range AMP_COASY.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_COASY_CODING
    plink --bfile AMP_COASY_CODING --recode 'vcf-fid' --out AMP_COASY_CODING
    
    bgzip AMP_COASY_CODING.vcf
    tabix -f -p vcf AMP_COASY_CODING.vcf.gz
    
    
    
    awk '$6=="exonic" {print}' AMP_CP.trimmed.annotation.txt > AMP_CP.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_CP.trimmed.annotation.txt > AMP_CP.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_CP.trimmed.annotation.txt > AMP_CP.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_CP.trimmed.annotation.exonic.variants.txt AMP_CP.trimmed.annotation.splicing.variants.txt AMP_CP.trimmed.annotation.exonicsplicing.variants.txt > AMP_CP.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_CP.trimmed.annotation.coding.variants.txt > AMP_CP.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_CP --extract range AMP_CP.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_CP_CODING
    plink --bfile AMP_CP_CODING --recode 'vcf-fid' --out AMP_CP_CODING
    
    bgzip AMP_CP_CODING.vcf
    tabix -f -p vcf AMP_CP_CODING.vcf.gz
    
    
    
    
    
    awk '$6=="exonic" {print}' AMP_DCAF17.trimmed.annotation.txt > AMP_DCAF17.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_DCAF17.trimmed.annotation.txt > AMP_DCAF17.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_DCAF17.trimmed.annotation.txt > AMP_DCAF17.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_DCAF17.trimmed.annotation.exonic.variants.txt AMP_DCAF17.trimmed.annotation.splicing.variants.txt AMP_DCAF17.trimmed.annotation.exonicsplicing.variants.txt > AMP_DCAF17.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_DCAF17.trimmed.annotation.coding.variants.txt > AMP_DCAF17.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_DCAF17 --extract range AMP_DCAF17.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_DCAF17_CODING
    plink --bfile AMP_DCAF17_CODING --recode 'vcf-fid' --out AMP_DCAF17_CODING
    
    bgzip AMP_DCAF17_CODING.vcf
    tabix -f -p vcf AMP_DCAF17_CODING.vcf.gz
    
    
    
    awk '$6=="exonic" {print}' AMP_FA2H.trimmed.annotation.txt > AMP_FA2H.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_FA2H.trimmed.annotation.txt > AMP_FA2H.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_FA2H.trimmed.annotation.txt > AMP_FA2H.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_FA2H.trimmed.annotation.exonic.variants.txt AMP_FA2H.trimmed.annotation.splicing.variants.txt AMP_FA2H.trimmed.annotation.exonicsplicing.variants.txt > AMP_FA2H.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_FA2H.trimmed.annotation.coding.variants.txt > AMP_FA2H.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_FA2H --extract range AMP_FA2H.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_FA2H_CODING
    plink --bfile AMP_FA2H_CODING --recode 'vcf-fid' --out AMP_FA2H_CODING
    
    bgzip AMP_FA2H_CODING.vcf
    tabix -f -p vcf AMP_FA2H_CODING.vcf.gz
    
    
    
    awk '$6=="exonic" {print}' AMP_FTL.trimmed.annotation.txt > AMP_FTL.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_FTL.trimmed.annotation.txt > AMP_FTL.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_FTL.trimmed.annotation.txt > AMP_FTL.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_FTL.trimmed.annotation.exonic.variants.txt AMP_FTL.trimmed.annotation.splicing.variants.txt AMP_FTL.trimmed.annotation.exonicsplicing.variants.txt > AMP_FTL.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_FTL.trimmed.annotation.coding.variants.txt > AMP_FTL.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_FTL --extract range AMP_FTL.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_FTL_CODING
    plink --bfile AMP_FTL_CODING --recode 'vcf-fid' --out AMP_FTL_CODING
    
    bgzip AMP_FTL_CODING.vcf
    tabix -f -p vcf AMP_FTL_CODING.vcf.gz
    
    
    
    
    awk '$6=="exonic" {print}' AMP_GTPBP2.trimmed.annotation.txt > AMP_GTPBP2.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_GTPBP2.trimmed.annotation.txt > AMP_GTPBP2.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_GTPBP2.trimmed.annotation.txt > AMP_GTPBP2.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_GTPBP2.trimmed.annotation.exonic.variants.txt AMP_GTPBP2.trimmed.annotation.splicing.variants.txt AMP_GTPBP2.trimmed.annotation.exonicsplicing.variants.txt > AMP_GTPBP2.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_GTPBP2.trimmed.annotation.coding.variants.txt > AMP_GTPBP2.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_GTPBP2 --extract range AMP_GTPBP2.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_GTPBP2_CODING
    plink --bfile AMP_GTPBP2_CODING --recode 'vcf-fid' --out AMP_GTPBP2_CODING
    
    bgzip AMP_GTPBP2_CODING.vcf
    tabix -f -p vcf AMP_GTPBP2_CODING.vcf.gz
    
    
    
    awk '$6=="exonic" {print}' AMP_PANK2.trimmed.annotation.txt > AMP_PANK2.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_PANK2.trimmed.annotation.txt > AMP_PANK2.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_PANK2.trimmed.annotation.txt > AMP_PANK2.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_PANK2.trimmed.annotation.exonic.variants.txt AMP_PANK2.trimmed.annotation.splicing.variants.txt AMP_PANK2.trimmed.annotation.exonicsplicing.variants.txt > AMP_PANK2.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_PANK2.trimmed.annotation.coding.variants.txt > AMP_PANK2.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_PANK2 --extract range AMP_PANK2.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_PANK2_CODING
    plink --bfile AMP_PANK2_CODING --recode 'vcf-fid' --out AMP_PANK2_CODING
    
    bgzip AMP_PANK2_CODING.vcf
    tabix -f -p vcf AMP_PANK2_CODING.vcf.gz
    
    
    
    awk '$6=="exonic" {print}' AMP_PLA2G6.trimmed.annotation.txt > AMP_PLA2G6.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_PLA2G6.trimmed.annotation.txt > AMP_PLA2G6.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_PLA2G6.trimmed.annotation.txt > AMP_PLA2G6.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_PLA2G6.trimmed.annotation.exonic.variants.txt AMP_PLA2G6.trimmed.annotation.splicing.variants.txt AMP_PLA2G6.trimmed.annotation.exonicsplicing.variants.txt > AMP_PLA2G6.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_PLA2G6.trimmed.annotation.coding.variants.txt > AMP_ATP13A2.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_PLA2G6 --extract range AMP_PLA2G6.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_PLA2G6_CODING
    plink --bfile AMP_PLA2G6_CODING --recode 'vcf-fid' --out AMP_PLA2G6_CODING
    
    bgzip AMP_PLA2G6_CODING.vcf
    tabix -f -p vcf AMP_PLA2G6_CODING.vcf.gz
    
    
    
    
    awk '$6=="exonic" {print}' AMP_VAC14.trimmed.annotation.txt > AMP_VAC14.trimmed.annotation.exonic.variants.txt
    awk '$6=="splicing" {print}' AMP_VAC14.trimmed.annotation.txt > AMP_VAC14.trimmed.annotation.splicing.variants.txt
    awk '$6=="exonic;splicing" {print}' AMP_VAC14.trimmed.annotation.txt > AMP_VAC14.trimmed.annotation.exonicsplicing.variants.txt
    cat AMP_VAC14.trimmed.annotation.exonic.variants.txt AMP_VAC14.trimmed.annotation.splicing.variants.txt AMP_VAC14.trimmed.annotation.exonicsplicing.variants.txt > AMP_VAC14.trimmed.annotation.coding.variants.txt
    awk '{print $1" "$2" "$2" "$7}' AMP_VAC14.trimmed.annotation.coding.variants.txt > AMP_VAC14.trimmed.annotation.coding.variants.SNPs.txt
    
    plink --bfile AMP_VAC14 --extract range AMP_VAC14.trimmed.annotation.coding.variants.txt --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out AMP_VAC14_CODING
    plink --bfile AMP_VAC14_CODING --recode 'vcf-fid' --out AMP_VAC14_CODING
    
    bgzip AMP_VAC14_CODING.vcf
    tabix -f -p vcf AMP_VAC14_CODING.vcf.gz

Need to make pheno file with rvtests requirements, MAT and PAT columns.

    cp /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt pheno_rvtests_AMP.txt
    # Added PAT and MAT col with 0 values
    head pheno_rvtests_AMP.txt
      # FID	IID	PAT	MAT	SEX	AGE_ANALYSIS	PC1	PC2	PC3	PC4	PC5	PD_PHENO
    # 001_10 001_10	0	0	1	62	-0.04517551	-0.006875955	0.01393995	0.04467011	-0.02931715	1

Run the rvtests analysis as a job. Running for all variant at MAF 0.01 and 0.03 and coding variants for MAF 0.01 and 0.03

    sbatch --mem=200g --cpus-per-task=10 --time=48:00:00 rvtests_batch_AMP.sh
   
   Example of what it looks like inside rvtests_batch_AMP.sh
   #!/bin/bash

module load rvtests

    rvtest --noweb --inVcf AMP_ATP13A2.vcf.gz  --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_ATP13A2_PD_BURDEN_maf003
    rvtest --noweb --inVcf AMP_ATP13A2.vcf.gz --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_ATP13A2_PD_BURDEN_maf001
    rvtest --noweb --inVcf AMP_ATP13A2_CODING.vcf.gz --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_ATP13A2_PD_BURDEN_maf003_CODING
    rvtest --noweb --inVcf AMP_ATP13A2_CODING.vcf.gz --pheno pheno_rvtests_AMP.txt --pheno-name PD_PHENO --covar pheno_rvtests_AMP.txt  --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_ATP13A2_PD_BURDEN_maf001_CODING
    # This example only has for ATP13A2, repeat commands for all other genes


Next running gene-set burden analysis for all variants and coding variants.

Like earlier, need to combine gene files. For all variants, combine the files from the very beginning (AMP_ATP13A2).

    vim gene_all_outputs.txt
    # Has AMP_$GENE list for all genes except ATP13A2
    plink --bfile AMP_ATP13A2 --merge-list gene_all_outputs.txt --make-bed --out GENE_SET_ALL
    plink --bfile GENE_SET_ALL --update-sex  /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out GENE_SET_ALL
    module load samtools
    plink --bfile GENE_SET_ALL --recode 'vcf-fid' --out GENE_SET_ALL
    
    bgzip GENE_SET_ALL.vcf
    tabix -f -p vcf GENE_SET_ALL.vcf.gz

Do same thing for coding files created above.

    vim coding_all_outputs_AMP.txt
    # Has AMP_$GENE_CODING list for all genes except ATP13A2.
    plink --bfile AMP_ATP13A2_CODING --merge-list coding_all_outputs_AMP.txt --make-bed --out GENE_SET_CODING_AMP
    plink --bfile GENE_SET_CODING_AMP --update-sex /data/CARD/projects/GBA_PILAR/NBIA/AMP_2.5_goodfiles/update_sex.txt --pheno /data/CARD/projects/GBA_PILAR/NBIA/pheno_ampv2.5.txt --make-bed --out GENE_SET_CODING_AMP

    plink --bfile GENE_SET_CODING_AMP --recode 'vcf-fid' --out GENE_SET_CODING_AMP
    
    bgzip GENE_SET_CODING_AMP.vcf
    tabix -f -p vcf GENE_SET_CODING_AMP.vcf.gz

   Run burden analysis
  

     sh run_burden_pathways_NBIA_AMP.sh GENE_SET_ALL.vcf.gz 0.05
     sh run_burden_pathways_NBIA_AMP.sh GENE_SET_CODING_AMP.vcf.gz 0.05
    # Look inside run_burden_pathways_NBIA_AMP.sh
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
Python script from Sara BC.

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


Compound heterozygotes were done for all variants/gene. Need to later go into results and only count those alt alleles that are within same gene. 
Another option is to do this compound heterozygotes analysis gene by gene using TEMP_$GENE files created earlier.

