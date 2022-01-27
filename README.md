
# **Neurodegeneration with brain iron accumulation (NBIA) and the risk for Parkinson’s disease (PD)**
 **Authors:** Sara Bandres Ciga, Pilar Alvarez Jerez, Raquel Duran etc
 
**Aim:** Neurodegeneration with brain iron accumulation (NBIA) represents a group of inherited heterogeneous neurodegenerative disorders characterized by iron accumulation and the presence of axonal spheroids in the basal ganglia and other brain areas. Given the rarity of the disease, it can often go unrecognized or misdiagnosed. In PD, iron accumulation is a cardinal feature of degenerating regions in the brain and seems to be a key player in mechanisms that precipitate cell death. The aim of this study was to further explore the relationship between NBIA associated genes and PD etiology. We examined whether a genetic burden of variants in NBIA genes could contribute to the risk of developing PD by performing gene-based and single variant association analyses.

## Structure of ReadMe
### 1) Locating and extracting data for our genes of interest 
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

# Fisher association and logistic regression

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

    > A <- fread("AMP_FA2H.trimmed.annotation.forR.txt", header = T)
    > B <- fread("/Users/alvarezjerezp2/Desktop/HGMD_info/HGMD_FA2H.txt", header=T)
    > B$V12 <- B$V13 <- B$V14 <- B$V15 <- B$V16 <- B$V17 <- B$V18 <- B$V19<- B$V20 <- B$V21 <- B$V22 <- B$V23 <-B$V24 <-  NULL
    > colnames(B) <- c("CHR", "Start","ID","REF","ALT","HGMD Description","HGVS nucleotide", "HGVS protein","Variant class","Variant type", "Reported phenotype")
    > total <- merge(A,B, by="Start")
    > write.table(total, file="AMP_FA2H_HGMD.txt",  quote=F, row.names = F, sep = "\t") ...

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
