########## Script check, 15 DEC #######################################################################################################
########## Comment: P-value for MR in Cardiogram/FG looks wrong, reported in MS as "OR 1.36, 95% CI (1.07 - 1.72), p = 0.023"    ######
#######################################################################################################################################

  ## upload FG and CardiogramC4D table (shared rsids)
cardiogram_fg <- read.table("~/Desktop/UKBB_project/post_revision_1_scripts/GWAS_associations_t2d_fg_snps/cardiogram_fg.txt")
cad_fg <- read.table("~/Desktop/UKBB_project/post_revision_1_scripts/GWAS_associations_t2d_fg_snps/cad_fg.txt")
  
  ## merge by rs-id
merge_fg <- merge(cardiogram_fg, cad_fg, by.x = "markername", by.y = "snp")

  #### Check 1: maybe fundamental problem with extracting wrong data?
  #### check by hand: take .txt files of downloaded SNP associations from MAGIC and Cardiogram, look up random SNPs
  ## CardiogramC4D: 
  ## rs11039182 11	47346723	T	C	.755238	.99145	FIXED	-.01216	.011054	.271307 ...
  ## rs3829109  9	139256766	G	A	.741192	.8335	FIXED	.020106	.0119773	.0932154 ...
  ## MAGIC, Scott supplement
  ## rs11039182  11	47303299	MADD	MADD	T	C	0.729911	0.023313	0.002414	4.82E-22 ...
  ## rs3829109  9	138376587	DNLZ	.	G	A	0.706634	0.017249	0.002673	1.13E-10 ...
  #### all correct, excactly as in merged data file in this script ####

  ## align effect alleles, i.e. Cardiogram-logOR becomes negative if the FG_effect_allele is Cardiogram's non-effect-allele
merge_fg[which(as.character(merge_fg$effect_allele) != as.character(merge_fg$effall)),]$beta <- 
  - merge_fg[which(as.character(merge_fg$effect_allele) != as.character(merge_fg$effall)),]$beta

BetaX <- merge_fg$Effect    # FG effect in original scale, mmol/l
BetaX_se <- merge_fg$SE     # FG-se
BetaY <- merge_fg$beta      # CAD log-OR effect, aligned
BetaY_se <- merge_fg$se_dgc # CAD-se

library(MendelianRandomization) ## Burgess R package  

mr_results <- mr_allmethods(mr_input(bx = BetaX, bxse = BetaX_se, by = BetaY, byse = BetaY_se))$Values
  
  ## We use: the penalized robust IVW-method (downweighs outlier effects, estimates robust bootstap SEs)

ivw <- mr_results[mr_results$Method == "Penalized robust IVW",]
    
    ## display these results 
###### >  ivw
######    Method                Estimate  Std Error   95% CI               P-value      
###### 7  Penalized robust IVW  0.3047398 0.1221126   0.06540341 0.5440761 0.02315709   

  #### Check 2: Calculate OR and 95% CIs - should be as in MS "OR 1.36, 95% CI (1.07 - 1.72), p = 0.023"
  
exp(ivw$Estimate)                                               # correct
exp(ivw$Estimate - 1.96*ivw[,which(names(ivw) == "Std Error")]) # correct
exp(ivw$Estimate + 1.96*ivw[,which(names(ivw) == "Std Error")]) # correct
ivw[,which(names(ivw) == "P-value")]                            # correct

#### bug in R package? highly unlikely, but we have the Suppl. from a previous Burgess paper that reports the
#### method used for IVW analysis (CAVE: not penalized robust IVW, which we use from the package as more robust estimate)
#### at least: compared R package "IVW" with original script "IVW" in Suppl. App.
#### literally from the appendix:

#   Web Appendix  (updated 31/7/15)
#   code and results for the lung function example
BetaYG    = BetaY # vector of gene-outcome associations 
BetaXG    = BetaX # vector of gene-exposure associations 
seBetaYG  = BetaY_se # vector of gene-outcome association standard errors 
seBetaXG  = BetaX_se # vector of gene-exposure association standard errors 
#############################################
# IVW approach (with MAF corrected weights) #  
#############################################  
IVWfit      = summary(lm(BetaYG ~ -1+BetaXG,weights=1/seBetaYG^2)) 
# Inference with correct standard errors
DF      = length(BetaYG)-1
IVWBeta = IVWfit$coef[1,1]
SE      = IVWfit$coef[1,2]/min(1,IVWfit$sigma)
IVW_p   = 2*(1-pt(abs(IVWBeta/SE),DF))
IVW_CI  = IVWBeta + c(-1,1)*qt(df=DF, 0.975)*SE
# IVWResults = (point estimate, corrected standard error, 
# 95% Confidence interval, t-statistic, p-value) 
IVWResults = c(IVWBeta,SE,IVW_CI,IVWBeta/SE,IVW_p)
IVWResults
mr_results[which(mr_results$Method == "IVW"),]

###   excactly similar
###   > IVWResults
###   [1]  0.33070548  0.18385601 -0.05719679  0.71860776  1.79872000  0.08984174
###   > mr_results[which(mr_results$Method == "IVW"),]
###   Method  Estimate Std Error     95% CI               P-value
###   4    IVW 0.3307055  0.183856 -0.02964568 0.6910566 0.08984174

#############################################################################################################################################
####### Conclusion (Chris, 15 DEC 2016): scripts are correct, methods correctly applied, data correctly extracted       #####################
####### the p-value for "OR 1.36, 95% CI (1.07 - 1.72), p = 0.023" may look large, but                                  #####################
####### a)  the CI is large                                                                                             #####################
####### b)  these are the exponentials to (logOR), the p-value looks less unusual when looking at the original scale    #####################
#######     causal logOR = 3047398, se = 0.1221126, 95% CI (0.06540341; 0.5440761)                                      #####################
#############################################################################################################################################