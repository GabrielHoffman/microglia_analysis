

library(data.table)
library(variancePartition)

# MAF
# AD_donors_eQTL (42 donors): 
df_maf_ad = fread("/sc/arion/projects/Microglia/genotyping/Combined_Phase1and2/preparing_eQTL_files/ADvsControl_eQTL/AD_donors_eQTL_MAF.frq")
colnames(df_maf_ad)[colnames(df_maf_ad) == "MAF"] = "MAF_AD"
colnames(df_maf_ad)[colnames(df_maf_ad) == "SNP"] = "rsid"

# Control_donors_eQTL (34 donors):
df_maf_ctrl = fread("/sc/arion/projects/Microglia/genotyping/Combined_Phase1and2/preparing_eQTL_files/ADvsControl_eQTL/Control_donors_eQTL_MAF.frq")
colnames(df_maf_ctrl)[colnames(df_maf_ctrl) == "MAF"] = "MAF_CT"
colnames(df_maf_ctrl)[colnames(df_maf_ctrl) == "SNP"] = "rsid"

WorkingDir = "/sc/arion/projects/Microglia/Combined_RNAseq_ATACseq/comparing_AD_control_QTL_results"

## GET THE WORKSPACE WITH ALL GENERATED OBJECTS IN THIS SCRIPT: The file is 2.5G:
load(file = paste0(WorkingDir, "/Workspace_preparingFilesForADvsCont_AND_metaQTL_Overlap_RK_1_12_22.RDa"))


# aggregate AD vs control eQTL results across chromosomes
df_ad_eqtl = do.call(rbind, lapply(ALL_AD_Cont_Combined_eQTL_list, function(x){
	x$Pair = rownames(x)
	x[,c("Pair", "pVal", "FDR", "beta_Cont", "beta_AD", 'se_AD', 'se_Cont')]
	}))
dim(df_ad_eqtl)

# merge with eQTL data
df_merge = merge(finemapped_eQTL[,c('PostProb', 'Pair')], df_ad_eqtl, by="Pair")

df_merge = data.table(df_merge)

# merge by rsid to add MAF infl
df_merge$rsid = gsub(".*__", "", df_merge$Pair)
df_merge = merge(df_merge, df_maf_ad, by="rsid")
df_merge = merge(df_merge, df_maf_ctrl, by="rsid")

# check that alleles are the same
isSame = with(df_merge, A1.x == A1.y)

head(df_merge[!isSame,])

par(mfrow=c(1,2))
with(df_merge[isSame,], plot(beta_Cont, beta_AD, main="Same alleles"))
fit = lm(beta_AD ~ beta_Cont, df_merge[isSame,])
abline(fit, col="red")

with(df_merge[!isSame,], plot(beta_Cont, beta_AD, main="Different alleles"))
fit = lm(beta_AD ~ beta_Cont, df_merge[!isSame,])
abline(fit, col="red")

# compute z-statistics for all SNPs
z = with(df_merge, abs(beta_AD - beta_Cont) / sqrt(se_AD^2 + se_Cont^2))

# for SNPs with flipped alleles
z[!isSame] = with(df_merge[!isSame,], abs(beta_AD + beta_Cont) / sqrt(se_AD^2 + se_Cont^2))

# compute p-values and FDR
p = 2*pnorm(z, lower.tail=FALSE)
f = p.adjust(p, "BH")
table(f < 0.05)






https://github.com/GabrielHoffman/microglia_analysis/blob/7ffe33c088d970828c8f889163c3bf6e11bf4a5e/AD_QTL_analysis.R#L45




table(finemapped_eQTL$PostProb > 0.01)


tab = with(df_merge, table(FDR < 0.05, PostProb > 0.01))
tab
fisher.test(tab)



head(df[df$FDR < 0.0005,])


dim(df_merge)


table(df_ad_eqtl$FDR < 0.05)


tab = with(df_merge, table(FDR < 0.05, PostProb > 0.9))
fisher.test(tab)


df_merge$FDR.perm = sample(df_merge$FDR, nrow(df_merge))

tab = with(df_merge, table(FDR.perm < 0.05, PostProb > 0.9))
fisher.test(tab)

df_merge$delta = with(df_merge, abs(beta_Cont - beta_AD))


fit = glm( I(PostProb > .9) ~ I(FDR < 0.05), data = df_merge, family="binomial")

coef(summary(fit))

calcVarPart(fit)



z = with(df, abs(beta_AD - beta_Cont) / sqrt(se_AD^2 + se_Cont^2))
p = 2*pnorm(z, lower.tail=FALSE)
f = p.adjust(p, "BH")


f = p.adjust(df_ad_eqtl$p)

# AD_donors_caQTL (42 donors): 
/sc/arion/projects/Microglia/genotyping/Combined_Phase1and2/preparing_caQTL_files/ADvsControl_caQTL/AD_donors_caQTL_MAF.frq
# Control_donors_caQTL (32 donors):
 /sc/arion/projects/Microglia/genotyping/Combined_Phase1and2/preparing_caQTL_files/ADvsControl_caQTL/Control_donors_eQTL_MAF.frq













