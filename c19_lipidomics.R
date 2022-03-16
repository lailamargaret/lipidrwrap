setwd("C:/Users/laila/OneDrive/Documents/rotations/reue/c19_obesity")
source("C:/Users/laila/OneDrive/Documents/rotations/reue/lipidr_wrap.R")

samples<-"c19_samples.csv"

## Note: using boxcox transformation for these since there are many zero values and cannot ltf.

liver_rd<-"c19_liver_raw_data.csv"
liver_exp_name<-"c19 liver"
d_c19_liver<-preprocess_lipidomics(liver_exp_name, liver_rd, samples, method = "boxcox")
d_c19_liver_results<-analyze_c19(d_c19_liver, liver_exp_name)



fat_rd<-"c19_fat_raw_data.csv"
fat_exp_name<-"c19 fat"
d_c19_fat<-preprocess_lipidomics(fat_exp_name, fat_rd, samples, method = "boxcox")
d_c19_fat_results<-analyze_c19(d_c19_fat, fat_exp_name)


serum_rd<-"c19_serum_raw_data.csv"
serum_exp_name<-"c19 serum"
d_c19_serum<-preprocess_lipidomics(serum_exp_name, serum_rd, samples, method = "boxcox")
d_c19_serum_results<-analyze_c19(d_c19_serum, serum_exp_name)



cliv<-significant_molecules(d_c19_liver_results[["deObesityInfected"]])
cser<-significant_molecules(d_c19_serum_results[["deObesityInfected"]])
cfat<-significant_molecules(d_c19_fat_results[["deObesityInfected"]])

csig<-unique(c(cliv[[1]], cser[[1]], cfat[[1]]))


d_c19_liver_sig<-d_c19_liver[rowData(d_c19_liver)$Molecule %in% csig,]
d_c19_liver_results<-analyze_c19(d_c19_liver_sig, "c19 liver sig")

d_c19_fat_sig<-d_c19_fat[rowData(d_c19_fat)$Molecule %in% csig,]
d_c19_fat_results<-analyze_c19(d_c19_fat_sig, "c19 fat sig")

d_c19_serum_sig<-d_c19_serum[rowData(d_c19_serum)$Molecule %in% csig,]
d_c19_serum_results<-analyze_c19(d_c19_serum_sig, "c19 serum sig")

hist(d_c19_liver_results[["deObesityInfected"]]$adj.P.Val, breaks=300)
hist(d_c19_fat_results[["deObesityInfected"]]$adj.P.Val, breaks=300)

