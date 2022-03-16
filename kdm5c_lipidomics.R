##kdm5c lipidomics
setwd("C:/Users/laila/OneDrive/Documents/rotations/reue/kdm5c_HFD")
source("C:/Users/laila/OneDrive/Documents/rotations/reue/lipidr_wrap.R")

gwat_exp_name<-"kdm5c gWAT HFD"
d_gwat<-preprocess_lipidomics(gwat_exp_name,
                              "kdm5c_HFD_gwat_raw_data.csv",
                              "kdm5c_HFD_gwat_samples.csv")
d_gwat_results<-analyze_kdm5c(d_gwat, gwat_exp_name)

liver_exp_name<-"kdm5c liver HFD"
d_liver<-preprocess_lipidomics(liver_exp_name, 
                               "kdm5c_HFD_liver_raw_data.csv", 
                               "kdm5c_HFD_liver_samples.csv" )
d_liver_results<-analyze_kdm5c(d_liver, liver_exp_name)


gwat_and_liver_exp_name<-"kdm5c gWAT and liver HFD"
d_gwat_liver<-preprocess_lipidomics(gwat_and_liver_exp_name,
                                    "kdm5c_HFD_gwat_liver_raw_data_combined.csv", 
                                    "kdm5c_HFD_liver_gwat_samples_combined.csv")
d_gwat_liver_results<-analyze_kdm5c_tissue(d_gwat_liver, gwat_and_liver_exp_name)





#####archive


########################## GWAT ONLY ###########################################










#read in data for just gwat
g_exp_df<-read.csv("kdm5c_HFD/kdm5c_HFD_gwat_raw_data.csv", header = T) #read raw data csv
g_removed<-g_exp_df[which(rowMeans(is.na(g_exp_df)) > 0.34),] #store removed species for later analysis
g_cleaned_exp_df<-g_exp_df[which(rowMeans(is.na(g_exp_df)) < 0.34),] #remove all species w > 34% NAs
write.csv(g_exp_df, "kdm5c_HFD/kdm5c_HFD_gwat_raw_data_33percna_removed.csv", row.names=FALSE)
g_old_mol_names=g_cleaned_exp_df[[1]] #save old molecule names for relabeling
g_cleaned_exp_df[[1]] = sub("^(PC|PE) ([OP])-", "\\1\\2 ", g_cleaned_exp_df[[1]]) #regex to make lipid names into standard nomenclature
g_cleaned_exp_df[[1]] = sub("^TG(.*)-FA.*", "TG\\1", g_cleaned_exp_df[[1]])#regex to make lipid names into standard nomenclature
write.csv(g_cleaned_exp_df, "kdm5c_HFD/kdm5c_HFD_gwat_raw_data_33percna_removed_cleaned_names.csv", row.names=FALSE)
g_sample_file="kdm5c_HFD/kdm5c_HFD_gwat_samples.csv" #which file is sample sheet
g_d<-as_lipidomics_experiment(g_cleaned_exp_df) #import to lipidr
rowData(g_d)$Molecule = g_old_mol_names #change rownames back
g_d=add_sample_annotation(g_d, g_sample_file)

g_d_knn = impute_na(g_d, measure = "Area", method = "knn", 10)
g_d_norm=normalize_pqn(g_d_knn, measure="Area", log=FALSE)
new_vals<-log2(g_d_norm@assays@data@listData[["Area"]])
g_d_norm@assays@data@listData[["Area"]]<-new_vals
g_d_norm<-set_logged(g_d_norm, measure="Area", T)
hist(g_d_norm@assays@data@listData[["Area"]])
drop<-which(apply(assay(g_d_norm), 1, function(r) any(r %in% c("Inf", "-Inf"))))
assay(g_d_norm)[is.infinite(assay(g_d_norm))]<-NA
qqnorm(g_d_norm@assays@data@listData[["Area"]], pch = 1, frame = FALSE)
qqline(g_d_norm@assays@data@listData[["Area"]], col = "steelblue", lwd = 2)
g_exp_name="kdm5c HFD gWAT"
g_d_norm_results<-analyze_kdm5c(g_d_norm, g_exp_name)

g_d<-preprocess_lipidomics(g_exp_name,"kdm5c_HFD/kdm5c_HFD_gwat_raw_data.csv","kdm5c_HFD/kdm5c_HFD_gwat_samples.csv")
qqnorm(g_d@assays@data@listData[["Area"]], pch = 1, frame = FALSE)
qqline(g_d@assays@data@listData[["Area"]], col = "steelblue", lwd = 2)

############################### LIVER ONLY #####################################

d_kdm5c_hfd_liver<-preprocess_lipidomics("kdm5c HFD liver", 
                                         "kdm5c_HFD/kdm5c_HFD_liver_raw_data.csv", 
                                         "kdm5c_HFD/kdm5c_HFD_liver_samples.csv" )


#read in data for jsut liver
l_exp_df<-read.csv("kdm5c_HFD/kdm5c_HFD_liver_raw_data.csv", header = T) #read raw data csv
l_removed<-l_exp_df[which(rowMeans(is.na(l_exp_df)) > 0.34),] #store removed species for later analysis
l_cleaned_exp_df<-l_exp_df[which(rowMeans(is.na(l_exp_df)) < 0.34),] #remove all species w > 34% NAs
write.csv(l_exp_df, "kdm5c_HFD/kdm5c_HFD_liver_raw_data_33percna_removed.csv", row.names=FALSE)
l_old_mol_names=l_cleaned_exp_df[[1]] #save old molecule names for relabeling
l_cleaned_exp_df[[1]] = sub("^(PC|PE) ([OP])-", "\\1\\2 ", l_cleaned_exp_df[[1]]) #regex to make lipid names into standard nomenclature
l_cleaned_exp_df[[1]] = sub("^TG(.*)-FA.*", "TG\\1", l_cleaned_exp_df[[1]])#regex to make lipid names into standard nomenclature
write.csv(l_cleaned_exp_df, "kdm5c_HFD/kdm5c_HFD_liver_raw_data_33percna_removed_cleaned_names.csv", row.names=FALSE)
l_sample_file="kdm5c_HFD/kdm5c_HFD_liver_samples.csv" #which file is sample sheet
l_d<-as_lipidomics_experiment(l_cleaned_exp_df) #import to lipidr
rowData(l_d)$Molecule = l_old_mol_names #change rownames back
l_d=add_sample_annotation(l_d, l_sample_file)

l_d_knn = impute_na(l_d, measure = "Area", method = "knn", 10)
l_d_norm=normalize_pqn(l_d_knn, measure="Area", log=FALSE)
new_vals<-log2(l_d_norm@assays@data@listData[["Area"]])
l_d_norm@assays@data@listData[["Area"]]<-new_vals
l_d_norm<-set_logged(l_d_norm, measure="Area", T)
hist(l_d_norm@assays@data@listData[["Area"]])
drop<-which(apply(assay(l_d_norm), 1, function(r) any(r %in% c("Inf", "-Inf"))))
assay(l_d_norm)[is.infinite(assay(l_d_norm))]<-NA
qqnorm(l_d_norm@assays@data@listData[["Area"]], pch = 1, frame = FALSE)
qqline(l_d_norm@assays@data@listData[["Area"]], col = "steelblue", lwd = 2)
l_exp_name="kdm5c HFD liver"
l_d_norm_results<-analyze_kdm5c(l_d_norm, l_exp_name)


############################# BOTH #############################################
#read in data for both 
gl_exp_df <- g_exp_df %>% full_join(l_exp_df)
write.csv(gl_exp_df, "kdm5c_HFD/kdm5c_HFD_gwat_liver_raw_data_combined.csv")

gl_removed<-gl_exp_df[which(rowMeans(is.na(gl_exp_df)) > 0.34),] #store removed species for later analysis
gl_cleaned_exp_df<-gl_exp_df[which(rowMeans(is.na(gl_exp_df)) < 0.34),] #remove all species w > 34% NAs
write.csv(gl_exp_df, "kdm5c_HFD/kdm5c_HFD_gwat_liver_raw_data_33percna_removed.csv", row.names=FALSE)
gl_old_mol_names=gl_cleaned_exp_df[[1]] #save old molecule names for relabeling
gl_cleaned_exp_df[[1]] = sub("^(PC|PE) ([OP])-", "\\1\\2 ", gl_cleaned_exp_df[[1]]) #regex to make lipid names into standard nomenclature
gl_cleaned_exp_df[[1]] = sub("^TG(.*)-FA.*", "TG\\1", gl_cleaned_exp_df[[1]])#regex to make lipid names into standard nomenclature
write.csv(gl_cleaned_exp_df, "kdm5c_HFD/kdm5c_HFD_gwat_liver_raw_data_33percna_removed_cleaned_names.csv", row.names=FALSE)

#save both as lipidr
gl_d<-as_lipidomics_experiment(gl_cleaned_exp_df) #import to lipidr
rowData(gl_d)$Molecule = gl_old_mol_names #change rownames back

#add annotations for both
liver_sample_sheet<-read.csv(l_sample_file)
gwat_sample_sheet<-read.csv(g_sample_file)
gl_sample_sheet<-rbind(gwat_sample_sheet, liver_sample_sheet)
write.csv(gl_sample_sheet, "kdm5c_HFD/kdm5c_HFD_liver_gwat_samples_combined.csv")
gl_d=add_sample_annotation(gl_d, gl_sample_sheet) #add sample annotations to lipidr experiment


gl_d_knn = impute_na(gl_d, measure = "Area", method = "knn", 10)
gl_d_norm=normalize_pqn(gl_d_knn, measure="Area", log=FALSE)
new_vals<-log2(gl_d_norm@assays@data@listData[["Area"]])
gl_d_norm@assays@data@listData[["Area"]]<-new_vals
gl_d_norm<-set_logged(gl_d_norm, measure="Area", T)
hist(gl_d_norm@assays@data@listData[["Area"]])
drop<-which(apply(assay(gl_d_norm), 1, function(r) any(r %in% c("Inf", "-Inf"))))
assay(gl_d_norm)[is.infinite(assay(gl_d_norm))]<-NA
qqnorm(gl_d_norm@assays@data@listData[["Area"]], pch = 1, frame = FALSE)
qqline(gl_d_norm@assays@data@listData[["Area"]], col = "steelblue", lwd = 2)
gl_exp_name="kdm5c HFD"
gl_d_norm_results<-analyze_kdm5c_tissue(gl_d_norm, gl_exp_name)


plot_samples(gl_d_knn, color="Group", type = "boxplot", log=TRUE)
plot_samples(gl_d_norm, color="Group", type = "boxplot", log=FALSE)


