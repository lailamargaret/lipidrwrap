library(lipidr)
library(ggpubr)
library(grid)
library(gridExtra)
library(gridGraphics)
library(dplyr)
library(plyr)
library(ggplotify)
library(cowplot)
library(car)
library(bestNormalize)
library(geoR)
library(forcats)
use_interactive_graphics(FALSE)

########## analyze_ functions: one liners to do ALL lipidomics analysis in one #########
#' 4-way analysis with gonSex and chrSex columns
analyze_fcg<-function(d, experiment_name){
  #PCAs - supervised and unsupervised
  mva_result<-get_mva(d, experiment_name)
  chr_sup_mva_result<-get_sup_mva(d, experiment_name, "chromSex")
  gon_sup_mva_result<-get_sup_mva(d, experiment_name, "gonSex")
  
  #do differential analysis by comparison
  chr_ov<-de_analysis(d, XX_ov - XY_ov)
  chr_test<-de_analysis(d, XX_test - XY_test)
  gon_XX<-de_analysis(d, XX_ov - XX_test)
  gon_XY<-de_analysis(d, XY_ov - XY_test)
  chr<-de_analysis(d, group_col="chromSex", XX - XY)
  gon<-de_analysis(d, group_col="gonSex", O - T)
  
  #lipidset enrichment for the above comparisons
  ls_chr_ov<-lipidset(chr_ov, experiment_name)
  ls_chr_test<-lipidset(chr_test, experiment_name)
  ls_gon_XX<-lipidset(gon_XX, experiment_name)
  ls_gon_XY<-lipidset(gon_XY, experiment_name)
  ls_chr<-lipidset(chr, experiment_name)
  ls_gon<-lipidset(gon, experiment_name)
  
  lsea_results<-rbind(ls_chr, ls_gon, ls_chr_ov, ls_chr_test, ls_gon_XX, ls_gon_XY)
  le<-lsea_results[["leadingEdge"]]
  lec<-c()
  for(i in 1:length(le)){
    lec<-append(lec,paste(le[[i]], collapse=" "))
  }
  lsea_results[["leadingEdge"]]<-lec
  
  filename=gsub(" ", "_", experiment_name)
  saveloc=paste(filename, "_lsea_results.csv", sep="")
  write.csv(lsea_results, saveloc, row.names=FALSE)
  
  #differential analysis plots
  fourway_list<-list(chr_ov, chr_test, gon_XY, gon_XX)
  composite_volcano(fourway_list, experiment_name, "4-way")
  twoway_list<-list(chr, gon)
  composite_volcano(twoway_list, experiment_name, "2-way")
  
  #anova
  chr_sex_anova<-de_design(d, ~chromSex)
  composite_volcano(chr_sex_anova, experiment_name, "chr_sex_ANOVA")
  ls_chr_sex_anova<-lipidset(chr_sex_anova, experiment_name)
  
  gon_sex_anova<-de_design(d, ~gonSex)
  composite_volcano(gon_sex_anova, experiment_name, "gon_sex_ANOVA")
  ls_gon_sex_anova<-lipidset(gon_sex_anova, experiment_name)
  

  de_results<-bind_rows(chr_ov, chr_test, gon_XX, gon_XY, chr, gon)
  de_save=paste(filename, "_de_results.csv", sep="")
  write.csv(de_results, de_save)
  
  
  lipidrresult<-list(mvaResult=mva_result,
                     chrSupMVAResult=chr_sup_mva_result,
                     gonSupMVAResult=gon_sup_mva_result,
                     deChrOv=chr_ov,
                     deChrTest=chr_test,
                     deGonXX=gon_XX,
                     deGonXY=gon_XY,
                     deChr=chr,
                     deGon=gon,
                     lsChrOv=ls_chr_ov,
                     lsChrTest=ls_chr_test,
                     lsGonXX=ls_gon_XX,
                     lsGonXY=ls_gon_XY,
                     lsChr=ls_chr,
                     lsGon=ls_gon,
                     deChrAnova=chr_sex_anova,
                     deGonAnova=gon_sex_anova,
                     lsChrAnova=ls_chr_sex_anova,
                     lsGonAnova=ls_gon_sex_anova)
  class(lipidrresult)<-"lipidrResults"
  return(lipidrresult)
  
}

#' 4-way analysis with kdm5c and tissue columns
analyze_kdm5c_tissue<-function(d, experiment_name){
  #PCAs - supervised and unsupervised
  mva_result<-get_mva(d, experiment_name)
  tissue_sup_mva_result<-get_sup_mva(d, experiment_name, "tissue")
  gt_sup_mva_result<-get_sup_mva(d, experiment_name, "kdm5c")
  
  #do differential analysis by comparison
  tissue_HET<-de_analysis(d, HET_liver - HET_gWAT)
  tissue_WT<-de_analysis(d, WT_liver - WT_gWAT)
  gt_liver<-de_analysis(d, WT_liver - HET_liver)
  gt_gwat<-de_analysis(d, WT_gWAT - HET_gWAT)
  tissue<-de_analysis(d, group_col="tissue", liver - gWAT)
  gt<-de_analysis(d, group_col="kdm5c", HET - WT)
  
  #lipidset enrichment for the above comparisons
  ls_tissue_HET<-lipidset(tissue_HET, experiment_name)
  ls_tissue_WT<-lipidset(tissue_WT, experiment_name)
  ls_gt_liver<-lipidset(gt_liver, experiment_name)
  ls_gt_gwat<-lipidset(gt_gwat, experiment_name)
  ls_tissue<-lipidset(tissue, experiment_name)
  ls_gt<-lipidset(gt, experiment_name)
  
  lsea_results<-rbind(ls_tissue_HET, ls_tissue_WT, ls_gt_liver, 
                      ls_gt_gwat, ls_tissue, ls_gt)
  le<-lsea_results[["leadingEdge"]]
  lec<-c()
  for(i in 1:length(le)){
    lec<-append(lec,paste(le[[i]], collapse=" "))
  }
  lsea_results[["leadingEdge"]]<-lec
  
  filename=gsub(" ", "_", experiment_name)
  saveloc=paste(filename, "_lsea_results.csv", sep="")
  write.csv(lsea_results, saveloc, row.names=FALSE)
  
  
  #differential analysis plots
  fourway_list<-list(tissue_HET, tissue_WT, gt_liver, gt_gwat)
  composite_volcano(fourway_list, experiment_name, "4-way")
  twoway_list<-list(tissue, gt)
  composite_volcano(twoway_list, experiment_name, "2-way")
  
  #anova
  tissue_anova<-de_design(d, ~tissue)
  composite_volcano(tissue_anova, experiment_name, "tissue_ANOVA")
  ls_tissue_anova<-lipidset(tissue_anova, experiment_name)
  
  gt_anova<-de_design(d, ~kdm5c)
  composite_volcano(gt_anova, experiment_name, "kdm5c_ANOVA")
  ls_gt_anova<-lipidset(gt_anova, experiment_name)
  
  lipidrresult<-list(mvaResult=mva_result,
                     tissueSupMVAResult=tissue_sup_mva_result,
                     gtSupMVAResult=gt_sup_mva_result,
                     deTissueHET=tissue_HET,
                     deTissueWT=tissue_WT,
                     deGtLiver=gt_liver,
                     deGtGwat=gt_gwat,
                     deTissue=tissue,
                     deGt=gt,
                     lsTissueHET=ls_tissue_HET,
                     lsTissueWT=ls_tissue_WT,
                     lsGtLiver=ls_gt_liver,
                     lsGtGwat=ls_gt_gwat,
                     lsTissue=ls_tissue,
                     lsGt=ls_gt,
                     deTissueAnova=tissue_anova,
                     deGtAnova=gt_anova,
                     lsTissueAnova=ls_tissue_anova,
                     lsGt=ls_gt_anova
  )
  class(lipidrresult)<-"lipidrResults"
  return(lipidrresult)
  
}

#' 4-way analysis with c19 and obesity columns 
analyze_c19<-function(d, experiment_name){
  #PCAs - supervised and unsupervised
  mva_result<-get_mva(d, experiment_name)
  c19_mva_result<-get_sup_mva(d, experiment_name, "c19")
  obesity_mva_result<-get_sup_mva(d, experiment_name, "obesity")
  
  #do differential analysis by comparison
  c19_obese<-de_analysis(d, obese_infected - obese_uninfected)
  c19_lean<-de_analysis(d, lean_infected - lean_uninfected)
  obesity_uninfected<-de_analysis(d, obese_uninfected - lean_uninfected)
  obesity_infected<-de_analysis(d, obese_infected - lean_infected)
  c19<-de_analysis(d, group_col="c19", infected - uninfected)
  obesity<-de_analysis(d, group_col="obesity", obese - lean)
  
  #lipidset enrichment for the above comparisons
  ls_c19_obese<-lipidset(c19_obese, experiment_name)
  ls_c19_lean<-lipidset(c19_lean, experiment_name)
  ls_obesity_uninfected<-lipidset(obesity_uninfected, experiment_name)
  ls_obesity_infected<-lipidset(obesity_infected, experiment_name)
  ls_c19<-lipidset(c19, experiment_name)
  ls_obesity<-lipidset(obesity, experiment_name)
  
  
  lsea_results<-rbind(ls_c19_obese, ls_c19_lean, ls_obesity_uninfected, 
                      ls_obesity_infected, ls_c19, ls_obesity)
  le<-lsea_results[["leadingEdge"]]
  lec<-c()
  for(i in 1:length(le)){
    lec<-append(lec,paste(le[[i]], collapse=" "))
  }
  lsea_results[["leadingEdge"]]<-lec
  
  filename=gsub(" ", "_", experiment_name)
  saveloc=paste(filename, "_lsea_results.csv", sep="")
  write.csv(lsea_results, saveloc, row.names=FALSE)
  
  
  #differential analysis plots
  fourway_list<-list(c19_obese, c19_lean, obesity_uninfected, obesity_infected)
  composite_volcano(fourway_list, experiment_name, "4-way")
  twoway_list<-list(c19, obesity)
  composite_volcano(twoway_list, experiment_name, "2-way")
  
  #anova
  c19_anova<-de_design(d, ~c19)
  composite_volcano(c19_anova, experiment_name, "c19_ANOVA")
  ls_c19_anova<-lipidset(c19_anova, experiment_name)
  
  obesity_anova<-de_design(d, ~obesity)
  composite_volcano(obesity_anova, experiment_name, "obesity_ANOVA")
  ls_obesity_anova<-lipidset(obesity_anova, experiment_name)
  
  filename=gsub(" ", "_", experiment_name)
  
  de_results<-bind_rows(c19_obese, c19_lean, obesity_uninfected, obesity_infected, c19, obesity)
  de_save=paste(filename, "_de_results.csv", sep="")
  write.csv(de_results, de_save)
  
  # lsea_results<-rbind(as.data.frame(ls_c19_obese), 
  #                     as.data.frame(ls_c19_lean), 
  #                     as.data.frame(ls_obesity_uninfected), 
  #                     as.data.frame(ls_obesity_infected), 
  #                     as.data.frame(ls_c19), 
  #                     as.data.frame(ls_obesity))
  # lsea_save=paste(filename, "_lsea_results.csv", sep="")
  # write.csv(lsea_results, lsea_save)
  
  lipidrresult<-list(mvaResult=mva_result,
                     c19SupMVAResult=c19_mva_result,
                     obesitySupMVAResult=obesity_mva_result,
                     dec19Obese=c19_obese,
                     dec19Lean=c19_lean,
                     deObesityUninfected=obesity_uninfected,
                     deObesityInfected=obesity_infected,
                     dec19=c19,
                     deObesity=obesity,
                     ls_c19Obese=ls_c19_obese,
                     ls_c19Lean=ls_c19_lean,
                     ls_ObesityUninfected=ls_obesity_uninfected,
                     ls_ObesityInfected=ls_obesity_infected,
                     ls_c19=ls_c19,
                     ls_Obesity=ls_obesity,
                     dec19Anova=c19_anova,
                     deObesityAnova=obesity_anova,
                     lsc19Anova=ls_c19_anova,
                     lsObesityAnova=ls_obesity_anova)
  class(lipidrresult)<-"lipidrResults"
  return(lipidrresult)
  
}

#' 2-way analysis with just kdm5c column
analyze_kdm5c<-function(d, experiment_name){
  #PCAs - supervised and unsupervised
  mva_result<-get_mva(d, experiment_name)
  gt_sup_mva_result<-get_sup_mva(d, experiment_name, "kdm5c")
  
  #do differential analysis by comparison
  gt<-de_analysis(d, group_col="kdm5c", WT-HET)
  
  composite_volcano(gt, experiment_name, "")
  
  #lipidset enrichment for the above comparisons
  ls_gt<-lipidset(gt, experiment_name)
  
  lipidrresult<-list(mvaResult=mva_result,
                     gtSupMVAResult=gt_sup_mva_result,
                     deGt=gt,
                     lsGt=ls_gt
  )
  class(lipidrresult)<-"lipidrResults"
  return(lipidrresult)
  
}



############ Helper functions used in analyze_ functions #######################

#' preprocess_lipidomics
#'
#' @param exp_name a string w experiment name, used for naming files
#' @param datafile csv file containing raw data, w lipid species as rows and samples as cols
#' @param samplesfile csv file containing sample data, including
#' @param knn how many neighbors to use for k-nearest neighbors imputation
#' @param percna of samples missing values for a lipid to be excluded
#'
#' @return lipidomicsExperiment class object with log transformed, normalized values
#' @export QCplot PNG with histograms of data dist and qq plot for LTF/norm data
#'
#' @examples
preprocess_lipidomics<-function(exp_name, datafile, samplesfile, knn=10, percna=0.33, method="log2"){
  exp_df<-read.csv(datafile, header = T) #read raw data csv
  removed<-exp_df[which(rowMeans(is.na(exp_df)) > (percna+0.01)),] #store removed species for later analysis
  cleaned_exp_df<-exp_df[which(rowMeans(is.na(exp_df)) < (percna+0.01)),] #remove all species w > 34% NAs
  old_mol_names=cleaned_exp_df[[1]] #save old molecule names for relabeling
  cleaned_exp_df[[1]] = sub("^(PC|PE) ([OP])-", "\\1\\2 ", cleaned_exp_df[[1]]) #regex to make lipid names into standard nomenclature
  cleaned_exp_df[[1]] = sub("^TG(.*)-FA.*", "TG\\1", cleaned_exp_df[[1]])#regex to make lipid names into standard nomenclature
  cleaned_exp_df[[1]] = sub("\\+H$", "", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("\\+NH4$", "", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("\\+Na$", "", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("\\+HCOO$", "", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("\\-H$", "", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("e)", ")", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("e_", "_", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("p_", "_", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("\\+O", "", cleaned_exp_df[[1]])
  cleaned_exp_df[[1]] = sub("\\+2O", "", cleaned_exp_df[[1]])
  
  
  
  d<-as_lipidomics_experiment(cleaned_exp_df) #import to lipidr
  rowData(d)$Molecule = old_mol_names #change rownames back
  d=add_sample_annotation(d, samplesfile) #add sample annotations to lipidr experiment
  d@elementMetadata@listData[["filename"]][d@elementMetadata@listData[["filename"]]=="dataframe"]<-exp_name
  
  
  par(cex=1, cex.axis=0.75, cex.main=1, cex.lab=1, mgp=c(1,0.4,0), mar=c(9,3,4,4))
  
  
  d_knn = impute_na(d, measure = "Area", method = "knn", knn)
  r<-ggplot()+
    (aes(d_knn@assays@data@listData[["Area"]]))+
    geom_histogram(bins=50) + 
    #scale_y_continuous(trans='log10') +
    labs(title="Raw intensity (imputed)", x= "Intensity")
  
  
  d_norm<-normalize_pqn(d_knn, measure="Area", log=FALSE)
  d_norm@elementMetadata@listData[["filename"]][d_norm@elementMetadata@listData[["filename"]]==exp_name]<-paste(exp_name, "PQN norm", sep=" ")
  n<-ggplot()+
    (aes(d_norm@assays@data@listData[["Area"]]))+
    geom_histogram(bins=50) + 
    #scale_y_continuous(trans='log10') +
    labs(title="Normalized intensity", x= "Intensity")
  
  if(method=="log2"){
    #d_norm@assays@data@listData[["Area"]]<-d_norm@assays@data@listData[["Area"]] + 0.1
    new_vals<-log2(d_norm@assays@data@listData[["Area"]])
    d_transform<-d_norm
    d_transform@assays@data@listData[["Area"]]<-new_vals
    d_transform<-set_logged(d_transform, measure="Area", T)
    assay(d_transform)[is.infinite(assay(d_transform))]<-NA
    
  }
  else if(method=="boxcox"){
    bc<-boxcoxfit(d_norm@assays@data@listData[["Area"]],lambda2=TRUE)
    l1<-rep(bc[["lambda"]][["lambda"]], 
            ncol(d_norm@assays@data@listData[["Area"]]))
    l2<-rep(bc[["lambda"]][["lambda2"]], 
            ncol(d_norm@assays@data@listData[["Area"]]))
    new_vals<-bcnPower(d_norm@assays@data@listData[["Area"]], l1, gamma=l2)
    d_transform<-d_norm
    d_transform@assays@data@listData[["Area"]]<-new_vals
    d_transform<-set_logged(d_transform, measure="Area", T)
    #assay(d_transform)[is.infinite(assay(d_transform))]<-NA
    
  }
  else{print("Error - method does not exist")
    return(NA)}
  
  l<-ggplot()+
    (aes(d_transform@assays@data@listData[["Area"]]))+
    geom_histogram(bins=50) +
    #scale_y_continuous(trans='log10') +
    labs(title="log2 Transformed Intensity", x= "Intensity")
  qqPlot(d_transform@assays@data@listData[["Area"]], 
         id=F,
         ylab="", 
         xlab="")
  title(xlab="norm quantiles", line=1.1)
  title(ylab="Intensity", line=1.1)
  title(main="Q-Q plot", line=0.5, font=1)
  grid.echo()
  lq<-grid.grab()
  
  title<-paste(exp_name, "QC Plots", sep=" ") 
  composite<-grid.arrange(r, n, l, lq, 
                          top = textGrob(title, gp=gpar(fontsize=14)),
                          nrow=4)
  filename=gsub(" ", "_", exp_name)
  saveloc=paste(filename, "_qc_plots.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 15, height = 10, dpi = 150, 
         units = "in", device='png')
  
  
  b<-plot_samples(d_knn, color="Group", type="tic", log=F)
  bn<-plot_samples(d_norm, color="Group", type="tic", log=F)
  s<-plot_samples(d_knn, color="Group", type = "boxplot", log=T)
  sn<-plot_samples(d_norm, color="Group", type="boxplot",log=T)
  title<-paste(exp_name, "total lipids (raw/pqn norm)")
  #if(norm){ title<-paste(exp_name, "total lipids (pqn norm)")}
  composite<-grid.arrange(b,s,bn,sn, top = textGrob(title, gp=gpar(fontsize=14)))
  filename=gsub(" ", "_", exp_name)
  saveloc=paste(filename, "_persample_total_lip.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
  
  
  return(d_transform)
}


#' samplePlot
#'
#' @param d - lipidomicsExperiment object
#' @param exp_name - str, experiment identifier
#' @param norm - whether or not to the d is pqn normalized, default F
#'
#' @return na
#' @export PNG with persample total lipids
#'
#' @examples
samplePlot<-function(d, exp_name, norm=F){
  file_name=gsub(" ", "_", exp_name)
  s<-plot_samples(d, color="Group", type = "tic", log=FALSE)
  #sn<-plot_samples(d_norm, color="Group", type="tic",log=FALSE)
  title<-paste(exp_name, "total lipids (raw)")
  if(norm){ title<-paste(exp_name, "total lipids (pqn norm)")}
  composite<-grid.arrange(s, top = textGrob(title, gp=gpar(fontsize=14)))
  saveloc=paste(file_name, "_persample_total_lip.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
}


#' get_mva
#'
#' @param d - lipidomicsExperiment object
#' @param exp_name - str, experiment identifier
#'
#' @return result of calling mva function in lipidr, DF results of PCA analysis
#' @export unsupervised PCA PNG
#'
#' @examples
get_mva<-function(d, exp_name){
  mva = mva(d, measure="Area", method="PCA")
  m<-plot_mva(mva, color_by = "Group")
  #l<-plot_mva_loadings(mva)
  title<-paste(exp_name, "PCA (unsupervised)") 
  composite<-grid.arrange(m,top = textGrob(title, gp=gpar(fontsize=14)))
  filename=gsub(" ", "_", exp_name)
  saveloc=paste(filename, "_unsup_pca.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
  return(mva)
}

#' get_sup_mva
#'
#' @param d - lipidomicsExperiment object
#' @param exp_name - str, experiment identifier
#' @param grouping - which group column on the samples sheet to use to group based on 
#'
#' @return result of mva function as in get_mva, but for OPLS-DA
#' @export PNG of OPLS-DA plot, plot with each lipid species' loadings, csv with tabular version of loadings
#'
#' @examples
get_sup_mva<-function(d, exp_name, grouping){
  gps<-unique(d[[grouping]])
  supmva = mva(d, measure="Area", method="OPLS-DA", group_col = grouping, groups=gps)
  
  m<-plot_mva(supmva, color_by = grouping)
  title<-paste(exp_name, " PCA (by ", grouping, ")", sep="") 
  composite<-grid.arrange(m, top = textGrob(title, gp=gpar(fontsize=14)))
  filename=gsub(" ", "_", exp_name)
  filename=paste(filename, grouping, sep="_")
  saveloc=paste(filename, "_sup_pca.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
  
  
  title2<-paste(exp_name, " OPLS-DA loadings (", grouping, ")", sep="")
  saveloc2=paste(filename, "_sup_pca_loadings.png", sep = "")
  l<-plot_mva_loadings(supmva, color_by="Class", top.n=20)
  composite2<-grid.arrange(l, top = textGrob(title2, gp=gpar(fontsize=14)))
  ggsave(saveloc2, plot=composite2, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
  
  csv_sl<-paste(filename, "_oplsda_top_lipids.csv", sep="")
  write.csv(top_lipids(supmva, top.n=25), csv_sl)
  
  
  return(supmva)
}


#' composite_volcano
#'
#' @param res_list - a list of de_analysis result data.tables to be used in one combined volcano plot
#' @param exp_name - str, experiment identifier
#' @param modifier - a str modifier to use as detail for filenames, defaults to empty string 
#'
#' @return na
#' @export png with multiple volcano plots saved as one image
#'
#' @examples
composite_volcano<-function(res_list, exp_name, modifier=""){
  composite_res<-bind_rows(res_list)
  p<-plot_results_volcano(composite_res)
  
  title<-paste(exp_name, modifier, "differential analysis", sep=" ") 
  composite<-grid.arrange(p, top = textGrob(title, gp=gpar(fontsize=14)))
  
  filename=gsub(" ", "_", exp_name)
  saveloc=paste(filename, "_", modifier, "_de_volcano.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
}


#' lipidset
#'
#' @param result a deanalysis result object
#' @param exp_name a string w experiment name, used for naming files 
#'
#' @return lsearesult 
#' @export lipidset_enrichment.png with plot of chain len, class, and saturation significance 
#'
#' @examples
lipidset<-function(result, exp_name, fold_change=NA){
  modifier<-deparse(substitute(result))
  lipset_res<-lsea(result)
  ss<-significant_lipidsets(lipset_res)
  
  csv_fn=gsub(" ", "_", exp_name)
  csv_fn=paste(csv_fn, modifier, sep="_")
  csv_sl=paste(csv_fn, "lipidset_enrichment.csv", sep="_")
  
  cl<-plot_enrichment_edit(result, ss, annotation="length")
  c<-plot_enrichment_edit(result, ss, annotation="class")
  s<-plot_enrichment_edit(result, ss, annotation="unsat")
  #clp<-plot_enrichment(result, ss, annotation="length", measure="adj.P.Val")
  #cp<-plot_enrichment(result, ss, annotation="class", measure="adj.P.Val")
  #sp<-plot_enrichment(result, ss, annotation="unsat", measure="adj.P.Val")
  
  title<-paste(exp_name, modifier, "set enrichment", sep = " ") 
  composite<-grid.arrange(c, cl, s,# clp, cp, sp, nrow=3, ncol=2,
                          #layout_matrix=rbind(c(1,4), c(2,5), c(3,6)),
                          #widths = c(1, 1, 1),
                          #layout_matrix = rbind(c(1, 2, 3),
                          #                     c(3, 3, 4)))
                          
                          top = textGrob(title, gp=gpar(fontsize=14)))
  csv_fn=gsub(" ", "_", exp_name)
  csv_fn=paste(csv_fn, modifier, sep="_")
  plot_sl=paste(csv_fn, "lipidset_enrichment.png", sep="_")
  ggsave(plot_sl, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
  
  return(lipset_res)
}



### A fix for lipidr's plot_enrichment function. Used in place of plot_enrichment in analyze functions
### Degree of unsat results are now in numerical order.
### can specify y-axis min/max (always centered on zero) if want to have same y axis scale.
### called only by lipidset function.
plot_enrichment_edit<-function(de.results, 
                               significant.sets, 
                               annotation=c("class", "length", "unsat"), 
                               measure = "logFC", ylim=NA) {
  annotation <- match.arg(annotation)
  collection <- c(class="Class", length="total_cl", unsat="total_cs")[[annotation]]
  x_label <- c(class="Lipid Class", length="Total Chain Length", unsat="Total Chain Unsaturation")[[annotation]]
  prefix = paste0("^", collection, "_")
  
  
  significant.sets <- lapply(
    significant.sets,
    function(c) sub(prefix, "", c[grep(prefix, c)])
  )
  #return(group_by(de.results, contrast))
  
  de_results <- de.results %>% dplyr::group_by(contrast) %>%
    dplyr::mutate(Significant = factor(
      de.results[[as.character(sym(collection))]] %in% significant.sets[[de.results[["contrast"]][[1]]]],
      levels = c(TRUE, FALSE)
    )) %>%
    dplyr::mutate(Enrichment = forcats::fct_recode(Significant, "Significant"="TRUE", "Not significant"="FALSE")) %>%
    ungroup() %>%
    dplyr::mutate(collection = as.character((!!sym(collection))), measure = (!!sym(measure)),
    )
  
  if(annotation=="class"){
    p <- ggplot(de_results, aes(x=(collection), y=measure, color = Enrichment)) +
      geom_boxplot() + geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(~contrast, scales = "free_x") +
      scale_color_manual(values = c(`Not significant`="black", `Significant`="red"), drop=FALSE) +
      labs(x=x_label, y=measure) +
      ylim(-ylim, ylim) + 
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
  } else {
    p <- ggplot(de_results, aes(x=fct_inseq(collection), y=measure, color = Enrichment)) +
      geom_boxplot() + geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(~contrast, scales = "free_x") +
      scale_color_manual(values = c(`Not significant`="black", `Significant`="red"), drop=FALSE) +
      labs(x=x_label, y=measure) +
      ylim(-ylim, ylim) +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
  }
  
  
  p
}
#' totallipPlots
#'
#' @param d - a lipidomicsExperiment object
#' @param exp_name - str, experiment identifier
#'
#' @return none
#' @export PNG with total lipids plots
#'
#' @examples
totallipPlots <- function(d, exp_name) {
  m<-plot_molecules(d, "boxplot", measure = "Area", log=TRUE)
  lc<-plot_lipidclass(d, "boxplot", log=TRUE)
  title<-paste(exp_name, "all lipids") 
  composite<-grid.arrange(lc, m,ncol=2,
                          top = textGrob(title, gp=gpar(fontsize=14)))
  
  filename=gsub(" ", "_", exp_name)
  saveloc=paste(filename, "_total_lipids.png", sep="")
  ggsave(saveloc, plot=composite, 
         width = 10, height = 8, dpi = 150, 
         units = "in", device='png')
}