##script set up##
setwd("C:/Users/laila/OneDrive/Documents/rotations/reue/FCG_HFD")
source("C:/Users/laila/OneDrive/Documents/rotations/reue/lipidr_wrap.R")



d_fcg_hfd<-preprocess_lipidomics("FCG HFD", "FCG_lipidomics_raw_data.csv", "FCG_samples.csv")
d_fcg_hfd_results<-analyze_fcg(d_fcg_hfd, "FCG HFD")

d_fcg_hfd_noTG<-d_fcg_hfd[!rowData(d_fcg_hfd)$Class=="TG",]
d_fcg_hfd_noTG_results<-analyze_fcg(d_fcg_hfd_noTG, "FCG HFD no TG")


cs_lips<-rowData(d_fcg_hfd)$Class %in% c("CE", "FFA", "PC", "PE", "PI", "PS", "SM", "TG")
d_fcg_top_classes<-d_fcg_hfd[cs_lips,]
d_fcg_top_classes_results<-analyze_fcg(d_fcg_top_classes, "FCG HFD top classes")




#### archive #####


#fourway<-function(d, experiment_name, var1, var2){
d<-d_fcg_hfd
experiment_name<-"test"
var1<-"chromSex"
var2<-"gonSex"

sv1<-sym(var1)
sv2<-sym(var2)
groups<-sort(unique(colData(d)$Group))
v1s<-sort(unique(colData(d)[[var1]]))
v2s<-sort(unique(colData(d)[[var2]]))


#mva_result<-get_mva(d, experiment_name)
#var1_sup_mva_result<-get_sup_mva(d, experiment_name, (sym(var1)))
#var2_sup_mva_result<-get_sup_mva(d, experiment_name, (sym(var2)))

#differential analysis by comparison
t<-(sym(groups[[1]]))
p<-(sym(groups[[2]]))

tg<-sym("XX_ov - XX_test")
tgq<-enquo(tg)

t2<-enquo(t)
p2<-enquo(p)


r<-sym("XX_ov - XX_test")
rq<-enquo(r)


#works
test<-de_analysis(d, !!rq)


#doesn't work
c12<-de_analysis(d, XX_ov - XX_test)

v1<-"XX_ov"
v2<-"XX_test"
testfun<-function(d, v1, v2){
  #v1<-enquo(v1)
  #v2<-enquo(v2)
  return(de_analysis(d, v1 - v2))
}

c12<-testfun(d, "XX_ov", "XX_test")


.quos_syms <- function(x) {
  print(x)
  if (rlang::is_syntactic_literal(x)) {
    print("test")
    return(lapply(sym(x), .quos_syms))
    
  } else if (is.symbol(x)) {
    print("test2")
    return(x)
    
  } else if (is.call(x)) {
    print("test3")
    print(x[TRUE][-1])
    return(unlist(lapply(x[TRUE][-1], .quos_syms)))
    
  } else if (is.pairlist(x)) {
    print("test4")
    return(unlist(lapply(x[TRUE], .quos_syms)))
    
  } else if (is.list(x)) {
    print("test5")
    return(unlist(lapply(x, .quos_syms)))
  }
}

de_analysis <- function(data, ..., measure = "Area", group_col = NULL) {
  if (is.null(group_col)) {
    if (ncol(colData(data)) > 0) {
      group_col <- names(colData(data))[[1]]
    } else {
      stop("Please add clinical data or specify a group column")
    }
  }
  
  #quos(...)
  quos(!!(...))
  
  symbols <- as.character(.quos_syms(quos(...)))
  if (length(symbols) == 0) {
    stop("No contrasts provided")
  }
  return(symbols)
  
  group <- colData(data)[[group_col]]
  if (!all(symbols %in% as.character(group))) {
    stop(
      "These constrast variables are not present in ", group_col, " column: ",
      paste(symbols[!symbols %in% as.character(group)], collapse=", ")
    )
  }
  data <- data[, group %in% symbols]
  
  group <- fct_drop(colData(data)[[group_col]])
  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  return(de_design(data = data, design = design, ..., measure = measure))
}


'''
analyze_fcg<-function(d, experiment_name){
  
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
'''


########################### Archive ############################################

ggplot(data=d_log_2_results[["deChr"]], aes(x=logFC, y=-log10(adj.P.Val), label=NA)) + 
  geom_point() + 
  theme_minimal()
ggplot(data=d_log_2_results[["deChr"]], aes(x=logFC, y=-log10(P.Value), label=NA)) + 
  geom_point() + 
  theme_minimal()
ggplot(data=d_log_2_results[["deChrTest"]], aes(x=logFC, y=-log10(P.Value), label=NA)) + 
  geom_point() + 
  theme_minimal()
ggplot(data=d_log_2_results[["deChrTest"]], aes(x=logFC, y=-log10(adj.P.Val), label=NA)) + 
  geom_point() + 
  theme_minimal()
ggplot(data=d_log_2_results[["deChrOv"]], aes(x=logFC, y=-log10(P.Value), label=NA)) + 
  geom_point() + 
  theme_minimal()
ggplot(data=d_log_2_results[["deChrOv"]], aes(x=logFC, y=-log10(adj.P.Val), label=NA)) + 
  geom_point() + 
  theme_minimal()

d_xxt = d_log_2[, d_log_2$Group == "XX_test"]
hist(d_xxt@assays@data@listData[["Area"]])
d_xxo = d_log_2[, d_log_2$Group == "XX_ov"]
hist(d_xxo@assays@data@listData[["Area"]])
d_xyt = d_log_2[, d_log_2$Group == "XY_test"]
hist(d_xyt@assays@data@listData[["Area"]])
d_xyo = d_log_2[, d_log_2$Group == "XY_ov"]
hist(d_xyo@assays@data@listData[["Area"]])

des<-c("deChrOv", "deChrTest")
for(i in 4:9){
  p<-hist(d_log_2_results[[i]]$P.Value)
  q<-hist(d_log_2_results[[i]]$adj.P.Val)
  title
  par(mfrow = c(1, 2))
  p
  q
  }

composite_volcano(list(d_log_2_results[["deChr"]], 
                       (d_log_2_results[["deGon"]])), "test", "4-way")
composite_volcano(list(d_log_2_results[["deChrOv"]], 
                       (d_log_2_results[["deChrTest"]]),
                       (d_log_2_results[["deGonXX"]]),
                       (d_log_2_results[["deGonXY"]])), "test", "4-way")

plot_samples(d_log_2, color="Group", type = "tic", log=FALSE)
plot_samples(d_log_2, color="Group", type = "boxplot", log=FALSE)
plot_samples(dl10_no_norm, color="Group", type = "boxplot", log=FALSE)

plot_mva_loadings(d_log_2_results[["chrSupMVAResult"]], color_by="Class", top.n=25)
plot_mva_loadings(d_log_2_results[["gonSupMVAResult"]], color_by="Class", top.n=25)


lsea_results<-rbind(as.data.frame(d_fcg_hfd_results[["lsChrOv"]]), 
                    as.data.frame(d_fcg_hfd_results[["lsChrTest"]]), 
                    as.data.frame(d_fcg_hfd_results[["lsGonXX"]]), 
                    as.data.frame(d_fcg_hfd_results[["lsGonXY"]]), 
                    as.data.frame(d_fcg_hfd_results[["lsChr"]]), 
                    as.data.frame(d_fcg_hfd_results[["lsGon"]]))


filename=gsub(" ", "_", "FCG HFD")
lsea_save=paste(filename, "_lsea_results.csv", sep="")

df<-ldply(lsea_results, data.frame)

le<-lsea_results[["leadingEdge"]]
lec<-c()
for(i in 1:length(le)){
  lec<-append(lec,paste(le[[i]], collapse=" "))
}
lsea_results[["leadingEdge"]]<-lec
write.csv(lsea_results, "FCG_HFD_lsea_results.csv")


test<-rbind(d_fcg_hfd_results[["lsChr"]], d_fcg_hfd_results[["lsGon"]], d_fcg_hfd_results[["lsGonXX"]])
