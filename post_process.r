library(dplyr)
library(tidyverse)
source("deepblue_query.r")

### annotate the snp variants
get_chr_loci_byRsidFile("candidate_rsid.txt") -> rsid_chr_loci_id_symbol
### load the query result file
query_result = read.table("final.txt",
                          colClasses = c("character", "integer", "integer","integer", "character", "integer", "character", "character"))
colnames(query_result) = c("rsid", "chr", "start", "end", "OBJECT_NAME", "n_peak", "peak_value", "signal_value")
save(query_result, file="query_result.RData")
load("query_result.RData")

generate_bigTable <- function(result_file, experiment_file, rsid_file, assays){
  load(result_file)
  #load("query_result.RData")
  load(assays)
  #load("assays.RData")
  query_rsid = read.table(rsid_file)
  num1 = length(unique(query_rsid$V1))
  ### load the query experiment file
  #experiment = get_assays(assays, experiment_file)
  experiment = get_assays(assays, experiment_file)
  ### inner join the file by narrowPeak file name
  df <- inner_join(query_result, experiment, by = c("OBJECT_NAME")) 
  df <- df %>% mutate(label1 = paste(EPIGENETIC_MARK, TISSUE, chr, start, end, sep = "_"),
                      label1 = as.factor(label1),
                      label2 = paste(EPIGENETIC_MARK, TISSUE, sep = "_"),
                      label2 = as.character(label2))
  df <- df[which(df$rsid %in% query_rsid$V1),]
  num2 = length(unique(df$rsid))
  print(paste("There are", num1, "rsids have been queried!"))
  print(paste("There are", num2, "rsids have been identified!"))
  ### calculate if there is a peak or peaks in the query region
  df$peak_or_not <- ifelse(df$n_peak == 0, 0, 1)
  ### calculate the peak_ave 
  df$peak_ave <- sapply(strsplit((df$peak_value), split=","), function(x) 
    mean(as.numeric(x)), simplify = "array")
  df$peak_ave <- ifelse(df$n_peak == 0, 0, df$peak_ave)
  ### calculate the mean signal_ave
  df$signal_ave <- sapply(strsplit((df$signal_value), split=","), function(x) 
    mean(as.numeric(x)), simplify = "array")
  df$signal_ave <- ifelse(df$n_peak == 0, 0, df$signal_ave)
  ### calculate the max signal_ave
  df$signal_max <- sapply(strsplit((df$signal_value), split=","), function(x) 
    mean(as.numeric(x)), simplify = "array")
  df$signal_max <- ifelse(df$n_peak == 0, 0, df$signal_max)
  return(df)
}

df = generate_bigTable("query_result.RData", "experiment.txt", "candidate_rsid.txt", "assays.RData")

### value_key = (1: peak_or_not, 2: n_peak, 3: peak_ave, 4: signal_ave)
generate_heatmap <- function(df, experiment){
  df_long <- df %>% dplyr::select(rsid, label1, label2, EPIGENETIC_MARK, signal_max) # need to change
  df_long <- df_long[order(df_long$rsid, df_long$label1),]
  assay_tissue <- unique(df_long$label1)
  df_long$label1 <- as.character(df_long$label1)
  rsid = c()
  value = c()
  label2 = c()
  epimark = c()
  for (i in c(1:length(assay_tissue)[1])){
    data = filter(df_long, df_long$label1 == assay_tissue[i])
    value[i] = mean(data$signal_max) # need to change
    rsid[i]=data$rsid[1]
    label2[i] = data$label2[1]
    epimark[i] = data$EPIGENETIC_MARK[1]
    }
  new_df = data.frame(rsid = rsid, label2 = label2, epimark = epimark, value = value)
  new_df <- new_df %>% dplyr::filter(epimark == experiment) %>%
                      dplyr::select(rsid, label2, value)
  row_names = unique(new_df$rsid)
  gene_names = rsid_chr_loci_id_symbol[order(match(rsid_chr_loci_id_symbol$rsid, row_names)),5, drop=FALSE]
  gene_names$hgnc_symbol <- ifelse(gene_names$hgnc_symbol == '', NA, gene_names$hgnc_symbol)
  row_names = paste(row_names, gene_names$hgnc_symbol, sep = '_')
  
  df_short <- new_df %>% spread(label2, value) %>%
             dplyr::select(-rsid)
  df_short <- data.matrix(df_short, rownames.force = NA)
  df_short <- t(df_short)
  #heatmap(df_short, scale = "column")
  #heatmap(df_short)
  png(paste0(experiment,"_40candidate_rsid_signal_max.png"), width = 800, height = 800, units = "px", pointsize = 12)
  par(oma=c(2,1,1,9))
  heatmap(df_short, labCol  = row_names)
  dev.off()
}

experiment = read.table("experiment.txt")
x = experiment$V1
x = x[!x %in% c("H3K23me2", "H2AK9ac","H4K12ac","H3T11ph")]
for (i in x){
    generate_heatmap(df, i)
    print(paste(i,"experiment is done!"))
}

22, 26,31
