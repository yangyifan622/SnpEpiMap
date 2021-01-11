library(DeepBlueR)
library(biomaRt)
library(dplyr)

######################### usage: ####################
# get_chr_loci_byRsidFile(file) -> rsid_chr_loci
# get chromosome_name and loci from biomaRt

get_chr_loci_byRsidFile <- function(file){
  
  rsid = read.table(file, header = FALSE, sep = "")
  rsid_uniq = unique(rsid$V1)
  ### annotate the chr, loci and ensembl_id of rsid
  snp_grch37 = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice")
  #datasets <- listDatasets(snp_grch37)
  #snp_grch37 = useMart("ENSEMBL_MART_SNP",dataset="hsapiens_snp")
  snp_grch37 = useDataset("hsapiens_snp",mart=snp_grch37)
  
  #filters = listFilters(snp_grch37)
  #attributes = listAttributes(snp_grch37)
  rsid_chr_loci <- getBM(attributes=c("refsnp_id", "chr_name", "chrom_start","ensembl_gene_stable_id"), 
                         filters = "snp_filter", 
                         values = rsid_uniq, 
                         mart = snp_grch37,
                         useCache = FALSE)
  
  ### remove the rsids if they are not from the regular chromosomes
  label = !grepl("_", rsid_chr_loci$chr_name)
  rsid_chr_loci[label,] -> rsid_chr_loci
  colnames(rsid_chr_loci) = c("rsid", "chr", "loci", "ensembl_gene_id")
  #save(rsid_chr_loci, file = "test.RData")
  
  ### annotate the gene_symbol by ensembl_id
  ensembl = useMart(biomart="ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  ensembl_id_gene_symbol <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                                  filters = "ensembl_gene_id", 
                                  values = rsid_chr_loci$ensembl_gene_id,
                                  mart = ensembl,
                                  useCache = FALSE)
  
  #save(ensembl_id_gene_symbol, file = "ensembl_id_gene_symbol.RData")
  #load("ensembl_id_gene_symbol.RData")
  
  ### merge the two dataframes
  rsid_chr_loci_id_symbol = merge(rsid_chr_loci, ensembl_id_gene_symbol, by = "ensembl_gene_id", all.x = T)
  ### remove duplicated rows based on rsid
  rsid_chr_loci_id_symbol <- rsid_chr_loci_id_symbol %>% distinct(rsid, .keep_all = TRUE) 
  
  return(rsid_chr_loci_id_symbol)
}

######################### usage: ####################
# get_assays(assays, file) -> experiment
# subset all experiments from assays.RData

get_assays <- function(assays, file){
  histone_marks = read.table(file, header = FALSE, sep = "", colClasses = c("character"))
  histone_marks = pull(histone_marks, V1) 
  experiment = subset(assays, assays$EPIGENETIC_MARK %in% histone_marks)
  return(experiment)
}

######################### usage: ####################
# rsid_exp_index(assays, start_index) -> c(rsid_index, exp_index)

rsid_exp_index <- function(assays, start_index){
  n_experiment = dim(assays)[1]
  rsid_index = (start_index - 1) %/% n_experiment + 1 
  exp_index = (start_index - 1) %% n_experiment + 1
  rsid_exp_index = c(rsid_index, exp_index)
}

######################### usage: ######################
# query(experiment_name, rsid, chr, loci, search_range) -> final.txt 

query <- function(experiment_name, rsid, chr, loci, search_range){
    ## start the upper bound of searching range
    ## end is the lower bound of searching range
    start = loci - 0.5 * search_range
    end = loci + 0.5 * search_range
    chromosome = paste0("chr", chr)
    query_id = deepblue_select_experiments(
      experiment_name = c(experiment_name),
      chromosome=chromosome, start=start, end=end
    )
    request_id = deepblue_get_regions(
      query_id=query_id, 
      output_format="CHROMOSOME,START,END,SIGNAL_VALUE,PEAK,@NAME,@BIOSOURCE"
    )
    check <- tryCatch({ 
      regions = deepblue_download_request_data(request_id=request_id)}, 
      warning = function(war){ 
        print(paste("MY WARNING: ", war)) 
        peak = 0}
      # error = function(err) { 
      #   print(paste("MY ERROR: ", err))
      #   peak = 1},  
      # finally = function(fun){ 
      #   print(paste("regions: ", regions, "has been processed!"))
      # }
    ) 
    # End of tryCatch 
    print(paste("Result: ", check)) 
    ### write the output for each SNP queried
    if (isS4(check)){
      num_regions = length(regions$PEAK)
      peak_value = paste0(regions$PEAK, collapse = ",")
      signal_value = paste0(regions$SIGNAL_VALUE, collapse = ",")
    }else{
      num_regions = 0
      peak_value = NA
      signal_value = NA
    }
    result_list = list()
    result_list <- append(result_list, rsid)
    result_list <- append(result_list, chr)
    result_list <- append(result_list, start)
    result_list <- append(result_list, end)
    result_list <- append(result_list, experiment_name)
    result_list <- append(result_list, num_regions)
    result_list <- append(result_list, peak_value)
    result_list <- append(result_list, signal_value)
    #remove(regions)
    write.table(result_list, file = "test.txt", append = TRUE, quote=FALSE, sep = "\t", row.names = F, col.names = F)
    Sys.sleep(1)
}
