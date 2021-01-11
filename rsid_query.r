### rsid.txt is a list from a combined GWAS study 
### (80 rsids are related to HF traits: 40 lead variants + 40 candidate variants)
### unique(rsid) to remove the repeated rsids, 62 left
### load all rsids to rsid

### package installation
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("biomaRt")

### start to query the chr and loci of rsids
library(biomaRt)
library(dplyr)
#listMarts()
#rsid = c("rs12541595", "rs35006907", "rs34866937")

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

get_chr_loci_byRsidFile("candidate_rsid.txt") -> rsid_chr_loci_id_symbol
#load('test.RData')

library(dplyr)

get_assays <- function(assays, file){
  histone_marks = read.table(file, header = FALSE, sep = "")
  histone_marks = pull(histone_marks, V1) 
  experiment = subset(assays, assays$EPIGENETIC_MARK %in% histone_marks)
  return(experiment)
}

x = get_assays(assays, "experiment.txt")
