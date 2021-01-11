### install DeepBlue 
#install.packages("BiocManager") 
#BiocManager::install("DeepBlueR")
library(DeepBlueR)

### Test the installation and connectivity.
deepblue_info("me")
#help(package="DeepBlueR")

### get familiar with some useful functions
#deepblue_list_genomes()
#grep("fetal heart", deepblue_list_biosources())
#grep("heart", deepblue_list_samples()) -> index
#deepblue_list_samples()[index]
#deepblue_list_epigenetic_marks()
#grep("DNaseI", deepblue_list_epigenetic_marks()$name) #73
#grep("narrowPeak", deepblue_list_experiments()$name) 
#deepblue_info("e16917")

### Search or List the database 
### Problems: 
### deepblue_search() can only return 50 items at most???
# user_key = "anonymous_key"
# experiments_found = deepblue_search(keyword="'H3K27ac', 'peaks'", type = "Experiments", user_key = "anonymous_key")
# 
# custom_table = do.call("rbind", apply(experiments_found, 1, function(experiment){
#   experiment_id = experiment[1]
#   # Obtain the information about the experiment_id
#   info = deepblue_info(experiment_id)
#   
#   # Print the experiment name, project, biosource, and epigenetic mark.
#   with(info, { data.frame(name = name, project = project, genome = genome,
#                           biosource = sample_info$biosource_name, epigenetic_mark = epigenetic_mark)
#   })
# }))

### List
# experiments = deepblue_list_experiments(genome="hg19", type="peaks", epigenetic_mark="DNaseI", biosource="heart", project="Roadmap Epigenomics")### list all assays with "narrowPeak"

### use a brutal way to search all "narrowPeak" files, returns 1032 files
narrowPeak.list = grep("narrowPeak", deepblue_list_experiments()$name)

## initialize vectors
id = c()
name = c()
type = c()
data_type = c()
description =c()
epigenetic_mark =c()
format =c()
genome =c()
object_name =c()
project = c()
sample_id = c()
technique =c()
mark = c()
Eid = c()
sex = c()
std_name = c()
age = c()
EDACC_NAME = c()
TYPE = c()
biosource_name = c()

### query the metadata information
for (i in 19:1032){
  j = narrowPeak.list[i]
  object = deepblue_list_experiments()[j]
  id[i] = object$id
  name[i] = object$name
  object = deepblue_info(id[i])
  type[i] = object$type
  data_type[i] = object$data_type
  description[i] = object$description
  epigenetic_mark[i] = object$epigenetic_mark
  format[i] = object$format
  genome[i] = object$genome
  object_name[i] = object$name
  project[i] = object$project
  sample_id[i] = object$sample_id
  technique[i] = object$technique
  mark[i] = object$extra_metadata$mark
  Eid[i] = object$extra_metadata$`roadmap epigenome`
  sex[i] = object$sample_info$SEX
  std_name[i] = object$sample_info$STD_NAME
  if(!is.null(object$sample_info$AGE)){
    age[i] = object$sample_info$AGE
  }else{
    age[i] = NA
  }
  EDACC_NAME[i] = object$sample_info$EDACC_NAME
  TYPE[i] = object$sample_info$TYPE
  biosource_name[i] = object$sample_info$biosource_name
}

### save the metadata information into assays.RData
assays = data.frame(
  "EID" = Eid,
  "SAMPLE_ID" = sample_id,
  "ID" = id,
  "OBJECT_NAME" = object_name,
  "PROJECT" = project,
  "TECHNIQUE" = technique,
  "TYPE" = TYPE,
  "DATA_TYPE" = data_type,
  "EDACC_NAME" = EDACC_NAME,
  "EPIGENETIC_MARK" = epigenetic_mark,
  "GENOME" = genome,
  "SEX" = sex,
  "AGE" = age,
  "TISSUE" = biosource_name
)
save(assays, file = "assays.RData")
