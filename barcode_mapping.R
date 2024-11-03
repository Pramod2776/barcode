setwd("/path/to/workdirectory/")

## load libraries
library("Biostrings")
library("tidyverse")

## read Fastq file from the  "/path/to/workdirectory/"
reference_R1 = readDNAStringSet("./R1_R1.fastq.gz", format="fastq", 
                                with.qualities=TRUE)

## save fastq file as a variable
dna_initial1 = reference_R1

## read barcodes (Genes, sequences) from the excel sheet in "/path/to/workdirectory/"
sequences = readxl::read_xlsx("./Barcodes.xlsx", sheet =1) %>%
  as.data.frame()


## change sequences to Upper case
sequences$Barcode = toupper(sequences$Barcode)
print(length(sequences$Gene))

## match each barcode sequence (sequences[[2]][1]) with the reads in the fastq file (dna_initial1)
## match barcode sequences to reads
list = lapply(1:length(sequences$Gene), function(num){
  num = 1
  print(num)
  print(sequences[[1]][num])
  print(sequences[[2]][num])
  
  ## find match hits
  mah1 = vmatchPattern(sequences[[2]][num], dna_initial1, max.mismatch=0, fixed=FALSE)
  hits.frame1 = as.data.frame(mah1)
  print(dim(hits.frame1))
  
  ##extract sequences for the matched hits
  
  dna_list1 <- mapply(subseq, dna_initial1[hits.frame1$group], hits.frame1$start, hits.frame1$end)
  dna_regions1 <- DNAStringSet(dna_list1)
  dna_region_df1 = data.frame(dna_regions1)
  names(dna_regions1)
  as.data.frame(dna_initial1[names(dna_regions1)])
  dat_r1 = data.frame(dna_initial1[names(dna_regions1)])
  
  
  ## column bind hits with corresponding sequences
  dat_r1 = cbind(dat_r1, dna_regions1, hits.frame1)
  dat_r1$group_name = NULL
  colnames(dat_r1)
  rownames(dat_r1)
  dat_r1$name = rownames(dat_r1)
  
  ## rearrange the columns
  dat_r1 = dat_r1 %>%
    dplyr::select("name", everything())
  ncol(dat_r1)
  
  ## rename the colnames
  colnames(dat_r1) = c("name", "sequence", "sequence_match", "line_number", "start", "end", "width")
  
  ## write the file as excel sheet
  
  file_name = paste0("R1_R1","_", sequences[[1]][num],"_", sequences[[2]][num], ".xlsx")
  print(file_name)
  writexl::write_xlsx(dat_r1, path = file_name)

  ## return output as a dataframe
  df = data.frame(gene = sequences[[1]][num],
                           sequence = sequences[[2]][num],
                         count = nrow(dat_r1),
                         max.mismatch = 0)
  print(df)
  
  
}) %>%
  bind_rows() %>%
  as.data.frame()

sum(list$count)

## write it to a file as pooled summary file (all barcodes match summary to the input fastq file)
file_name = paste0("R1_R1","_summary", ".xlsx")
writexl::write_xlsx(list, path = file_name)


