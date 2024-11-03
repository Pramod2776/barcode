setwd("~/path/to/workdir/")
library(tidyverse)  
library(Biostrings)
library(writexl)
reference_combined = readDNAStringSet("./pass/co.fastq", format="fastq", with.qualities=TRUE)
dna_initial1 = reference_combined
indexes = readxl::read_xlsx("./pass/illumina_indexing.xlsx", sheet =1) %>%
  as.data.frame()
indexes$i7sequence = toupper(indexes$i7sequence)

sequences = indexes$i7sequence

## Generate Reverse complement of the indexes
list_rc = lapply(1:length(sequences), function(l){
  
  print(sequences[l])
  rc = as.character(reverseComplement(DNAString(sequences[l])))
  data.frame(reversecomp = rc)
  
}
) %>%
  bind_rows()


Gene_barcodes = readxl::read_xlsx("./pass/Barcodes.xlsx", sheet =1) %>%
  as.data.frame()

Gene_barcodes$Barcode = toupper(Gene_barcodes$Barcode)

## map reverse complement of indexes
list_reverse_comp_indexes = lapply(1:length(list_rc$reversecomp), function(num){
 
  mah1 = vmatchPattern(list_rc$reversecomp[num], dna_initial1, max.mismatch=0, fixed=FALSE)
  hits.frame1 = as.data.frame(mah1)

  dna_list1 <- mapply(subseq, dna_initial1[hits.frame1$group], hits.frame1$start, hits.frame1$end)
  dna_regions1 <- DNAStringSet(dna_list1)
  dna_region_df1 = data.frame(dna_regions1)
  as.data.frame(dna_initial1[names(dna_regions1)])
  dat_r1 = data.frame(dna_initial1[names(dna_regions1)])
  dat_r1 = cbind(dat_r1, dna_regions1, hits.frame1)
  dat_r1$group_name = NULL
  colnames(dat_r1)
  rownames(dat_r1)
  dat_r1$name = rownames(dat_r1)
  dat_r1 = dat_r1 %>%
    dplyr::select("name", everything())

  colnames(dat_r1) = c("name", "sequence", "sequence_match", "line_number", "start", "end", "width")
 
  file_name = paste0("combined.fastq","_", indexes$`Sample ID`[num],"_", list_rc$reversecomp[num], ".xlsx")
  print(file_name)
  
 
  
  writexl::write_xlsx(dat_r1, path = paste0(file.path("./pass_analysis/"), file_name))
  
  index_mapped = dna_initial1[names(dna_regions1)]
  dna_initial1
  index_mapped
  
  list_reverse_comp_indexes_genes = lapply(1:length(Gene_barcodes$Barcode), function(g){

    
    mah1_gene = vmatchPattern(Gene_barcodes$Barcode[g], index_mapped, max.mismatch=0, fixed=FALSE)
    hits.frame1_gene = as.data.frame(mah1_gene)
  
    dna_list1_gene <- mapply(subseq, index_mapped[hits.frame1_gene$group], hits.frame1_gene$start, hits.frame1_gene$end)
    dna_regions1_gene <- DNAStringSet(dna_list1_gene)
    dna_region_df1_gene = data.frame(dna_regions1_gene)
    names(dna_regions1_gene)
    head(as.data.frame(index_mapped[names(dna_regions1_gene)]))
    dat_r1_gene = data.frame(index_mapped[names(dna_regions1_gene)])
    dat_r1_gene
    dat_r1_gene = cbind(dat_r1_gene, dna_regions1_gene, hits.frame1_gene)
    dat_r1_gene$group_name = NULL
    colnames(dat_r1_gene)
    rownames(dat_r1_gene)
    dat_r1_gene$name = rownames(dat_r1_gene)
    dat_r1_gene = dat_r1_gene %>%
      dplyr::select("name", everything())
    ncol(dat_r1_gene)
    colnames(dat_r1_gene) = c("name", "sequence", "sequence_match", "line_number", "start", "end", "width")
    
    file_name_gene = paste0("combined.fastq","_", indexes$`Sample ID`[num],"_", list_rc$reversecomp[num],"_",Gene_barcodes$Barcode[g],".xlsx")
    print(file_name_gene)
    writexl::write_xlsx(dat_r1_gene, path = paste0(file.path("./pass_analysis/"), file_name_gene))
    
    
    df_gene = data.frame(Sample = Gene_barcodes$Gene[g],
                    Gene_barcode = Gene_barcodes$Barcode[g],
                    count = nrow(dat_r1_gene),
                    max.mismatch = 0)
    print(df_gene)
    
    
    
  }
    
  ) %>%
    bind_rows() %>%
    as.data.frame()

  file_name_gene_sum = paste0("Comb_Indexes_reverse_complement_mapping_genes","_summary", ".xlsx")
  writexl::write_xlsx(list_reverse_comp_indexes_genes, path = paste0(file.path("./pass_analysis/"),file_name_gene_sum))
  
    
  
  
  df = data.frame(Sample = indexes$`Sample ID`[num],
                  index_reverse_complement = list_rc$reversecomp[num],
                  count = nrow(dat_r1),
                  max.mismatch = 0)
  print(df)
  
  
}) %>%
  bind_rows() %>%
  as.data.frame()

sum(list_reverse_comp_indexes$count)
file_name = paste0("Combineed_Indexes_reverse_complement_mapping","_summary", ".xlsx")
writexl::write_xlsx(list_reverse_comp_indexes, path = paste0(file.path("./pass_analysis/"), file_name))




## map indexes 

list_indexes = lapply(1:length(sequences), function(num){

  
  mah1 = vmatchPattern(sequences[num], dna_initial1, max.mismatch=0, fixed=FALSE)
  hits.frame1 = as.data.frame(mah1)
  print(dim(hits.frame1))

  dna_list1 <- mapply(subseq, dna_initial1[hits.frame1$group], hits.frame1$start, hits.frame1$end)
  dna_regions1 <- DNAStringSet(dna_list1)
  dna_region_df1 = data.frame(dna_regions1)
  names(dna_regions1)
  as.data.frame(dna_initial1[names(dna_regions1)])
  dat_r1 = data.frame(dna_initial1[names(dna_regions1)])
  dat_r1 = cbind(dat_r1, dna_regions1, hits.frame1)
  dat_r1$group_name = NULL
  colnames(dat_r1)
  rownames(dat_r1)
  dat_r1$name = rownames(dat_r1)
  dat_r1 = dat_r1 %>%
    dplyr::select("name", everything())
  ncol(dat_r1)
  colnames(dat_r1) = c("name", "sequence", "sequence_match", "line_number", "start", "end", "width")
 
  file_name = paste0("combined.fastq","_", indexes$`Sample ID`[num],"_", sequences[num], ".xlsx")
  print(file_name)
  
  writexl::write_xlsx(dat_r1, path = paste0(file.path("./pass/"), file_name))
  
  df = data.frame(Sample = indexes$`Sample ID`[num],
                  index_sequence = sequences[num],
                  count = nrow(dat_r1),
                  max.mismatch = 0)
  print(df)
  
  
}) %>%
  bind_rows() %>%
  as.data.frame()

sum(list_indexes$count)
file_name = paste0("Combineed_Indexes_mapping","_summary", ".xlsx")
writexl::write_xlsx(list_reverse_comp_indexes, path = file_name)



