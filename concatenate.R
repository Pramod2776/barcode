setwd("~/path/to/workdir/")
library(tidyverse)
forward_summary_files = list.files(pattern = "_summary.xlsx")

list = lapply(1:length(forward_summary_files), function(l){
  print(forward_summary_files[l])
  
  name = gsub("_forward_barcode_summary.xlsx", "", forward_summary_files[l])
  summary_file = readxl::read_xlsx(forward_summary_files[l]) %>%
    as.data.frame()
  
  df = data.frame(gene = summary_file$gene,
                  name = summary_file$count)
  colnames(df)[2] = name
  print(df)
  
}

)  %>%
  bind_cols() %>%
  as.data.frame()

list_sort = list %>% select(-contains('gene'))
list_sort$gene = rownames(list_sort)
rownames(list_sort) = NULL

list_sort = list_sort %>%
  dplyr::select("gene", everything())
writexl::write_xlsx(list_sort, "forward_barcode_summary.xlsx")


reverse_summary_files = list.files(pattern = "_summary.xlsx")
reverse_summary_files[1]

list = lapply(1:length(reverse_summary_files), function(l){
  print(reverse_summary_files[l])
  
  name = gsub("_reverse_barcode_summary.xlsx", "", reverse_summary_files[l])
  summary_file = readxl::read_xlsx(reverse_summary_files[l]) %>%
    as.data.frame()
  
  df = data.frame(gene = summary_file$gene,
                  name = summary_file$count)
  colnames(df)[2] = name
  print(df)
  
}

)  %>%
  bind_cols() %>%
  as.data.frame()

list_sort = list %>% select(-contains('gene'))
list_sort$gene = rownames(list_sort)
rownames(list_sort) = NULL
list_sort = list_sort %>%
  dplyr::select("gene", everything())
writexl::write_xlsx(list_sort, "reverse_barcode_summary.xlsx")

data = readxl::read_xlsx("./forward_barcode_summary_9.xlsx", sheet = 2) %>%
  as.data.frame()
df = data[1:44,]

df[, 2:10] = apply(df[, 2:10],2,function(x){x/sum(x) * 100})
writexl::write_xlsx(df, "Barcode_summary_percent_100_9_samples.xlsx")

