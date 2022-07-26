# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165252/matrix/GSE165252_series_matrix.txt.gz -O ~/Documents/ICBCuration/data_source/VanDenEnde/GSE165252_series_matrix.txt.gz
# wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE165252&format=file&file=GSE165252%5Fnorm%2Ecnt%5FPERFECT%2Etxt%2Egz -O ~/Documents/ICBCuration/data_source/VanDenEnde/GSE165252_norm.cnt_PERFECT.txt.gz

library(GEOquery)
library(Biobase)
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

# CLIN.txt
gunzip(file.path(work_dir, "GSE165252_series_matrix.txt.gz"))
clin <- getGEO(filename=file.path(work_dir, 'GSE165252_series_matrix.txt'), destdir=work_dir)
clin <- pData(clin)
colnames(clin) <- str_replace_all(colnames(clin), '\\W', '_')
clin = clin[ grep( "baseline" , clin$title ) , ]
write.table(clin, file=file.path(work_dir, 'CLIN.txt'), sep = "\t" , quote = FALSE , row.names = TRUE)

file.remove(file.path(work_dir, 'GPL20301.soft'))
file.remove(file.path(work_dir, 'GSE165252_series_matrix.txt'))

# EXPR.txt.gz
expr <- read.table(file.path(work_dir, 'GSE165252_norm.cnt_PERFECT.txt.gz'), sep='\t', header=TRUE)
rownames(expr) <- expr[, 1]
expr <- expr[, -1]
expr <- expr[, colnames(expr)[colnames(expr) %in% str_replace_all(clin$description, '-', '.')]]
colnames(expr) <- unlist(lapply(colnames(expr), function(col){
  return(clin[str_replace_all(clin$description, '-', '.') == col, 'geo_accession'])
}))

gz <- gzfile(file.path(work_dir, 'EXPR.txt.gz'), "w")
write.table( expr , file=gz , quote=FALSE , sep="\t" , col.names=TRUE, row.names = TRUE )
close(gz)
