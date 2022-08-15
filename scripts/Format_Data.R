library(data.table)
library(stringr)
library(tibble)

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

#############################################################################
#############################################################################
## Get Clinical data

clin = read.table( file.path(input_dir, "CLIN.txt") , sep="\t" , header=TRUE,  stringsAsFactors = FALSE )
clin = clin[ grep( "baseline" , clin$title ) , ]

clin_original <- clin
selected_cols <- c( "geo_accession" , "response_ch1" ) 
clin = as.data.frame( cbind( clin[ , selected_cols ] , "PD-1/PD-L1" , "Esophageal" , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA ) )
colnames(clin) = c( "patient" , "response" , "drug_type" , "primary" , "recist" , "age" , "histo" , "response" , "pfs" ,"os" , "t.pfs" , "t.os" , "stage" , "sex" , "response.other.info" , "dna" , "rna" )

rownames(clin) = clin$patient

clin$response = ifelse( clin$response %in% "responder" , "R" , ifelse( clin$response %in% "non-responder" , "NR" , NA ) )

clin$rna = "tpm"
clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

clin <- format_clin_data(clin_original, 'geo_accession', selected_cols, clin)

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='VanDenEnde', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, unique_drugid=annotate_drug('VanDenEnde', str_extract(clin$characteristics_ch1_1, '(?<=treatment: ).*'), annotation_drug), .after='unique_tissueid')

#############################################################################
#############################################################################

# Expression data
load(file.path(annot_dir, 'Gencode.v19.annotation.RData'))
expr <- read.table(file.path(input_dir, 'EXPR.txt.gz'), sep='\t', header=TRUE, stringsAsFactors = FALSE)
expr <- expr[rownames(expr) %in% rownames(features_gene), ]

## Compute TPM data
genes <- features_gene[rownames(features_gene) %in% rownames(expr), c('start', 'end', 'gene_id')]
size <- genes$end - genes$start
names(size) <- rownames(genes)

GetTPM <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

expr_tpm = log2( GetTPM(expr, size) + 0.001 )

# patient = intersect( colnames(expr) , rownames(clin) )
# clin = clin[ patient , ]
# expr =  expr[ , patient ]

# cased_sequenced
case = cbind( clin$patient , 0 , 0 , 1 )
colnames(case ) = c( "patient" , "snv" , "cna" , "expr" )

write.table( case , file = file.path(output_dir, "cased_sequenced.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( clin , file = file.path(output_dir, "CLIN.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( expr_tpm , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
