library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

load( file.path(input_dir, "Vandenende.RData") )

expr = tpm

#############################################################################
#############################################################################
## Get Clinical data

clin = as.data.frame( t( read.table( file.path(input_dir, "GSE165252_series_matrix.txt") , sep="\t" , header=FALSE , stringsAsFactors = FALSE ) ) )
colnames(clin) = clin[ 1 , ]
clin = clin[ grep( "baseline" , clin$sample.1 ) , ]

clin = as.data.frame( cbind( clin[ , c( "sample.2" , "response" ) ] , "PD-1/PD-L1" , "Esophageal" , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA ) )
colnames(clin) = c( "patient" , "response" , "drug_type" , "primary" , "recist" , "age" , "histo" , "response" , "pfs" ,"os" , "t.pfs" , "t.os" , "stage" , "sex" , "response.other.info" , "dna" , "rna" )

clin$patient = sapply( clin$patient , function( x ){ paste( unlist( strsplit( x , "-" , fixed = TRUE )) , collapse = "." ) } ) 

rownames(clin) = clin$patient

clin$response = ifelse( clin$response %in% "responder" , "R" , ifelse( clin$response %in% "non-responder" , "NR" , NA ) )

clin$rna = "tpm"
clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

#############################################################################
#############################################################################

patient = intersect( colnames(expr) , rownames(clin) )
clin = clin[ patient , ]
expr =  expr[ , patient ]

case = cbind( patient , 0 , 0 , 1 )
colnames(case ) = c( "patient" , "snv" , "cna" , "expr" )

write.table( case , file = file.path(output_dir, "cased_sequenced.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( clin , file = file.path(output_dir, "CLIN.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( expr , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )