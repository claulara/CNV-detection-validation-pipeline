#!/usr/bin/Rscript
library(cn.mops)
library(optparse)


###################################################################################################
# Argument definition
option_list <- list(
  make_option(c("-b", "--bedfile"), type = "character", help = "Ruta al archivo BED"),
  make_option(c("-o", "--outputfile"), type = "character", default = "output.csv", help = "Nombre del archivo de salida"),
  make_option(c("-t", "--bam_test"), type = "character", help = "Rutas al bam de interés"),
  make_option(c("-c", "--bam_control"), type = "character", help = "Rutas a los archivos de los controles")
)

# Argument parsing
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Argument to variables
exon_bed <- opt$bedfile
output_file <- opt$outputfile

bam_test <- opt$bam_test
sample_name <- basename(bam_test)

bam_control <- strsplit(opt$bam_control, ",")[[1]]
bam_files <- c(bam_test, bam_control)



###################################################################################################
segments <- read.table(exon_bed,  
                       sep = "\t", 
                       as.is = TRUE)
gr <- GRanges(segments[,1], IRanges(segments[,2], segments[,3]))

# Read counts
X <- getSegmentReadCountsFromBAM(bam_files, GR=gr)

# Call CNVS
resCNMOPS <- exomecn.mops(X, normType = "poisson", qu = 0.99, minReadCount = 1)
resultExomeData <- calcIntegerCopyNumbers(resCNMOPS)
cnvs <- as.data.frame(cnvs(resultExomeData))

# Function to write a particular sample in cn.mops result GRanges to bed file.
write_sample_cnvs <- function(cnv_result, bam_test, output_file) {
  filtered_cnv <- cnv_result[cnv_result$sampleName %in% sample_name, ]
  write.csv(filtered_cnv, file = output_file, row.names = FALSE)
}

write_sample_cnvs(cnv_result = cnvs, bam_test = bam_test, output_file = output_file)










