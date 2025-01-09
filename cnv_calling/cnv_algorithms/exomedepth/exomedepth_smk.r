#!/usr/bin/Rscript

library(ExomeDepth)
library(optparse)

###################################################################################################
# Argument definition
option_list <- list(
  make_option(c("-b", "--bedfile"), type = "character", help = "Ruta al archivo BED"),
  make_option(c("-r", "--referencefasta"), type = "character", help = "Ruta al archivo FASTA de referencia"),
  make_option(c("-o", "--outputfile"), type = "character", default = "output.csv", help = "Nombre del archivo de salida"),
  make_option(c("-f", "--bam_files"), type = "character", help = "Rutas a los archivos BAM"),
  make_option(c('--debug'), action='store_true', type='logical', default=FALSE, help='Debug mode')
)

# Argument parsing
opt <- parse_args(opt_parser <- OptionParser(option_list = option_list))
# Argument to variables
bed_file <- opt$bedfile
referenceFasta <- opt$referencefasta
output_file <- opt$outputfile
bam_files <- strsplit(opt$bam_files, ",")[[1]]

sample_names <- basename(bam_files)


###################################################################################################
# Data load
# Bed exons file
col_names = c("chromosome", 
              "start", 
              "end", 
              "gene")
bed_frame <- read.table(bed_file, 
                        header = FALSE, 
                        sep = "\t")

exome_df <- bed_frame[, 1:4]
colnames(exome_df) <- col_names
exome_df <- as.data.frame(exome_df)


###################################################################################################
# Read count data
# the function getBamCounts in ExomeDepth is set up to parse the BAM files
ExomeCount <- getBamCounts(bed.frame = exome_df,
                           bam.files =  bam_files,
                           include.chr = FALSE,
                           referenceFasta = referenceFasta)
ExomeCount.df <- as.data.frame(ExomeCount)
colnames(ExomeCount.df)[6:(5 + length(sample_names))] <- sample_names



###################################################################################################
# build the most appropriate reference set
my.test <- ExomeCount.df[, sample_names[1]]
my.ref.samples <- sample_names[-1]
my.reference.set <- as.matrix(ExomeCount.df[, my.ref.samples])
my.choice <- select.reference.set (test.counts = my.test,
                                   reference.counts = my.reference.set,
                                   bin.length = (ExomeCount.df$end - ExomeCount.df$start),
                                   n.bins.reduced = 10000)
my.matrix <- as.matrix( ExomeCount.df[, my.choice$reference.choice, drop = FALSE])
my.reference.selected <- apply(X = my.matrix,
                               MAR = 1,
                               FUN = sum)


###################################################################################################
# CNV calling
all.exons <- new('ExomeDepth',
                 test = my.test,
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.df$chromosome,
                      start = ExomeCount.df$start,
                      end = ExomeCount.df$end,
                      name = ExomeCount.df$exon)

write.csv(all.exons@CNV.calls, output_file)

