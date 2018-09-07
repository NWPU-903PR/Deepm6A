# Deepm6A
A CNN based R function to identify single base m6A site from MeRIP-Seq peak region.

The input files of this example and the files of CNN model can be downloaded from  



####################################################################################################################
#################################################### Annotation ####################################################
####################################################################################################################
## package requiered: keras, exomePeak, Guitar, BSgenome.Hsapiens.UCSC.hg19
##
## parameters:
##
## filepath: the filepath where you save all the code and data downloaded from our GitHub link
##
## input: the input of the CNN model you generated using 'getinput' function
##
## ip_bam: the IP bam file of MeRIP-Seq data
##
## peak_bed: the bed file of peaks gets from result of exomePeak
##
## outfilepath: the filepath where you want to save the result
##
## inputfile: the input file of the CNN model you generated using 'getinput' function
##
## modelfile: the file path where you save all the CNNmodel files 
##
## The output of the function is a excel file. Each row is a candidate single base m6A site extracted from exomePeak
## peak region and each colum denoting chrome, chrome start, chrome end, Entrez gene ID, the probability to be a m6A 
## site, strand and the motif at this location separately.
#####################################################################################################################
#####################################################################################################################


## make input of the CNN model ---------------------------------------------------------------------------------------

filepath <- "F:/Deepm6A/Function"
outfilepath <- "F:/Deepm6A/example"

ip_bam <- "F:/Deepm6A/example/p001_HEK293T_S1_SYSY_ip/SRR494614_accept_hits.bam"
peak_bed <- "F:/Deepm6A/example/p001_HEK293T_S1_SYSY_ip/exomePeak_output/peak.bed"

source(paste(filepath, "makeInputFromBamBed.R", sep = "/"))
input <- getinput(filepath, ip_bam, peak_bed, outfilepath)


## get single base m6A sites -----------------------------------------------------------------------------------------

filepath <- "F:/Deepm6A/Function"
outfilepath <- "F:/Deepm6A/example"

inputfile <- "F:/Deepm6A/example/input.RData"
modelfile <- "F:/Deepm6A/CNNmodel"
load(inputfile)

library(keras)
library(exomePeak)
source(paste(filepath, "getprobsFun.R", sep = "/"))

peak <- getpredict(input, modelfile, outfilepath)
