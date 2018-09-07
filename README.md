# Deepm6A
A CNN based R function to identify single base m6A site from MeRIP-Seq peak region.

See example.R for the usage of Deepm6A. The input of Deepm6A are MeRIP-Seq IP sample bam file and bed file which annotates the peaks identified by exomePeak R package from MeRIP-Seq data. The output is a excel file. Each row is a candidate single base m6A site extracted from exomePeak peak region and each colum denoting chrome, chrome start, chrome end, Entrez gene ID, the probability to be a m6A site, strand and the motif at this location separately.

The input files of the example and the files of CNN model can be downloaded from
https://drive.google.com/drive/folders/16Jnv-LDyFq64BmMcOvdT6fTI8GaS-5mv?usp=sharing
