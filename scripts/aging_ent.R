## calculate H3K27ac entropy for each sample in
## /media/Data/zhangz/chip/analysis/summary2/aging/ihec/bigBed.txt
## these samples are from IHEC project, human brains, some are from young, some are from old
## the yound ones were probably illed
## the old ones were probably healthy

## the entropy is calculated as follows:
## 1. read the bigBed.bed file
## 2. discretize the H3K27ac signal for each bin
## 3. calculate the entropy for each sample

# load necessary libraries
library(entropy)
library(plyr)
library(doMC)
library(magrittr)
library(stringr)
doMC::registerDoMC(cores = 4)

# read the bigBed.bed file
bed_files <- read.table("/media/Data/zhangz/chip/analysis/summary2/aging/ihec/bedlist", header = FALSE, sep = "\t")
setwd("/media/Data/zhangz/chip/analysis/summary2/aging/ihec/")
entropies <- adply(bed_files$V1, 1, function(x) {
    # read the bigBed file
    bed <- read.table(paste0(getwd(), "/", x), header = FALSE, sep = "\t")
    name <- str_split(x, "[.]")[[1]][3]
    print(name)
    # discretize the H3K27ac signal
    bed_bin <- discretize(bed$V7, numBins = 25)
    # calculate the entropy
    ent <- entropy(bed_bin, unit = "log2")
    return(data.frame(acc = name, entropy = ent))
})
# young ones dose not have the signal value, failed

# try mouse data
setwd("/media/Data/zhangz/chip/analysis/summary2/aging/mouse_aging")
## these data are from GEO, GSE122867, article: https://onlinelibrary.wiley.com/doi/10.1111/acel.12996
## mouse muscles, young and old, 2,10,20 months
bed_files <- read.table("bedlist", header = FALSE, sep = "\t")
#   17884 GSM3487511_TA_2months.bed
#   35221 GSM3487512_TA_10months.bed
#   44752 GSM3487513_TA_20months.bed
# top 15000 peaks
entropies <- adply(bed_files$V1, 1, function(x) {
    # read the bigBed file
    bed <- read.table(paste0(getwd(), "/", x), header = FALSE, sep = "\t")
    name <- str_split(x, "_")[[1]][3] %>%
        gsub(".bed", "", .)
    print(name)
    # discretize the H3K27ac signal
    bed <- bed[order(bed$V7, decreasing = TRUE), ][1:15000, ]
    bed_bin <- discretize(bed$V7, numBins = 15000)
    # calculate the entropy
    ent <- entropy(bed_bin, unit = "log2")
    return(data.frame(age = name, entropy = ent))
})
write.csv(entropies, "mouse_muscle_entropy.csv", row.names = FALSE)
