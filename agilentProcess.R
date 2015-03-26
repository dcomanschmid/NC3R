##################################################################################
# Analyze Agilent dual Cy3/Cy5 channel microarrays as Cy3 and Cy5 single channel #
#    >> personal communication & use scripts from                                #
#                https://github.com/MWSchmid/microarray                          #
#                "agilentDualChannelExample.R"                                   #
# Diana Coman Schmid                                                             #
# Eawag 2015                                                                     #
# diana.comanschmid@eawag.ch                                                     #
##################################################################################

library("limma")

workDir <- "..." 

# read an annotation table
annotation <- read.table(file.path(workDir, "annotation.txt"), header = TRUE, sep = '\t', quote = "", stringsAsFactors = FALSE)
annotation$addName <- paste(file.path(workDir, annotation$File), annotation$Channel, sep = '_')
if (length(unique(annotation$Name)) != nrow(annotation)) {stop("The names in the annotation file must be unique!")}


# read raw data and create a raw data list which mimics a one-channel raw data list
# the column ProbeName must exist in the raw microarray data file

twoChannelRawData <- read.maimages(unique(file.path(workDir, annotation$File)), source = "agilent", annotation = c("ProbeName"))

colnames(twoChannelRawData$G) <- paste(colnames(twoChannelRawData$G), ".txt_Cy3", sep = '')
colnames(twoChannelRawData$Gb) <- paste(colnames(twoChannelRawData$Gb), ".txt_Cy3", sep = '')
colnames(twoChannelRawData$R) <- paste(colnames(twoChannelRawData$R), ".txt_Cy5", sep = '')
colnames(twoChannelRawData$Rb) <- paste(colnames(twoChannelRawData$Rb), ".txt_Cy5", sep = '')
summary(twoChannelRawData)

greenTargets <- twoChannelRawData$targets; rownames(greenTargets) <- paste(rownames(twoChannelRawData$targets), ".txt_Cy3", sep = '')
redTargets <- twoChannelRawData$targets; rownames(redTargets) <- paste(rownames(twoChannelRawData$targets), ".txt_Cy5", sep = '')

rawDataE <- matrix(0, nrow = nrow(twoChannelRawData$G), ncol = nrow(annotation), dimnames = list(rownames(twoChannelRawData$G), annotation$Name))
rawDataEb  <- matrix(0, nrow = nrow(twoChannelRawData$Gb), ncol = nrow(annotation), dimnames = list(rownames(twoChannelRawData$G), annotation$Name))

rawDataE[,annotation$Name[annotation$Channel == "Cy3"]] <- twoChannelRawData$G[,annotation$addName[annotation$Channel == "Cy3"]]
rawDataE[,annotation$Name[annotation$Channel == "Cy5"]] <- twoChannelRawData$R[,annotation$addName[annotation$Channel == "Cy5"]]
rawDataEb[,annotation$Name[annotation$Channel == "Cy3"]] <- twoChannelRawData$Gb[,annotation$addName[annotation$Channel == "Cy3"]]
rawDataEb[,annotation$Name[annotation$Channel == "Cy5"]] <- twoChannelRawData$Rb[,annotation$addName[annotation$Channel == "Cy5"]]
rawData <- list(E = rawDataE, Eb = rawDataEb, targets = rbind(redTargets, greenTargets), source = "agilent", genes = twoChannelRawData$genes)


# background correction and normalization

bgData <- backgroundCorrect(rawData,method="normexp")
normData <- normalizeBetweenArrays(bgData, method = "quantile")
normEset <- normData$E; rownames(normEset) <- normData$genes$ProbeName

# load the locus ID to probe name mappings (based on the latest genome annotation Zv9)

LOCUStoPNtable <- read.table(file.path(workDir, "IDtoProbeName.txt"), header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE, row.names = 1)
colnames(LOCUStoPNtable) <- c("probeNames")
LOCUStoPN <- sapply(LOCUStoPNtable$probeNames, function(x) strsplit(x, ";", fixed = TRUE)); names(LOCUStoPN) <- rownames(LOCUStoPNtable)

# summarize the probes per locus ID (if there are several, take the mean)

locusData <- list()
for (sample in colnames(normEset)) {
	cat(paste(sample, '\n', sep = ''))
	sub <- normEset[,sample]
	locusData[[sample]] <- sapply(LOCUStoPN, function(x) ifelse(length(x)>1, mean(sub[x]), sub[x])) 
}

# normalized expression data
eset <- do.call("cbind", locusData)
