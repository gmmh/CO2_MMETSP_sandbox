# norm_RNAseq_data.R
# Normalize RNA seq data to combine with Microarray data
# Input: raw counts (rows = gene IDs, columns = conditions), which condition(s) = Control
# Output: normalized log2 ratios (rows= gene IDs, columns = conditions)

# Steps:
# 1. Normalize reads to eliminate batch effects (TMM or other method)
# 2. compute log2 ratios to the control sample(s)
# 3. Remove batch effects with sva? or save this step for RNAseq + microarray data


#Load libraries
library(NOISeq) # simple TMM normalization

#Load merged RNASEq Data
dir <- "~/Desktop/CO2_MMETSP/Coscinodiscophyceae/counts/"
m <- read.table(paste0(dir,"merged_counts_core.txt"), header=T)
control.name =  "MMETSP0088" #"MMETSP0143" #"MMETSP1363"

for(col in colnames(m)){
	replace <- which(is.na(m[,col]))
	m[replace,col] <- 0
}
##Exploratory analyses
boxplot(m, ylim= c(0,15000))
	
## Normalization Step
TMM <- tmm(m, long =1000, lc=0, k=0.4) #TMM norm
	#Visualize
boxplot(TMM, ylim=c(0,10000))

## Calculate log2 ratios with Control
control.mean <- TMM[,grep(control.name, colnames(TMM))]
#control.sd <- apply(TMM[,grep("seastar", colnames(TMM))], 1, sd)
norm.data <- log(sweep(TMM, 1, control.mean, "/"), 2) #Divide by control.mean and log2 transform
norm.data <- norm.data[ ,-grep(control.name, colnames(norm.data))] # Remove control samples
boxplot(norm.data)

write.table(norm.data, paste0(dir,"norm_RNAseq_log2ratios.tsv"), sep="\t")
