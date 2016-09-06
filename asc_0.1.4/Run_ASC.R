#Run_ASC.R
#load ASC code
source("./asc_0.1.4.R")

#Load libraries
library(NOISeq) # simple TMM normalization
library(plyr)

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Chaetoceros-affinis-CCMP159_cds_counts_orthoID.txt"

#load the counts_orthoID.txt file
t= read.table(file, header=T) #open the file and read in data
t2=t[-which(is.na(t$Orthogroup)),] #get rid of genes that were not placed in Orthogroup
t.final=ddply(t2,"Orthogroup",numcolwise(sum)) #sum counts for contigs in same orthogroup

#make orthogroups into row names and delete Orthogroup column
rownames(t.final)<-t.final$Orthogroup
t.final <- t.final[,-which(colnames(t.final)=="Orthogroup")]

##Exploratory analyses
boxplot(t.final, ylim= c(0,15000))
plot(log(t.final[,4],base=2), log(t.final[,4]/t.final[,1], base=2))
#TMM normalize and recheck boxplot distributution
t.norm = tmm(t.final)
boxplot(t.norm, ylim=c(0,10000))
plot(log(t.norm[,4],base=2), log(t.norm[,4]/t.norm[,1], base=2))

X1 = t.final$MMETSP0088
names(X1) = rownames(t.final)
X2 = t.final$MMETSP0092
names(X2) = rownames(t.final)
boxplot(X1,X2, ylim=c(0,10000))

#This step does not work, not sure why, works on test data!
D <- DGE(X1,X2)