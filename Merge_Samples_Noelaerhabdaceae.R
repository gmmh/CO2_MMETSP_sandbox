#Merge_Samples.R
# This script will take the orthoID count files from a directory and output a file merged by orthoID where counts are summed by orthogroup
library(plyr)

#find all filenames from count directory with "orthoID"
dir="~/Desktop/CO2_MMETSP/Haptophyta/counts/"
files=list.files(dir, "orthoID")
length(files)

#loop to open files, discard unassigned contigs, sum reads for othrogroups and merge
m=NA
for(file in files){
	t= read.table(paste0(dir,file), header=T) #open the file and read in data
	t2=t[-which(is.na(t$Orthogroup)),] #get rid of genes that were not placed in Orthogroup
	t.final=ddply(t2,"Orthogroup",numcolwise(sum)) #sum counts for contigs in same orthogroup
	if(is.na(m)){m=t.final}else{m= merge(t.final,m, all=T)} #merge samples into one table
}

#make orthogroups into row names and delete Orthogroup column
rownames(m)<-m$Orthogroup
m <- m[,-which(colnames(m)=="Orthogroup")]

#keep only Noelaerhabdaceae samples
calc <- c("MMETSP1006", "MMETSP1007", "MMETSP1008", "MMETSP1009", "MMETSP0994", "MMETSP0995", "MMETSP0996", "MMETSP0997", "MMETSP1154", "MMETSP1155", "MMETSP1156", "MMETSP1157", "MMETSP1150", "MMETSP1151", "MMETSP1152", "MMETSP1153", "MMETSP1363", "MMETSP1364", "MMETSP1365", "MMETSP1366")
m <- m[,which(colnames(m) %in%  calc)]

#trim orthogroups found in <50% of samples
num.na =apply(m,1,function(x){length(which(is.na(x)))}) 
cut = which(num.na > ncol(m)/2)
m = m[-cut,]
dim(m)
#write new file
write.table(m,paste0(dir,"merged_counts_calcifiers.txt"))

#Cut table down to shared "core" orthologous groups
cut <- which(is.na(apply(m,1,mean)))
m.core = m[-cut,]
dim(m.core)
write.table(m.core,paste0(dir,"merged_counts_core_calcifiers.txt"))


