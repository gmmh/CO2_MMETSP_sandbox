#Merge_Samples.R
# This script will take the orthoID count files from a directory and output a file merged by orthoID where counts are summed by orthogroup
library(plyr)

#find all filenames from count directory with "orthoID"
dir="~/Desktop/CO2_MMETSP/Ochrophyta/counts/"
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

#keep only raphidophyte samples
raph <- c("MMETSP0947", "MMETSP0948", "MMETSP0949", "MMETSP0950", "MMETSP0292", "MMETSP0294", "MMETSP0295", "MMETSP0296", "MMETSP0409", "MMETSP0410", "MMETSP0411", "MMETSP0894", "MMETSP0895", "MMETSP0896", "MMETSP0897", "MMETSP0414", "MMETSP0415", "MMETSP0416")
m <- m[,which(colnames(m) %in%  raph)]

#trim orthogroups found in <50% of samples
num.na =apply(m,1,function(x){length(which(is.na(x)))}) 
cut = which(num.na > ncol(m)/2)
m = m[-cut,]
dim(m)

#write new file
write.table(m,paste0(dir,"merged_counts_Raph.txt"))

#Cut table down to shared "core" orthologous groups
cut <- which(is.na(apply(m,1,mean)))
m.core = m[-cut,]
dim(m.core)
write.table(m.core,paste0(dir,"merged_counts_core_Raph.txt"))


