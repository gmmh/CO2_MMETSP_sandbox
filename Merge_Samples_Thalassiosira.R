#Merge_Samples_Thalassiosira.R
#Trim and merge files to make orthogroup x sample, log2FC tables
library(plyr)

#find all filenames from count directory with "orthoID"
dir="~/Desktop/CO2_MMETSP/Thalassiosira/counts/"
files=list.files(dir, "orthoID")

#loop to open files, discard unassigned contigs, sum reads for othrogroups and merge
m=NA
for(file in files){
	t= read.table(paste0(dir,file), header=T)
	t2=t[-which(is.na(t$Orthogroup)),]
	t.final=ddply(t2,"Orthogroup",numcolwise(sum))
	if(is.na(m)){m=t.final}else{m= merge(t.final,m, all=T)}
}

#Get Reference counts with orthoIDs from count/Ref/ directory
ref.dir=paste0(dir,"Ref/")
ref.files=list.files(ref.dir, "orthoID")
for(file in ref.files){
	t= read.table(paste0(ref.dir,file), header=T)
	t2=t[-which(is.na(t$Orthogroup)),]
	t.final=ddply(t2,"Orthogroup",numcolwise(sum))
	m= merge(t.final,m, all=T)
}

#cut any remaining rows with NA in ref sample
m<-m[-which(is.na(m$A314.antisense.seastar.tab)),] 

#make orthogroups into row names and delete Orthogroup column
rownames(m)<-m$Orthogroup
m <- m[,-which(colnames(m)=="Orthogroup")]

#write new file
write.table(m,paste0(dir,"merged_counts.txt"))
