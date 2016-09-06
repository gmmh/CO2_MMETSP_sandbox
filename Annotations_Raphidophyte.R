#Annotations_Raphidophyte


#load annotations from Chaetoceros affinis
counts = read.table("~/Desktop/CO2_MMETSP/Ochrophyta/counts/Heterosigma-akashiwo-CCMP2393_cds_counts_orthoID.txt", header =T)
counts = counts[-which(is.na(counts$Orthogroup)),]
annot = read.table("~/Desktop/CO2_MMETSP/Ochrophyta/annot/Heterosigma-akashiwo-CCMP2393_pfam.txt", sep="\t")
annot = annot[,c(1,9)]
colnames(annot)= c("Contig", "pfam")

#merge OG IDs and annotations
annotations = merge(counts, annot, by="Contig")

ogid <- unique(annotations[,"Orthogroup"])
keep <- as.data.frame(matrix(nrow=0, ncol=length(colnames(annotations))))
for(id in ogid){
	location <- grep(id, annotations[,"Orthogroup"])[1]
	keep <- rbind(keep,annotations[location,])
}

keep <- keep[order(keep[,"Orthogroup"]),c("Contig", "Orthogroup", "pfam")]
write.table(keep, file="~/Desktop/CO2_MMETSP/Ochrophyta/annot/Hak_ortho_pfam.txt", sep=",", quote =F, row.names=F)

