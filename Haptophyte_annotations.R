#Haptophyte_annotations

#fix Chrysochomulina annotations
annot = read.table("~/Desktop/CO2_MMETSP/DiffExp_CO2/annotations/Chrysochromulina-polylepis-CCMP1757_pfams.txt", sep ="\t", header =T)
key = read.table("~/Desktop/CO2_MMETSP/DiffExp_CO2/annotations/Cpo_conversion.tab", sep ="\t", header =T)
rownames(annot) = annot$Gene_ID
rownames(key) = key$New.Seq.Name
out = merge(annot, key, by="row.names")
out2 = out[,c("Old.Seq.Name","PFAM")]
Contig = tapply(out2$Old.Seq.Name, 1:length(out2$Old.Seq.Name), function(x){paste0(x,"_1")})
out3 <-cbind(Contig,out2)
out3 <- out3[,-2]
write.table(out3, file="~/Desktop/CO2_MMETSP/DiffExp_CO2/annotations/Chrysochromulina-polylepis-CCMP1757_fixed_pfams.txt", sep="\t", row.names =F, quote=F)
write.table(out3, file="~/Desktop/CO2_MMETSP/Haptophyta/annot/Chrysochromulina-polylepis-CCMP1757_pfams.txt", sep="\t", row.names =F, quote=F)
cpo.counts <- read.table("~/Desktop/CO2_MMETSP/Haptophyta/counts/Chrysochromulina-polylepis-CCMP1757_cds_counts_orthoID.txt", header=T)
cpo.counts = cpo.counts[-which(is.na(cpo.counts$Orthogroup)),]
cpo.annotations =merge(cpo.counts, out3, by="Contig")
cpo.annot.ortho = cpo.annotations[,c("Contig","Orthogroup", "PFAM")]

#load annotations from Gephyrocapsa
counts = read.table("~/Desktop/CO2_MMETSP/Haptophyta/counts/Gephyrocapsa-oceanica-RCC1303_cds_counts_orthoID.txt", header =T)
counts = counts[-which(is.na(counts$Orthogroup)),]
annot = read.table("~/Desktop/CO2_MMETSP/Haptophyta/annot/Gephyrocapsa-oceanica-RCC1303_pfam.txt", sep="\t")
annot = annot[,c(1,9)]
colnames(annot)= c("Contig", "pfam")
annotations = merge(counts, annot, by="Contig")
annotations = annotations[,c("Contig", "Orthogroup", "pfam")]

ogid <- unique(annotations[,"Orthogroup"])
keep <- as.data.frame(matrix(nrow=0, ncol=length(colnames(annotations))))
for(id in ogid){
	location <- grep(id, annotations[,"Orthogroup"])[1]
	keep <- rbind(keep,annotations[location,])
}
keep <- keep[order(keep[,"Orthogroup"]),]
write.table(keep, file="~/Desktop/CO2_MMETSP/Haptophyta/annot/Goc_ortho_pfam.txt", sep=",", quote =F, row.names=F)

consensus <- merge(cpo.annot.ortho, annotations, by="Orthogroup")
