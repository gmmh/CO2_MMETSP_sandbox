#Calc_FC_DiffCO2.R
#Load libraries
library(NOISeq) # simple TMM normalization
library(plyr)
library(gplots)

#find all filenames from count directory with "orthoID"
# dir="~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/"
# files=list.files(dir, "orthoID")
# length(files)

calc_FC_CO2 <- function(file, control, co2){	
	#load the counts_orthoID.txt file
	t= read.table(file, header=T) #open the file and read in data
	t2=t[-which(is.na(t$Orthogroup)),] #get rid of genes that were not placed in Orthogroup
	t.final=ddply(t2,"Orthogroup",numcolwise(sum)) #sum counts for contigs in same orthogroup
		
	#make orthogroups into row names and delete Orthogroup column
	rownames(t.final)<-t.final$Orthogroup
	t.final <- t.final[,-which(colnames(t.final)=="Orthogroup")]
	
	##Exploratory analyses
	#boxplot(t.final, ylim= c(0,15000))
	#plot(log(t.final[,co2],base=2), log(t.final[,co2]/t.final[,control], base=2))
	
	#TMM normalize and recheck boxplot distributution
	TMM = tmm(t.final)
	#boxplot(TMM, ylim=c(0,10000))
	#plot(log(TMM[,co2],base=2), log(TMM[,co2]/TMM[,control], base=2))
	
	## Calculate log2 ratios with Control
	control.vals <- TMM[,control]
	norm.data <- log(sweep(TMM, 1, control.vals, "/"), 2) #Divide by control.mean and log2 transform
	norm.data <- norm.data[ ,-which(colnames(TMM) ==control)] # Remove control samples
	#boxplot(norm.data)
	
	bfc <- norm.data[which(abs(norm.data[,co2]) > 1), co2]
	print(paste(length(bfc), "= genes with >2FC under CO2"))
	norm.data
}

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Alexandrium-monilatum-CCMP3105_cds_counts_orthoID.txt"
control = "MMETSP0093"
co2 = "MMETSP0097"
norm.data.Amo = calc_FC_CO2(file, control, co2)
up.Amo <- norm.data.Amo[which(norm.data.Amo[,co2] > 1), co2]
down.Amo <- norm.data.Amo[which(norm.data.Amo[,co2] < -1), co2]

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Chaetoceros-affinis-CCMP159_cds_counts_orthoID.txt"
control = "MMETSP0088"
co2 = "MMETSP0092"
norm.data.Caf = calc_FC_CO2(file, control, co2)
up.Caf <- norm.data.Caf[which(norm.data.Caf[,co2] > 1), co2]
down.Caf <- norm.data.Caf[which(norm.data.Caf[,co2] < -1), co2]

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Chrysochromulina-polylepis-CCMP1757_cds_counts_orthoID.txt"
control = "MMETSP0143"
co2 = "MMETSP0147"
norm.data.Cpo = calc_FC_CO2(file, control, co2)
up.Cpo <- norm.data.Cpo[which(norm.data.Cpo[,co2] > 1), co2]
down.Cpo <- norm.data.Cpo[which(norm.data.Cpo[,co2] < -1), co2]

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Gephyrocapsa-oceanica-RCC1303_cds_counts_orthoID.txt"
control = "MMETSP1363"
co2 = "MMETSP1364"
norm.data.Goc = calc_FC_CO2(file, control, co2)
up.Goc <- norm.data.Goc[which(norm.data.Goc[,co2] > 1), co2]
down.Goc <- norm.data.Goc[which(norm.data.Goc[,co2] < -1), co2]

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Heterosigma-akashiwo-CCMP2393_cds_counts_orthoID.txt"
control = "MMETSP0292"
co2 = "MMETSP0296"
norm.data.Hak = calc_FC_CO2(file, control, co2)
up.Hak <- norm.data.Hak[which(norm.data.Hak[,co2] > 1), co2]
down.Hak <- norm.data.Hak[which(norm.data.Hak[,co2] < -1), co2]

file = "~/Desktop/CO2_MMETSP/DiffExp_CO2/counts/Prorocentrum-minimum-CCMP1329_cds_counts_orthoID.txt"
control = "MMETSP0053"
co2 = "MMETSP0057"
norm.data.Pmi = calc_FC_CO2(file, control, co2)
up.Pmi <- norm.data.Pmi[which(norm.data.Pmi[,co2] > 1), co2]
down.Pmi <- norm.data.Pmi[which(norm.data.Pmi[,co2] < -1), co2]

#venn(list(rownames(norm.data.Cpo), rownames(norm.data.Goc), rownames(norm.data.Amo), rownames(norm.data.Pmi), rownames(norm.data.Caf)))

par(mfrow=c(1,2))
# Venn of up regulated genes:
#venn(list(names(up.Cpo), names(up.Goc), names(up.Amo), names(up.Pmi), names(up.Caf)))
Cpo = names(up.Cpo)
Goc = names(up.Goc)
Pmi = names(up.Pmi)
Caf = names(up.Caf)
Hak = names(up.Hak)
data.up =list(Cpo=Cpo, Goc=Goc, Pmi=Pmi, Caf=Caf, Hak=Hak)
venn(data.up)
#Venn of down regulated genes:
Cpo = names(down.Cpo)
Goc = names(down.Goc)
Pmi = names(down.Pmi)
Caf = names(down.Caf)
Hak = names(down.Hak)
data.down =list(Cpo=Cpo, Goc=Goc, Pmi=Pmi, Caf=Caf, Hak=Hak)
venn(data.down)


venn(list(names(down.Cpo), names(down.Goc), names(down.Amo), names(down.Pmi), names(down.Caf)))
venn(list(names(down.Cpo), names(down.Goc), names(down.Pmi), names(down.Caf), names(down.Hak)))


names(up.Cpo)[which(names(up.Cpo) %in% names(up.Goc))]
names(down.Cpo)[which(names(down.Cpo) %in% names(down.Goc))]
data.frame()

names(up.Amo)[which(names(up.Amo) %in% names(up.Pmi))]
names(down.Amo)[which(names(down.Amo) %in% names(down.Pmi))]