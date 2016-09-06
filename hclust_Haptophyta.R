#hclust_Haptophyta

#Load libraries
library(NOISeq) # simple TMM normalization
library(fastcluster) #this just speeds up hclust function by implementing in c++
library(hyperSpec) #need this for pearson distance calculation function
library(impute) #impute missing values


#Load merged RNASEq Data
dir <- "~/Desktop/CO2_MMETSP/Haptophyta/counts/"
m <- read.table(paste0(dir,"merged_counts_all.txt"), header=T)
control.name =  "MMETSP0143" #Chrysochromulina polylepsis exponential
co2.cpo =  "MMETSP0147" #Chrysochromulina polylepsis high CO2
control.goc =  "MMETSP1363" #Gephyrocapsa oceanica exponential
co2.goc = "MMETSP1364" #Gephyrocapsa oceanica CO2

#Impute missing data
	#impute.knn (k=10), uses 10 nearest neighbors determined by hclust to infer what the expression should be of missing genes
	data.imp <- impute.knn(as.matrix(m))$data 


#Replace NA with 0
for(col in colnames(m)){
	replace <- which(is.na(m[,col]))
	m[replace,col] <- 0
}


##Exploratory analyses, compare NA replaced with 0 to NA imputed
par(mfrow=c(2,1))
boxplot(m, outline =F)
boxplot(data.imp, outline=F)

## Normalization Step, compare NA replaced with 0 to NA imputed
TMM <- tmm(m, long =1000, lc=0, k=0.4) #TMM norm
TMM.imp <- tmm(data.imp, long =1000, lc=0, k=0.4)
	#Visualize
par(mfrow=c(2,1))
boxplot(TMM, outline=F)
boxplot(TMM.imp, outline=F)

## Calculate log2 ratios with Control
control.mean <- TMM.imp[,grep(control.name, colnames(TMM))]
#control.sd <- apply(TMM[,grep("seastar", colnames(TMM))], 1, sd)
norm.data <- log(sweep(TMM.imp, 1, control.mean, "/"), 2) #Divide by control.mean and log2 transform
norm.data <- norm.data[ ,-grep(control.name, colnames(norm.data))] # Remove control samples
boxplot(norm.data, outline=F)

#Calculate distances
distances <- pearson.dist(norm.data) #pearson correlation distances
#distances <- pearson.dist(data.imp) 

#Perform hierarchical clustering 
hclusters <- hclust(distances, method="ward.D") #other method choices are possible

#Cut the tree to yield cluster assignments
k = 500 #arbitrary number of clusters
clusters <- cutree(hclusters, k)

#load annotations:
annot= read.table(file="~/Desktop/CO2_MMETSP/Haptophyta/annot/Goc_ortho_pfam.txt", sep=",", header=T)

#Write file with new formatted cluster membership, residuals and date stamp
out="clust.id: residual, mean.co2.Cpo, se.co2.Cpo, mean.co2.Goc, se.co2.Goc , ortho.ids"
#out.CO2 = "clust.id: residual, mean.co2, sd.co2, ortho.ids"
for(n in 1:k){
	clust.memb <- names(which(clusters == n))
	clust.mat=t(norm.data[clust.memb,])
	mean=apply(clust.mat,1,mean)
	mean.co2.cpo <- mean[co2.cpo]
	mean.co2.goc <- mean[co2.goc]/mean[control.goc]
	sd=apply(clust.mat,1,sd)/sqrt(length(clust.memb))
	sd.co2.cpo = sd[co2.cpo]
	sd.co2.goc = sd[co2.goc]/mean[control.goc]
	#boxplot(data.sub)
	#lines(mean,col="red")
	#lines(mean+sd,col="orange")
	#lines(mean-sd,col="orange")
	Rss=sum(apply(sweep(clust.mat,1,mean),2,function(x){sum(x^2)}))
	Rtot=sum(apply(sweep(clust.mat,1,mean(mean)),2,function(x){sum(x^2)}))
	residual= Rss/Rtot
	out[n+1]=paste(c(paste0(n, ":"), paste0(round(residual, digits=3), ", ", round(mean.co2.cpo, digits=3),", ", round(sd.co2.cpo, digits =3), ", ", round(mean.co2.goc, digits=3), ", ", round(sd.co2.goc, digits=3), ", "), clust.memb), collapse=" ")
	if(abs(mean.co2.cpo) >1){ print(paste("Diff in Cpo",n))}
	if(abs(mean.co2.goc) >1){print(paste("Diff in Goc", n))}
}

control.goc =  "MMETSP1363" #Gephyrocapsa oceanica exponential
co2.goc = "MMETSP1364" #Gephyrocapsa oceanica CO2


