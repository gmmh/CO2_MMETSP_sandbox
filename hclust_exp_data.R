# hclust_exp_data.R
# Hierarchical Clustering of Gene Expression Data
# Inputs: file with normalized log2 ratios (rows= gene IDs, columns = conditions), number of clusters (k), filename for output
# Output: cluster assignments

library(fastcluster) #this just speeds up hclust function by implementing in c++
library(hyperSpec) #need this for pearson distance calculation function
library(impute) #impute missing values


#load data for clustering
data.dir <- "~/Desktop/CO2_MMETSP/Haptophyta/counts/"
data <- read.table(paste0(data.dir, "norm_RNAseq_log2ratios.tsv"), header=T)
#control = "MMETSP0088"
#co2 = "MMETSP0092"

control = "MMETSP0143"
co2 = "MMETSP0147"

#Subset data
	#If desired can cut experiments with fewer genes (like TM.microarrays)
 data.sub <- data
	#keep only genes with no NA values in any experiments
# na.genes <- which(is.na(apply(data.sub, 1,mean))) #should find a gene if NA in any condition
# data.sub <- data.sub[-na.genes,]

#Impute missing data
	#impute.knn (k=10), uses 10 nearest neighbors determined by hclust to infer what the expression should be of missing genes
	#data.imp <- impute.knn(as.matrix(data))$data 


#Calculate distances
distances <- pearson.dist(data.sub) #pearson correlation distances
#distances <- pearson.dist(data.imp) 

#Perform hierarchical clustering 
hclusters <- hclust(distances, method="ward.D") #other method choices are possible

#Cut the tree to yield cluster assignments
k = 100 #arbitrary number of clusters
clusters <- cutree(hclusters, k)

#make heatmap plot
pdf(file= paste0(data.dir,"heatmap.pdf"))
heatmap(as.matrix(data.sub), Rowv=as.dendrogram(hclusters))
dev.off()

#Visualizing select clusters across all conditions
#clust.id <- clusters["OG0000253"] #this is delta CA
#clust.memb <- names(which(clusters == clust.id))
clust.id <- 100 #clusters["OG0001486"]
clust.memb <- names(which(clusters == clust.id))
pdf(file = paste0(data.dir,"cluster100_log2FC.pdf"))
par(mfrow=c(2,1))
boxplot(data.sub, ylab="Log2FC", xlab="Sample", cex.axis=0.5,outline=F)
abline(v=which(colnames(data.sub)=="MMETSP0147")+0.5, lty=3, col="red")
abline(v=which(colnames(data.sub)=="MMETSP0147")-0.5, lty=3, col="red")
abline(v=which(colnames(data.sub)=="MMETSP1364")+0.5, lty=3, col="red")
abline(v=which(colnames(data.sub)=="MMETSP1364")-0.5, lty=3, col="red")
for(i in 1:length(clust.memb)){
	lines(t(data.sub[clust.memb[i],]), col="red")
}
dev.off()

# clust.id <- clusters["OG0005831"] #this is PGP
# clust.memb <- names(which(clusters == clust.id))
# boxplot(data.sub)
# for(i in 1:length(clust.memb)){
	# lines(t(data.sub[clust.memb[i],]), col="red")
# }

#Write file with cluster membership and date stamp
#write.table(clusters, file=paste0("~/Desktop/CO2_MMETSP/Thalassiosira/all_clusters_", Sys.Date(),".txt"))

#Write file with new formatted cluster membership, residuals and date stamp
out="clust.id: residual, mean.co2, sd.co2, ortho.ids"
for(n in 1:k){
	clust.memb <- names(which(clusters == n))
	clust.mat=t(data.sub[clust.memb,])
	mean=apply(clust.mat,1,mean)
	mean.co2 <- mean[co2]
	sd=apply(clust.mat,1,sd)
	sd.co2 = sd[co2]
	#boxplot(data.sub)
	#lines(mean,col="red")
	#lines(mean+sd,col="orange")
	#lines(mean-sd,col="orange")
	Rss=sum(apply(sweep(clust.mat,1,mean),2,function(x){sum(x^2)}))
	Rtot=sum(apply(sweep(clust.mat,1,mean(mean)),2,function(x){sum(x^2)}))
	residual= Rss/Rtot
	out[n+1]=paste(c(paste0(n, ":"), paste0(round(residual, digits=3), ", ", round(mean.co2, digits=3),", ", round(sd.co2, digits =3), ", "), clust.memb), collapse=" ")
}

write.table(out, file=paste0(data.dir,"core_clusters", Sys.Date(),".txt"), row.names=F, col.names=F, quote=F)
