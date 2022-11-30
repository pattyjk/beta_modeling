```
setwd("./Github/")



#need below packages, can skip if installed
install.packages('calibrate')
install.packages('gplots')
install.packages('indicspecies')
install.packages('limma 3.26.9')
install.packages('MASS') 
install.packages('outliers')
install.packages('reshape2')
install.packages('reldist')
install.packages('bipartite')
install.packages('GUniFrac')
install.packages('ape')
install.packages('phangorn')
install.packages('vegan')
install.packages("ggplot2")

#load packages
library(ggplot2)
library(vegan)

#read in meta data, ASV table, taxonomic assignments
meta<-read.delim("becker_map.txt", header=T)
comm<-read.delim("becker_table.txt", header=T, row.names=1)

#transpose matrix
comm.t=t(comm)

#rarefy data to same depth
min(rowSums(comm.t))
#19619
com.rare<-as.data.frame(rrarefy(comm.t, sample=19000))



#null model

#Source for model fits is from Burns et al. ISMEJ 2015, downloaded R code from their supporting materials
#Source code requires:  minpack.lm, Hmisc, stats4 packages - make sure they are installed (and their dependencies)
source("sncm.fit_function.r")

#fit null model to full dataset
total.model<-sncm.fit(com.rare, pool = NULL, stats=T, taxon=NULL)


#subset data to only inculde infected/non infected
infected.rows<-meta[which(meta$Infection=="infected"),]
com.infect<-com.rare[infected.rows$SampleID,]
infect.model<-sncm.fit(com.infect, pool = NULL, stats=T, taxon=NULL)

uninfected.rows<-meta[which(meta$Infection=="uninfected"),]
com.uninfect<-com.rare[uninfected.rows$SampleID,]
uninfect.model<-sncm.fit(com.uninfect, pool = NULL, stats=T, taxon=NULL)




#beta null models
source("MetacommunityDynamicsFctsOikos.R")
source("PANullDevFctsOikos.R")


### Prepare and calculate abundance beta-null deviation metric
## Adjusted from Stegen et al 2012 GEB

#load packages needed
library(GUniFrac)
library(ape)
library(phangorn)

#read in rooted tree
tree<-read.tree("tree.nwk")

bbs.sp.site <- comm.t
patches=nrow(bbs.sp.site)
rand <- 999

null.alphas <- matrix(NA, ncol(comm.t), rand)
null.alpha <- matrix(NA, ncol(comm.t), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)
bucket_wuf_res <- matrix(NA, patches, rand) #als add

bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
gamma <- ncol(bbs.sp.site) #gamma
obs_beta <- 1-mean.alpha/gamma
obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma

##Generate null patches, shit takes a while on a laptop :()
for (randomize in 1:rand) {  
  null.dist = comm.t
  for (species in 1:ncol(null.dist)) {
    tot.abund = sum(null.dist[,species])
    null.dist[,species] = 0
    for (individual in 1:tot.abund) {
      sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
      null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
    }
  }
  
  ##Calculate null deviation for null patches and store
  null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
  null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
  expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
  null.alpha <- mean(null.alphas[,randomize])
  null.alpha.comp <- c(null.alpha.comp, null.alpha)
  
  bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
  wuf<-(GUniFrac(null.dist, tree, alpha=1)) #als add
  #wuf<-(GUniFrac(comm.t, tree, alpha=1)) #als add test that comparable  values are calculated as with QIIME
  bucket_wuf <- as.matrix(wuf$unifracs[,,"d_1"]) #als add
  diag(bucket_bray) <- NA
  diag(bucket_wuf) <- NA #als add
  bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
  bucket_wuf_res[,randomize] <- apply(bucket_wuf, 2, FUN="mean", na.rm=TRUE) #als add
} ## end randomize loop

## Calculate beta-diversity for obs metacommunity
beta_comm_abund <- vegdist(comm.t, "bray")
wuf_comm_abund <- GUniFrac(comm.t, tree, alpha=1) #als add
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
res_wuf_comm_abund <- as.matrix(as.dist(wuf_comm_abund$unifracs[,,"d_1"])) #als add
diag(res_beta_comm_abund) <- NA
diag(res_wuf_comm_abund) <- NA #als add

# output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
wuf_div_abund_stoch <- apply(res_wuf_comm_abund, 2, FUN="mean", na.rm=TRUE) #als add

# output abundance beta-null deviation
bray_abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)
wuf_abund_null_dev <- wuf_div_abund_stoch - mean(bucket_wuf_res) #als add

### Outputs:
#beta_div_stoch  - Jaccard beta-diversity for the metacommunity, average value (of all pairwise comparisons) for each patch
#beta_div_abund_stoch - Bray-Curtis beta-diversity for the metacommunity, average value (of all pairwise comparisons) for each patch
#PA_null_dev - presence-absence null deviation values or the metacommunity, average value (of all pairwise comparisons) for each patch
#abund_null_dev - abundance null deviation values or the metacommunity, average value (of all pairwise comparisons) for each patch
###
#END script by Tucker et al.


####
#examining change in the taxa between the infected/uninfected
library(ggplot2)
library(plyr)
library(reshape2)
library(vegan)

#read in data
meta<-read.delim("becker_map.txt", header=T)
comm<-read.delim("becker_table.txt", header=T, row.names=1)

#install.packages("BiodiversityR")
library(BiodiversityR)

#transpose matrix
comm.t=t(comm)

#rarefy data to same depth
min(rowSums(comm.t))
#19619
com.rare<-as.data.frame(rrarefy(comm.t, sample=19000))

#calculate for total community
total_com_rank<-rankabundance(com.rare)

#calculate for infected/not infected
infected.rows<-meta[which(meta$Infection=="infected"),]
com.infect<-com.rare[infected.rows$SampleID,]


uninfected.rows<-meta[which(meta$Infection=="uninfected"),]
com.uninfect<-com.rare[uninfected.rows$SampleID,]

uninfected_rank<-rankabundance(com.uninfect)
infected_rank<-rankabundance(com.infect)

#add ASV names as a column
uninfected_rank$ASV<-row.names(uninfected_rank)
infected_rank$ASV<-row.names(infected_rank)

#rename a few columns for merging
names(infected_rank)<-c("Infected_rank",  "Infected_abun",  "Infected_plower",  "Infected_pupper",  "Infected_accumfreq",  "Infected_logabun",  "Infected_rankfreq")


```
