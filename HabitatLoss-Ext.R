install.packages(c("vegan", "dplyr", "ade4", "funrar", "PMCMR")) #Install required packages
#Load required packages
library(vegan)
library(dplyr)
library(ade4)
library(funrar)
library(PMCMR)

#### INDEXES CALCULATION #### ----

abun<-read.table('clipboard', h=T,dec = '.'); abun #Load abundance matrix, where samples are rows and species are collumns
spp<-decostand(abun, "pa"); spp #Transform the abundance matrix in a presence-absence matrix. Needed for regional taxonomic index
traits<-read.table('clipboard', h=T); traits #Load trait matrix, where species are in the rows and traits in collumns
env<-read.table('clipboard', h=T,dec = '.'); env #Load minimal habitat integrity needed for each species to occur, where species are in the rows and the integrity in collumns

decostand(env, "standardize")->env.pad; env.pad

Maxlength<-traits[,1]
tabDiet<-prep.fuzzy(traits[,2:8], col.blocks=7) #Identify that the trait is a fuzzy type. If your traits are of other type, check what function to use in its preparation using help(dist.ktab)
VertStrat<-prep.binary(traits[,9:10], col.blocks = 2)
Protect<-prep.binary(traits[,11:18], col.blocks = 8)
Active<-prep.circular(traits[, 19:21])
kTab<-ktab.list.df(list(as.data.frame(Maxlength), tabDiet, VertStrat, Protect, Active)) #Prepares the fuzzy matrix to the calculation of Gower's distance 
dist.fun<-dist.ktab(kTab, type=c("Q", "F", "B", "B", "C"), option=c("scaledBYrange")) #Gower's distance

abun.rel<-make_relative(as.matrix(abun)) #Distinctiveness index is suggested to be calculated with relative abundance of species by the authors that described the index (see Grenié, Denelle, Tucker, Munoz, & Violle, 2017 for more) 
abun.rel<-as.matrix(abun.rel) #Turns relative abundance in a quadratic matrix object
dist.fun<-as.matrix(dist.fun) #Turns Gower's distance in a quadratic matrix object
rownames(dist.fun)=colnames(abun)
colnames(dist.fun)=colnames(abun)
spp<-as.matrix(spp) #Turns species occurrences in a quadratic matrix object
nSamp<-nrow(abun) #Establish the number of Samples in your dataset.
nSpp<-ncol(abun) #Establish the number of species in your dataset

loc.di = distinctiveness(abun.rel, dist.fun); loc.di #Calculates local distinctiveness index
write.table(loc.di, file="dist.txt", sep="\t") #Exports local distinctiveness

loc.scar<-scarcity(abun.rel); loc.scar #Calculates local scarcity index
write.table(loc.scar, file="scar.txt", sep="\t")  #Exports local scarcity index

reg.uni<-uniqueness(spp, dist.fun); reg.uni  #Calculates regional uniqueness index
write.table(reg.uni, file="uni.txt", sep="\t") #Exports regional uniqueness index

reg.res<-restrictedness(spp); reg.res #Calculates regional restrictedness index
write.table(reg.res, file="restr", sep="\t") #Exports regional restrictedness index

pca1 <- dudi.pca(env.pad, scan = FALSE)
nic1 <- niche(pca1, as.data.frame(abun.rel), scan = FALSE)
amp.hab<-niche.param(nic1)
as.data.frame(amp.hab)->amp.hab
as.data.frame(amp.hab$Tol)->tol; tol
rownames(tol)<-rownames(traits); tol
as.numeric(as.matrix(tol/max(tol)))->tol
(tol-1)*-1->vul; vul
write.table(vul, file="vul.txt", sep="\t") #Exports regional restrictedness index

t(loc.scar)->loc.scar.t
matrix(NA, nrow = nSpp, ncol = nSamp)->newmatrix
FinalRI=lapply(seq(1, nSamp, by = 1), function(x) {
  loc.scar.t[,x]+vul->sum
  sum/2
})
unlist(FinalRI)->newmatrix[]
rownames(newmatrix)<-rownames(loc.scar.t);newmatrix
t(newmatrix)->ExtRiskLoc; ExtRiskLoc
write.table(ExtRiskLoc, file="FinalRI_loc.txt", sep="\t")  #Exports local scarcity index

as.data.frame((vul+reg.res[,2])/2)->ExtRiskReg
write.table(ExtRiskReg, file="FinalRI_reg.txt", sep="\t")  #Exports local scarcity index
row.names(ExtRiskReg)=colnames(abun)

### HOW THE UNIQUENESS AND DISTINCTIVENESS OF THE COMMUNITY CHANGES AS SPECIES GOES EXTINCT? ###
#############################
#REGIONAL ----
Uni<-as.data.frame(reg.uni[,2]) #Turns the uniqueness matrix in a data frame object
t(spp)->t_spp
as.data.frame(rowSums(t_spp))->oc_spp
data<-cbind(Uni, ExtRiskReg, oc_spp)
colnames(data)=c("Uniqueness", "Extinction_Risk", "Occurrence"); data

unimatrix=lapply(seq(1, nrow(data), by = 1), function(u) {
  data[u,3]->ab
  esp=print(Uni[rep(u, ab), 1],)
  esp
})
restmatrix=lapply(seq(1, nrow(data), by = 1), function(u) {
  data[u,3]->ab
  esp=print(ExtRiskReg[rep(u, ab), 1],)
  esp
})
unlist(unimatrix)->unl.unimtx
unlist(restmatrix)->unl.restmtx
as.data.frame(unlist(restmatrix))->indlist
fac<-nrow(indlist)
ind_matrix=matrix(NA, nrow = fac, ncol = 2)
unl.unimtx->ind_matrix[,1]
unl.restmtx->ind_matrix[,2]; ind_matrix
#first Sensible -----

uni_ext=matrix(NA, nrow = ind_matrix, ncol = 1) #Creates an empty matrix to input our results
data_Rst<-ind_matrix[order(ind_matrix[,2]),]; data_Rst #Orders species from the most restrict to the widespread ones.
uni_rst<-data_Rst[, 1]; uni_rst #Extract only uniqueness from the table
as.data.frame(uni_rst)->uni_rst #Turns uniqueness to a data frame object
output_rst=lapply(seq(1, nrow(ind_matrix), by = 1), function(i) {
  mean(uni_rst[1:i,])
}) #Extincts species from the most restrict to the widespread species and calculates mean uniqueness at each step
unlist(output_rst)->output_rst #Unlist mean uniqueness at each extinction step
as.data.frame(output_rst)->output_rst; output_rst ##Turns mean uniqueness at each extinction step to a data frame object

#first Tolerant -----

data_Wid<-ind_matrix[order(ind_matrix[,2], decreasing=T),]; data_Wid #Orders species for the most widespread to the restrict ones.
uni_wid<-data_Wid[, 1]; uni_wid #Extract only uniqueness from the table
as.data.frame(uni_wid)->uni_wid #Turns uniqueness to a data frame object
output_wid=lapply(seq(1, nrow(ind_matrix), by = 1), function(i) {
  mean(uni_wid[1:i,])
}) #Extincts species from the most widespread to the restrict species and calculates mean uniqueness at each step
unlist(output_wid)->output_wid #Unlist mean uniqueness at each extinction step
as.data.frame(output_wid)->output_wid; output_wid ##Turns mean uniqueness at each extinction step to a data frame object

#Random ----
dataperm=ind_matrix
randSPE=matrix(NA,nrow=nrow(ind_matrix),ncol=1000) #Creates an empty matrix to input our results. ncol=number of randomizations in the random model
perm_mtx=matrix(NA,nrow=nrow(ind_matrix),ncol=1)##Creates an empty matrix to input our results
for	(j	in	1:1000)
{
  perm_mtx=as.data.frame(sample(dataperm[,1]))
  Null_Simulat<-perm_mtx[order(dataperm[,2]),]		
  for	(i	in	1:nrow(ind_matrix))
  {
    randSPE[i,j]<-mean(Null_Simulat)
    Null_Simulat=Null_Simulat[-1]
  }
} #Shuffles species order 1000 times while preserving species uniqueness. Results are inserted in randSPE object
##Calculating	mean	and	quantiles	from	the	null	models. 
Null_Uni=matrix(NA,nrow=nrow(ind_matrix),ncol=3)
for	(i	in	1:nrow(ind_matrix)){
  Null_Uni[i,1]=mean(randSPE[i,])		
  Null_Uni[i,2]=quantile(randSPE[i,],probs=0.025)
  Null_Uni[i,3]=quantile(randSPE[i,],probs=0.975)		
}; Null_Uni #First collumn in mean, second is lower quantile and third is higher quantile

#Plot regional pattern ----
Spp_erosion<-c(1:nrow(ind_matrix))
plot(Spp_erosion,Null_Uni[,1],type="l",xlim=c(0,nrow(ind_matrix)),xlab="Regional	 Species	 loss	 (%)",ylab="Uniqueness",font.lab=2,cex.lab=1.4,cex.axis=1.2,col="gray",pch=16,cex=1.5,xaxt="n
     ",ylim=c(min(Null_Uni),max(Null_Uni)))
axis(side=1,	at=c(0,nrow(ind_matrix)/4,nrow(ind_matrix)/2,nrow(ind_matrix)/1.333,nrow(ind_matrix)),	
     labels=c(0,25,50,75,100),line=F,tick=-0.3,cex.axis=1.2,mgp=c(3,1,0))
polygon(c(1:nrow(ind_matrix),rev(1:nrow(ind_matrix))),c(Null_Uni[,2],rev(Null_Uni[,3])),col="grey88",border=
          F)
points(Null_Uni[,1],col="gray",cex=1.6,type="l",lty=1,lwd=4) #Plots random scenario means
apply(output_rst, 2, rev)->FSpe_Simulat_rare
points(FSpe_Simulat_rare,col="black",type="l",lty=3,lwd=4) #Plots scenario where restrict species are extinct first
apply(output_wid, 2, rev)->FSpe_Simulat_comm
points(FSpe_Simulat_comm,col="black",type="l",lty=1,lwd=4)  #Plots scenario where widespread species are extinct first

#Save Output ----
cbind(Null_Uni, apply(output_wid, 2, rev), apply(output_rst, 2, rev))->reg_uni_pat
colnames(reg_uni_pat)<-c("Random (Mean)", "Random (0.025)", "Random (0.975)", "first Widespread", "first Restrict")
rownames(reg_uni_pat)<-c(1:nrow(ind_matrix)); reg_uni_pat
write.table(reg_uni_pat, file="Uniqueness_regional_extinction.txt", sep="\t")


#LOCAL ----

traits_localfilter=t(spp)
traits_localfilter[traits_localfilter == 0] <- "Absent" #Replace all 0s in the presence-absence matrix to "Absent"
noquote(traits_localfilter)->traits_localfilter
replace(traits_localfilter, traits_localfilter == 1, "Present")->traits_localfilter #Replace all 1s in the presence-absence matrix to "Present"
dist_localfilter=traits_localfilter; dist_localfilter
dist_local<-as.data.frame(t(loc.di)); dist_local #Transposes and converts the matrix of distinctiveness
scar_local<-as.data.frame(t(ExtRiskLoc)); scar_local #Transposes and converts the matrix of scarcity
as.data.frame(t(abun))->abun_t

#first Sensible ----

dist.scar.perc=lapply(seq(1, nSamp, by = 1), function(x) #Extincts species from the most scarce to the abundant ones in each sample. Then, it calculates the mean distinctiveness of the community in each sample in the nine steps (from 10% to 90%, by ten).
{
  dist_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(UAfilter[,x], Scarfilter[,x], Abunfilter[,x])->dist_scar_abun
  distmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 1], ab)
  })
  scarmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 2], ab)
  })
  unlist(distmatrix)->unl.distmtx
  unlist(scarmatrix)->unl.scarmtx
  as.data.frame(unlist(scarmatrix))->indlist
  fac<-nrow(indlist)
  ind_matrix=matrix(NA, nrow = fac, ncol = 2)
  unl.distmtx->ind_matrix[,1]
  unl.scarmtx->ind_matrix[,2]
  UAfilter<-ind_matrix[order(ind_matrix[,2]),]
  nrow(UAfilter)->numsp
  UAfilter[,1]->UAfiltered
  output_scar=lapply(seq(1, numsp, by = 1), function(i) {
    mean(UAfiltered[1:i])
  })
  unlist(output_scar)->output_scar
  ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
    numsp*i/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    output_scar[ua.perc[w]]
  })
  as.numeric(output.ext)->output.ext
})
unlist(dist.scar.perc)->dist.scar.perc #Unlist output
perc.dist.first_scar=matrix(NA, nrow = 10, ncol = nSamp) #Creates a empty matrix to insert our results
dist.scar.perc->perc.dist.first_scar[] #Input our results in the empty matrix
rownames(perc.dist.first_scar)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); perc.dist.first_scar
write.table(perc.dist.first_scar, file="Dist_percentage_first_scar_local.txt", sep="\t") #Saves output from the nine steps method.

#first Tolerant ----

dist.abun.perc=lapply(seq(1, nSamp, by = 1), function(x) #Extincts species from the most scarce to the abundant ones in each sample. Then, it calculates the mean distinctiveness of the community in each sample in the nine steps (from 10% to 90%, by ten).
{
  dist_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(UAfilter[,x], Scarfilter[,x], Abunfilter[,x])->dist_scar_abun
  distmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 1], ab)
  })
  scarmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 2], ab)
  })
  unlist(distmatrix)->unl.distmtx
  unlist(scarmatrix)->unl.scarmtx
  as.data.frame(unlist(scarmatrix))->indlist
  fac<-nrow(indlist)
  ind_matrix=matrix(NA, nrow = fac, ncol = 2)
  unl.distmtx->ind_matrix[,1]
  unl.scarmtx->ind_matrix[,2]
  UAfilter<-ind_matrix[order(ind_matrix[,2], decreasing = T),]
  nrow(UAfilter)->numsp
  UAfilter[,1]->UAfiltered
  output_scar=lapply(seq(1, numsp, by = 1), function(i) {
    mean(UAfiltered[1:i])
  })
  unlist(output_scar)->output_scar
  ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
    numsp*i/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    output_scar[ua.perc[w]]
  })
  as.numeric(output.ext)->output.ext
})

unlist(dist.abun.perc)->dist.abun.perc #Unlist output
perc.dist.first_abun=matrix(NA, nrow = 10, ncol = nSamp) #Creates a empty matrix to insert our results
dist.abun.perc->perc.dist.first_abun[] #Input the results in the empty matrix
rownames(perc.dist.first_abun)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); perc.dist.first_abun
write.table(perc.dist.first_abun, file="Dist_percentage_first_abun_local.txt", sep="\t") #Saves output from the nine steps method.

#Random ----

RandPerm_mtx=matrix(NA,nrow=10,ncol=nSamp) #Creates a empty matrix to input our results
for (k in 1:nSamp) {
  randDIS=matrix(NA,nrow=10,ncol=1000)#ncol=randomizations
  dist_local %>% filter(dist_localfilter[, k] == "Present") -> UAfilter
  scar_local %>% filter(dist_localfilter[, k] == "Present") -> Scarfilter
  abun_t %>% filter(dist_localfilter[,k] == "Present") -> Abunfilter
  cbind(UAfilter[,k], Scarfilter[,k], Abunfilter[,k])->dist_scar_abun
  distmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 1], ab)
  })
  scarmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 2], ab)
  })
  unlist(distmatrix)->unl.distmtx
  unlist(scarmatrix)->unl.scarmtx
  as.data.frame(unlist(scarmatrix))->indlist
  fac<-nrow(indlist)
  ind_matrix=matrix(NA, nrow = fac, ncol = 2)
  unl.distmtx->ind_matrix[,1]
  unl.scarmtx->ind_matrix[,2]
  UAfilter<-ind_matrix
  nrow(UAfilter)->numsp
  as.data.frame(UAfilter[,1])->UAfiltered
  for (j in 1:1000) {
    sample(UAfiltered[,1])->UAfilter_sampled
    as.data.frame(UAfilter_sampled)->Null_Simulat
    means_dist=lapply(seq(1, nrow(Null_Simulat), by = 1), function(y) {
      mean(Null_Simulat[1:y,])
    })
    unlist(means_dist)->output_rand
    ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
      nrow(Null_Simulat)*i/100->perc
      round(perc, digits = 0)->perc
    })
    as.numeric(ua.perc)->ua.perc
    output.ext=lapply(seq(1, 10, by = 1), function(w){
      output_rand[ua.perc[w]]
    })
    as.numeric(output.ext)->output.ext
    randDIS[,j]<-output.ext
    print(paste("Wait for it... We are randomizing sample", k, "at the moment. This is permutation", j, "of 1000."))
  }
  RandPerm_mtx[,k]=rowMeans(randDIS)
} #Shuffles the order of extinction 1000 times for each sample while preserving species distinctiveness values. Nine levels of extraction for mean distiviness values are defined (10%-90%).
rownames(RandPerm_mtx)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); RandPerm_mtx

write.table(RandPerm_mtx, file="Dist_percentage_random_local.txt", sep="\t") #Exports the output of the random extinction scenario.

#Friedman test ----

Friedman.X2=lapply(seq(1, 9, by = 1), function(i)
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  cbind(RandPerm_mtx[i,],perc.dist.first_abun[i,], perc.dist.first_scar[i,])->perc10[]
  friedman.test(perc10)->fried.test
  print(fried.test$statistic)
}) #Extract chi-squared values of Friedman's paired test for the comparisons between the distinctiveness of the three scenarios in all nine levels of extinction

Friedman.p=lapply(seq(1, 9, by = 1), function(i)
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(RandPerm_mtx[i,],perc.dist.first_abun[i,], perc.dist.first_scar[i,])->perc10[]
  friedman.test(perc10)->fried.test
  print(fried.test$p.value)
}) #Extract p values of Friedman's paired test for the comparisons between the distinctiveness of the three scenarios in all nine levels of extinction

Post.Hoc=lapply(seq(1, 9, by = 1), function(i) 
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(RandPerm_mtx[i,],perc.dist.first_abun[i,], perc.dist.first_scar[i,])->perc10[]
  posthoc.friedman.nemenyi.test(as.matrix(perc10))->post.fried
  print(post.fried$p.value)
}) #Returns a post-hoc pairwise test of mean distinctviness in the three scenarios for all nine levels of extinction 

capture.output(c(Friedman.X2, Friedman.p, Post.Hoc), file = "Friedman_Test_Index.txt") #Exports Friedman's test results

#Plot local pattern----
dist.rand=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  mean(RandPerm_mtx[x,], na.rm=T)->dist.rand[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 
dist.scar=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  mean(perc.dist.first_scar[x,], na.rm=T)->dist.scar[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 
dist.abun=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  mean(perc.dist.first_abun[x,], na.rm=T)->dist.abun[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 

#Inverts matrixes so results can be represented in incresing order of extinction (from 0-90%).
rev(dist.abun)->dist.abun
rev(dist.rand)->dist.rand
rev(dist.scar)->dist.scar

#Plotting the graph
erosion<-seq(0,90,10)
plot(erosion,dist.rand,type="l",ylim=c(min(dist.scar),max(dist.abun)),xlim=c(0,90),col="gray",pch=16,cex=1.5,xlab="Local	 Species	 Loss (%)",ylab="Distinctiveness",font.lab=2,cex.lab=1.4,cex.axis=1)
points(erosion[-1],dist.rand[-1],type="p",col="gray",pch=21,bg="gray",cex=1)
points(erosion,dist.scar,type="l",ylim=c(0,1),col="black",lty=2)
points(erosion[-1],dist.scar[-1],type="p",col="black",pch=21,bg="white",cex=1)
points(erosion,dist.abun,type="l",ylim=c(0,1),col="black")
points(erosion[-1],dist.abun[-1],type="p",col="black",pch=21,bg="black",cex=1)
legend(0,0.4060,legend=c("scarce extinct first","random extinction","abundant extinct first"),col=c("black","gray","black"),cex=1,lty=c(3,1,1),lwd=1.3,xjust=-0.01,yjust=0.5,horiz=F,text.font=0.5,x.intersp=0,bty="n")


#########################################################################
### EFFECTS OF EXTINCTION IN THE NUMBER OF FUNCTIONS IN THE COMMUNITY ### ----
#########################################################################

##### Maxlength ----
###############
#REGIONAL -----

BodySize<-Maxlength #Transforms the trait matrix in a binary table, where if a species plays that function it will be assigned with an 1, if it does not plays that function, then a 0 is assigned
data_trait<-cbind(ExtRiskReg, BodySize); data_trait #Combine rariy index and traits to order
as.data.frame(rowSums(t(spp)))->oc_spp
data<-cbind(oc_spp, data_trait); data #Combines uniqueness and restrictedness
traitmatrix=lapply(seq(1, nrow(data), by = 1), function(u) {
  data[u,1]->ab
  esp=print(data_trait[rep(u, ab), 1:ncol(data_trait)],)
  esp
})

do.call(rbind.data.frame, traitmatrix)->ind_matrix

#first Restrict -----

traits_Rst<-ind_matrix[order(ind_matrix[, 1]),]; traits_Rst #Order species from the most restrict to the less (widespread).
write.table(traits_Rst, file="Restrict_classified_Maxlength_BAT_Temp.csv", sep="\t") #Saves species functions ordered.
Length<-traits_Rst[, 2:ncol(ind_matrix)]; Length #Extracts only species' functions as previously ordered, without the restrictiveness index.
as.data.frame(Length)->Length #Turns uniqueness to a data frame object
Length_rst=lapply(seq(1, nrow(ind_matrix), by = 1), function(i) {
  mean(Length[1:i,])
}) #Extincts species from the most restrict to the widespread species and calculates mean uniqueness at each step
unlist(Length_rst)->Length_rst #Unlist mean uniqueness at each extinction step
as.data.frame(Length_rst)->Length_rst; Length_rst ##Turns mean uniqueness at each extinction step to a data frame object

#first Widespread -----

traits_Wid<-ind_matrix[order(ind_matrix[,1], decreasing=T),]; traits_Wid #Order species from the most widespread to the less (restrict).
write.table(traits_Wid, file="widespread_classified_Maxlength_BAT_Temp.csv", sep="\t")#Saves species functions ordered.
Length<-traits_Wid[, 2:ncol(ind_matrix)]; Length #Extracts only species' functions as previously ordered, without the restrictiveness index.
as.data.frame(Length)->Length #Turns uniqueness to a data frame object
Length_wid=lapply(seq(1, nrow(ind_matrix), by = 1), function(i) {
  mean(Length[1:i,])
}) #Extincts species from the most restrict to the widespread species and calculates mean uniqueness at each step
unlist(Length_wid)->Length_wid #Unlist mean uniqueness at each extinction step
as.data.frame(Length_wid)->Length_wid; Length_wid ##Turns mean uniqueness at each extinction step to a data frame object

#Random ----
Lengthperm=ind_matrix
randLth=matrix(NA,nrow=nrow(ind_matrix),ncol=1000) #Creates an empty matrix to input our results. ncol=number of randomizations in the random model
for	(j	in	1:1000) {
  Null_Simulat=sample(Lengthperm[,2])
  for	(i	in	1:nrow(ind_matrix))
  {
    randLth[i,j]<-mean(Null_Simulat)
    Null_Simulat=Null_Simulat[-1]
  }
} #Shuffles species order 1000 times while preserving species uniqueness. Results are inserted in randSPE object
##Calculating	mean	and	quantiles	from	the	null	models. 
Null_Lth=matrix(NA,nrow=nrow(ind_matrix),ncol=3)
for	(i	in	1:nrow(ind_matrix)){
  Null_Lth[i,1]=mean(randLth[i,])		
  Null_Lth[i,2]=quantile(randLth[i,],probs=0.025)
  Null_Lth[i,3]=quantile(randLth[i,],probs=0.975)		
}; Null_Lth #First collumn in mean, second is lower quantile and third is higher quantile


#Plot regional extinction pattern ----
Spp_erosion<-c(nrow(ind_matrix):1)
plot(Spp_erosion,Null_Lth[,1],type="l",xlim=c(0,nrow(ind_matrix)),xlab="Regional	 Species 	loss (%)",ylab="Mean Length (cm)",font.lab=2,cex.lab=1.4,cex.axis=1.2,col="gray",pch=16,cex=1.5,xaxt="n
     ",ylim=c(min(Null_Lth),max(Null_Lth)))
axis(side=1,	at=c(0,nrow(ind_matrix)/4,nrow(ind_matrix)/2,nrow(ind_matrix)/1.333,nrow(ind_matrix)),	
     labels=rev(c(0,25,50,75,100)),line=F,tick=-0.3,cex.axis=1.2,mgp=c(3,1,0))
polygon(c(nrow(ind_matrix):1,rev(nrow(ind_matrix):1)),c(Null_Lth[,2],rev(Null_Lth[,3])),col="grey88",border=
          F)
points(Null_Lth[nrow(Null_Lth):1,1],col="gray",cex=1.6,type="l",lty=1,lwd=4)
points(rev(Length_rst),col="black",type="l",lty=3,lwd=4)
points(rev(Length_wid),col="black",type="l",lty=1,lwd=4)

#Save Output ----

cbind(Null_Lth, Length_wid, Length_rst)->ext.regional
colnames(ext.regional)<-c("Random (Mean)", "Random (0.025)", "Random (0.975)", "first Widespread", "first Restrict"); ext.regional
write.table(ext.regional, file="Regional_ext_MaxLength.txt", sep="\t")


#LOCAL ----

traits_localfilter=t(spp)
traits_localfilter[traits_localfilter == 0] <- "Absent" #Replace all 0s in the presence-absence matrix to "Absent"
noquote(traits_localfilter)->traits_localfilter
replace(traits_localfilter, traits_localfilter == 1, "Present")->traits_localfilter #Replace all 1s in the presence-absence matrix to "Present"
dist_localfilter=traits_localfilter; dist_localfilter
dist_local<-as.data.frame(traits[,1]); dist_local #Transposes and converts the matrix of distinctiveness
scar_local<-as.data.frame(t(ExtRiskLoc)); scar_local #Transposes and converts the matrix of scarcity
as.data.frame(t(abun))->abun_t

#first Sensible ----

length.scar.perc=lapply(seq(1, nSamp, by = 1), function(x) #Extincts species from the most scarce to the abundant ones in each sample. Then, it calculates the mean distinctiveness of the community in each sample in the nine steps (from 10% to 90%, by ten).
{
  dist_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(UAfilter, Scarfilter[,x], Abunfilter[,x])->dist_scar_abun
  distmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 1], ab)
  })
  scarmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 2], ab)
  })
  unlist(distmatrix)->unl.distmtx
  unlist(scarmatrix)->unl.scarmtx
  as.data.frame(unlist(scarmatrix))->indlist
  fac<-nrow(indlist)
  ind_matrix=matrix(NA, nrow = fac, ncol = 2)
  unl.distmtx->ind_matrix[,1]
  unl.scarmtx->ind_matrix[,2]
  UAfilter<-ind_matrix[order(ind_matrix[,2]),]
  nrow(UAfilter)->numsp
  UAfilter[,1]->UAfiltered
  output_scar=lapply(seq(1, numsp, by = 1), function(i) {
    mean(UAfiltered[1:i])
  })
  unlist(output_scar)->output_scar
  ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
    numsp*i/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    output_scar[ua.perc[w]]
  })
  as.numeric(output.ext)->output.ext
})
unlist(length.scar.perc)->length.scar.perc #Unlist output
perc.dist.first_scar=matrix(NA, nrow = 10, ncol = nSamp) #Creates a empty matrix to insert our results
length.scar.perc->perc.dist.first_scar[] #Input our results in the empty matrix
rownames(perc.dist.first_scar)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); perc.dist.first_scar
write.table(perc.dist.first_scar, file="Length_percentage_first_scar_local.txt", sep="\t") #Saves output from the nine steps method.

#first Tolerant ----

length.abun.perc=lapply(seq(1, nSamp, by = 1), function(x) #Extincts species from the most scarce to the abundant ones in each sample. Then, it calculates the mean distinctiveness of the community in each sample in the nine steps (from 10% to 90%, by ten).
{
  dist_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(UAfilter, Scarfilter[,x], Abunfilter[,x])->dist_scar_abun
  distmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 1], ab)
  })
  scarmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 2], ab)
  })
  unlist(distmatrix)->unl.distmtx
  unlist(scarmatrix)->unl.scarmtx
  as.data.frame(unlist(scarmatrix))->indlist
  fac<-nrow(indlist)
  ind_matrix=matrix(NA, nrow = fac, ncol = 2)
  unl.distmtx->ind_matrix[,1]
  unl.scarmtx->ind_matrix[,2]
  UAfilter<-ind_matrix[order(ind_matrix[,2], decreasing = T),]
  nrow(UAfilter)->numsp
  UAfilter[,1]->UAfiltered
  output_scar=lapply(seq(1, numsp, by = 1), function(i) {
    mean(UAfiltered[1:i])
  })
  unlist(output_scar)->output_scar
  ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
    numsp*i/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    output_scar[ua.perc[w]]
  })
  as.numeric(output.ext)->output.ext
})

unlist(length.abun.perc)->length.abun.perc #Unlist output
perc.dist.first_abun=matrix(NA, nrow = 10, ncol = nSamp) #Creates a empty matrix to insert our results
length.abun.perc->perc.dist.first_abun[] #Input the results in the empty matrix
rownames(perc.dist.first_abun)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); perc.dist.first_abun
write.table(perc.dist.first_abun, file="length_percentage_first_abun_local.txt", sep="\t") #Saves output from the nine steps method.

#Random ----

RandPerm_mtx=matrix(NA,nrow=10,ncol=nSamp) #Creates a empty matrix to input our results
for (k in 1:nSamp) {
  randDIS=matrix(NA,nrow=10,ncol=1000)#ncol=randomizations
  dist_local %>% filter(dist_localfilter[, k] == "Present") -> UAfilter
  scar_local %>% filter(dist_localfilter[, k] == "Present") -> Scarfilter
  abun_t %>% filter(dist_localfilter[,k] == "Present") -> Abunfilter
  cbind(UAfilter, Scarfilter[,k], Abunfilter[,k])->dist_scar_abun
  distmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 1], ab)
  })
  scarmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,3]->ab
    esp=rep(dist_scar_abun[u, 2], ab)
  })
  unlist(distmatrix)->unl.distmtx
  unlist(scarmatrix)->unl.scarmtx
  as.data.frame(unlist(scarmatrix))->indlist
  fac<-nrow(indlist)
  ind_matrix=matrix(NA, nrow = fac, ncol = 2)
  unl.distmtx->ind_matrix[,1]
  unl.scarmtx->ind_matrix[,2]
  UAfilter<-ind_matrix
  nrow(UAfilter)->numsp
  as.data.frame(UAfilter[,1])->UAfiltered
  for (j in 1:1000) {
    sample(UAfiltered[,1])->UAfilter_sampled
    as.data.frame(UAfilter_sampled)->Null_Simulat
    means_dist=lapply(seq(1, nrow(Null_Simulat), by = 1), function(y) {
      mean(Null_Simulat[1:y,])
    })
    unlist(means_dist)->output_rand
    ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
      nrow(Null_Simulat)*i/100->perc
      round(perc, digits = 0)->perc
    })
    as.numeric(ua.perc)->ua.perc
    output.ext=lapply(seq(1, 10, by = 1), function(w){
      output_rand[ua.perc[w]]
    })
    as.numeric(output.ext)->output.ext
    randDIS[,j]<-output.ext
    print(paste("Wait for it... We are randomizing sample", k, "at the moment. This is permutation", j, "of 1000."))
  }
  RandPerm_mtx[,k]=rowMeans(randDIS)
} #Shuffles the order of extinction 1000 times for each sample while preserving species distinctiveness values. Nine levels of extraction for mean distiviness values are defined (10%-90%).
rownames(RandPerm_mtx)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); RandPerm_mtx

write.table(RandPerm_mtx, file="Length_percentage_random_local.txt", sep="\t") #Exports the output of the random extinction scenario.

#Friedman test ----

Friedman.X2=lapply(seq(1, 9, by = 1), function(i)
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  cbind(RandPerm_mtx[i,],perc.dist.first_abun[i,], perc.dist.first_scar[i,])->perc10[]
  friedman.test(perc10)->fried.test
  print(fried.test$statistic)
}) #Extract chi-squared values of Friedman's paired test for the comparisons between the distinctiveness of the three scenarios in all nine levels of extinction

Friedman.p=lapply(seq(1, 9, by = 1), function(i)
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(RandPerm_mtx[i,],perc.dist.first_abun[i,], perc.dist.first_scar[i,])->perc10[]
  friedman.test(perc10)->fried.test
  print(fried.test$p.value)
}) #Extract p values of Friedman's paired test for the comparisons between the distinctiveness of the three scenarios in all nine levels of extinction

Post.Hoc=lapply(seq(1, 9, by = 1), function(i) 
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(RandPerm_mtx[i,],perc.dist.first_abun[i,], perc.dist.first_scar[i,])->perc10[]
  posthoc.friedman.nemenyi.test(as.matrix(perc10))->post.fried
  print(post.fried$p.value)
}) #Returns a post-hoc pairwise test of mean distinctviness in the three scenarios for all nine levels of extinction 

capture.output(c(Friedman.X2, Friedman.p, Post.Hoc), file = "Friedman_Test_Length.txt") #Exports Friedman's test results

#Plot local pattern----
dist.rand=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  mean(RandPerm_mtx[x,], na.rm=T)->dist.rand[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 
dist.scar=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  mean(perc.dist.first_scar[x,], na.rm=T)->dist.scar[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 
dist.abun=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  mean(perc.dist.first_abun[x,], na.rm=T)->dist.abun[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 

#Inverts matrixes so results can be represented in incresing order of extinction (from 0-90%).
rev(dist.abun)->dist.abun
rev(dist.rand)->dist.rand
rev(dist.scar)->dist.scar

#Plotting the graph
erosion<-seq(0,90,10)
plot(erosion,dist.rand,type="l",ylim=c(min(dist.scar),max(dist.abun)),xlim=c(0,90),col="gray",pch=16,cex=1.5,xlab="Local	 Species	 Loss (%)",ylab="Body size",font.lab=2,cex.lab=1.4,cex.axis=1)
points(erosion[-1],dist.rand[-1],type="p",col="gray",pch=21,bg="gray",cex=1)
points(erosion,dist.scar,type="l",ylim=c(0,1),col="black",lty=2)
points(erosion[-1],dist.scar[-1],type="p",col="black",pch=21,bg="white",cex=1)
points(erosion,dist.abun,type="l",ylim=c(0,1),col="black")
points(erosion[-1],dist.abun[-1],type="p",col="black",pch=21,bg="black",cex=1)
legend(0,0.4060,legend=c("scarce extinct first","random extinction","abundant extinct first"),col=c("black","gray","black"),cex=1,lty=c(3,1,1),lwd=1.3,xjust=-0.01,yjust=0.5,horiz=F,text.font=0.5,x.intersp=0,bty="n")


##### N of functions ----
####################

#REGIONAL -----

traits.PA<-decostand(traits[,2:ncol(traits)], method = "pa", na.rm = TRUE) #Transforms the trait matrix in a binary table, where if a species plays that function it will be assigned with an 1, if it does not plays that function, then a 0 is assigned
data_trait<-cbind(ExtRiskReg, traits.PA); data_trait #Combine rariy index and traits to order
as.data.frame(rowSums(t(spp)))->oc_spp
data<-cbind(oc_spp, data_trait); data #Combines uniqueness and restrictedness
traitmatrix=lapply(seq(1, nrow(data), by = 1), function(u) {
  data[u,1]->ab
  esp=print(data_trait[rep(u, ab), 1:ncol(data_trait)],)
  esp
})

do.call(rbind.data.frame, traitmatrix)->ind_matrix
#first Sensible -----

traits_Rst<-ind_matrix[order(ind_matrix[, 1]),]; traits_Rst #Order species from the most restrict to the less (widespread).
write.table(traits_Rst, file="Restrict_classified_BAT_Temp.csv", sep="\t") #Saves species functions ordered.
traits.PA<-traits_Rst[, 2:ncol(ind_matrix)]; traits.PA #Extracts only species' functions as previously ordered, without the restrictiveness index.
ext <- specaccum(na.omit(traits.PA), "collector", permutations = 1000); ext #Extincts species from the most restrict to the less (widespread).
ext$richness->ext.Rst #Saves extinction pattern to a new object

#first Tolerant -----

traits_Wid<-ind_matrix[order(ind_matrix[,1], decreasing=T),]; traits_Wid #Order species from the most widespread to the less (restrict).
write.table(traits_Wid, file="widespread_classified_BAT_Temp.csv", sep="\t")#Saves species functions ordered.
traits.PA<-traits_Wid[, 2:ncol(ind_matrix)]; traits.PA #Extracts only species' functions as previously ordered, without the restrictiveness index.
ext <- specaccum(na.omit(traits.PA), "collector", permutations = 1000); ext #Extincts species from the most widespread to the less (restrict).
ext$richness->ext.Wid #Saves extinction pattern to a new object

#Random -----
ext <- specaccum(na.omit(traits.PA), "random", permutations = 1000); ext #Extincts species at random
ext$perm->ext_perm #Saves extinction pattern to a new object
Null_Ext=matrix(NA,nrow=nrow(na.omit(traits.PA)),ncol=3) 
for	(i	in	1:nrow(na.omit(traits.PA))) #Extracts mean values and confidence intervals of the random pattern 
{
  Null_Ext[i,1]=mean(ext_perm[i,], na.rm = TRUE)		
  Null_Ext[i,2]=quantile(ext_perm[i,],probs=0.025)
  Null_Ext[i,3]=quantile(ext_perm[i,],probs=0.975)		
}

#Plot regional extinction pattern ----
Spp_erosion<-c(nrow(na.omit(traits.PA)):1)
plot(Spp_erosion,Null_Ext[,1],type="l",xlim=c(0,nrow(na.omit(traits.PA))),xlab="Regional	 Species 	loss (%)",ylab="Number of functions",font.lab=2,cex.lab=1.4,cex.axis=1.2,col="gray",pch=16,cex=1.5,xaxt="n
     ",ylim=c(min(Null_Ext),max(Null_Ext, na.rm = TRUE)))
axis(side=1,	at=c(0,nrow(na.omit(traits.PA))/4,nrow(na.omit(traits.PA))/2,nrow(na.omit(traits.PA))/1.333,nrow(na.omit(traits.PA))),	
     labels=c(0,25,50,75,100),line=F,tick=-0.3,cex.axis=1.2,mgp=c(3,1,0))
polygon(c(nrow(na.omit(traits.PA)):1,rev(nrow(na.omit(traits.PA)):1)),c(Null_Ext[,2],rev(Null_Ext[,3])),col="grey88",border=
          F)
points(Null_Ext[nrow(Null_Ext):1,1],col="gray",cex=1.6,type="l",lty=1,lwd=4)
points(rev(ext.Rst),col="black",type="l",lty=3,lwd=4)
points(rev(ext.Wid),col="black",type="l",lty=1,lwd=4)

#Save Output ----

cbind(Null_Ext, ext.Wid, ext.Rst)->ext.regional
colnames(ext.regional)<-c("Random (Mean)", "Random (0.025)", "Random (0.975)", "first Widespread", "first Restrict"); ext.regional
write.table(ext.regional, file="Regional_ext_NFunc.txt", sep="\t")

#LOCAL -----
traits_local<-decostand(traits[,2:ncol(traits)], method = "pa", na.rm = TRUE) #Transforms the trait matrix in a binary table, where if a species plays that function it will be assigned with an 1, if it does not plays that function, then a 0 is assigned
traits_localfilter=t(spp)
traits_localfilter[traits_localfilter == 0] <- "Absent" #Replace all 0s in the presence-absence matrix to "Absent"
noquote(traits_localfilter)->traits_localfilter
replace(traits_localfilter, traits_localfilter == 1, "Present")->traits_localfilter #Replace all 1s in the presence-absence matrix to "Present"
dist_localfilter=traits_localfilter; dist_localfilter
scar_local<-as.data.frame(t(ExtRiskLoc)); scar_local #Transposes and converts the matrix of scarcity

#first Scarce ----
percentage.scar=lapply(seq(1, nSamp, by = 1), function(x) #Calculates the number of functions remaining in the community, removing mean values at the interval of 0-90% of species extinction in each sample 
{
  traits_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(Abunfilter[,x], Scarfilter[,x], UAfilter)->dist_scar_abun
  traitmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,1]->ab
    esp=print(dist_scar_abun[,2:ncol(dist_scar_abun)][rep(u, ab), ],)
    esp
  })
  do.call(rbind.data.frame, traitmatrix)->ind_matrix
  UAfilter<-ind_matrix[order(ind_matrix[,1]),]
  row.names(UAfilter)<-c(1:nrow(UAfilter))
  ext.ua = specaccum(na.omit(UAfilter[,2:ncol(UAfilter)]), method = "collector", permutations = 1000); ext.ua
  as.matrix(ext.ua$richness)->ex
  nrow(na.omit(UAfilter))->numsp
  ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
    numsp*i/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    ex[ua.perc[w],1]
  })
  as.numeric(output.ext)->output.ext
  print(output.ext)
})

ext.perc.scar = matrix(NA,nrow=10,ncol=nSamp) #Creates an empty matrix to insert results.
unlist(percentage.scar)->ext.perc.scar[] #Unlist results and insert results in the empty matrix
rownames(ext.perc.scar)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); ext.perc.scar
write.table(ext.perc.scar, file="Ext_percentage_scar_local.txt", sep="\t") #Saves results that consider intervals (mean values)

## Which is the first lost function in each sample? ## -----

FstLostFunction.Scar=lapply(seq(1, nSamp, by = 1), function(x) #Calculates the number of functions remaining in the community, removing mean values at the interval of 0-90% of species extinction in each sample 
{
  traits_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(Abunfilter[,x], Scarfilter[,x], UAfilter)->dist_scar_abun
  na.omit(dist_scar_abun)->dist_scar_abun
  traitmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,1]->ab
    esp=print(dist_scar_abun[,2:ncol(dist_scar_abun)][rep(u, ab), ],)
    esp
  })
  do.call(rbind.data.frame, traitmatrix)->ind_matrix
  UAfilter<-ind_matrix[order(ind_matrix[,1]),]
  UAfilter[, colSums(UAfilter != 0, na.rm = TRUE) > 0]->UAfilter
  IdentifyCat<-lapply(seq(nrow(UAfilter), 1, by=(-1)), function(y){
    colSums(UAfilter[1:y,2:ncol(UAfilter)])->AllCat
    NamesFunc<-lapply(seq(1:length(AllCat)), function(k) {
      if (AllCat[k] == 0) {
        print(names(AllCat[k]))
      } 
    })
    print(unlist(NamesFunc))
  }
  )
  print(unlist(IdentifyCat)[!duplicated(unlist(IdentifyCat))])
})

capture.output(FstLostFunction.Scar, file = "FstLostFunction.Scar.txt") #Extracts results

#first Abundant ----
percentage.abun=lapply(seq(1, nSamp, by = 1), function(x) #Calculates the number of functions remaining in the community, removing mean values at the interval of 0-90% of species extinction in each sample 
{
  traits_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(Abunfilter[,x], Scarfilter[,x], UAfilter)->dist_scar_abun
  traitmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,1]->ab
    esp=print(dist_scar_abun[,2:ncol(dist_scar_abun)][rep(u, ab), ],)
    esp
  })
  do.call(rbind.data.frame, traitmatrix)->ind_matrix
  UAfilter<-ind_matrix[order(ind_matrix[,1], decreasing = TRUE),]
  row.names(UAfilter)<-c(1:nrow(UAfilter))
  ext.ua = specaccum(na.omit(UAfilter[,2:ncol(UAfilter)]), method = "collector", permutations = 1000); ext.ua
  as.matrix(ext.ua$richness)->ex
  nrow(na.omit(UAfilter))->numsp
  ua.perc=lapply(seq(10, 100, by = 10), function(i){ 
    numsp*i/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    ex[ua.perc[w],1]
  })
  as.numeric(output.ext)->output.ext
  print(output.ext)
})

ext.perc.abun = matrix(NA,nrow=10,ncol=nSamp) #Creates an empty matrix to insert results.
unlist(percentage.abun)->ext.perc.abun[] #Unlist results and insert results in the empty matrix
rownames(ext.perc.abun)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); ext.perc.abun
write.table(ext.perc.abun, file="Ext_percentage_abun_local.txt", sep="\t") #Saves results that consider intervals (mean values)

## Which is the first lost function in each sample? ## -----

FstLostFunction.Abun=lapply(seq(1, nSamp, by = 1), function(x) #Calculates the number of functions remaining in the community, removing mean values at the interval of 0-90% of species extinction in each sample 
{
  traits_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(Abunfilter[,x], Scarfilter[,x], UAfilter)->dist_scar_abun
  na.omit(dist_scar_abun)->dist_scar_abun
  traitmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,1]->ab
    esp=print(dist_scar_abun[,2:ncol(dist_scar_abun)][rep(u, ab), ],)
    esp
  })
  do.call(rbind.data.frame, traitmatrix)->ind_matrix
  UAfilter<-ind_matrix[order(ind_matrix[,1], decreasing = TRUE),]
  UAfilter[, colSums(UAfilter != 0, na.rm = TRUE) > 0]->UAfilter
  IdentifyCat<-lapply(seq(nrow(UAfilter), 1, by=(-1)), function(y){
    colSums(UAfilter[1:y,2:ncol(UAfilter)])->AllCat
    NamesFunc<-lapply(seq(1:length(AllCat)), function(k) {
      if (AllCat[k] == 0) {
        print(names(AllCat[k]))
      } 
    })
    print(unlist(NamesFunc))
  }
  )
  print(unlist(IdentifyCat)[!duplicated(unlist(IdentifyCat))])
})


capture.output(FstLostFunction.Abun, file = "FstLostFunction.Abun.txt") #Extracts results

#Random ----
percentage=lapply(seq(1, nSamp, by = 1), function(x) #Calculates the number of functions remaining in the community, removing mean values at the interval of 0-90% of species extinction in each sample. Beyond, creates 1000 random scenarios for each sample, returning mean values of the number of functions in each step of extinction for each sample.
{
  traits_local %>% filter(dist_localfilter[,x] == "Present") -> UAfilter
  abun_t %>% filter(dist_localfilter[,x] == "Present") -> Abunfilter
  scar_local %>% filter(dist_localfilter[,x] == "Present") -> Scarfilter
  cbind(Abunfilter[,x], Scarfilter[,x], UAfilter)->dist_scar_abun
  traitmatrix=lapply(seq(1, nrow(dist_scar_abun), by = 1), function(u) {
    dist_scar_abun[u,1]->ab
    esp=dist_scar_abun[,2:ncol(dist_scar_abun)][rep(u, ab), ]
  })
  do.call(rbind.data.frame, traitmatrix)->ind_matrix
  UAfilter<-ind_matrix[,2:ncol(ind_matrix)]
  row.names(UAfilter)<-c(1:nrow(UAfilter))
  ext.ua = specaccum(na.omit(UAfilter), method = "random", permutations = 1000); ext.ua
  ext.ua$richness->ex
  as.matrix(ex)->ex
  nrow(na.omit(UAfilter))->numsp
  ua.perc=lapply(seq(10, 100, by = 10), function(x){ 
    numsp*x/100->perc
    round(perc, digits = 0)->perc
  })
  as.numeric(ua.perc)->ua.perc
  output.ext=lapply(seq(1, 10, by = 1), function(w){
    ex[ua.perc[w], 1]
  })
  as.numeric(output.ext)->output.ext
  print(output.ext)
})
ext.perc = matrix(NA,nrow=10,ncol=nSamp) #Creates a new matrix
unlist(percentage)->ext.perc[] #Unlist results and insert results in the empty matrix
rownames(ext.perc)=c("90%ext", "80%ext", "70%ext", "60%ext", "50%ext", "40%ext", "30%ext", "20%ext", "10%ext", "0%ext"); ext.perc

write.table(ext.perc, file="Ext_percentage_rand_local.txt", sep="\t") #Saves results that consider intervals (mean values)

#Friedman test ----

Friedman.X2=lapply(seq(1, 9, by = 1), function(i) #Extracts chi-shared values of the paired Friedman's test for each extinction step.
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(ext.perc[i,],ext.perc.abun[i,], ext.perc.scar[i,])->perc10[]
  friedman.test(perc10)->fried.test
  print(fried.test$statistic)
})

Friedman.p=lapply(seq(1, 9, by = 1), function(i) #Extracts p-values of the paired Friedman's test for each extinction step.
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(ext.perc[i,],ext.perc.abun[i,], ext.perc.scar[i,])->perc10[]
  friedman.test(perc10)->fried.test
  print(fried.test$p.value)
})

Post.Hoc=lapply(seq(1, 9, by = 1), function(i)  #Extracts results of the post-hoc test for each extinction step.
{ 
  perc10=matrix(NA,nrow=nSamp,ncol=3)
  colnames(perc10)=c("Rand", "Abun", "Scar")
  c(ext.perc[i,],ext.perc.abun[i,], ext.perc.scar[i,])->perc10[]
  posthoc.friedman.nemenyi.test(as.matrix(perc10))->post.fried
  print(post.fried$p.value)
})

capture.output(c(Friedman.X2, Friedman.p, Post.Hoc), file = "Friedman_TestNFunc.txt") #Combines e saves Friedman's test results

#Plot local pattern----
nfunc.rand=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  median(ext.perc[x,], na.rm=T)->nfunc.rand[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 
nfunc.scar=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  median(ext.perc.scar[x,], na.rm=T)->nfunc.scar[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 
nfunc.abun=matrix(NA,nrow=10,ncol=1) #Creates a new matrix to insert mean values of distinctiveness to each step of extinction
for (x in 10:1) {
  median(ext.perc.abun[x,], na.rm=T)->nfunc.abun[x,1]} #Calculates mean values of distinctiveness based in the values of all samples 

#Inverts matrixes so results can be represented in incresing order of extinction (from 0-90%).
rev(nfunc.abun)->nfunc.abun
rev(nfunc.rand)->nfunc.rand
rev(nfunc.scar)->nfunc.scar

#Plotting the graph
erosion<-seq(0,90,10)
plot(erosion,nfunc.rand,type="l",ylim=c(min(nfunc.scar),max(nfunc.rand)),xlim=c(0,90),col="gray",pch=16,cex=1.5,xlab="Local	 Species	 Loss (%)",ylab="Number  of  functions",font.lab=2,cex.lab=1.4,cex.axis=1)
points(erosion[-1],nfunc.rand[-1],type="p",col="gray",pch=21,bg="gray",cex=1)
points(erosion,nfunc.scar,type="l",ylim=c(0,1),col="black",lty=2)
points(erosion[-1],nfunc.scar[-1],type="p",col="black",pch=21,bg="white",cex=1)
points(erosion,nfunc.abun,type="l",ylim=c(0,1),col="black")
points(erosion[-1],nfunc.abun[-1],type="p",col="black",pch=21,bg="black",cex=1)
legend(-5,3.5,legend=c("scarce extinct first","random extinction","abundant extinct first"),col=c("black","gray","black"),cex=1,lty=c(3,1,1),lwd=1.3,xjust=-0.01,yjust=0.5,horiz=F,text.font=0.5,x.intersp=0,bty="n")
