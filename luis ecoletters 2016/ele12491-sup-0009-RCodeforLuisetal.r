################################################################### 
#    R Code for analyses in Luis, et al. 2015 "Network analysis of
#    host-virus communities in bats and rodents reveals determinants
#    of cross-species transmission"
###################################################################




################################################################### 
#### import the data
###################################################################
batvirusdata=read.csv("batvirusdata.csv")

battraits=read.csv("battraitdata.csv") 
#####change some of the traits to factors (categorical rather than numerical)
battraits$cave=factor(battraits$cave)
battraits$migration=factor(battraits$migration)
battraits$torpor.category=factor(battraits$torpor.category)



###################################################################
#  organize bat virus data into lists then make the adjacency matrix
###################################################################
batnames=sort(unique(batvirusdata$binomial))

#create a list of each bat species and their viruses
bat.host.virus.list=list(NA) 
for(i in 1:length(batnames)){
bat.host.virus.list[[i]]= batvirusdata $virus.species[which(batvirusdata $binomial==batnames[i])]
	}
names(bat.host.virus.list)=batnames	

#create the adjacency matrix
batsharedvirusesmat=matrix(NA,ncol=length(batnames),nrow=length(batnames))
colnames(batsharedvirusesmat)= batnames
rownames(batsharedvirusesmat)= batnames

for (i in 1:length(batnames)){
	for (j in 1:length(batnames)){	 
	batsharedvirusesmat[i,j]=length(which(is.finite(match(unique(bat.host.virus.list[[i]]),unique(bat.host.virus.list[[j]])))))
	}
}


###################################################################
# using the igraph package, create the network
###################################################################

library(igraph)

batsharedvirusesgraph=graph.adjacency(batsharedvirusesmat,mode="undirected",weighted=TRUE,diag=FALSE)

summary(batsharedvirusesgraph) #gives number of nodes and links

###### some basic network statistics
#connectance=links/(nodes^2)
degree(batsharedvirusesgraph)
degree.distribution(batsharedvirusesgraph)
betweenness(batsharedvirusesgraph,directed=FALSE)
transitivity(batsharedvirusesgraph)


###################################################################
# Phylogenetic Matrices
###################################################################

library(ape)

tree=read.tree("Mammaliansupertree.txt")
del=1:length(tree$tip.label)
batnames2=sub(" ","_",batnames)
which.tip.label.allbats=match(batnames2,tree$tip.label)
which.tip.label.allbats=which.tip.label.allbats[which(is.finite(which.tip.label.allbats))] #remove R. horsfeldi which prob doesn't exist
bat.tree=drop.tip(tree,del[-which.tip.label.allbats])
plot.phylo(bat.tree)

#phylgenetic correlation matrix
batphylo=vcv.phylo(bat.tree, model = "Brownian", cor = TRUE)
batphylo = batphylo[sort(rownames(batphylo)),sort(colnames(batphylo))]


###################################################################
# Sympatry Matrices
###################################################################

library(sp)
library(maptools)
library(rgeos)

#shape files for the distributions of every terrestrial mammal in the IUCN database, can be downloaded from http://www.iucnredlist.org/technical-documents/spatial-data
mammterr=readShapeSpatial("MAMMTERR.shp") 

summary(mammterr)
class(mammterr)

#all of the species names in the shape files
all.binomial=as.character(mammterr$BINOMIAL)

#### index of which records to keep (bat species in database)
keep=list()
for(i in 1:length(batnames)){
keep[[i]]=which(all.binomial== batnames[i])
}
x=keep[[1]]
for(i in 2:length(batnames)){
x=c(x,keep[[i]])
}
keep=x


batspecies.distr=mammterr[keep,]


##create a separate file for for each species (many have multiple polygons)
for(i in 1:length(batnames)){
nam=paste("batpolygon.sp",i,sep="")

if(length(which(as.character(batspecies.distr$BINOMIAL)== batnames[i])>0)){
assign(nam,
batspecies.distr[which(as.character(batspecies.distr$BINOMIAL)== batnames[i]),])
}
}


#### create matrix with a 1 if the two species distributions' overlap, and 0 if not
batsympatry=matrix(NA,length(batnames),length(batnames),dimnames=list(batnames, batnames))

index=which(lower.tri(batsympatry),arr.ind=TRUE)

for(i in 1:dim(index)[1]){
if(exists(paste("batpolygon.sp",index[i,1],sep=""))&exists(paste("batpolygon.sp",index[i,2],sep=""))){
batsympatry[index[i,1],index[i,2]]=
ifelse(sum(over(get(paste("batpolygon.sp",index[i,1],sep="")),
    get(paste("batpolygon.sp",index[i,2],sep=""))),na.rm=TRUE)>0,1,0)
    
    }
}


#make symmetric
diag(batsympatry)=1

indl=which(lower.tri(batsympatry),arr.ind=TRUE)
indu=which(lower.tri(batsympatry),arr.ind=TRUE)[,c(2,1)]
batsympatry[indu]=batsympatry[indl]

batsympatry= batsympatry[-113,-113] ###########remove R. horsfeldi which probably doens't exist


###################################################################
# Multiple regression on distance matrices
###################################################################

library(ecodist)


#function to run models and put in table ordered by R squared

MRM.table=
function(mods){
	n=length(mods)
	table=data.frame(Model=rep(NA,n),R.squared=rep(NA,n),p=rep(NA,n))
	
	for(i in 1:n){
		mod=MRM(formula(mods[i]),nperm=10000)
		table[i,1]=mods[i]
		table[i,2]=mod$r.squared[1]
		table[i,3]=mod$r.squared[2]		
	}
	table=table[order(table$R.squared,decreasing = TRUE),]
	table
}


## create a matrix for sampling effort with the product of the 2 species' log citations
batcitationsmat=matrix(NA,dim(batsympatry),dim(batsympatry))
colnames(batcitationsmat)=colnames(batsympatry)
rownames(batcitationsmat)=rownames(batsympatry)
for(i in 1:dim(batcitationsmat)[1]){
	for(j in 1:dim(batcitationsmat)[1]){
		batcitationsmat[i,j]=battraits$log.citations[which(battraits$binomial==rownames(batcitationsmat)[i])]* battraits$log.citations[which(battraits$binomial==rownames(batcitationsmat)[j])]
		}
}


y=as.dist(batsharedvirusesmat[match(rownames(batsympatry),rownames(batsharedvirusesmat)),match(colnames(batsympatry),colnames(batsharedvirusesmat))])
sympatry=as.dist(batsympatry)
phylo=as.dist(batphylo)
cit=as.dist(batcitationsmat)

models=c("y~phylo","y~overlap","y~phylo+overlap","y~cit2","y~phylo+cit2","y~overlap+cit2","y~phylo+overlap+cit2","y~overlap+phylo+cit2")
MRM.table(models)



###################################################################
# PGLS Models
###################################################################


### function to run gls model and optimize Pagel's lambda
pgls.lambda=function (model, #model to be run
					  VCV, #variance-covariance matrix, here the phylogenetic correlation matrix
					  lam)  # a starting value for optimization
{    tree.weights <- varIdent(value = diag(VCV), fixed = T)
    tree.cor <- corSymm(VCV[lower.tri(VCV)], fixed = T)
    lam = logit(lam)
    lmax <- optim(lam, fn = calc.lam.l, model = model, tree.cor = tree.cor, 
        tree.weights = tree.weights, method="Brent",hessian = TRUE,lower=-100,upper=100)
    model.lmax <<- gls(model = model, cor = revlogit(lmax[[1]]) * 
        tree.cor, weights = tree.weights, method = "ML", control = glsControl(maxIter = 1e+05))
    n.obs <- as.numeric(model.lmax$dims[1])
    loglik <- logLik(model.lmax)[1]
    n.par <- as.numeric(model.lmax$dims[2]) + 1
    se <- sqrt(diag(model.lmax$varBeta))
    ucl <- coef(model.lmax) + qt(0.95, n.obs - n.par) * se
    lcl <- coef(model.lmax) - qt(0.95, n.obs - n.par) * se
    aic <- -2 * loglik + 2 * n.par
    aicc <- -2 * loglik + 2 * n.par + (2 * n.par * (n.par + 1))/(n.obs - 
        n.par - 1)
    pgls.lam.out = list(model, model.lmax, coef(model.lmax), 
        se, ucl, lcl, n.obs, n.par, loglik, aic, aicc, revlogit(lmax$par), 
        lmax$convergence)
    names(pgls.lam.out) = c("model.name", "model", "summary", 
        "se", "ucl", "lcl", "n.obs", "n.par", "logLik", "AIC", 
        "AICc", "lambda", "converge")
    pgls.lam.out
}

######functions used in above
logit=function(x){
	log(x/(1-x))}
revlogit=function(x){
	exp(x)/(1+exp(x))}

calc.lam.l=function(lam,model,tree.weights,tree.cor) {
    m1 <- gls(model=model,cor=revlogit(lam)*tree.cor,weights=tree.weights,method="ML",control = glsControl(opt="optim",msMaxIter=200,maxIter = 100000))
    return(-logLik(m1))
  }

###function to produce an output table
LAICc.table=function (mods) #list of models to be run
{   n = length(mods)
    table = data.frame(model.name = rep(NA, n), AICc = rep(NA, 
        n), npar = rep(NA, n), weight = rep(NA, n), lambda = rep(NA, n))
    for (i in 1:n) {
        table[i, 1] = as.character(mods[[i]]$model.name[3])
        table[i, 2] = mods[[i]]$AICc
        table[i, 3] = mods[[i]]$n.par
        table[i, 5] = mods[[i]]$lambda
    }
    table = table[order(table$AIC), ]
    n = dim(table)[1]
    deltaAIC = rep(NA, n)
    dvec = rep(NA, n)
    for (i in 1:n) {
        deltaAIC[i] = table[i, 2] - min(table[, 2])
        dvec[i] = exp(-0.5*deltaAIC[i])
        }
     dsum = sum(dvec)
    for (i in 1:n) {
        table[i, 4] = exp(-0.5*deltaAIC[i])/dsum
    }
    table
}


reduced.battraits= battraits [which(!is.na(battraits$PC1)),]

#y= reduced.battraits $num.viruses
#y= reduced.battraits $betweenness
y= reduced.battraits $degree
cit= reduced.battraits $log.citations
PC1= reduced.battraits $PC1
symp= reduced.battraits $num.sympatric
area= reduced.battraits $distr.area
lat= reduced.battraits $Latitude.abs
PC2= reduced.battraits $PC2
PC3= reduced.battraits $PC3
IUCN= reduced.battraits $IUCN
torpor= reduced.battraits $torpor.category
migration= reduced.battraits $migration ###problem here
greg= reduced.battraits $gregariousness
cave= reduced.battraits $cave
diet= reduced.battraits $food

reduced.batphylo=batphylo[match(reduced.battraits$binomial,sub("_"," ",rownames(batphylo))),match(reduced.battraits$binomial,sub("_"," ",rownames(batphylo)))]



library(nlme)

GLSbat1=pgls.lambda(y~cit+greg+symp+diet, VCV=reduced.batphylo, lam=0)
GLSbat2=pgls.lambda(y~cit+greg, VCV=reduced.batphylo, lam=0)
GLSbat3=pgls.lambda(y~cit+greg+symp, VCV=reduced.batphylo, lam=0)
GLSbat4=pgls.lambda(y~cit+symp+PC2, VCV=reduced.batphylo, lam=0)
GLSbat5=pgls.lambda(y~cit+PC1+greg, VCV=reduced.batphylo, lam=0)
GLSbat6=pgls.lambda(y~cit+symp, VCV=reduced.batphylo, lam=0)
GLSbat7=pgls.lambda(y~cit, VCV=reduced.batphylo, lam=0)
GLSbat8=pgls.lambda(y~greg+symp, VCV=reduced.batphylo, lam=0)
GLSbat9=pgls.lambda(y~cit+PC1, VCV=reduced.batphylo, lam=0)
GLSbat10=pgls.lambda(y~cit+symp+migration, VCV=reduced.batphylo, lam=0)
GLSbat11=pgls.lambda(y~cit+lat, VCV=reduced.batphylo, lam=0)
GLSbat12=pgls.lambda(y~cit+symp+diet, VCV=reduced.batphylo, lam=0)
GLSbat13=pgls.lambda(y~cit+diet, VCV=reduced.batphylo, lam=0)
GLSbat14=pgls.lambda(y~cit+symp+PC1, VCV=reduced.batphylo, lam=0)
GLSbat15=pgls.lambda(y~cit+greg+cave, VCV=reduced.batphylo, lam=0)
GLSbat16=pgls.lambda(y~cit+migration, VCV=reduced.batphylo, lam=0)
GLSbat17=pgls.lambda(y~cit+greg+torpor, VCV=reduced.batphylo, lam=0)
GLSbat18=pgls.lambda(y~cit+greg+symp+torpor, VCV=reduced.batphylo, lam=0)
GLSbat19=pgls.lambda(y~cit+area, VCV=reduced.batphylo, lam=0)
GLSbatn=pgls.lambda(y~1,VCV=reduced.batphylo, lam=0)


mods=list(GLSbat1,GLSbat2,GLSbat3,GLSbat4,GLSbat5,GLSbat6, GLSbat7,GLSbat8,GLSbat9,GLSbat10,GLSbat11,GLSbat12,GLSbat13,
 GLSbat14,GLSbat15,GLSbat16,GLSbat17,GLSbat18,GLSbat19,GLSbatn)


LAICc.table(mods)

###################################################################
# Community Detection 
###################################################################

## We used Matlab code from http://perso.uclouvain.be/vincent.blondel/research/louvain.html
#reference: Blondel et al. 2008. Fast unfolding of communities in large networks. Journal of Statistical Mechanics: Theory and Experiment. 10:P10008.


## save the communities as a csv file 'batcomsMatlab.csv'

## import and make into a list of species names in each community
batcomsMatlab=read.csv("batcomsMatlab.csv")
bc=sort(unique(unlist(batcomsMatlab)))
batcoms=list()
for(i in bc){
	batcoms[[i]]=rownames(batsharedvirusesmat)[which(batcoms0matlab==i)]
}



###################################################################
# Arc Diagram
###################################################################

# install devtools
#install.packages("devtools")

# load devtools
#library(devtools)

# install arcdiagram
#install_github('arcdiagram', username='gastonstat')


library(arcdiagram)
library(RColorBrewer)


# get edgelist
edgelist = get.edgelist(batsharedvirusesgraph)

# get vertex labels
#arcdiagram orders by edgelist not by igraph vertex labels- so doesn't include those species that don't have edges.
nam=character()
for(i in 1:dim(edgelist)[1]){
nam=c(nam,edgelist[i,])	
}
vlabels = unique(nam)


# vertex groups (communities)
b=character()
com=numeric()
for(i in 1:length(batcoms)){
	b=c(b,as.character(batcoms[[i]]))
	com=c(com,rep(i,length(batcoms[[i]])))
}
vgroups = com[match(vlabels,b)]
#rename groups so that largest communities first
x=sort(unique(vgroups))
lx=numeric()
for(i in 1:length(x)){
	lx [i]=length(which(vgroups==x[i]))
}
y=x[order(lx,decreasing=TRUE)]
vgroups2=numeric()
for(i in 1:length(vgroups)){
	vgroups2[i]=which(y==vgroups[i])
}


# get vertex fill color
col1=brewer.pal(12,"Paired")
cols=rep(col1[c(1,3,5,7,9,11)],3)
vfill = rep("gray",length(vgroups2))
for(i in 1:length(y)){
vfill[which(vgroups2==i)]=cols[i]
}
# get vertex border color
cols=rep(col1[c(2,4,6,8,10,12)],3)
vborders = rep("gray",length(vgroups2))
for(i in 1:length(y)){
vborders[which(vgroups2==i)]=cols[i]
}


# get vertex degree
degrees=batdata.network.njab$degree[match(vlabels,batdata.network.njab$binomial)]

# get edges value
values = get.edge.attribute(batsharedvirusesgraph.nr.njab, "weight")

#sort by community, then degree
ord=order(vgroups2,max(degrees)-degrees)



#with species labels 
par(mar=c(8,0,2,0))
arcplot(edgelist, ordering=ord, labels=vlabels, cex.labels=0.6,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = degrees/max(degrees)+0.6, pch.nodes=21,col.labels="grey30",
        lwd.nodes = 1.5, line=-0.6, lwd.arcs = 2*values,
        col.arcs = "#5998ff30")


#without species labels 
par(mar=c(2,0,2,0))
arcplot(edgelist, ordering=ord, labels=vlabels, show.labels=FALSE,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = degrees/max(degrees)+0.8, pch.nodes=21,
        lwd.nodes = 1.65, line=-0.5, lwd.arcs = 2*values,
        col.arcs = "#5998ff30")



###################################################################
# Maps
###################################################################

library(sp)
library(maptools)
library(maps)


## Colored Community maps

x=batcoms[[1]]
tiff("batCommBlue.tiff",width=6,height=2.6,units="in",res=200);map("world",col="gray",mar=c(1,1,1,1),xlim=c(-180,180),ylim=c(-60,90),lwd=2)
 for(i in 1:length(x)){
 plot(batspecies.distr[which(batspecies.distr $BINOMIAL==x[i]),],add=TRUE,col="#1F78B450",border=FALSE)}
dev.off()


x=batcoms[[2]]
tiff("batCommGreen.tiff",width=6,height=2.6,units="in",res=200);map("world",col="gray",mar=c(1,1,1,1),xlim=c(-180,180),ylim=c(-60,90),lwd=2)
 for(i in ord){
 plot(batspecies.distr[which(batspecies.distr $BINOMIAL==x[i]),],add=TRUE,col="#33A02C50",border=FALSE)}
dev.off()


x=batcoms[[3]]
tiff("batCommRed.tiff",width=6,height=2.6,units="in",res=200);map("world",col="gray",mar=c(1,1,1,1),xlim=c(-180,180),ylim=c(-60,90),lwd=2)
 for(i in ord){
 plot(batspecies.distr[which(batspecies.distr $BINOMIAL==x[i]),],add=TRUE,col="#E31A1C70",border=FALSE)}
dev.off()


x=batcoms[[4]]
tiff("batCommOrange.tiff",width=6,height=2.6,units="in",res=200);map("world",col="gray",mar=c(1,1,1,1),xlim=c(-180,180),ylim=c(-60,90),lwd=2)
 for(i in ord){
 plot(batspecies.distr[which(batspecies.distr $BINOMIAL==x[i]),],add=TRUE,col="#FF7F0070",border=FALSE)}
dev.off()



###################################################################
# Quantitative linkage density function
###################################################################

qLD=function(mat){
  Hr=numeric()
  h=numeric()
  for(i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[1]){
      h[j]=mat[i,j]/sum(mat[i,])*log(mat[i,j]/sum(mat[i,]))
    }
    Hr[i]=-sum(h,na.rm=TRUE)
    
  }
  
  Nr=exp(Hr)
  x=numeric()
  for(i in 1:dim(mat)[1]){
    x[i]=sum(mat[i,])*Nr[i]
  }
  
  sum(x)/sum(mat)
}

###################
bat.qLD=qLD(batsharedvirusesmat) #quantitative linkage density of the bat network
bat.qConnectance=bat.qLD/dim(batsharedvirusesmat)[1] #quantitative connectance of the bat network

rodent.qLD=qLD(rodentsharedvirusesmat.njab)
rodent.qConnectance=rodent.qLD/dim(rodentsharedvirusesmat.njab)[1]



###################################################################
# Create networks adjusted by sampling effort
###################################################################

batedgelist=data.frame(edge=get.edgelist(batsharedvirusesgraph),weight=get.edge.attribute(g, name="weight"))
log.cit1=numeric()
log.cit2=numeric()
lower.cit=numeric()
for(i in 1:dim(batedgelist)[1]){
  log.cit1[i]=battraits$log.citations[which(battraits$binomial==batedgelist$edge.1[i])]
  log.cit2[i]=battraits$log.citations[which(battraits$binomial==batedgelist$edge.2[i])]
  lower.cit[i]=min(c(log.cit1[i],log.cit2[i]))
}
batedgelist$log.cit1=log.cit1
batedgelist$log.cit2=log.cit2
batedgelist$lower.cit = lower.cit

##############################

rodentedgelist=data.frame(edge=get.edgelist(rodentsharedvirusesgraph),weight=get.edge.attribute(g, name="weight"))
log.cit1=numeric()
log.cit2=numeric()
lower.cit=numeric()
for(i in 1:dim(rodentedgelist)[1]){
  log.cit1[i]=rodenttraits$log.citations[which(rodenttraits$binomial==rodentedgelist$edge.1[i])]
  log.cit2[i]=rodenttraits$log.citations[which(rodenttraits$binomial==rodentedgelist$edge.2[i])]
  lower.cit[i]=min(c(log.cit1[i],log.cit2[i]))
}
rodentedgelist$log.cit1=log.cit1
rodentedgelist$log.cit2=log.cit2
rodentedgelist$lower.cit = lower.cit


#################

bothedgelist=rbind(batedgelist,rodentedgelist)
plot(bothedgelist$lower.cit,bothedgelist$weight)
bothweightLM=lm(bothedgelist$weight~bothedgelist$lower.cit-1)
abline(bothweightLM)

bothedgelist$residuals.lowercit =residuals(bothweightLM)

bothedgelist$residuals.lowercit.std=bothedgelist$residuals.lowercit-min(bothedgelist$residuals.lowercit)+1

batedgelist$adjusted.weights=bothedgelist$residuals.lowercit.std[1:dim(batedgelist)[1]]
rodentedgelist$adjusted.weights=bothedgelist$residuals.lowercit.std[(dim(batedgelist)[1]+1):dim(bothedgelist)[1]]


##################################### make the new viral sharing matrix
#####bats
batsharedvirusesmat.adjustedweights1=matrix(0,nrow=dim(batsharedvirusesmat)[1],ncol=dim(batsharedvirusesmat)[2])
colnames(batsharedvirusesmat.adjustedweights1)=colnames(batsharedvirusesmat)
rownames(batsharedvirusesmat.adjustedweights1)=rownames(batsharedvirusesmat)

for(i in 1:dim(batedgelist)[1]){
  sp1=which(colnames(batsharedvirusesmat.adjustedweights1)==batedgelist$edge.1[i])
  sp2=which(colnames(batsharedvirusesmat.adjustedweights1)==batedgelist$edge.2[i])
  batsharedvirusesmat.adjustedweights1[sp1,sp2]= batsharedvirusesmat.adjustedweights1[sp2,sp1]= batedgelist$adjusted.weights[i]
}

####rodents
rodentsharedvirusesmat.adjustedweights1=matrix(0,nrow=dim(rodentsharedvirusesmat)[1],ncol=dim(rodentsharedvirusesmat)[2])
colnames(rodentsharedvirusesmat.adjustedweights1)=colnames(rodentsharedvirusesmat)
rownames(rodentsharedvirusesmat.adjustedweights1)=rownames(rodentsharedvirusesmat)

for(i in 1:dim(rodentedgelist)[1]){
  sp1=which(colnames(rodentsharedvirusesmat.adjustedweights1)==rodentedgelist$edge.1[i])
  sp2=which(colnames(rodentsharedvirusesmat.adjustedweights1)==rodentedgelist$edge.2[i])
  rodentsharedvirusesmat.adjustedweights1[sp1,sp2]= rodentsharedvirusesmat.adjustedweights1[sp2,sp1]= rodentedgelist$adjusted.weights[i]
}


bat.qLD.adj=qLD(batsharedvirusesmat.adjustedweights1) #quantitative linkage density of the sampling effort-adjusted bat network
bat.qConnectance.adj=bat.qLD.adj/dim(batsharedvirusesmat)[1] #quantitative connectance of the sampling effort-adjusted bat network
rodent.qLD.adj=qLD(rodentsharedvirusesmat.adjustedweights1)
rodent.qConnectance.adj=rodent.qLD.adj/dim(rodentsharedvirusesmat.njab)[1]


###################################################################
# Calculate p values for Table 1 
# by creating 10,000 permutations of the bat and rodent networks
# and comparing the differences in the observed versus permuted statistics
###################################################################

bnames=rownames(batsharedvirusesmat)
bat.edges=t(combn(bnames,2))
bat.edges=as.data.frame(bat.edges)
weight=numeric()
for(i in 1:dim(bat.edges)[1]){
  sp1=which(bnames==bat.edges[i,1])
  sp2=which(bnames==bat.edges[i,2])  
  weight[i]=batsharedvirusesmat[sp1,sp2]
}
bat.edges$weight=weight

rnames=rownames(rodentsharedvirusesmat.njab)
rodent.edges=t(combn(rnames,2))
rodent.edges=as.data.frame(rodent.edges)
weight=numeric()
for(i in 1:dim(rodent.edges)[1]){
  sp1=which(rnames==rodent.edges[i,1])
  sp2=which(rnames==rodent.edges[i,2])  
  weight[i]=rodentsharedvirusesmat.njab[sp1,sp2]
}
rodent.edges$weight=weight


both.edges=rbind(bat.edges,rodent.edges)
bothnames=c(bnames,rnames)

n=10000 #number of permutations
both.edges.permuteweight=matrix(NA,nrow=dim(both.edges)[1],ncol=n)
for(i in 1:n){
  both.edges.permuteweight[,i]=sample(both.edges$weight,length(both.edges$weight))
}


both.permute.data=data.frame(bat.assortativity=rep(NA,n),bat.transitivity=rep(NA,n),bat.meandegree=rep(NA,n),bat.mean.wt.degree=rep(NA,n),bat.links=rep(NA,n),bat.qLdens=rep(NA,n),
                             rodent.assortativity=rep(NA,n),rodent.transitivity=rep(NA,n),rodent.meandegree=rep(NA,n),rodent.mean.wt.degree=rep(NA,n),rodent.links=rep(NA,n),rodent.qLdens=rep(NA,n))

for(i in 1:n){
  EL=both.edges
  EL$weight=both.edges.permuteweight[,i]
  adjmat=matrix(0,ncol=length(bothnames),nrow=length(bothnames),dimnames=list(bothnames,bothnames))
  for(j in 1:dim(both.edges)[1]){
    sp1=which(bothnames==both.edges[j,1])
    sp2=which(bothnames==both.edges[j,2])
    adjmat[sp1,sp2]=adjmat[sp2,sp1]=EL$weight[j]
  }
  graph=graph.adjacency(adjmat[1:143,1:143],mode="undirected",weighted=TRUE,diag=FALSE) #just the part of adjmat that relates to bats
  both.permute.data$bat.assortativity[i]=assortativity.degree(graph,directed=FALSE)
  both.permute.data$bat.transitivity[i]=transitivity(graph)
  both.permute.data$bat.meandegree[i]=mean(degree(graph))
  both.permute.data$bat.mean.wt.degree[i]=mean(apply(adjmat[1:143,1:143],1,sum,na.rm=TRUE))
  both.permute.data$bat.links[i]=length(E(graph))
  both.permute.data$bat.qLdens[i]=qLD(adjmat[1:143,1:143])
  
  graph=graph.adjacency(adjmat[144:339,144:339],mode="undirected",weighted=TRUE,diag=FALSE) #just the part of adjmat that relates to rodents
  both.permute.data$rodent.assortativity[i]=assortativity.degree(graph,directed=FALSE)
  both.permute.data$rodent.transitivity[i]=transitivity(graph)
  both.permute.data$rodent.meandegree[i]=mean(degree(graph))
  both.permute.data$rodent.mean.wt.degree[i]=mean(apply(adjmat[144:339,144:339],1,sum,na.rm=TRUE))
  both.permute.data$rodent.links[i]=length(E(graph))
  both.permute.data$rodent.qLdens[i]=qLD(adjmat[144:339,144:339])
  
  cat(i)
}


# for each permutation compare the difference in values for the bat vs rodent network
# how many times out of 10,000 permutations was the result more extreme than the observed difference
both.permute.data$meandegree.diff=both.permute.data$bat.meandegree-both.permute.data$rodent.meandegree
both.permute.data$mean.wt.degree.diff=both.permute.data$bat.mean.wt.degree-both.permute.data$rodent.mean.wt.degree
both.permute.data$connect.diff=both.permute.data$bat.links/143^2-both.permute.data$rodent.links/196^2
both.permute.data$links.diff=both.permute.data$bat.links-both.permute.data$rodent.links
both.permute.data$quantconnect.diff=both.permute.data$bat.qLdens/143-both.permute.data$rodent.qLdens/196

# (also do this for the sampling effort- adjusted network for mean weighted degree and quantitative connectance)


###################################################################
# Clique Percolation
###################################################################


clique.mat.fun=function(clique.list){
  mat=matrix(NA,length(clique.list),length(clique.list))
  for(i in 1:length(clique.list)){
    for(j in 1:length(clique.list)){
      mat[i,j]=length(which(is.finite(match(clique.list[[i]], clique.list[[j]]))))
    }
  }
  mat
}


k.clique.mat.fun=function(matrix,k){
  mat=matrix(NA,dim(matrix)[1],dim(matrix)[2])
  for(i in 1:dim(matrix)[1]){
    for(j in 1:dim(matrix)[2]){
      if(i==j){
        mat[i,j]=ifelse(matrix[i,j]>=(k-1),1,0)
      }
      if(i!=j){
        mat[i,j]=ifelse(matrix[i,j]>=k,1,0)
        
      }
    }
  }
  mat
}


k.clique.comm.fun=function(mat,k,w=0){
  if(w>0){
    mat=replace(mat,mat<w,0)
  }
  
  graph=graph.adjacency(mat,mode="undirected",weighted=TRUE,diag=FALSE)
  
  max.cliques=maximal.cliques(graph)
  clique.mat=clique.mat.fun(max.cliques)
  k.clique.mat=k.clique.mat.fun(clique.mat,k)
  
  comms=list()
  comms[[1]]=which(k.clique.mat[1,]==1)	
  
  for(i in 2:dim(k.clique.mat)[1]){
    x=which(k.clique.mat[i,]==1)	
    x2=numeric()
    for(j in (1:length(comms))){
      x2=c(x2,sum(is.finite(match(x,comms[[j]])))>0) 
    }
    if(sum(x2)==0){
      comms[[length(comms)+1]]=x
    }
    if(sum(x2)>0){
      if(sum(x2)==1){
        comms[[which(x2>0)]]=unique(c(comms[[which(x2>0)]],x))
      }
      
      if(sum(x2>1)){
        x3=which(x2>0)
        comms[[x3[1]]]=unique(c(comms[[x3[1]]],x)) 
        for(k in 2:length(x3)){
          comms[[x3[1]]]=unique(c(comms[[x3[1]]],comms[[x3[k]]]))
          comms[[x3[k]]]=NA
        }	
      }
    }
  }
  
  keep=(1:length(comms))[-which(lapply(comms,length)==0)]
  new.comms=list()
  for(i in 1:length(keep)){
    new.comms[[i]]=comms[[keep[i]]]
  }
  
    k.clique.communities=list()
  for(i in 1:length(new.comms)){
    x=new.comms[[i]]
    k.clique.communities[[i]]=character()
    for(j in 1:length(x)){
      
      k.clique.communities[[i]]=unique(c(k.clique.communities[[i]],V(graph)$name[max.cliques[[x[j]]]]))
    }
  }
  
  k.clique.communities
}

###############################
 
bat.k3.clique.communities=k.clique.comm.fun(batsharedvirusesmat,3)
bat.k4.clique.communities=k.clique.comm.fun(batsharedvirusesmat,4)
bat.k5.clique.communities=k.clique.comm.fun(batsharedvirusesmat,5)
bat.k6.clique.communities=k.clique.comm.fun(batsharedvirusesmat,6)
bat.k7.clique.communities=k.clique.comm.fun(batsharedvirusesmat,7)

rodent.k3.clique.communities=k.clique.comm.fun(rodentsharedvirusesmat,3)
rodent.k4.clique.communities=k.clique.comm.fun(rodentsharedvirusesmat,4)
rodent.k5.clique.communities=k.clique.comm.fun(rodentsharedvirusesmat,5)
rodent.k6.clique.communities=k.clique.comm.fun(rodentsharedvirusesmat,6)
rodent.k7.clique.communities=k.clique.comm.fun(rodentsharedvirusesmat,7)
##################################################################################

