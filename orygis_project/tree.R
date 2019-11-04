  
library("ape")
library("phytools")
setwd('~/gitprojects/mtbdiff/orygis_project/')
anim <- read.table('ani_matrix.csv',sep=',',header=TRUE, row.names=1)
samples <- read.table('orygis_samples.csv',sep=',',header=TRUE,row.names=1)
#rownames(anim) <- anim$query

distances <- as.matrix(1 - (anim/100))
dm <- dist(anim)
njtree  <- nj(distances)
rooted <- root(njtree,"MTB_beijing")
rooted$edge.length[rooted$edge.length<0] <- 0
labels <- samples[rooted$tip.label,]$host
cols<-setNames(c('red','blue','green2','white'),levels(as.factor(labels)))
png('orygis_tree.png',width=800,height=800)
plot(rooted, cex=2.5,label.offset=.00005)
tiplabels(pie=to.matrix(labels, levels(labels)),cex=0.3,piecol=cols)
dev.off()
