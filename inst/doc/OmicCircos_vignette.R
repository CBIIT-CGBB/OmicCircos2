### R code from vignette source 'OmicCircos_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: OmicCircos_vignette.Rnw:132-137
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos); 
## input hg19 cytogenetic band data
data(UCSC.hg19.chr);
head(UCSC.hg19.chr);


###################################################
### code chunk number 2: OmicCircos_vignette.Rnw:156-161
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos); 
## TCGA gene expression data
data(TCGA.BC.gene.exp.2k.60);
head(TCGA.BC.gene.exp.2k.60[,c(1:5)]);


###################################################
### code chunk number 3: OmicCircos_vignette.Rnw:179-184
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos); 
## TCGA fusion gene data
data(TCGA.BC.fus);
head(TCGA.BC.fus[,c(1:6)]);


###################################################
### code chunk number 4: OmicCircos_vignette.Rnw:208-216
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 10;
sim.out     <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, link.pg=link.pg.num);


###################################################
### code chunk number 5: OmicCircos_vignette.Rnw:237-238
###################################################
head(sim.out$seg.frame[,c(1:3)])


###################################################
### code chunk number 6: OmicCircos_vignette.Rnw:246-247
###################################################
head(sim.out$seg.mapping[,c(1:5)])


###################################################
### code chunk number 7: OmicCircos_vignette.Rnw:255-256
###################################################
head(sim.out$seg.link)


###################################################
### code chunk number 8: OmicCircos_vignette.Rnw:264-265
###################################################
head(sim.out$seg.link.pg)


###################################################
### code chunk number 9: OmicCircos_vignette.Rnw:296-319
###################################################
library(OmicCircos);
options(stringsAsFactors = FALSE);
set.seed(1234);
## initial values for simulation data 
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 4;
## output simulation data
sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
  link.pg=link.pg.num);
seg.f     <- sim.out$seg.frame;
seg.v     <- sim.out$seg.mapping;
link.v    <- sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num   <- length(unique(seg.f[,1]));
## select segments
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);
db[,2]   <- round(as.numeric(db[,2]), 3);
db[,3]   <- round(as.numeric(db[,3]), 3);
db


###################################################
### code chunk number 10: OmicCircos_vignette.Rnw:372-390
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);
set.seed(1234);
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 4;
sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
  link.pg=link.pg.num);
seg.f     <- sim.out$seg.frame;
seg.v     <- sim.out$seg.mapping;
link.v    <- sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num   <- length(unique(seg.f[,1]));
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);
colors   <- rainbow(seg.num, alpha=0.5);


###################################################
### code chunk number 11: OmicCircos4vignette1
###################################################
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(R=400, cir=db, type="chr",  col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=360, cir=db, W=40, mapping=seg.v, col.v=3, type="l",   B=TRUE, col=colors[1], lwd=2, scale=TRUE);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=3, type="ls",  B=FALSE, col=colors[9], lwd=2, scale=TRUE);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=3, type="lh",  B=TRUE, col=colors[7], lwd=2, scale=TRUE);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=19, type="ml",  B=FALSE, col=colors, lwd=2, scale=TRUE);
circos(R=200, cir=db, W=40, mapping=seg.v, col.v=19, type="ml2", B=TRUE, col=colors, lwd=2);
circos(R=160, cir=db, W=40, mapping=seg.v, col.v=19, type="ml3", B=FALSE, cutoff=5, lwd=2);
circos(R=150, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors[c(1,7)]);
circos(R=150, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(colors,link.pg.num));


###################################################
### code chunk number 12: OmicCircos_vignette.Rnw:471-496
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);
set.seed(1234);

## initial values for simulation data 
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 4;
## output simulation data
sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
  link.pg=link.pg.num);

seg.f     <- sim.out$seg.frame;
seg.v     <- sim.out$seg.mapping;
link.v    <- sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num   <- length(unique(seg.f[,1]));

## select segments
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);

colors   <- rainbow(seg.num, alpha=0.5);


###################################################
### code chunk number 13: OmicCircos4vignette2
###################################################
par(mar=c(0, 0, 0, 0));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, type="chr", cir=db, col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=360, cir=db, W=40, mapping=seg.v, col.v=8, type="box",   B=TRUE, col=colors[1], lwd=0.1, scale=TRUE);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=8, type="hist",  B=TRUE, col=colors[3], lwd=0.1, scale=TRUE);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=8, type="ms",  B=TRUE, col=colors[7], lwd=0.1, scale=TRUE);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=3, type="h",  B=FALSE,  col=colors[2], lwd=0.1);
circos(R=200, cir=db, W=40, mapping=seg.v, col.v=3, type="s", B=TRUE, col=colors, lwd=0.1);
circos(R=160, cir=db, W=40, mapping=seg.v, col.v=3, type="b", B=FALSE, col=colors, lwd=0.1);
circos(R=150, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors[c(1,7)]);
circos(R=150, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(colors,link.pg.num));


###################################################
### code chunk number 14: OmicCircos_vignette.Rnw:572-592
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);
set.seed(1234);
## initial values for simulation data 
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 4;
## output simulation data
sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
  link.pg=link.pg.num);
seg.f     <- sim.out$seg.frame;
seg.v     <- sim.out$seg.mapping;
link.v    <- sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num   <- length(unique(seg.f[,1]));
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);
colors   <- rainbow(seg.num, alpha=0.5);


###################################################
### code chunk number 15: OmicCircos4vignette3
###################################################
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, type="chr", cir=db, col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=360, cir=db, W=40, mapping=seg.v, col.v=8, type="quant90", B=FALSE, col=colors, lwd=2, scale=TRUE);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=3, type="sv", B=TRUE, col=colors[7],  scale=TRUE);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=3, type="ss", B=FALSE, col=colors[3],  scale=TRUE);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=20, type="heatmap", lwd=3);
circos(R=200, cir=db, W=40, mapping=seg.v, col.v=3, type="s.sd", B=FALSE, col=colors[4]);
circos(R=160, cir=db, W=40, mapping=seg.v, col.v=3, type="ci95", B=TRUE, col=colors[4], lwd=2);
circos(R=150, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors[c(1,7)]);
circos(R=150, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(colors,link.pg.num));

the.col1=rainbow(10, alpha=0.5)[3];
highlight <- c(160, 410, 6, 2, 6, 10, the.col1, the.col1);
circos(R=110, cir=db, W=40, mapping=highlight, type="hl", lwd=1);

the.col1=rainbow(10, alpha=0.1)[3];
the.col2=rainbow(10, alpha=0.5)[1];
highlight <- c(160, 410, 3, 12, 3, 20, the.col1, the.col2);
circos(R=110, cir=db, W=40, mapping=highlight, type="hl", lwd=2);



###################################################
### code chunk number 16: OmicCircos_vignette.Rnw:712-758
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);
set.seed(1234);

data("UCSC.mm10.chr", package="OmicCircos");
ref     <- UCSC.mm10.chr;
ref[,1] <- gsub("chr", "", ref[,1]);

colors <- rainbow(10, alpha=0.8);
lab.n  <- 50;
cnv.n  <- 200;
arc.n  <- 30;
fus.n  <- 10;

## make arc data
arc.d <- c();
for (i in 1:arc.n){
  chr     <- sample(1:19, 1);
  chr.i   <- which(ref[,1]==chr);
  chr.arc <- ref[chr.i,];
  arc.i   <- sample(1:nrow(chr.arc), 2);
  arc.d   <- rbind(arc.d, c(chr.arc[arc.i[1],c(1,2)], chr.arc[arc.i[2],c(2,4)]));
}
colnames(arc.d) <- c("chr", "start", "end", "value");

## make fusion
fus.d <- c();
for (i in 1:fus.n){
  chr1    <- sample(1:19, 1);
  chr2    <- sample(1:19, 1);
  chr1.i  <- which(ref[,1]==chr1);
  chr2.i  <- which(ref[,1]==chr2);
  chr1.f  <- ref[chr1.i,];
  chr2.f  <- ref[chr2.i,];
  fus1.i  <- sample(1:nrow(chr1.f), 1);
  fus2.i  <- sample(1:nrow(chr2.f), 1);
  n1      <- paste0("geneA", i);
  n2      <- paste0("geneB", i);
  fus.d   <- rbind(fus.d, c(chr1.f[fus1.i,c(1,2)], n1, chr2.f[fus2.i,c(1,2)], n2));
}
colnames(fus.d) <- c("chr1", "po1", "gene1", "chr2", "po2", "gene2");

cnv.i <- sample(1:nrow(ref), cnv.n);
vale  <- rnorm(cnv.n);
cnv.d <- data.frame(ref[cnv.i,c(1,2)], value=vale);



###################################################
### code chunk number 17: OmicCircos4vignette4
###################################################
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");

circos(R=400, type="chr", cir="mm10", print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=340, cir="mm10", W=60, mapping=cnv.d, type="b3", B=TRUE, col=colors[7]);
circos(R=340, cir="mm10", W=60, mapping=cnv.d, type="s2", B=FALSE, col=colors[1], cex=0.5);
circos(R=280, cir="mm10", W=60, mapping=arc.d, type="arc2", B=FALSE, col=colors, lwd=10, cutoff=0);
circos(R=220, cir="mm10", W=60, mapping=cnv.d, col.v=3, type="b2", B=TRUE, cutoff=-0.2, col=colors[c(7,9)], lwd=2);
circos(R=160, cir="mm10", W=60, mapping=arc.d, col.v=4, type="arc",  B=FALSE, col=colors[c(1,7)], lwd=4, scale=TRUE);
circos(R=150, cir="mm10", W=10, mapping=fus.d, type="link",  lwd=2, col=colors[c(1,7,9)]);



###################################################
### code chunk number 18: OmicCircos_vignette.Rnw:836-862
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);

data("TCGA.PAM50_genefu_hg18");
data("TCGA.BC.fus");
data("TCGA.BC.cnv.2k.60");
data("TCGA.BC.gene.exp.2k.60");
data("TCGA.BC.sample60");
data("TCGA.BC_Her2_cnv_exp");

pvalue <- -1 * log10(TCGA.BC_Her2_cnv_exp[,5]);
pvalue <- cbind(TCGA.BC_Her2_cnv_exp[,c(1:3)], pvalue);

Her2.i <- which(TCGA.BC.sample60[,2] == "Her2");
Her2.n <- TCGA.BC.sample60[Her2.i,1];

Her2.j <- which(colnames(TCGA.BC.cnv.2k.60) %in% Her2.n);
cnv    <- TCGA.BC.cnv.2k.60[,c(1:3,Her2.j)]; 
cnv.m  <- cnv[,c(4:ncol(cnv))];
cnv.m[cnv.m >  2] <- 2;
cnv.m[cnv.m < -2] <- -2;
cnv <- cbind(cnv[,1:3], cnv.m);

Her2.j   <- which(colnames(TCGA.BC.gene.exp.2k.60) %in% Her2.n);
gene.exp <- TCGA.BC.gene.exp.2k.60[,c(1:3,Her2.j)]; 
colors <- rainbow(10, alpha=0.5);


###################################################
### code chunk number 19: OmicCircos4vignette5
###################################################
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");

circos(R=300, type="chr", cir="hg18", print.chr.lab=FALSE, W=4);
circos(R=310, cir="hg18", W=20, mapping=TCGA.PAM50_genefu_hg18, type="label", 
       side="out", col=c("black", "blue","red"), cex=0.4);
circos(R=250, cir="hg18", W=50, mapping=cnv, col.v=11, type="ml3", B=FALSE, col=colors[7], cutoff=0, scale=TRUE);
circos(R=200, cir="hg18", W=50, mapping=gene.exp, col.v=11, type="ml3", B=TRUE, col=colors[3], cutoff=0, scale=TRUE);
circos(R=140, cir="hg18", W=50, mapping=pvalue, col.v=4, type="l", B=FALSE, col=colors[1], scale=TRUE);
## set fusion gene colors
cols  <- rep(colors[7], nrow(TCGA.BC.fus));
col.i <- which(TCGA.BC.fus[,1]==TCGA.BC.fus[,4]);
cols[col.i] <- colors[1];
circos(R=132, cir="hg18", W=50, mapping=TCGA.BC.fus, type="link", col=cols, lwd=2);


###################################################
### code chunk number 20: OmicCircos4vignette6
###################################################
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");
circos(R=300, type="chr", cir="hg18", col=TRUE, print.chr.lab=FALSE, W=4);
circos(R=290, cir="hg18", W=20, mapping=TCGA.PAM50_genefu_hg18, type="label", side="in", col=c("black", "blue"), cex=0.4);
circos(R=310, cir="hg18", W=50, mapping=cnv, col.v=11, type="ml3", B=TRUE, col=colors[7], cutoff=0, scale=TRUE);
circos(R=150, cir="hg18", W=50, mapping=gene.exp, col.v=11, type="ml3", B=TRUE, col=colors[3], cutoff=0, scale=TRUE);
circos(R=90,  cir="hg18", W=50, mapping=pvalue, col.v=4, type="l", B=FALSE, col=colors[1], scale=TRUE);
circos(R=82, cir="hg18", W=50, mapping=TCGA.BC.fus, type="link", col=cols, lwd=2);


###################################################
### code chunk number 21: OmicCircos_vignette.Rnw:977-1004
###################################################
options(stringsAsFactors = FALSE);
library(OmicCircos);

data("TCGA.PAM50_genefu_hg18");
data("TCGA.BC.fus");
data("TCGA.BC.cnv.2k.60");
data("TCGA.BC.gene.exp.2k.60");
data("TCGA.BC.sample60");
data("TCGA.BC_Her2_cnv_exp");

pvalue <- -1 * log10(TCGA.BC_Her2_cnv_exp[,5]);
pvalue <- cbind(TCGA.BC_Her2_cnv_exp[,c(1:3)], pvalue);

Her2.i <- which(TCGA.BC.sample60[,2] == "Her2");
Her2.n <- TCGA.BC.sample60[Her2.i,1];

Her2.j <- which(colnames(TCGA.BC.cnv.2k.60) %in% Her2.n);
cnv    <- TCGA.BC.cnv.2k.60[,c(1:3,Her2.j)]; 
cnv.m  <- cnv[,c(4:ncol(cnv))];
cnv.m[cnv.m >  2] <- 2;
cnv.m[cnv.m < -2] <- -2;
cnv <- cbind(cnv[,1:3], cnv.m);

Her2.j   <- which(colnames(TCGA.BC.gene.exp.2k.60) %in% Her2.n);
gene.exp <- TCGA.BC.gene.exp.2k.60[,c(1:3,Her2.j)]; 

colors <- rainbow(10, alpha=0.5);


###################################################
### code chunk number 22: OmicCircos4vignette7
###################################################
par(mar=c(0, 0, 0, 0));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, cir="hg18", W=4,   type="chr", print.chr.lab=TRUE, scale=TRUE);
circos(R=300, cir="hg18", W=100, mapping=gene.exp,  col.v=8,  type="heatmap2", 
       cluster=TRUE, col.bar=TRUE, lwd=0.1, col="blue");
circos(R=220, cir="hg18", W=80,  mapping=cnv,   col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0);
circos(R=140, cir="hg18", W=80,  mapping=pvalue,  col.v=4,    type="l",   B=TRUE, lwd=1, col=colors[1]);

cols        <- rep(colors[7], nrow(TCGA.BC.fus));
col.i       <- which(TCGA.BC.fus[,1]==TCGA.BC.fus[,4]);
cols[col.i] <- colors[1];
circos(R=130, cir="hg18", W=10,  mapping=TCGA.BC.fus, type="link2", lwd=2, col=cols);



