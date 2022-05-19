rm(list=ls());
options(stringsAsFactors = FALSE);
library(OmicCircos2);

data(UCSC.hg18);
data("TCGA.BC.fus");
data("TCGA.BC.cnv.2k.60");
data("TCGA.BC.gene.exp.2k.60");
data("TCGA.BC.sample60");

## gene expression data for PCA
exp.m <- TCGA.BC.gene.exp.2k.60[,c(4:ncol(TCGA.BC.gene.exp.2k.60))];

cnv   <- TCGA.BC.cnv.2k.60;
cnv[,1]      <- paste0("chr", cnv[,1]);

## PCA
type.n  <- unique(TCGA.BC.sample60[,2]);
colors  <- rainbow(length(type.n), alpha=0.5);

pca.col <- rep(NA, nrow(TCGA.BC.sample60));
for (i in 1:length(type.n)){
  n   <- type.n[i];
  n.i <- which(TCGA.BC.sample60[,2] == n);
  n.n <- TCGA.BC.sample60[n.i,1];
  g.i <- which(colnames(exp.m) %in% n.n);
  pca.col[g.i] <- colors[i];
}

exp.m   <- na.omit(exp.m);
pca.out <- prcomp(t(exp.m), scale = TRUE);

## subtype cnv
cnv.i <- c();
for (i in 1:length(type.n)){
  n     <- type.n[i];
  n.i   <- which(TCGA.BC.sample60[,2] == n);
  n.n   <- TCGA.BC.sample60[n.i,1];
  cnv.i <- which(colnames(cnv) %in% n.n);
}

## main
pdf("OmicCircos4vignette7.pdf", 8,8);
par(mar=c(5, 5, 5, 5));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

legend(680,800, c("Basal","Her2","LumA","LumB"), pch=19, col=colors[c(2,4,1,3)], cex=0.7, 
     title ="Gene Expression (PCA)", box.col="white");

legend(5,800, c("1 Basal", "2 Her2", "3 LumA", "4 LumB", "(center)"), cex=0.7, 
     title ="CNV (OmicCircos)", box.col="white");

circos(R=385, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.label.h", cex=1.3);
circos(R=350, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.scale", lwd=0.001, cex=0.3);
circos(R=345, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr2", W=10)
R.v <- 280;
for (i in 1:length(type.n)){
  n     <- type.n[i];
  n.i   <- which(TCGA.BC.sample60[,2] == n);
  n.n   <- TCGA.BC.sample60[n.i,1];
  cnv.i <- which(colnames(cnv) %in% n.n);
  cnv.v <- cnv[,cnv.i];
  cnv.v[cnv.v > 2]  <- 2;
  cnv.v[cnv.v < -2] <- -2;
  cnv.m <- cbind(cnv[,c(1:3)], cnv.v);
  circos(xc=400, yc=400, R=R.v, cir=UCSC.hg18, W=40, mapping=cnv.m, col.v=4,  type="ml3", B=FALSE, lwd=1, cutoff=0, scale=TRUE, cex=0.4);
  R.v <- R.v - 40;
}

points(pca.out$x[,1]*8+400, pca.out$x[,2]*8+400, pch=19, col=pca.col, cex=2);

dev.off() 



