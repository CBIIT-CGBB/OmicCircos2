rm(list=ls());
options(stringsAsFactors = FALSE);
library(OmicCircos3);

data(UCSC.hg18);
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

cols <- rainbow(10, alpha=0.5);
cnv[,1]      <- paste0("chr", cnv[,1]);
gene.exp[,1] <- paste0("chr", gene.exp[,1]);
pvalue[,1]   <- paste0("chr", pvalue[,1]);
TCGA.BC.fus[,1]   <- paste0("chr", TCGA.BC.fus[,1]);
TCGA.BC.fus[,4]   <- paste0("chr", TCGA.BC.fus[,4]);

pdf("OmicCircos4vignette5.pdf", 8,8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=385, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.label.h", cex=1.3);
circos(R=350, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.scale", lwd=0.001, cex=0.3);
circos(R=345, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr2", W=10)
circos(R=240, cir=UCSC.hg18, W=100, mapping=gene.exp,  col.v=4,  type="heatmap2", cluster=TRUE, col.bar=TRUE, lwd=0.01);
circos(R=180, cir=UCSC.hg18, W=60,  mapping=cnv,   col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0);
circos(R=120, cir=UCSC.hg18, W=60,  mapping=pvalue,  col.v=4,    type="l",   B=TRUE, lwd=1, col=cols[1]);
circos(R=120, cir=UCSC.hg18, W=10,  mapping=TCGA.BC.fus, type="link", lwd=2);

dev.off()



