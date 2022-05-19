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
gene.exp[,1] <- paste0("chr", gene.exp[,1]);
cnv[,1]      <- paste0("chr", cnv[,1]);
pvalue[,1]   <- paste0("chr", pvalue[,1]);
TCGA.BC.fus[,1]   <- paste0("chr", TCGA.BC.fus[,1]);
TCGA.BC.fus[,4]   <- paste0("chr", TCGA.BC.fus[,4]);

pdf("OmicCircos4vignette10.pdf", 8,8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

zoom <- c("chr1", "chr22", 939245.5, 154143883, 0, 180);
circos(R=385, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.label.h", cex=1.3, zoom=zoom);
circos(R=350, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.scale", lwd=0.001, cex=0.3, zoom=zoom);
circos(R=345, cir=UCSC.hg18, W=4, mapping=UCSC.hg18.chr,  type="chr", zoom=zoom);
circos(R=240, cir=UCSC.hg18, W=100, mapping=gene.exp, col.v=4,  type="heatmap2", cluster=TRUE, col.bar=TRUE, lwd=0.01, zoom=zoom);
circos(R=180, cir=UCSC.hg18, W=60,  mapping=cnv,      col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0, zoom=zoom);
circos(R=120, cir=UCSC.hg18, W=60,  mapping=pvalue,   col.v=4,    type="l",   B=TRUE, lwd=1, col=cols[1], zoom=zoom);
circos(R=120, cir=UCSC.hg18, W=10,  mapping=TCGA.BC.fus, type="link", lwd=2, zoom=zoom);

## zoom in links by using the hightlight functions 
## highlight
the.col1=rainbow(10, alpha=0.5)[1];

highlight <- c(140, 400, "chr11", 282412.5, "chr11", 133770314.5, the.col1, the.col1);
circos(R=110, cir=UCSC.hg18, W=40, mapping=highlight, type="hl", lwd=2, zoom=zoom);
the.col2=rainbow(10, alpha=0.5)[6];
highlight <- c(140, 400, "chr17", 739525, "chr17", 78385909, the.col2, the.col2);
circos(R=110, cir=UCSC.hg18, W=40, mapping=highlight, type="hl", lwd=2, zoom=zoom);
## highlight link
highlight.link1 <- c(400, 400, 140, 376.8544, 384.0021, 450, 540.5);
circos(cir=UCSC.hg18, mapping=highlight.link1, type="highlight.link", col=the.col1, lwd=1);
highlight.link2 <- c(400, 400, 140, 419.1154, 423.3032, 543, 627);
circos(cir=UCSC.hg18, mapping=highlight.link2, type="highlight.link", col=the.col2, lwd=1);

## zoom in chromosome 11
zoom <- c("chr11", "chr11", 282412.5, 133770314.5, 180, 270);
circos(R=385, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.label.h", cex=1.3, zoom=zoom);
circos(R=350, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.scale", lwd=0.001, cex=0.3, zoom=zoom);
circos(R=345, cir=UCSC.hg18, W=4, mapping=UCSC.hg18.chr,  type="chr", zoom=zoom);
circos(R=240, cir=UCSC.hg18, W=100, mapping=gene.exp, col.v=4,  type="heatmap2", cluster=TRUE, lwd=0.01, zoom=zoom);
circos(R=180, cir=UCSC.hg18, W=60,  mapping=cnv,      col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0, zoom=zoom);
circos(R=120, cir=UCSC.hg18, W=60,  mapping=pvalue,   col.v=4,    type="l",   B=TRUE, lwd=1, col=cols[1], zoom=zoom);

## zoom in chromosome 17
zoom <- c("chr17", "chr17", 739525, 78385909, 274, 356);
circos(R=385, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.label.h", cex=1.3, zoom=zoom);
circos(R=350, cir=UCSC.hg18, mapping=UCSC.hg18.chr, type="chr.scale", lwd=0.001, cex=0.3, zoom=zoom);
circos(R=345, cir=UCSC.hg18, W=4, mapping=UCSC.hg18.chr,  type="chr", zoom=zoom);
circos(R=240, cir=UCSC.hg18, W=100, mapping=gene.exp, col.v=4,  type="heatmap2", cluster=TRUE, lwd=0.01, zoom=zoom);
circos(R=180, cir=UCSC.hg18, W=60,  mapping=cnv,      col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0, zoom=zoom);
circos(R=120, cir=UCSC.hg18, W=60,  mapping=pvalue,   col.v=4,    type="l",   B=TRUE, lwd=1, col=cols[1], zoom=zoom);

dev.off()



