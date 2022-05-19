rm(list=ls());
options(stringsAsFactors = FALSE);
library(OmicCircos3);

set.seed(1234);

## initial
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 10;

sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
  link.pg=link.pg.num);

seg.f     <- sim.out$seg.frame;
seg.v     <- sim.out$seg.mapping;
link.v    <- sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num   <- length(unique(seg.f[,1]));

## 
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);

colors   <- rainbow(seg.num, alpha=0.5);
pdffile  <- "OmicCircos4vignette3.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, cir=db, W= 5, mapping=db, type="chr.label.h", col=colors, cex=1.5);
circos(R=365, cir=db, W= 5, mapping=seg.f, type="chr.scale", col=colors, cex=0.4);
circos(R=365, cir=db, W= 5, mapping=db, type="chr.noB", col=colors);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=8, type="quant90", B=FALSE, col=colors[c(1,7,9)], lwd=0.1, scale=TRUE, cex=0.4);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=3, type="sv", B=TRUE, col=colors[7], lwd=0.1, scale=TRUE, cex=0.4);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=3, type="ss", B=FALSE, col=colors[3], lwd=0.1, scale=TRUE, cex=0.4);
circos(R=200, cir=db, W=40, mapping=seg.v, col.v=8, type="heatmap", lwd=3);
circos(R=160, cir=db, W=40, mapping=seg.v, col.v=3, type="s.sd", B=FALSE, col=colors[4], lwd=0.1);
circos(R=120, cir=db, W=40, mapping=seg.v, col.v=3, type="ci95", B=TRUE, col=colors[c(1,7,9)], lwd=0.1);
circos(R=115, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors);
circos(R=115, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=colors);

the.col1=rainbow(10, alpha=0.3)[3];
highlight <- c(120, 360, "chr6", 2, "chr6", 10, the.col1, the.col1);
circos(R=110, cir=db, W=40, mapping=highlight, type="hl", lwd=2);

the.col1=rainbow(10, alpha=0.01)[3];
the.col2=rainbow(10, alpha=0.8)[1];
highlight <- c(120, 360, "chr6", 12, "chr7", 10, the.col1, the.col2);
circos(R=110, cir=db, W=40, mapping=highlight, type="hl", lwd=2);

dev.off()
