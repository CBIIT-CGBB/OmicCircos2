rm(list=ls());
options(stringsAsFactors = FALSE);
library(OmicCircos2);

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

## select segments
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);

colors   <- rainbow(seg.num, alpha=0.5);

pdffile  <- "OmicCircos4vignette2.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(R=400, cir=db, W= 5, mapping=db, type="chr.label.h", col=colors, cex=1.5);
circos(R=365, cir=db, W= 5, mapping=seg.f, type="chr.scale", col=colors, cex=0.4);
circos(R=365, cir=db, W= 5, mapping=db, type="chr.noB", col=colors);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=8, type="box",   B=TRUE, col=colors[1], lwd=0.1, scale=TRUE, cex=0.2);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=8, type="hist",  B=TRUE, col=colors[3], lwd=0.1, scale=TRUE, cex=0.2);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=8, type="ms",  B=TRUE, col=colors[7], lwd=0.1, scale=TRUE, cex=0.2);
circos(R=200, cir=db, W=40, mapping=seg.v, col.v=3, type="h",  B=FALSE,  col=colors[2], lwd=0.1);
circos(R=160, cir=db, W=40, mapping=seg.v, col.v=3, type="s", B=TRUE, col=colors, lwd=0.1);
circos(R=120, cir=db, W=40, mapping=seg.v, col.v=3, type="b", B=FALSE, col=colors, lwd=0.1);
circos(R=115, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors);
circos(R=115, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=colors);

dev.off()

