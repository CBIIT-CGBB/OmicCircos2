rm(list=ls());

library(OmicCircos3);
options(stringsAsFactors = FALSE);

set.seed(1234);

## load mm cytogenetic band data
data("UCSC.mm10");
data("UCSC.mm10.chr");
ref     <- UCSC.mm10.chr;
ref[,1] <- gsub("chr", "", ref[,1]);
## initial values for simulation data 
colors <- rainbow(10, alpha=0.8);
lab.n  <- 50;
cnv.n  <- 200;
arc.n  <- 100;
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

## make fusion data
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
cnv.d[,1] <- paste0("chr", cnv.d[,1]);
arc.d[,1] <- paste0("chr", arc.d[,1]);
fus.d[,1] <- paste0("chr", fus.d[,1]);
fus.d[,4] <- paste0("chr", fus.d[,4]);

pdffile  <- "OmicCircos4vignette8.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="");

circos(R=385, cir=UCSC.mm10, mapping=UCSC.mm10.chr, type="chr.label.h", cex=1.4);
circos(R=350, cir=UCSC.mm10, mapping=UCSC.mm10.chr, type="chr.scale", lwd=0.001, cex=0.2);
circos(R=345, cir=UCSC.mm10, mapping=UCSC.mm10.chr, type="chr2", W=10)

circos(R=300, cir=UCSC.mm10, W=40, mapping=cnv.d, type="b3", B=T, col=colors[7]);
circos(R=300, cir=UCSC.mm10, W=40, mapping=cnv.d, type="s2", B=F, col=colors[1], cex=0.5);

circos(R=260, cir=UCSC.mm10, W=40, mapping=arc.d, type="arc2", B=F, col=colors, lwd=10);
circos(R=220, cir=UCSC.mm10, W=40, mapping=cnv.d, col.v=3, type="b2", B=T, cutoff=-0.2, col=colors[c(7,9)], lwd=2);
circos(R=180, cir=UCSC.mm10, W=40, mapping=arc.d, col.v=4, type="arc",  B=F, col=colors[c(1,7)], lwd=4, scale=T, cex=0.4);
circos(R=150, cir=UCSC.mm10, W=10, mapping=fus.d, type="link",  lwd=2, col=colors[c(1,7,9)]);

dev.off();
