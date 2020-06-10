%% sim circos
\name{sim.circos}
\alias{sim.circos}
\title{
 circular data simulation
}
\description{
 This function generates data for user to test the circos functions
}
\usage{
  sim.circos(seg=10, po=c(20,50), ind=10, link=10, link.pg=10);
}
\arguments{
   \item{seg}{integer, the segment number. The default is 10.}
   \item{po}{vector, the segment positions. The default is c(20:50)}
   \item{ind}{integer, the number of samples. The default is 10.}
   \item{link}{integer, the number of links. The default is 10.}
   \item{link.pg}{integer, the number of link ploygons. The default is 10.}
}
\value{
%% 
 sim.circos returns a list containing at least the following components: 
 \item{seg.frame}{data.frame, segment data}
 \item{seg.mapping}{data.frame, mapping data}
 \item{seg.link}{data.fame, link data}
 \item{seg.link.pg}{data.frame, link polygon data}
}
\examples{

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

colors   <- rainbow(seg.num, alpha=0.5);
pdffile  <- "OmicCircos4vignette3.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, type="chr", cir=db, col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=360, cir=db, W=40, mapping=seg.v, col.v=8, type="quant90", B=FALSE, col=colors, lwd=2, scale=TRUE);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=3, type="sv", B=TRUE, col=colors[7],  scale=TRUE);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=3, type="ss", B=FALSE, col=colors[3],  scale=TRUE);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=8, type="heatmap", lwd=3);
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

dev.off()

}
\author{
 Ying Hu <yhu@mail.nih.gov>
 Chunhua Yan <yanch@mail.nih.gov>
}
\references{
  OmicCircos: an R package for simple and circular visualization of omics data. Cancer Inform. 2014 Jan 16;13:13-20. doi: 10.4137/CIN.S13495. eCollection 2014. PMID: 24526832 [PubMed] PMCID: PMC3921174
}
