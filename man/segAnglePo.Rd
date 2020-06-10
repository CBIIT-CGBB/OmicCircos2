%% segAnglePo
\name{segAnglePo}
\alias{segAnglePo}

\title{
 generate circular skeleton data from user's mapping data
}
\description{
 This function creates a data frame and converts the segment 
 pointer positions (linear coordinates) into angle values (the angle based 
 coordinates along circumference). In the data frame, column 1 is unique 
 segment or chromosome names; column 2 is the start angle; column 3 is the end angle; 
 column 4 is the accumulative start position; column 5 is the accumulative 
 end position; column 6 is the start position and column 7 is the end position for each segment or chromosome.
}

\usage{
  segAnglePo(seg.dat=seg.dat, seg=seg, angle.start=angle.start, angle.end=angle.end);
}
\arguments{
   \item{seg.dat}{the segment data should be a matrix or a data frame: column 1 is the segment name or chromosome name;
     column 2 is the segment start; column 3 is the segment end; column 4 is segment name2 (optional);
     and column 5 is segment description (optional).}
   \item{seg}{vector: segment names (optional)}
   \item{angle.start}{numeric: plot start angle, angle.start=0 by default (optional)}
   \item{angle.end}{numeric: plot end angle, angle.end=360 by default (optional)}
}

\author{
 Ying Hu <yhu@mail.nih.gov>
 Chunhua Yan <yanch@mail.nih.gov>
}
\references{
OmicCircos: an R package for simple and circular visualization of omics data. Cancer Inform. 2014 Jan 16;13:13-20. doi: 10.4137/CIN.S13495. eCollection 2014. PMID: 24526832 [PubMed] PMCID: PMC3921174
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
pdffile  <- "OmicCircos4vignette2.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
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

dev.off()

}

