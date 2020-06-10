rm(list=ls());

library(OmicCircos);
options(stringsAsFactors = FALSE);

perm_list <- function (n, r, v = 1:n){ 
    if (r == 1) 
       X <- matrix(v, n, 1) 
    else if (n == 1) 
       X <- matrix(v, 1, r) 
    else { 
       X <- NULL 
       for (i in 1:n){ 
            X <- rbind(X, cbind(v[i], perm_list(n-1 , r-1 , v[-i]))) 
       } 
    } 
    return(X);
} 

## initial
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 10;

## select segments
seg.name <- paste("chr", 1:seg.num, sep="");

center   <- c(200, 400, 600);
center.i <- perm_list(3, 2);
for (i in 1:3){
  center.i <- rbind(center.i, rep(i,2));
}

colors   <- rainbow(seg.num, alpha=0.8);
color2   <- rainbow(10,      alpha=0.3);
pdffile  <- "OmicCircos4vignette9.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

for (i in 1:nrow(center.i)){
 xc <- center[center.i[i,1]];
 yc <- center[center.i[i,2]];

 sim.out   <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
                   link.pg=link.pg.num);
 db        <- segAnglePo(sim.out$seg.f, seg=seg.name);
 link.pg.v <- sim.out$seg.link.pg
 circos(xc=xc, yc=yc, R=90, type="chr", cir=db, col=colors, print.chr.lab=FALSE, W=4);
 cols <- sample(color2, nrow(link.pg.v), replace=TRUE);
 circos(xc=xc, yc=yc, R=86, cir=db, mapping=link.pg.v, type="link.pg", lwd=2, col=cols);

}
dev.off()

## detach(package:OmicCircos, unload=TRUE)

