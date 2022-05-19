%% circos
\name{circos}
\alias{circos}

\title{
 draw circular
}
\description{
 This is the main function of OmicCircos to draw circular plots.
}
\usage{
  circos(mapping=mapping, xc=400, yc=400, R=400, W=W,
         cir="", type="n", col.v=3, B=FALSE, print.chr.lab=TRUE,
         col.bar=FALSE, col.bar.po = "topleft", cluster=FALSE, order=NULL,
         scale=FALSE, cutoff = "n", zoom="", col=rainbow(10, alpha=0.8)[7], side="", ...)
}

\arguments{
	\item{mapping}{data frame or matrix containing mapping information
  and values. Column 1: the segment or chromosome ID; column 2: the position;
  column 3: the position label (optional) or the value and additional columns are the values.
such as gene expression and copy number.  Missing values are allowed and will be ignored.}

   \item{xc}{integer, the circle center x coordinate}
   \item{yc}{integer, the circle center y coordinate}
   \item{R}{integer, the circle radius}
   \item{W}{integer, the circle width}
   \item{cir}{genome reference name (hg19, mm10 ...) or data frame from
   segAnglePo function or data frame from user's mapping data}.
   \item{type}{the type is one of
      \itemize{
         \item "arc": arcs with variable radii
         \item "arc2": arcs with the fixed radius
         \item "b": bar charts
         \item "b2": bar charts (opposite side of cutoff value)
         \item "b3": bar charts with the same height
         \item "box": box plots
         \item "chr": plots of banded chromosome
         \item "chr.label.h": horizontal label of chromosome or segment
         \item "chr.label.v": vertical label of chromosome or segment
         \item "chr.noB": plots of chromosome or segment without banding 
         \item "chr.scale": draw scale of chromosome or segment 
         \item "ci95": 95\% confidence interval lines
         \item "h": histograms
         \item "heatmap": heatmaps
         \item "heatmap2": heatmaps with genomic coordinates
         \item "hightlight.link": link lines for zoom
         \item "hist": polygons for multiple samples
         \item "hl": highlight
         \item "l": lines
         \item "label": gene labels or text annotations
         \item "label2": gene labels or text annotations with the same circumference coordinate
         \item "lh": horizontal lines
         \item "link.pg": link polygons based on Bezier curve
         \item "link": link lines based on Bezier curve
         \item "link2": link lines with smaller intra-chromosome arcs
         \item "ls": lines in stair steps
         \item "ml": multiple lines (for more than 1 samples)
         \item "ml2": multiple horizontal lines
         \item "ml3": multiple lines in stair steps
         \item "ms":  multiple points
         \item "quant75": 75\% quantile lines
         \item "quant90": 90\% quantile lines
         \item "s": dots
         \item "s2": dots with the fixed radius
         \item "s.sd": dots proportional to standard deviation
         \item "ss": dot sizes proportional to the values
         \item "sv": dot sizes proportional to the variances
       }
   }
   \item{col.v}{column number. The column value will be displayed. If type=heatmap, the
                number is as the first column.}
   \item{B}{logical: draw background?}
   \item{print.chr.lab}{logical: draw chromosomes or segment labels?}
   \item{col.bar}{logical: draw col.bar? It is for type=heatmap.}
   \item{col.bar.po}{draw col.bar position, e.g. topleft, bottomright.}
   \item{cluster}{logical: cluster and draw Dendrogram at left coner? It is for type=heatmap only.}
   \item{order}{vector: chromosome or segment order}
   \item{scale}{logical: draw scale?}
   \item{cutoff}{numeric: for multiple samples}
   \item{zoom}{vector containing six values: start chromosome, end chromosome, start position, end
              position, start angle and end angle}
   \item{col}{character or vector: color names}
   \item{side}{character (in or out): for type=label(2) only}
   \item{...}{...}
}


\author{
 Ying Hu <yhu@mail.nih.gov>
 Chunhua Yan <yanch@mail.nih.gov>
}
\references{
  OmicCircos: an R package for simple and circular visualization of omics data. Cancer Inform. 2014 Jan 16;13:13-20. doi: 10.4137/CIN.S13495. eCollection 2014. PMID: 24526832 [PubMed] PMCID: PMC3921174
}
