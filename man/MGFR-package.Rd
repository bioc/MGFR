\name{MGFR-package}
\alias{MGFR-package}
\alias{MGFR}
\docType{package}
\title{
Marker Gene Finder in RNA-seq data
}
\description{
The package is designed to detect marker genes from RNA-seq data
}
\details{
\tabular{ll}{
Package: \tab MGFR\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr 
License: \tab GPL-3 \cr
}

}
\author{
Khadija El Amrani
Maintainer: Khadija El Amrani <khadija.el-amrani@charite.de>
}

\keyword{ package }

\examples{
data(ref.mat)
res.list <- getMarkerGenes.rnaseq(ref.mat, samples2compare="all", annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
names(res.list)
## show the first 20 markers of liver
res.list[["liver_markers"]][1:20]
}
