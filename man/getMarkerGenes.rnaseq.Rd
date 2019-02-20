\name{getMarkerGenes.rnaseq}
\alias{getMarkerGenes.rnaseq}
\title{Marker Gene Detection

}
\description{
Function to detect marker genes using normalized RNA-seq data
}
\usage{
getMarkerGenes.rnaseq(data.mat, class.vec=colnames(data.mat), samples2compare="all", annotate=FALSE, gene.ids.type="ensembl", score.cutoff=1)
}

\arguments{
  \item{data.mat}{
RNA-seq gene expression matrix with genes corresponding to rows
and samples corresponding to columns.
}

 \item{class.vec}{
A character vector containing the classes of samples (columns) of data.mat in the same order as provided in the matrix.
}
  \item{samples2compare}{
A character vector with the sample names to be compared (e.g.
c("liver", "lung", "brain")). By default all samples in the reference matrix are used.
}
  \item{annotate}{
A boolean value. If TRUE the gene symbol and the entrez gene id are shown.
}

  \item{gene.ids.type}{
Type of the used gene identifiers, the following gene identifiers are supported: ensembl, refseq and ucsc gene ids. default is ensembl.
}
  \item{score.cutoff}{
A value in the interval [0,1] to filter the markers according to the  specificity score. The default value is 1 (no filtering).}

}
\details{
For each marker in the output list, the gene id and
the corresponding score are shown. If annotate is TRUE, the gene symbol and the entrez gene id are shown. The score is used to
rank the markers according to their specificity. A lower value means a higher specificity.
}
\value{
A list with marker genes associated with each sample type.
}

\author{
Khadija El Amrani <a.khadija@gmx.de>
}

\examples{
data(ref.mat)
res.list <- getMarkerGenes.rnaseq(ref.mat, class.vec = colnames(ref.mat), samples2compare="all", annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
names(res.list)
## show the first 20 markers of liver
res.list[["liver_markers"]][1:20]
}

\keyword{ marker genes }
\keyword{ RNA-seq data }
