\name{getMarkerGenes.rnaseq.html}
\alias{getMarkerGenes.rnaseq.html}
\title{Marker Gene Detection

}
\description{
Function to detect marker genes using normalized RNA-seq data and show the marker genes in HTML tables with links to various online annotation sources (Ensembl, GenBank and EntrezGene repositories)
}
\usage{
getMarkerGenes.rnaseq.html(data.mat, samples2compare="all", gene.ids.type="ensembl",score.cutoff=1, directory = getwd())
}

\arguments{
  \item{data.mat}{
RNA-seq gene expression matrix with genes corresponding to rows
and samples corresponding to columns.
}
  \item{samples2compare}{
A character vector with the sample names to be compared (e.g.
c("liver", "lung", "brain")). By default all samples in the reference matrix are used.
}
 
  \item{gene.ids.type}{
Type of the used gene identifiers, the following gene identifiers are supported: ensembl, refseq and ucsc gene ids. default is ensembl.
}
  \item{score.cutoff}{
A value in the interval [0,1] to filter the markers according to the  specificity score. The default value is 1 (no filtering).}

  \item{directory}{
Path to the directory where to save the html pages, default is the current working directory.}
}

\details{
  This function is based on the function \link[annotate]{htmlpage} from the R-package
'annotate'.
}
\value{
This function is used only for the side effect of creating HTML tables.
}

\author{
Khadija El Amrani <khadija.el-amrani@charite.de>
}

\examples{
data(ref.mat)
getMarkerGenes.rnaseq.html(ref.mat, samples2compare="all", gene.ids.type="ensembl", score.cutoff=1,directory = getwd())
}
\keyword{ marker genes }
\keyword{ RNA-seq data }
