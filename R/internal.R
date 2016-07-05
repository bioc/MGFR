# Function to test if a given gene is a marker
.isMarker.rnaseq <- function(named.vec, rep.vec) {
    sort.vec <- sort(named.vec, decreasing = TRUE)
    sv.len <- length(sort.vec)
    names.vec <- names(sort.vec)
    rep.fe <- unname(rep.vec[names.vec[1]])
    poss.marker <- (length(unique(names.vec[1:rep.fe])) == 1 && sum((as.numeric(sort.vec)[1:rep.fe]) > 
        1) == rep.fe)
    if (poss.marker) {
        mean.vec <- c()
        sort.num <- as.numeric(sort.vec)
        mean.vec <- c(mean.vec, mean(sort.num[1:rep.fe]))
        start.p <- rep.fe + 1
        rep.se <- unname(rep.vec[names.vec[start.p]])
        end.pos <- start.p + rep.se - 1
        cp.found <- FALSE
        while (end.pos <= sv.len && !cp.found) {
            len.sa <- length(unique(names.vec[start.p:end.pos]))
            if (sum(unique(names.vec[start.p:end.pos]) %in% names.vec[(end.pos + 1):sv.len]) == 
                0 | end.pos == sv.len) {
                mean.vec <- c(mean.vec, mean(sort.num[start.p:end.pos]))
                cp.found <- TRUE
            } else {
                end.pos <- start.p + sum(rep.vec[unique(names.vec[start.p:end.pos])]) - 1
            }
        }
        return(c(names.vec[1], (mean.vec[2]/mean.vec[1])))
    }
}

## Function to map gene identifiers (ensembl, refseq or ucsc ids) to gene symbols and entrez
## gene ids using biomaRt.
.get.genes.rnaseq <- function(IDs = "", gene.identifiers = "ensembl") {
    
    # if(!require('biomaRt')) { stop('Please install the R package biomaRt to map the gene
    # identifiers to gene symbols!') } require('biomaRt')
    if (gene.identifiers == "ensembl") {
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))
        ann.df <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol", 
            "entrezgene", "description"), values = IDs, mart = mart)
        return(ann.df)
    }
    if (gene.identifiers == "ucsc") {
        # Convert from UCSC IDs to Gene Name and Description
        ann.df <- getBM(attributes = c("ucsc", "hgnc_symbol", "entrezgene", "description"), 
            filters = "ucsc", values = IDs, mart = mart)  #, uniqueRows=T)
        return(ann.df)
    }
    if (gene.identifiers == "refseq") {
        # Convert RefSEq IDs
        ann.df <- getBM(attributes = c("refseq_mrna", "hgnc_symbol", "entrezgene", "description"), 
            filters = "refseq_mrna", values = IDs, mart = mart)  #, uniqueRows=T)
        return(ann.df)
    }
}


