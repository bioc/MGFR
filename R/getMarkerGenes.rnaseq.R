## Function to get marker genes parameters: data.mat: RNA-Seq gene expression matrix with
## genes corresponding to rows and samples corresponding to columns.  samples2compare: a
## character vector with the sample names to be compared, default is to compare all annotate:
## a boolean value indicating if an annotation (mapping gene identifiers to gene symbols) is
## desired, default is TRUE gene.ids.type: type of the used gene identifiers, the following
## gene identifiers are supported: ensembl, refseq and ucsc gene ids score.cutoff: an integer
## value to filter the resulted markers output: a list with the corresponding markers to each
## sample type

getMarkerGenes.rnaseq <- function(data.mat, samples2compare = "all", annotate = TRUE, gene.ids.type = "ensembl", 
    score.cutoff = 1) {
    
    if (length(samples2compare) > 1) {
        if (any(!samples2compare %in% colnames(data.mat))) {
            stop("The samples to compare should be included in the data matrix!!")
        } else {
            ii <- which(colnames(data.mat) %in% samples2compare)
            data.mat <- data.mat[, ii]
        }
    }
    ## remove all 0 rows
    x <- data.mat
    data.mat <- x[!apply(x == 0, 1, all), , drop = FALSE]
    markers.list <- list()
    rep.vec <- table(colnames(data.mat))
    cat("Detecting marker genes...\n")
    res.list <- apply(data.mat, 1, .isMarker.rnaseq, rep.vec = rep.vec)
    res.list[sapply(res.list, is.null)] <- NULL
    mar.len <- length(unlist(res.list))
    samples.vec <- unname(unlist(res.list))[seq(1, mar.len, 2)]
    scores.vec <- round(as.numeric(unname(unlist(res.list))[seq(2, mar.len, 2)]), digits = 4)
    markers.vec <- names(res.list)
    names(scores.vec) <- markers.vec
    u.snames <- unique(colnames(data.mat))
    
    if (annotate) {
        if (!gene.ids.type %in% c("ensembl", "refseq", "ucsc")) {
            stop("To map gene identifiers to gene symbols, only the following gene identifiers are supported: ensembl, refseq and ucsc gene ids!")
        }
        cat("Mapping gene IDs to gene symbols...\n")
        annot.df <- .get.genes.rnaseq(IDs = rownames(data.mat), gene.identifiers = gene.ids.type)
        for (i in seq_along(u.snames)) {
            inds <- which(samples.vec == u.snames[i])
            sort.scores <- sort(scores.vec[inds])
            sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
            genes.inds <- match(names(sort.scores), annot.df[, "ensembl_gene_id"])
            markers.list[[i]] <- paste(names(sort.scores), annot.df[genes.inds, "hgnc_symbol"], 
                annot.df[genes.inds, "entrezgene"], unname(sort.scores), sep = " : ")
        }
        
    } else {
        
        
        for (i in seq_along(u.snames)) {
            inds <- which(samples.vec == u.snames[i])
            sort.scores <- sort(scores.vec[inds])
            sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
            markers.list[[i]] <- paste(names(sort.scores), unname(sort.scores), sep = " : ")
        }
    }
    names(markers.list) <- paste(u.snames, "markers", sep = "_")
    cat("Done! \n")
    return(markers.list)
}
