## Function to get marker genes and show the results in html pages parameters: data.mat:
## RNA-Seq gene expression matrix with genes corresponding to rows and samples corresponding
## to columns.  samples2compare: a character vector with the sample names to be compared,
## default is to compare all gene.ids.type: type of the used gene identifiers, only ensembl
## gene ids are supported score.cutoff: an integer value to filter the resulted markers,
## default is 1 (no filtering) directory: path to a directory to save the html files, default
## is the current working directory output: HTML tables containing marker genes for each
## sample type, as well as links to various online annotation sources (Ensembl, GenBank and
## EntrezGene repositories)
getMarkerGenes.rnaseq.html <- function(data.mat, samples2compare = "all", gene.ids.type = "ensembl", 
    score.cutoff = 1, directory = getwd()) {
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
    
    
    if (!gene.ids.type == "ensembl") {
        stop("Only ensembl ids are supported!")
    }
    cat("Mapping gene IDs to gene symbols...\n")
    annot.df <- .get.genes.rnaseq(IDs = rownames(data.mat), gene.identifiers = gene.ids.type)
    cat("Creating html pages...\n")
    for (i in seq_along(u.snames)) {
        inds <- which(samples.vec == u.snames[i])
        sort.scores <- sort(scores.vec[inds])
        sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
        
        genes.inds <- match(names(sort.scores), annot.df[, "ensembl_gene_id"])
        
        res.df <- as.data.frame(cbind(Gene_ID = names(sort.scores), Gene_Symbol = annot.df[genes.inds, 
            "hgnc_symbol"], Entrez_Id = annot.df[genes.inds, "entrezgene"], Description = annot.df[genes.inds, 
            "description"], Score = as.numeric(sort.scores)))
        
        suppressWarnings(htmlpage(res.df[, c(1, 2, 3)], paste(directory, paste(u.snames[i], 
            "html", sep = "."), sep = "/"), " ", table.head = c("Gene ID", "Gene Symbol", "Entrez Id", 
            "Description", "Score"), othernames = res.df[, c(4, 5)], repository = list("ens", 
            "gb", "en"), species = "Homo_sapiens"))
        
    }
    cat("The html files are stored in the following directory: ", directory, "\n")
    cat("Done!")
    
}
