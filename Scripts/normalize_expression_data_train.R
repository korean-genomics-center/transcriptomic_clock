{
    suppressMessages(library("stringr"))
    suppressMessages(library("DESeq2"))
}
{
    expData <- "../Data/GEO/GSE134080/GSE134080_count_preprocessed.txt"
    metaData <- "../Data/GEO/GSE134080/GSE134080_metadata_preprocessed.txt"
    fileout <- "../Data/GEO/GSE134080/GSE134080_count_preprocessed_normalized.txt"
    design <- ~1
    id_column <- "Project.ID"
    gene_column <- "Gene_ID"
    file_geoMeans <- "../Data/GEO/GSE134080/GSE134080_geoMeans.txt"
    file_sizeFactors <- "../Data/GEO/GSE134080/GSE134080_sizeFactors.txt"
}
{
    metainfo <- read.csv(metaData, sep="\t")
    metasamplenames <- metainfo[[id_column]]
    counts <- read.csv(expData, sep="\t")
    countsamplenamesraw <- names(counts)
    countsamplenamesnew <- sapply(strsplit(countsamplenamesraw, "\\."), function(x) paste(x, collapse = "-"))
    names(counts) <- countsamplenamesnew
    countsfilt <- counts[, names(counts) %in% metasamplenames]
    roundcountsfilt <- as.data.frame(round(countsfilt))
    idcol <- as.data.frame(counts[["Gene_ID"]])
    names(idcol) <- "Gene_ID"
    countData <- cbind(idcol, roundcountsfilt)
    rownames(countData) <- countData[[gene_column]]
    countData <- countData[, -1]
}
{
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=metainfo, design=as.formula(design))
}
{   
    ddssize <- estimateSizeFactors(dds)
    ddsdisp <- estimateDispersions(ddssize)
    ddsnorm <- nbinomWaldTest(ddsdisp, maxit=1000)
}
{
    normcount <- counts(ddsnorm, normalized=TRUE)
    rawcount <- counts(ddsnorm, normalized=FALSE)
    geoMeans <- exp(rowMeans(log(rawcount + 1)))
}
{
    write.table(as.data.frame(round(normcount)), file=paste0(fileout), row.names=TRUE, sep="\t", quote=FALSE)
}
{
    write.table(as.data.frame(geoMeans), file=paste0(file_geoMeans), row.names=TRUE, sep="\t", quote=FALSE)
}
{
    write.table(as.data.frame(ddssize$sizeFactor), file=paste0(file_sizeFactors), row.names=TRUE, sep="\t", quote=FALSE)     
}
