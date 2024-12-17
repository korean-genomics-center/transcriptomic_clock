{
    suppressMessages(library("stringr"))
    suppressMessages(library("DESeq2"))
    suppressMessages(library("argparse"))
}
parse_arguments <- function(){
    parser <- ArgumentParser()
    parser$add_argument("--expData", type="character", required=TRUE)
    parser$add_argument("--metaData", type="character", required=TRUE)
    parser$add_argument("--fileout", type="character", required=TRUE)
    parser$add_argument("--design", type="character", required=TRUE, default="~1")
    parser$add_argument("--id_column", type="character", required=TRUE, default="Project.ID")
    parser$add_argument("--gene_column", type="character", required=TRUE, default="Gene_ID")
    parser$add_argument("--geoMean_column", type="character", required=TRUE, default="geoMeans")
    parser$add_argument("--file_geoMeans", type="character", required=TRUE)
    parser$add_argument("--file_sizeFactors", type="character", required=TRUE)
    inputargs <- parser$parse_args()
    return (inputargs)
}
{
    inputargs <- parse_arguments()
}
{
    metainfo <- read.csv(inputargs$metaData, sep="\t")
    metasamplenames <- metainfo[[inputargs$id_column]]
    counts <- read.csv(inputargs$expData, sep="\t")
    countsamplenamesraw <- names(counts)
    # DO NOT LET SAMPLE ID TO BE NUMERIC!!!
    countsamplenamesnew <- sapply(strsplit(countsamplenamesraw, "\\."), function(x) paste(x, collapse = "-"))
    names(counts) <- countsamplenamesnew
    countsfilt <- counts[, names(counts) %in% metasamplenames]
    roundcountsfilt <- as.data.frame(round(countsfilt))
    idcol <- as.data.frame(counts[[inputargs$gene_column]])
    names(idcol) <- inputargs$gene_column
    countData <- cbind(idcol, roundcountsfilt)
    rownames(countData) <- countData[[inputargs$gene_column]]
    countData <- countData[, -1]
}
{
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=metainfo, design=as.formula(inputargs$design))
}
{   
    geoMeans <- read.csv(inputargs$file_geoMeans, sep="\t")[[inputargs$geoMean_column]]
    ddssize <- estimateSizeFactors(dds, geoMeans=geoMeans)
    ddsdisp <- estimateDispersions(ddssize)
    ddsnorm <- nbinomWaldTest(ddsdisp, maxit=1000)
}
{
    normcount <- counts(ddsnorm, normalized=TRUE)
}
{
    write.table(as.data.frame(round(normcount)), file=inputargs$fileout, row.names=TRUE, sep="\t", quote=FALSE)
}
{
    write.table(as.data.frame(ddssize$sizeFactor), file=inputargs$file_sizeFactors, row.names=TRUE, sep="\t", quote=FALSE)     
}
