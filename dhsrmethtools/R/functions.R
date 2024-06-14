# dhsrmethtools/R/functions.R

#' modkit2bsseq
#'
#' @param file An ONT bedmethyl file
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' modkit2bsseq(file="test.bedmethyl.gz",samplename="test")
#' @export
modkit2bsseq <- function(file,samplename=basename(file)){
  require(data.table)
  require(HDF5Array)
  require(bsseq)
  samplename <- gsub(".basemods.bedmethyl(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- fread(cmd = paste("zcat", file, "| awk -F'\t' '$4 == \"m\"'"), header = FALSE)
  } else {
    # Use awk to filter regular file
    data <- fread(cmd = paste("awk -F'\t' '$4 == \"m\"' ", file), header = FALSE)
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[V1 %in% chromosomes]
  data$V3 <- data$V3 + 1
  hdf5_M <- writeHDF5Array(as.matrix(data$V12))
  hdf5_C <- writeHDF5Array(as.matrix(data$V10))
  bsOut <- BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$V1,IRanges(data$V2,data$V3)),sampleNames=samplename)
  return(bsOut)
}

#' pbcpg2bsseq
#'
#' @param file A pb-CpG tools bedmethyl file
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' pbcpg2bsseq(file="test.bedmethyl.gz",samplename="test")
#' @export
pbcpg2bsseq <- function(file,samplename=basename(file)){
  require(data.table)
  require(HDF5Array)
  require(bsseq)
  samplename <- gsub(".basemods.bedmethyl(.gz)?|.bed|.bed.gz|.5mC.bed.gz", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- fread(cmd = paste("zcat", file), header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
  } else {
    # Use awk to filter regular file
    data <- fread(file, header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  dat <- read.table(file,sep="\t",)
  dat$end <- dat$end + 1
  dat <- dat[which(dat$chr %in% paste0("chr",c(1:22,"X","Y"))),]
  hdf5_M <- writeHDF5Array(data.matrix(dat[,"mod"]))
  hdf5_C <- writeHDF5Array(data.matrix(dat[,"coverage"]))
  bsOut <- BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(dat$chr,IRanges(dat$start,dat$end)),sampleNames=samplename)
  return(bsOut)
}
