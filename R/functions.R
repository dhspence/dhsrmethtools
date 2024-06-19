# dhsrmethtools/R/functions.R

#' bed2bsseq
#'
#' @param file A simplfied bed file with methylation: chr\tstart\tend\tmethylation\coverage
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' bed2bsseq(file=system.file("extdata", "test.bed.gz", package = "dhsrmethtools"),samplename="test")
#' @export
bed2bsseq <- function(file,samplename=basename(file)){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }

  samplename <- gsub(".bed(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- fread(cmd = paste("gunzip -c ", file, " | awk -F'\t' '$1~/^chr[0-9XYM]+$/'"), header = FALSE, col.names=c("chr","start","end","meth","cov"))
  } else {
    # Use awk to filter regular file
    data <- fread(cmd = paste("awk -F'\t' '$1~/^chr[0-9XYM]+$/' ", file), header = FALSE, col.names=c("chr","start","end","meth","cov"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  data$meth <- data$meth * data$cov
  hdf5_M <- writeHDF5Array(as.matrix(data$meth))
  hdf5_C <- writeHDF5Array(as.matrix(data$cov))
  bsOut <- BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$V1,IRanges(data$V2,data$V3)),sampleNames=samplename)
  return(bsOut)
}

#' modkit2bsseq
#'
#' @param file An ONT bedmethyl file
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' modkit2bsseq(file=system.file("extdata", "test.basemods.bedmethyl.gz", package = "dhsrmethtools"),samplename="test")
#' @export
modkit2bsseq <- function(file,samplename=basename(file)){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }

  samplename <- gsub(".basemods.bedmethyl(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- fread(cmd = paste("gunzip -c ", file, "| awk -F'\t' '$4 == \"m\"'"), header = FALSE)
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
#' pbcpg2bsseq(file=system.file("extdata", "test.pbcpg.bedmethyl.gz", package = "dhsrmethtools"),samplename="test")
#' @export
pbcpg2bsseq <- function(file,samplename=basename(file)){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }

  samplename <- gsub(".basemods.bedmethyl(.gz)?|.bed|.bed.gz|.5mC.bed.gz", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- fread(cmd = paste("gunzip -c ", file), header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
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
