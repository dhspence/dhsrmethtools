# dhsrmethtools/R/functions.R

# Import necessary functions
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom HDF5Array writeHDF5Array
NULL

#' bed2bsseq
#'
#' @param file A simplfied bed file with methylation: chr\\tstart\\tend\\tmethylation\\tcoverage
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
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }

  samplename <- gsub(".bed(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- data.table::fread(cmd = paste("gunzip -c ", file, " | awk -F'\t' '$1~/^chr[0-9XYM]+$/'"), header = FALSE, col.names=c("chr","start","end","meth","cov"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(cmd = paste("awk -F'\t' '$1~/^chr[0-9XYM]+$/' ", file), header = FALSE, col.names=c("chr","start","end","meth","cov"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  data$meth <- data$meth * data$cov
  hdf5_M <- HDF5Array::writeHDF5Array(as.matrix(data$meth))
  hdf5_C <- HDF5Array::writeHDF5Array(as.matrix(data$cov))
  bsOut <- bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$start,data$end)),sampleNames=samplename)
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
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }

  samplename <- gsub(".basemods.bedmethyl(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    # include column names to match the ONT bedmethyl file from modkit
    data <- data.table::fread(cmd = paste("gunzip -c ", file, "| awk -F'\t' '$4 == \"m\"'"), header = FALSE, col.names=c("chr","start","end","type","score","strand","start2","end2","color","valid_cov","fraction_mod","n_mod","n_canonical","n_othermod","n_deleted","n_fail","n_diff","n_nocall"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(cmd = paste("awk -F'\t' '$4 == \"m\"' ", file), header = FALSE, col.names=c("chr","start","end","type","score","strand","start2","end2","color","valid_cov","fraction_mod","n_mod","n_canonical","n_othermod","n_deleted","n_fail","n_diff","n_nocall"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  data$end <- data$end + 1
  hdf5_M <- HDF5Array::writeHDF5Array(as.matrix(data$n_mod))
  hdf5_C <- HDF5Array::writeHDF5Array(as.matrix(data$valid_cov))
  bsOut <- bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$start,data$end)),sampleNames=samplename)
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
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }

  samplename <- gsub(".basemods.bedmethyl(.gz)?|.bed|.bed.gz|.5mC.bed.gz", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- data.table::fread(cmd = paste("gunzip -c ", file), header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(file, header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data$end <- data$end + 1
  data <- data[chr %in% chromosomes]
  hdf5_M <- HDF5Array::writeHDF5Array(data.matrix(data$mod))
  hdf5_C <- HDF5Array::writeHDF5Array(data.matrix(data$coverage))
  bsOut <- bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$start,data$end)),sampleNames=samplename)
  return(bsOut)
}
