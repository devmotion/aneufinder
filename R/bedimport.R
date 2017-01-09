#' Import BED file into GRanges
#'
#' This is a simple convenience function to read a BED file into a \code{\link{GRanges}} object. The BED file is expected to contain at least the following fields \code{chromosome, start, end}.
#'
#' @param bedfile Filename of a BED file.
#' @param skip Number of lines to skip at the beginning.
#' @return A \code{\link{GRanges}} object with the contents of the BED file.
#' @author Aaron Taudt, David Widmann
#' @importFrom utils count.fields read.table
#' @export
#'
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Import the file and skip the first 10 lines
#'data <- importBed(bedfile, skip=10)
#'
importBed <- function(bedfile, skip=0) {
  ncols <- max(utils::count.fields(bedfile, skip=skip))

  if ( ncols < 3 )
    stop("File ", bedfile, " is not a correct BED file. Please specify at least the three fields 'chromosome', 'start', 'end'.")

  classes <- c("character", rep("numeric", 2), "character", "numeric", "character", rep("NULL", max(0, ncols-6)))[1:ncols]
  data <- utils::read.table(bedfile, colClasses=classes, skip=skip)
  data <- cbind(data, data.frame('.', 1000, '*')[,0:(6-ncol(data))])
  names(data) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

  ## adjust strand information
  data$strand <- sub("^[^+-]$", "*", data$strand)

  ## convert to GRanges object
  ranges <- GenomicRanges::GRanges(seqnames=data$chrom,
                                   ranges=IRanges(start=data$chromStart+1,     # Convert from 0-based half open to 1-based closed
                                                  end=data$chromEnd),
                                   name=data$name,
                                   mapq=data$score,
                                   strand=data$strand)

  return(ranges)
}
