#' Import reads into GRanges
#'
#' Import aligned reads from a BAM or bed(.gz) file into a \code{\link{GRanges}} object.
#'
#' @param file Sorted BAM file or BED file with aligned reads.
#' @param bamindex BAM index file. Only necessary when importing BAM files. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @param assembly Please see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies. Only necessary when importing BED files. BAM files are handled automatically. Alternatively a data.frame with columns 'chromosome' and 'length'.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your BAM files (not implemented for BED files).
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments.
#' @param blacklist A \code{\link{GRanges}} or a BED file with blacklisted regions. Reads falling into those regions will be discarded.
#' @param what A character vector of fields that are returned (not implemented for BED files). Type \code{\link[Rsamtools]{scanBamWhat}} to see what is available.
#' @return A \code{\link{GRanges}} object containing the reads.
#' @export
#'
importReads <- function(file, bamindex=file, assembly, chromosomes=NULL, pairedEndReads=FALSE, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000, blacklist=NULL, what='mapq') {
  ## Determine format
  format <- tolower(sub('\\(.bam|bed)((?<=\\.bed)\\.gz)?', '\\1', file, ignore.case=TRUE, perl=TRUE))

  ## Determine chromosome length information
  chrom.lengths <- chromsomeLengths(bamfile=if(format == 'bam') file else NA, assembly=assembly, chromosomes=chromosomes)

  ## Import reads
  if (format == 'bam') {
    reads <- bam2GRanges(file, bamindex=bamindex, chromosomes=names(chrom.lengths), pairedEndReads=pairedEndReads, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist, what=what)
  } else if (format == 'bed') {
    reads <- bed2GRanges(file, chromosomes=names(chrom.lengths), remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
  }

  return(reads)
}

#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param bamfile A sorted BAM file.
#' @param bamindex BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your BAM files (not implemented for BED files).
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.
#' @param blacklist A \code{\link{GRanges}} or a bed(.gz) file with blacklisted regions. Reads falling into those regions will be discarded.
#' @param what A character vector of fields that are returned. Type \code{\link[Rsamtools]{scanBamWhat}} to see what is available.
#' @return A \code{\link{GRanges}} object containing the reads.
#' @author Aaron Taudt, David Widmann
#' @importFrom Rsamtools indexBam BamFile ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first
#' @importFrom S4Vectors queryHits
#'
#'@examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata", "BB150803_IV_074.bam", package="AneuFinderData")
#'## Read the file into a GRanges object
#'reads <- bam2GRanges(bamfile, chromosomes=c(1:19,'X','Y'), pairedEndReads=FALSE,
#'                     min.mapq=10, remove.duplicate.reads=TRUE)
#'print(reads)
#'
bam2GRanges <- function(bamfile, bamindex=bamfile, chromosomes=NULL, pairedEndReads=FALSE, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000, blacklist=NULL, what='mapq') {
  ## Check if bamindex exists and create it (if necessary)
  bamindex <- createBamIndex(bamfile, bamindex=bamindex)

  ## Read chromosome lengths
  chroms2use <- chromosomeLengths(bamfile, chromosomes)

  ## Import the file into GRanges
  ptm <- startTimedMessage("Reading file ",basename(bamfile)," ...")

  ## Define import options
  gr <- GenomicRanges::GRanges(seqnames=chroms2use$chromosome, ranges=IRanges(start=1, end=chroms2use$length))
  gr <- blacklistGRanges(gr, blacklist) # Exclude blacklisted regions
  isDuplicate <- if(remove.duplicate.reads) FALSE else NA
  mapqFilter <- if(is.null(min.mapq)) NA_integer_ else as.integer(min.mapq)

  if (pairedEndReads) {
    data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter=mapqFilter, flag=scanBamFlag(isDuplicate=isDuplicate)))
  } else {
    data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter=mapqFilter, flag=scanBamFlag(isDuplicate=isDuplicate)))
  }

  stopTimedMessage(ptm)

  ## Convert to GRanges
  ptm <- startTimedMessage("Converting to GRanges ...")
  data <- as(data.raw, 'GRanges')
  stopTimedMessage(ptm)

  ## Filter out too long fragments
  ptm <- startTimedMessage("Filtering reads ...")
  data <- data[width(data)<=max.fragment.width]
  stopTimedMessage(ptm)

  ## Check if any reads are imported
  if (length(data) == 0)
    stop("No reads imported. Try with 'pairedEndReads=FALSE', 'min.mapq=NULL', and 'blacklist=NULL' to keep all reads.")

  return(data)
}


#' Import BED file into GRanges
#'
#' Import aligned reads from a BED file into a \code{\link{GRanges}} object.
#'
#' @param bedfile A file with aligned reads in BED format.
#' @param assembly Please see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies. Alternatively a data.frame with columns 'chromosome' and 'length'.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments.
#' @param blacklist A \code{\link{GRanges}} or a BED file with blacklisted regions. Reads falling into those regions will be discarded.
#' @return A \code{\link{GRanges}} object containing the reads.
#' @author Aaron Taudt, David Widmann
#' @importFrom utils read.table
#' @importFrom S4Vectors queryHits
#'
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Read the file into a GRanges object
#'reads <- bed2GRanges(bedfile, assembly='mm10', chromosomes=c(1:19,'X','Y'),
#'                     min.mapq=10, remove.duplicate.reads=TRUE)
#'print(reads)
#'
bed2GRanges <- function(bedfile, assembly, chromosomes=NULL, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000, blacklist=NULL) {
  ## Import BED file
  data <- importBed(bedfile, skip=skip)

  ## Read chromosome length information
  chroms2use <- chromosomeLengths(assembly, chromosomes)

  ## Select only desired chromosomes
  ptm <- startTimedMessage("Subsetting chromosomes ...")
  seqlevelsStyle(data) <- seqlevelsStyle(chroms2use$chromosome)
  data <- keepSeqlevels(data, chroms2use$chromosome)
  seqlengths(data) <- chroms2use$seqlengths
  data <- cleanSeqlevels(data)
  stopTimedMessage(ptm)

  ptm <- startTimedMessage("Filtering reads ...")
  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(data)$mapq))) {
      warning(paste0(bedfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
    }
    data <- data[mcols(data)$mapq >= min.mapq]
  }
  ## Filter out too long fragments
  data <- data[width(data)<=max.fragment.width]
  stopTimedMessage(ptm)

  ## Exclude reads falling into blacklisted regions
  blacklistGRanges(data, blacklist)

  ## Check if any reads are imported
  if (length(data) == 0)
    stop("No reads imported. Try with 'min.mapq=NULL' and 'blacklist=NULL' to keep all reads.")

  return(data)
}

#' Remove blacklisted regions from GRanges object.
#'
#' This is a simple convenience function to remove blacklisted regions, given either as \code{\link{GRanges}} object or BED file, from a \code{\link{GRanges}} object. The BED file is expected to contain at least the three fields \code{chromosome, start, end}.
#'
#' @param data A \code{\link{GRanges}} object.
#' @param blacklist A \code{\link{GRanges}} object or the filename of a BED file which contains the blacklisted regions.
#' @return A \code{\link{GRanges}} object with the regions not overlapping any regions of the blacklist.
#' @author Aaron Taudt, David Widmann
#' @importFrom utils write.table
#'
blacklistGRanges <- function(data, blacklist=NULL) {
  # Input checks
  if (is.null(blacklist))
    return(data)

  if ( !is.character(blacklist) & class(blacklist)!='GRanges' )
    stop("'blacklist' has to be either a BED file or a GRanges object")

  ptm <- startTimedMessage("Filtering blacklisted regions ...")

  # Import blacklist file as GRanges object
  if (is.character(blacklist))
    blacklist <- importBed(blacklist, skip=0)

  # Convert 'blacklist' to the chromosome format of 'data'
  seqlevelsStyle(blacklist) <- GenomeInfoDb::seqlevelsStyle(data)

  # Remove overlaps between 'data' and 'blacklist'
  overlaps <- findOverlaps(data, blacklist)
  idx <- setdiff(1:length(data), S4Vectors::queryHits(overlaps))
  data <- data[idx]

  stopTimedMessage(ptm)

  return(data)
}

#' Get chromosome length information from BAM file or assembly.
#'
#' This is a simple convenience function to obtain chromosome length information from either a BAM file or an assembly.
#'
#' @param bamfile A BAM file.
#' @param assembly Please see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies. Not used if \code{\link{bamfile}} specified. Alternatively a data.frame or the filename of a table with columns either 'chromosome' and 'length' or 'UCSC_seqlevel' and 'UCSC_seqlength'.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @return A \code{\link{data.frame}} with a column 'chromosome' of chromosome names and a column 'length' of corresponding lengths.
#' @author Aaron Taudt, David Widmann
#'
chromosomeLengths <- function(bamfile=NA, assembly, chromosomes=NULL) {
  ## Input check
  if (is.na(bamfile) & is.null(assembly))
    stop("Please specify a 'bamfile' or an 'assembly'")

  ## Define helper function to obtain chromosome length information from data.frame
  extractChromosomeLengths  <- function(df) {
    if (all(c('UCSC_seqlevel', 'UCSC_seqlength') %in% names(df))) {
      chrom.lengths <- df[,'UCSC_seqlength']
      names(chrom.lengths) <- df[,'UCSC_seqlevel']
    } else if (all(c('chromosome', 'length') %in% names(df))) {
      chrom.lengths  <- df[,'length']
      names(chrom.lengths)  <- df[,'chromosome']
    } else {
      stop("Can't read chromosome length information. Table has to contain either columns 'UCSC_seqlevel' and 'UCSC_seqlength' or columns 'chromosome' and 'length'.")
    }

    return(chrom.lengths)
  }

  if (!is.na(bamfile)) {
    ## Read chromosome length information from BAM file
    ptm <- startTimedMessage("Obtaining chromosome length information from file ", bamfile, " ...")
    chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
    stopTimedMessage(ptm)
  } else {
    if (is.character(assembly)) {
      ## Read chromosome length information from assembly
      if (file.exists(assembly)) {
        ptm <- startTimedMessage("Obtaining chromosome length information from file ", assembly, " ...")
        df <- utils::read.table(assembly, sep='\t', header=TRUE)
        stopTimedMessage(ptm)
      } else {
        ptm <- startTimedMessage("Obtaining chromosome length information from UCSC ...")
        df <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(assembly)[,c('UCSC_seqlevel','UCSC_seqlength')]
        stopTimedMessage(ptm)
      }
    } else if (is.data.frame(assembly)) {
      df <- assembly
    } else {
      stop("'assembly' must be either a data.frame or a character specifying the assembly.")
    }

    chrom.lengths <- extractChromosomeLengths(df)
  }

  chrom.names <- names(chrom.lengths)

  ## Check if chromosomes are unique
  if (anyDuplicated(chrom.names))
    stop("Names of chromosomes must be unique.")

  chroms2use <- filterChromosomes(chrom.names, chromosomes=chromosomes)

  return(data.frame(chromosome=chroms2use, length=chrom.lengths[chroms2use]))
}

filterChromosomes <- function(data, chromosomes=NULL) {
  if (is.null(chromosomes))
    return(data)

  ## Input check
  if ( !is.character(data) & class(data)!='GRanges' )
    stop("'data' must be chromosome names or a GRanges object.")

  if ( is.character(data) ) {
    chroms.names <- data
  } else {
    chroms.names <- seqlevels(data)
  }

  seqlevelsStyle(chroms.names) <- seqlevelsStyle(chromosomes) # Convert to same chromosome style
  chroms2use <- intersect(chromosomes, chrom.names)

  ## Stop if non of the specified chromosomes exist
  if (length(chroms2use)==0)
    stop('The specified chromosomes ', toString(chromosomes), ' do not exist in the data.')

  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chrom.names)
  if (length(diff)>0)
    warning(paste0('Not using chromosomes ', toString(diff), ' because they are not in the data.'))

  return(chroms2use)
}

createBamIndex <- function(bamfile, bamindex=bamfile) {
  ## Sanitize bamindex filename
  bamindex <- sub('^(.*)(\\.bai)?$', '\\1\\.bai', bamindex, ignore.case=TRUE)
  ## Check if bamindex exists
  if (!file.exists(bamindex)) {
    ptm <- startTimedMessage("Making bam-index file ...")
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
    stopTimedMessage(ptm)
  }

  return(bamindex)
}

cleanSeqlevels <- function(data) {
  ## Input check
  if (class(data)!='GRanges' & !(is.data.frame(data) & c("chromosome", "length") %in% names(data)))
    stop("Please specify a GRanges object or a data frame with columns 'chromosome' and 'length'.")

  ## Get chromosome length information
  chrom.lengths <- seqlengths(data)
  chrom.names <- names(chrom.lengths)

  ## Filter chromosomes without length information
  na.seqlevels <- chrom.names[is.na(chrom.lengths) | is.na(chrom.names)]
  keep.seqlevels <- chrom.names[!is.na(chrom.lengths) & !is.na(chrom.names)]
  data <- keepSeqlevels(data, keep.seqlevels)

  if (length(na.seqlevels) > 0)
    warning("Dropped seqlevels because no length information was available: ", toString(na.seqlevels))

  return(data)
}

