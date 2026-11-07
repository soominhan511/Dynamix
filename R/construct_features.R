# Script : construct_features.R 
#   This script contains functions to extract genomic features :
#   GC content and mappability scores around restriction enzyme cutsites.
#   Adapted and modified for the Dynamix package pipeline.
#
# Original Source / Repository:
#   - HiCDCPlus : https://bitbucket.org/leslielab/hicdcplus/src/master/R/construct_features.R
#
# Reference :
#   Sahin, M. et al. (2021)
#   HiC-DC+ enables systematic 3D interaction calls and differential analysis for Hi-C and HiChIP. Nat Commun 12, 3366.

#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data

get_chrs <- function(gen = "Hsapiens", gen_ver = "hg19") {
    genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
    library(genome, character.only = TRUE)
    chrs <- GenomeInfoDb::seqnames(get(gen))[vapply(GenomeInfoDb::seqnames(get(gen)), nchar, FUN.VALUE = 1) <= 5]
    chrs <- chrs[!chrs %in% c("chrMT","chrY","chrM")]
    return(chrs)
}


get_enzyme_cutsites <- function(sig, gen = "Hsapiens", gen_ver = "hg19", chrs = NULL) {
    genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
    library(genome, character.only = TRUE)
    if (is.null(chrs)) {
        # get list of chromosomes
        chrs <- get_chrs(gen, gen_ver)
        chrs <- chrs[!chrs %in% c("chrY", "chrM")]
    }
    genome.chromSizes <- get_chr_sizes(gen, gen_ver, chrs)
    while (any(grepl("N", sig))) {
        if (length(sig) <= 1) {
            sig <- unique(c(sub("N", "A", sig), sub("N", "C", sig), sub("N", "G", sig), sub("N", "T", sig)))
        } else {
            rng <- length(sig)
            sig_new <- NULL
            for (i in seq(rng)) {
                if (grepl("N", sig[i])) {
                sig_new_i <- unique(c(sub("N", "A", sig[i]), sub("N", "C", sig[i]), sub("N", "G", sig[i]), sub("N", "T", sig[i])))
                sig_new <- c(sig_new, sig_new_i)
                } else {
                sig_new <- unique(c(sig_new, sig[i]))
                }
            }
            sig <- unique(sig_new)
        }
    }
    if (length(sig) <= 1) {
        enzymeCuts <- lapply(chrs, function(x) {
            range_chr <- IRanges::ranges(Biostrings::matchPattern(sig, get(gen)[[x]]))
            if (length(IRanges::start(range_chr)) <= 0) {
                warning(paste0("chromosome ", x, " does not have any cutsites with
                pattern ", sig, " .Using full chromosome instead"))
                range_chr <- IRanges::IRanges(start = 1, end = genome.chromSizes[x])
            }
            sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
        })
    } else {
        patternobj <- Biostrings::DNAStringSet(sig)
        enzymeCuts <- lapply(chrs, function(x) {
            range_chr <- IRanges::ranges(Biostrings::extractAllMatches(get(gen)[[x]], Biostrings::matchPDict(patternobj, get(gen)[[x]])))
            if (length(IRanges::start(range_chr)) <= 0) {
                warning(paste0("chromosome ", x, " does not have any cutsites with
                pattern set ", paste(sig, collapse = ","), " .Using full chromosome instead"))
                range_chr <- IRanges::IRanges(start = 1, end = genome.chromSizes[x])
            }
            sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
        })
    }
    enzymeCuts <- suppressWarnings(sort(
        unlist(methods::as(enzymeCuts, "GRangesList"))))
    return(enzymeCuts)
}

get_chr_sizes <- function(gen = "Hsapiens", gen_ver = "hg19", chrs = NULL) {
    genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
    library(genome, character.only = TRUE)
    if (is.null(chrs)) {
        # get list of chromosomes
        chrs <- get_chrs(gen, gen_ver)
    }
    sizes <- GenomeInfoDb::seqlengths(get(gen))[chrs]
    return(sizes)
}


construct_features_chr <- function(chrom, gen = "Hsapiens", gen_ver = "hg19",
                               sig = "GATC", bin_type = "Bins-uniform", binsize = 5000, 
                               wg_file = NULL, feature_type = "RE-based") {
    
    getGC <- function(regions) {
        genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
        library(genome, character.only = TRUE)
        gc <- rowSums(Biostrings::alphabetFrequency(BSgenome::Views(BiocGenerics::get(genome), regions), as.prob = TRUE)[, 
                                                                                                                         seq(2,3,1)])
        return(gc)
    }
    
    getMap <- function(regions, data) {
        hits <- GenomicRanges::findOverlaps(data, regions, type = "within")
        DT <- data.frame(queryHits = as.data.frame(hits)[[1]], subjectHits = as.data.frame(hits)[[2]])
        DT <- DT %>% dplyr::mutate(map = data$score[.data$queryHits], len = GenomicRanges::end(data)[.data$queryHits] - GenomicRanges::start(data)[.data$queryHits] + 
                                       1) %>% dplyr::group_by(.data$subjectHits) %>% dplyr::summarize(avgmap = stats::weighted.mean(map, w = len))
        map <- rep(0, length(regions))
        map[DT$subjectHits] <- DT$avgmap
        return(map)
    }
    
    # get chromosome size
    genome.chromSize <- get_chr_sizes(gen, gen_ver, chrom)[chrom]
    genome.chromGR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = 1, end = genome.chromSize))
    
    if (feature_type=="RE-based"){
    msg<-paste0("Using ", paste(chrom, collapse = " "), "and cut patterns ", paste(sig, collapse = " "))
    message(msg)
    
    
    enzymeCutsites <- get_enzyme_cutsites(sig, gen, gen_ver, chrom)
    RE_sites <- enzymeCutsites
    newends <- dplyr::lead(BiocGenerics::start(RE_sites)) - 1
    newends <- c(newends[-length(newends)], genome.chromSize)
    BiocGenerics::end(RE_sites) <- newends
    minstart <- min(BiocGenerics::start(RE_sites))
    if (minstart > 1) {
        RE_sites <- GenomeInfoDb::sortSeqlevels(c(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = 1, 
                                                                                                                       end = minstart)), RE_sites))
        RE_sites <- sort(RE_sites)
    }
    } else {
        if (bin_type!="Bins-uniform") stop("feature_type=RE-agnostic requires a fixed binsize.")
        msg<-paste0("Using ", paste(chrom, collapse = " "), " RE agnostic")
        message(msg)
    }
    
    
    if (!(is.null(wg_file))) {
        wgdata <- rtracklayer::import.bw(wg_file, which = genome.chromGR, as = "GRanges")
    } else {
        wgdata <- NULL
    }
    
    if (bin_type == "Bins-RE-sites") {
        endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
                                                                                                                    width = 200))
        BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
        BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)

        endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(end = GenomicRanges::end(RE_sites), 
                                                                                                                    width = 200))

        BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
        BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
        gcL <- getGC(regions = endL)
        gcR <- getGC(regions = endR)
        gc <- (gcL + gcR)/2
        RE_sites$gc <- gc
        
        if (!is.null(wgdata)) {
            endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
                                                                                                                        width = 500))
                BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
                BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)
            endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), IRanges::IRanges(end = GenomicRanges::end(RE_sites), 
                                                                                                                        width = 500))
                BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
                BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
            mapL <- getMap(regions = endL, data = wgdata)
            mapR <- getMap(regions = endR, data = wgdata)
            map <- (mapL + mapR)/2
            RE_sites$map <- map
        } else {
            RE_sites$map <- 0
        }
        # generate bintolen from RE_sites
        bintolen <- as.data.frame(RE_sites, stringsAsFactors = FALSE)
        bintolen <- bintolen %>% dplyr::mutate(RE_id = dplyr::row_number())
        if (binsize > 1) {
            # We cluster enzyme cutting sites every 'binsize' consecutive RE sites per each chromosome
            bintolen <- bintolen %>% dplyr::mutate(binNumb = floor((dplyr::row_number() - 1)/binsize) + 
                                                                                     1)
            bintolen <- bintolen %>% dplyr::group_by(.data$seqnames, .data$binNumb) %>% dplyr::summarize(start = min(start), 
                                                                                                         end = max(end), gc = mean(gc), map = mean(map), RE_id = min(.data$RE_id)) %>% dplyr::mutate(width = .data$end - 
                                                                                                                                                                                                         .data$start + 1) %>% dplyr::rename(chr = "seqnames") %>% dplyr::mutate(bins = paste(.data$chr, .data$start, 
                                                                                                                                                                                                                                                                                             .data$end, sep = "-"))
        } else {
            bintolen <- bintolen %>% dplyr::rename(chr = "seqnames") %>% dplyr::mutate(bins = paste(.data$chr, .data$start, 
                                                                                                    .data$end, sep = "-"))
        }
        # construct a similar bintolen object
        bintolen <- bintolen %>% dplyr::select(.data$bins, .data$gc, .data$map, .data$width, .data$RE_id)
        bintolen <- as.data.frame(bintolen, stringsAsFactors = FALSE)
        rownames(bintolen) <- bintolen$bins
        return(bintolen)
    }
    
    if (bin_type == "Bins-uniform") {
        bins.chrom <- function(chrom, binsize) {
            cuts_start <- seq(1, genome.chromSize, by = binsize)
            cuts_end <- seq(binsize, genome.chromSize - (genome.chromSize%%binsize) + binsize, by = binsize)
            cuts_end <- c(cuts_end[-length(cuts_end)], genome.chromSize)
            bins <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = cuts_start, end = cuts_end))
            return(bins)
        }
        binsGR <- bins.chrom(chrom, binsize)
        names(binsGR) <- paste(as.character(GenomicRanges::seqnames(binsGR)), GenomicRanges::start(binsGR), GenomicRanges::end(binsGR), 
                               sep = "-")
        
        #Read mappability data if exists
        if (!is.null(wgdata)) {
            wgdata <- wgdata[GenomicRanges::seqnames(wgdata) == chrom]
        } else {
            wgdata <- NULL
        }        
        
        if (feature_type=="RE-based"){
        medians <- (GenomicRanges::start(enzymeCutsites) + GenomicRanges::end(enzymeCutsites))/2

        # 500bp from ends for mappability and len features
        FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::restrict(IRanges::IRanges(end = medians - 
                                                                                                                               1, width = 500), start = 1, end = genome.chromSize))
        FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::restrict(IRanges::IRanges(start = medians, 
                                                                                                                           width = 500), start = 1, end = genome.chromSize))
        ends <- c(FragmentendsL, FragmentendsR)
        hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
        LR <- data.frame(bins = names(binsGR[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), end = GenomicRanges::end(ends[hits[[1]]]), 
                        stringsAsFactors = FALSE)
        } else {
        LR <- data.frame(bins=names(binsGR), start=GenomicRanges::start(binsGR), end=GenomicRanges::end(binsGR),stringsAsFactors = FALSE)
        }
        LRgr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR$start, end = LR$end))
        if (!is.null(wgdata)) {
            LR$map <- getMap(LRgr, wgdata)
        } else {
            LR$map <- 0
        }
        map <- LR %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(map = mean(map))
        if (feature_type=="RE-based"){
        # effective length
        ir <- IRanges::IRanges(LR$start, LR$end)
        LR$group <- as.data.frame(GenomicRanges::findOverlaps(ir, IRanges::reduce(ir)))[[2]]
        len <- LR %>% dplyr::group_by(.data$group, .data$bins) %>% dplyr::summarize(start = min(start), end = max(end)) %>% 
            dplyr::mutate(width = .data$end - .data$start + 1) %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(len = sum(.data$width))
        } else {
        len <- NULL
        }
        
        if (feature_type=="RE-based"){    
        # 200bp from ends for gc feature
        FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(end = medians - 1, 
                                                                                                        width = 200))
        FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = medians, 
                                                                                                         width = 200))
        ends <- c(FragmentendsL, FragmentendsR)
        hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
        LR2 <- data.frame(bins = names(binsGR[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), end = GenomicRanges::end(ends[hits[[1]]]), 
                              stringsAsFactors = FALSE)
        } else {
        LR2<-LR%>%dplyr::select(.data$bins,.data$start,.data$end)
        }
        LR2gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR2$start, end = LR2$end))
        LR2$gc <- getGC(LR2gr)
        gc <- LR2 %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(gc = mean(gc))
        if (!is.null(len)){
        LR <- dplyr::left_join(dplyr::left_join(gc, map), len)
        } else {
        LR <- dplyr::left_join(gc, map)
        }
        bintolen <- as.data.frame(LR)
        if (sum(bintolen$map)==0){
            bintolen<-bintolen%>%dplyr::select(-.data$map)
        }
        rownames(bintolen) <- bintolen$bins
        allbins <- data.frame(bins = names(binsGR), stringsAsFactors = FALSE)
        bintolen <- suppressWarnings(dplyr::left_join(allbins, bintolen) %>% 
                            tidyr::replace_na(list(gc = 0, map = 0, len = 0)))
        return(bintolen)


    }
}
    



#' @export
construct_features <- function(gen = "Hsapiens", gen_ver = "hg19",
                      sig = "GATC", bin_type = "Bins-uniform", binsize = 5000, 
                      wg_file = NULL, chrs = NULL, feature_type = "RE-based") {
    
    if (is.null(chrs)) {
        chrs <- get_chrs(gen, gen_ver)
    }
    bintolen<-lapply(chrs,function(x) data.frame(stringsAsFactors = FALSE))
    names(bintolen)<-chrs
    for (chrom in names(bintolen)){
    bintolen[[chrom]] <- construct_features_chr(chrom=chrom,
                                          gen = gen, 
                                          gen_ver = gen_ver,
                                          sig = sig, 
                                          bin_type = bin_type, 
                                          binsize = binsize, 
                                          wg_file = wg_file,
                                          feature_type = feature_type)
    }
    bintolen<-suppressWarnings(dplyr::bind_rows(bintolen))
    
    return(bintolen)
}
