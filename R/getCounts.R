#setDTthreads(8)
.importFragments <- function(fragData, chunkSize=10)
  {
    if(is.list(fragData))
    {
      if(grepl(".bam", fragData[[1]], fixed=TRUE))
      {
        print("Importing .bam files")
        param <- Rsamtools::ScanBamParam(what=c('pos', 'qwidth', 'isize'))
        fragData <- lapply(fragData, function(bamPath){
          readPairs <- GenomicAlignments::readGAlignmentPairs(bamPath, 
                                                              param=param)
          # get fragment coordinates from read pairs
          frags <- GRanges(seqnames(GenomicAlignments::first(readPairs)), 
            IRanges(start=pmin(
              GenomicAlignments::start(GenomicAlignments::first(readPairs)), 
              GenomicAlignments::start(GenomicAlignments::second(readPairs))), 
                    end=pmax(
                      GenomicAlignments::end(GenomicAlignments::first(readPairs)), 
                      GenomicAlignments::end(GenomicAlignments::second(readPairs)))))
          frags <- granges(frags, use.mcols=TRUE)
          
          # ATAC shift
          start(frags) <- ifelse(strand(frags) == "+", 
            start(frags) + 4, start(frags))
          end(frags) <- ifelse(strand(frags) == "-", 
            end(frags) - 5, end(frags))
          frags <- as.data.table(frags)
          #setnames(frags, c("seqnames"), c("chr"))
          frags$count <- 1
          frags
        })
      }
      else if(grepl(".bed", fragData[[1]], fixed=TRUE))
      {
        fragData <- lapply(fragData, data.table::fread, select=1:3, 
                           col.names=c("seqnames", "start", "end"))
        fragData <- lapply(fragData, function(dt){dt$count <- 1
        dt <- dt[,c("seqnames", "start", "end", "count"), with=FALSE]
        dt})
      }
      else if(grepl(".rds", fragData[[1]], fixed=TRUE))
      {
        fragData <- lapply(fragData, function(gr){
          gr <- readRDS(gr)
          dt <- as.data.table(gr)
          #setnames(dt, c("seqnames"), c("chr"))
          dt$count <- 1
          dt <- dt[,c("seqnames", "start", "end", "count")]
          dt
        })
      }
      else if(is.data.table(fragData[[1]]))
      {
        if(!(grepl(paste(c("seqnames", "start", "end"), collapse=";"),
                   paste(colnames(fragData[[1]]),collapse=";"))))
        {stop("data.table list elements of data need to contain columns: seqnames, start, end")}
      }
      else stop("List elements of data need to be either data.tables, .bams-, .beds- .rds-files")
    }
    else
    {
      stop("Data needs to be a list")
    }
    
    return(fragData)
}


#' @description
#' Sanity check to ensure the input arguments have the correct classes
#' 
#' @param atacFrag: a list of data.tables that contain the ranges of fragments
#' @param peakRanges: a GRange object that contains the ranges of peaks
#' @param motifRanges:  a GRange object that contains the ranges of motifs, 
#' check metadata columns
#' check seqnames to factor in datatable
.sanityCheck <- function(atacFrag, 
  ranges, 
  type = c("peaks", "motifs")) {
  lapply(atacFrag, function(frag) {
      if (!is.data.table(frag)) {
        stop("Each element in atacFrag should be a data.table")
      }
  })
  
  if (!class(ranges)=="GRanges") {
    stop("ranges should be a GRanges object")
  }
  
  type <- match.arg(type, choices = c("peaks", "motifs"))
  if (type=="motifs") {
    if (!("motif" %in% names(mcols(ranges)))) {
      stop("There is no motif names in metadata columns")
    }
  }
    
}




#' @description
#' Resize the peaks
#' 
#' @param peakRanges: a GRanges object of peak ranges
#' @param width: the re-defined size of each peak
#' @return a GRange object with resized ranges

.resizeRanges <- function(peakRanges, 
  width = 300, 
  fix = c("center", "start", "end", "summit"),
  ...) {
  
      fix <- match.arg(fix, choices = c("center", "start", "end", "summit"))
      # Sanity check
      if (!class(peakRanges) == "GRanges") {
        stop("peakRanges must be a GRanges object")
      }
      
      if (fix == "summit") {
        start(peakRanges) <- round(peakRanges$summit-width/2)
        end(peakRanges) <- start(peakRanges)+width-1
      } else {
        peakRanges <- resize(peakRanges, width = width, fix = fix)
      }
      
      return(peakRanges)
      
}

#' @param peakRanges: a GRanges object of peak ranges
#' @param genome: a BSgenome object, the corresponding genome 
#' @return a GRanges objects with an additional metadata column gc that contains
#' GC content
.getGCContent <- function(peakRanges, genome) {
    # Sanity check
    if (!class(peakRanges) == "GRanges") {
        stop("peakRanges must be a GRanges object")
    }
        
    peakSeqs <- Biostrings::getSeq(x = genome, peakRanges)
    mcols(peakRanges)$gc <- letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
    peakRanges
} 


.filterFrags <- function(atacFrag, min = 30, max = 2000) {
    res <- lapply(atacFrag, function(frag) {
        frag[,width:=end-start+1]
        frag[width>=min & width<=max,]
        frag
    })
    res
}

.matchSeqlevels <- function(atacFrag, ranges) {
    frags <- rbindlist(atacFrag)
    fragSeq <- unique(frags$seqnames)
    rangeSeq <- GenomicRanges::seqnames(ranges) 
    common <- intersect(fragSeq,rangeSeq)
    atacFrag <- lapply(atacFrag, function(frag) {
     frag <- frag[seqnames %in% common,]
     frag$seqnames <- factor(frag$seqnames)
     frag})
    ranges <- ranges[seqnames(ranges) %in% common,]
    # turn seqnames to factor
    list(atacFrag=atacFrag, ranges=ranges)
}

#TODO: Functionality
# - filter max and min length of fragments (30 & 2000 as defaults)
# - match seqLevels of ranges and atac frags in the beginning
# - write helper functions for sanity checks, i.e. for data.table / granges conversions etc
# - write a helper function for ATAC shifting
# - conversion seqLevels to integer

#TODO: Speed-up & memory usage
# - use dtToGr helper for data.table to GenomicRanges conversion
# - overlap with findOverlaps and use indices 
# - convert samples (!), chr & motif ids to integers / numerical factors wherever possible
# - test chunking across chromosomes: Different chunk sizes

#TODO:
# 1. unify the arguments 
# 2. filter functions
# 3. seqLevels

#' @Author: Emanuel Sonder
dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end"){
  setnames(dt, seqCol, "seqnames", skip_absent = TRUE)
  gr <- GRanges(seqnames=dt[["seqnames"]], ranges=IRanges(start=dt[[startCol]], 
                                                          end=dt[[endCol]]))
  mcols(gr)$count <- dt$count
  if (!is.null(dt$sample)) {
    mcols(gr)$sample <- dt$sample
  }
  
  if (!is.null(dt$motif)) {
    mcols(gr)$motif <- dt$motif
  }
  
  if (!is.null(dt$barcode)) {
    mcols(gr)$barcode <- dt$barcode
  }
  
  if(startCol==endCol)
  {
    gr <- GPos(seqnames=dt[["seqnames"]], pos=dt[[startCol]])
  }
  return(gr)
}

.standardChromosomes <- function(gr, species) {
  gr <- keepStandardChromosomes(gr,
                                species=species,
                                pruning.mode="coarse")
  seqlevelsStyle(gr) <- "UCSC"
  gr
}


#' @Author: Emanuel Sonder
#' @Author: Emanuel Sonder
.getInsertionProfiles <- function(atacFrag, 
    motifRanges,
    margin=200,
    #minWidth=30,
    #maxWidth=2000,
    nullModel=FALSE,
    returnAll=FALSE,
    symmetric=TRUE,
    libNorm=FALSE,
    chunk=TRUE){
    
    # prep motif data
    motifData <- as.data.table(motifRanges)
    motifMarginRanges <- as.data.table(resize(motifRanges, width=2*margin, fix="center"))
    setnames(motifMarginRanges, c("start", "end"), c("start_margin", "end_margin"))
    motifData <- cbind(motifData, motifMarginRanges[,c("start_margin", "end_margin"), with=FALSE])
    
    setnames(motifData, "seqnames", "chr")
    chrLevels <- unique(motifData$chr)
    motifLevels <- unique(motifData$motif_id)
    
    # convert to factors (memory usage)
    motifData[,chr:=as.integer(factor(chr, 
        levels=chrLevels, ordered=TRUE))]
    motifData[,motif_id:=factor(motif_id, 
        levels=motifLevels, ordered=TRUE)]
    
    # determine motif center
    motifData[,motif_center:=floor((end-start)/2)+start]
    
    # determine margins
    #motifData[,ml:=end-start+1, by=motif_id]
    #motifData[,start_margin:=start-(margin-ceiling(ml/2))]
    #motifData[,end_margin:=end+(margin-ceiling(ml/2))]
    
    # determine motif center
    #motifData[,motif_center:=floor((end_margin-start_margin)/2)+start_margin]
    
    # convert to factors (memory usage)
    if("seqnames" %in% colnames(atacFrag)){
        setnames(atacFrag, "seqnames", "chr")
    }
    atacFrag <- copy(atacFrag) #TODO: take out that copy 
    atacFrag[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
    
    
    medZero <- function(x, len){median(c(rep(0,len-length(x)),x))}
    
    nSamples <- length(unique(atacFrag$sample))
    
    setorder(motifData, chr)
    setorder(atacFrag, chr)
    motifData[,motif_match_id:=1:nrow(motifData)]
    motifData <- split(motifData, by="chr")
    atacFrag <- split(atacFrag, by="chr")
    
    ptm <- proc.time()
    atacProfiles <- mapply(function(md,af){
        
        #md[,motif_match_id:=1:nrow(md)]
        
        # convert to granges for faster overlapping
        motifMarginRanges <- dtToGr(md, startCol="start_margin", endCol="end_margin", seqCol="chr")
        atacStartRanges <- dtToGr(af, startCol="start", endCol="start", seqCol="chr")
        atacEndRanges <- dtToGr(af, startCol="end", endCol="end", seqCol="chr")
        
        startHits <- findOverlaps(atacStartRanges, 
            motifMarginRanges, type="within") # check if type within faster or slower
        endHits <- findOverlaps(atacEndRanges, motifMarginRanges, type="within") 
        
        # get overlapping insertion sites
        atacStartInserts <- af[queryHits(startHits), c("sample", "start"), with=FALSE]
        atacEndInserts <- af[queryHits(endHits), c("sample", "end"), with=FALSE]
        setnames(atacStartInserts, "start", "insert")
        setnames(atacEndInserts, "end", "insert")
        
        ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
            rbindlist(list(
                md[subjectHits(startHits), c("motif_center", "start", 
                    "end", "motif_id", "motif_match_id")],
                md[subjectHits(endHits), c("motif_center", "start", "end", "motif_id", 
                    "motif_match_id")])))
        
        # count insertions around motif
        ai[,rel_pos:=insert-motif_center]
        ai[,type:=fifelse(insert>=start & insert<=end, 1,0)]
        
        # calculate motif length
        ai[,ml:=end-start+1, by=motif_id]
        
        ap <- ai[,.(pos_count_global=.N), by=.(ml, rel_pos, sample, motif_id, type)]
        gc()
        ap
    }, 
        motifData, 
        atacFrag, 
        SIMPLIFY=FALSE)
    print(proc.time()-ptm)
    atacProfiles <- rbindlist(atacProfiles, idcol="seqnames")
    
    if(!nullModel){
        
        atacProfiles <- atacProfiles[,.(pos_count_global=sum(pos_count_global)),
            by=.(rel_pos, motif_id, type, ml)]
        
        atacProfilesMotif <- subset(atacProfiles, type==1)
        #atacProfilesMotif$rel_pos_m <-0
        atacProfilesMotif[,med_pos_count_global:=medZero(pos_count_global, first(ml)), by=.(motif_id)]
        atacProfilesMotif <- atacProfilesMotif[,.(pos_count_global=(pos_count_global+first(med_pos_count_global))/2), 
            by=.(rel_pos, motif_id)]
        
        atacProfilesMargin <- subset(atacProfiles, type==0)
        atacProfiles <- rbind(atacProfilesMargin, atacProfilesMotif, fill=TRUE)
        
        # calculate weights
        setorder(atacProfiles, motif_id, rel_pos)
        atacProfiles[,w:=smooth(pos_count_global, twiceit=TRUE), by=motif_id]
        if(symmetric) atacProfiles[,w:=rev(w)+w, by=motif_id]
        
        atacProfiles[,w:=length(w)*w/sum(w), by=motif_id]
        #atacProfiles[,w:=w/sum(w), by=motif_id]
        
        #atacProfiles[,motif_name:=motifLevels[motif_id]]
    }
    else
    {
        # uniform weighting
        atacProfiles <- atacProfiles[,.(w=1), by=.(motif_id, rel_pos)]
        atacProfiles[,w:=w/sum(w), by=motif_id]
        atacProfiles[,w:=length(w)*w/sum(w), by=motif_id]
        # got the error atacInserts does not exist
        #atacProfiles <- atacProfiles[,.(w=1/.N), by=.(motif_id, rel_pos_m)]
        #atacProfiles[,motif_name:=motifLevels[motif_id]]
    }
    
    # get match scores
    motifScores <- mapply(function(md,af){
        
        # convert to granges for faster overlapping
        motifMarginRanges <- dtToGr(md, startCol="start_margin", endCol="end_margin", seqCol="chr")
        atacStartRanges <- dtToGr(af, startCol="start", endCol="start", seqCol="chr")
        atacEndRanges <- dtToGr(af, startCol="end", endCol="end", seqCol="chr")
        
        startHits <- findOverlaps(atacStartRanges, 
            motifMarginRanges, type="within") # check if type within faster or slower
        endHits <- findOverlaps(atacEndRanges, motifMarginRanges, type="within") 
        
        # get overlapping insertion sites
        atacStartInserts <- af[queryHits(startHits), c("sample", "start"), with=FALSE]
        atacEndInserts <- af[queryHits(endHits), c("sample", "end"), with=FALSE]
        setnames(atacStartInserts, "start", "insert")
        setnames(atacEndInserts, "end", "insert")
        
        ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
            rbindlist(list(
                md[subjectHits(startHits), c("motif_center", "seqnames", "start", 
                    "end", "motif_id", "motif_match_id")],
                md[subjectHits(endHits), c("motif_center","seqnames", "start", "end", "motif_id", 
                    "motif_match_id")])))
        
        # count insertions around motif
        ai[,rel_pos:=insert-motif_center]
        ai <- ai[,.(pos_count=.N), 
            by=.(motif_match_id, motif_id, sample, rel_pos)]
        
        ai <- merge(ai, atacProfiles[,c("rel_pos", "motif_id", "w"), with=FALSE], 
            by.x=c("motif_id","rel_pos"), 
            by.y=c("motif_id","rel_pos"), all.x=TRUE)
        ai[,score:=w*pos_count]
        as <- ai[,.(score=sum(score), 
            tot_count=sum(pos_count)), 
            by=.(motif_match_id, motif_id, sample)]
        gc()
        # get ai[,rel_pos_m:=fifelse(insert<start, insert-start, 0)]
        #     ai[,rel_pos_m:=fifelse(insert>end, insert-end, rel_pos_m)]
        as
    }, 
        motifData, 
        atacFrag, 
        SIMPLIFY=FALSE)
    
    motifScores <- rbindlist(motifScores)
    # motifScores <- motifScores[,.(score=sum(score),
    #     tot_inserts=sum(tot_count)),
    #     by=.(seqnames, start, end, motif_match_id, motif_id, sample)]
    
    if(libNorm)
    {
        motifScores[,tot_count:=sum(tot_count), by=.(sample)]
        motifScores[,score:=score/tot_count]
    }
    
    motifData <- rbindlist(motifData)
    motifData[,seqnames:=chrLevels[seqnames]]
    motifScores <- cbind(motifScores, 
                         motifData[motifScores$motif_match_id,
                                   c("start", "end", "seqnames"), with=FALSE])
    #chrLevels
    motifScores[,seqnames:=chrLevels[seqnames]]
    #motifScores[,nmotif_name:=motifLevels[motif_id]]
    #return(list(ms=motifScores, ap=atacProfiles))
    #return(list(motifScores, atacProfiles))
    return(motifScores)
}


# Comment: 
# - I removed the sanity checks as we either write a seperate function for that 
# or we do it once in the main function at the beginning
# - same for the atac shifts

.atacShift <- function(atacFrag, shifts=c(4L, 5L)){
  atacFrag[, start := ifelse(strand == "+", start + shifts[1], start)]
  atacFrag[, end   := ifelse(strand == "-", end   - shifts[2], end)]
  atacFrag
}

# something like this: might add more checks
.checkRanges <- function(ranges){
  # Sanity check
  if (!is.data.table(ranges)) {
    dt <- data.table::as.data.table(ranges)
  } else {
    dt <- data.table::copy(ranges)
  }
  
  # if ("seqnames" %in% colnames(dt))
  #   colnames(dt)[which(colnames(dt)=="seqnames")] <- "chr"
  
  return(dt)
}

#' @param atacFrag a GRange or data.table object that contains atac fragment ranges
#' @param motifRanges a GRange or data.table object that contains motif ranges
#' @param flankSize integer, the number of nucleotides to define the buffer/flanking 
#' region near the motif instance
#' @param shiftATAC logic if shifting the ATAC-seq fragment data table
#' @Author: Emanuel Sonder
.getInsertionCounts <- function(atacFrag, 
                                motifRanges,
                                mode=c("total", "weight"),
                                flankSize=30,
                                shiftATAC=FALSE,
                                weightCol="weight", 
                                #addProfile=FALSE,
                                ...) {
  
  mode <- match.arg(mode, choices=c("total", "weight"))
  
  atacFrag <- .checkRanges(atacFrag)
  motifData <- .checkRanges(motifRanges)
  
  # ATAC insertion shifts
  if (shiftATAC == TRUE) {
    atacFrag <- .atacShift(atacFrag, ...)
  }
  
  # set up flanking region
  motifData[, start_margin := start - flankSize]
  motifData[, end_margin   := end + flankSize]
  motifData$match_id <- seq_len(nrow(motifData))
  
  # Either count weights or total number of fragments
  if(mode=="total"){
    atacFrag$weight <- 1
  } else {
    setnames(atacFrag, weightCol, "weight")
  }
  
  # get number of samples and motif matches
  nSample <- length(unique(atacFrag$sample))
  nMatch <- length(unique(motifData$match_id))
  
  # maybe its faster to chunk that part across samples -------------------------
  # convert to granges for faster overlaps
  fragStarts <- dtToGr(atacFrag, startCol="start", endCol="start")
  fragEnds <- dtToGr(atacFrag, startCol="end", endCol="end")
  motifMarginRanges <- dtToGr(motifData, 
                              startCol="start_margin", endCol="end_margin")
  
  # refactor this: Try 
  # - with sparse matrix depending on zero fraction
  # - findOverlaps invert=TRUE
  # - (!) if this is run in a loop by sample also countByOverlaps could be used
  startInsertsMargin <- as.data.table(findOverlaps(motifMarginRanges, fragStarts))
  startInsertsMargin <- cbind(startInsertsMargin, 
                              atacFrag[startInsertsMargin$subjectHits, 
                                c("sample", "weight")])
  endInsertsMargin <- as.data.table(findOverlaps(motifMarginRanges, fragEnds))
  endInsertsMargin <- cbind(endInsertsMargin, 
                            atacFrag[endInsertsMargin$subjectHits, 
                              c("sample", "weight")])
  
  # get matrix with inserts within the motif and the margin
  insertsMargin <- rbind(startInsertsMargin, endInsertsMargin)
  # get zero inserts
  insertsMargin <- rbind(insertsMargin, data.table(
    queryHits=setdiff(seq_len(length(motifRanges)), 
      unique(insertsMargin$queryHits)),
    weight=0, 
    sample=insertsMargin$sample[1]), # give some dummy sample
    use.names=TRUE, fill=TRUE)
  
  insertsMargin <- dcast(insertsMargin, queryHits~sample, 
                         value.var="weight",
                         fill=0,
                         fun.aggregate=sum, drop=FALSE)
  insertsMargin$queryHits <- NULL
  
  # get the inserts within the motif
  startInsertsWithin <- as.data.table(findOverlaps(motifRanges, fragStarts))
  startInsertsWithin <- cbind(startInsertsWithin, 
                              atacFrag[startInsertsWithin$subjectHits, 
                                c("sample", "weight")])
  endInsertsWithin <- as.data.table(findOverlaps(motifRanges, fragEnds))
  endInsertsWithin <- cbind(endInsertsWithin, 
                            atacFrag[endInsertsWithin$subjectHits, 
                              c("sample", "weight")])
  
  insertsWithin <- rbind(startInsertsWithin, endInsertsWithin)
  # get missing combinations 
  insertsWithin <- rbind(insertsWithin, data.table(
    queryHits=setdiff(seq_len(length(motifRanges)), 
      unique(insertsWithin$queryHits)),
    weight=0, 
    sample=insertsWithin$sample[1]), # give some dummy sample
    use.names=TRUE, fill=TRUE)
  
  insertsWithin <- dcast(insertsWithin, queryHits~sample, 
                         value.var="weight",
                         fill=0,
                         fun.aggregate=sum, drop=FALSE)
  insertsWithin$queryHits <- NULL
  
  # ----------------------------------------------------------------------------
  
  flankingCounts <- insertsMargin - insertsWithin
  
  return(list(total_counts=insertsMargin, 
              flanking_counts=flankingCounts, 
              within_counts=insertsWithin))
}

#' @description
#' count the number of fragments in each peak region
#' 
#' @param peakRanges a GRange or data.table object that contains the peak ranges
#' @param atacFrag a list of data.table that contains the fragment ranges, 
#' each data.table represents a sample
#' return a list of count table; by all, nucleosome-free, mono...
.getOverlapCounts <- function(peakRanges, 
    atacFrag,
    mode = c("total", "weight"),
    cuts,
    genome,
    overlap = c("any", "start", "end", "within", "equal"),
    smooth, 
    aRange,
    nWidthBins,
    nGCBins,
    peakWeight = c("none", "loess", "lm", "wlm", "tlm", "wtlm"),
    moderating = FALSE,
    singleCell,
    ...
  ) {
  
    by <- match.arg(mode, choices = c("total", "weight"))
    peakWeight <- match.arg(peakWeight, 
        choices = c("none", "loess", "lm", "wlm", "tlm", "wtlm"))
      
    peaks <- data.table::as.data.table(peakRanges)
    peakGR <- peakRanges
    
    frags <- data.table::copy(atacFrag)
    peaks$peakID <- seq_len(nrow(peaks))
    
    if (by == "weight") {
      frags <- .weightFragments(frags, 
        genome = genome, 
        smooth = smooth,
        nGCBins = nGCBins,
        nWidthBins = nWidthBins,
        aRange = aRange,
        moderating = moderating,
        singleCell = singleCell)
    }
    
    frags <- .getType(frags, cuts = cuts)
    
    fragCounts <- lapply(frags, function(frag) {
      types <- names(frag)[grepl("^type_", names(frag))]
      if (by == "weight") {
        frag[,count:=weight*count]
        frag[,(types) := lapply(.SD, function(x) x*weight), 
          .SDcols = types]
      } 
      fragGR <- dtToGr(frag)
      hits <- findOverlaps(fragGR, peakGR, type = overlap)
      overlaps <- cbind(frag[queryHits(hits),],
               peaks[subjectHits(hits), c("peakID")])
      tmp <- overlaps[, c(list(counts = sum(count, na.rm = TRUE)), 
          list(mean_width = mean(width, na.rm = TRUE)),
          list(median_width = median(width, na.rm = TRUE)),
          lapply(.SD, sum, na.rm = TRUE)), 
          by = peakID, .SDcols = c(types)]
      
      res <- data.table(peakID = seq_len(nrow(peaks)))
      res <- merge(res, tmp, all =TRUE)
      res[is.na(res)] <- 0
      res
      
    })
    
    cols <- names(fragCounts[[1]])[grepl("^type_|counts|_width", 
      names(fragCounts[[1]]))]
    allCounts <- lapply(cols, function(x) {
      lst <- lapply(fragCounts, function(.) data.frame(.)[,x])
      mat <- do.call(cbind, lst)
      colnames(mat) <- names(fragCounts)
      mat
    })
    names(allCounts) <- cols
    if (peakWeight != "none") {
        allCounts[["counts"]] <- .weightPeaks(allCounts[["counts"]], 
          method=peakWeight)
        allCounts[["type_1"]] <- .weightPeaks(allCounts[["type_1"]], 
          method=peakWeight)
    }
    
    
    allCounts
    
    
}

#' change cuts into bin size; given the parameter of nBin;
#' first calculate intervals
.getType <- function(atacFrag, cuts=c(0,120,300,500)) {
    # check the min of cuts
    if (cuts[1] != 0) cuts <- c(0,cuts)
    res <- lapply(names(atacFrag), function(sample) {
        dt <- atacFrag[[sample]]
        dt[,width:=end-start+1]
        # check if max(cuts) covers max(width)
        if (max(dt$width) > cuts[length(cuts)]) 
            cuts[length(cuts)+1] <- max(dt$width)
        
        dt[,type:=as.numeric(cut(width, breaks = cuts))]
        if (length(unique(dt$type))==1) {
            warnings("All fragments fell into 1 type")
        } else if(length(unique(dt$type))>=6) {
            stop("Too many types!")
        }
        for (i in unique(dt$type)) {
          dt[, paste0("type_",i) := as.integer(type == i)]
        }
        dt

    })
    names(res) <- names(atacFrag)
    res
}

## bin by both fragment length and GC content
.getBins <- function(atacFrag, 
  nWidthBins = 30, 
  nGCBins = 10, 
  genome) {
    
    fragDts <- lapply(names(atacFrag), function(x){
      dt <- atacFrag[[x]]
      dt[,width:=end-start+1]
      gr <- dtToGr(dt)
      gr <- .getGCContent(gr, genome = genome)
      dt <- as.data.table(gr)
      dt[,sample:=x]
      dt
    })
    fragDt <- rbindlist(fragDts)
    
    widthIntervals <- unique(quantile(fragDt$width, 
                                      probs = seq(0,1,by=1/nWidthBins)))
    GCIntervals <-  unique(quantile(fragDt$gc, 
                                    probs = seq(0,1,by=1/nGCBins)))
    
    
    fragDt[,widthBin:=cut(width, 
        breaks=widthIntervals, 
        include.lowest=TRUE)]
    fragDt[,GCBin:=as.numeric(cut(gc, 
        breaks=GCIntervals, 
        include.lowest=TRUE))]
    fragDt
    
}
  

#' moderateBinFrequencies
#'
#' @param bins A vector of FL/GC bin labels for each element of `counts`
#' @param samples A vector of sample labels for each element of `counts`
#' @param counts A vector of counts per sample/bin
#'
#' @return A vector of moderated frequencies for each element of `counts`
#' @author Pierre-Luc
moderateBinFrequencies <- function (bins, samples, counts) {
  totPerSamp <- tapply(counts, samples, sum)
  binFreq <- counts/totPerSamp[samples]
  mv <- aggregate(binFreq, by=list(bin=bins), FUN=function(x){
    c(mu=mean(x,na.rm=TRUE), v=var(x, na.rm=TRUE))
  })
  mv$mu <- mv$x[,1]
  mv$v <- mv$x[,2]
  mv$alpha <- ((1-mv$mu)/mv$v - 1/mv$mu)*mv$mu^2
  mv$beta <- mv$alpha*(1/mv$mu-1)
  (counts + mv$alpha[bins])/(totPerSamp[samples] + mv$alpha[bins] + mv$beta[bins])
}  




.weightFragments <- function (atacFrag, 
    genome,
    smooth = c("none", "smooth.2d"),
    nWidthBins,
    nGCBins,
    aRange,
    moderating,
    ...) {
    smooth <- match.arg(smooth, choices = c("none", "smooth.2d"))
    #' TODO: check if .getBins can be improved
    fragDt <- .getBins(atacFrag, genome = genome, 
      nWidthBins = nWidthBins, nGCBins = nGCBins)
    fragDt[, bin:=paste0(widthBin, GCBin)]
    fragDt[, bin:=as.integer(as.factor(bin))]
    fragDt[,count_bin:=sum(count), by=c("sample", "bin")]
    # add estimateBetaParams
    if (moderating) {
      dt <- unique(fragDt, by=c("bin","sample"))
      dt <- dt[, c("sample", "count_bin", "bin"), with = FALSE]
      dt$freq_bin <- moderateBinFrequencies(dt$bin, dt$sample, dt$count_bin)
      fragDt <- merge(fragDt, dt, by = c("bin","sample"))
    } else {
      fragDt[,freq_bin:=(count_bin+1L)/(sum(count_bin)+1L),by=sample]
    }
    tmp <- fragDt[,.(mean_freq_bin=mean(freq_bin, na.rm=TRUE)), 
      by=c("bin")]
    fragDt <- merge(fragDt, tmp, by = "bin")
    fragDt[,weight:=mean_freq_bin/freq_bin]
    if (smooth=="none") {
      return(split(fragDt, fragDt$sample))
    } else if (smooth=="smooth.2d") {
      fb <- fragDt[,c("sample","GCBin","widthBin","weight")]
      fb <- fb[!duplicated(fb),]
      dts <- lapply(split(fb, fb$sample), function(dt) {
        dt[,logWeight:=log2(weight)]
        sm <- smooth.2d(dt$logWeight, x=cbind(dt$widthBin, dt$GCBin), 
          surface=FALSE, 
          nrow=length(unique(fragDt$widthBin)), 
          ncol=length(unique(fragDt$GCBin)),
          aRange=aRange)
        dimnames(sm) <- list(levels(dt$widthBin), levels(factor(dt$GCBin)))
        sm
      })
      tbl <- melt(dts)
      names(tbl) <- c("widthBin", "GCBin", "smooth", "sample")
      fragDt <- merge(fragDt, 
        tbl, 
        by = c("GCBin", "widthBin", "sample"),
        allow.cartesian=TRUE)
      fragDt[,weight:=2^smooth]
      return(split(fragDt, fragDt$sample))
      
    }
}




#loess, lm, weighted lm, trimmed lm, and weighted trimmed lm
.weightPeaks <- function(counts, 
    method,
    #group_id = NULL,
    #trimM = 0.15,
    #trimA = 0.05,
    ...) {
    # lfc <- sapply(seq_len(ncol(counts)), \(i) log2((counts[,i]+1L)/(rowMeans(counts+1L)))) 
    # colnames(lfc) <- colnames(counts)
    # avg <- rowMeans(log2(counts+1L))
    
    if (method == "loess") {
      if(nrow(counts) <= 1e4) {
          res <- affy::normalize.loess(counts+1L, span=0.3, 
              subset=1:nrow(counts), family.loess="symmetric")
      } else {
        nBins <- 100
        nSample <- 1e4
        allAvg <- log2(rowMeans(counts+1L))
        bins <- cut(allAvg,
          breaks = seq(min(allAvg), max(allAvg), (max(allAvg)-min(allAvg))/nBins),
          include.lowest = TRUE)
        idx <- unlist(sapply(seq_len(length(levels(bins))), \(i) {
          ids <- which(bins==levels(bins)[i])
          if (length(ids) > 1)
            sample(ids,
              size = min(length(ids), round(nSample/nBins)))
          else if (length(ids)==1)
            ids
          else NULL
        }))
        res <- affy::loess.normalize(counts+1L, span=0.3, 
            subset=idx, family.loess="symmetric")
      }
    } else if (method == "lm") {
        nBins <- 20
        nSample <- 1e4
        logAvg <- log1p(avg)
        bins <- cut(logAvg,
            breaks = seq(min(logAvg), max(logAvg), (max(logAvg)-min(logAvg))/nBins),
            include.lowest = TRUE)
        idx <- unlist(sapply(seq_len(length(levels(bins))), \(i) {
            ids <- which(bins==levels(bins)[i])
            if (length(ids) > 1)
                sample(ids,
                    size = min(length(ids), round(nSample/nBins)))
            else if (length(ids)==1)
                ids
            else NULL
        }))
        subcounts <- counts[idx,]
        sublfc <- lfc[idx,]
        subavg <- avg[idx]
        models <- lapply(colnames(counts), \(x) {
          df <- data.frame(lfc=sublfc[,x], avg=subavg)
          lm(lfc ~ 0+avg, df)
      })
    } else if (method=="wlm") {
        w <- rowMeans(sqrt(cpm(counts)))
        models <- lapply(colnames(counts), \(x) {
          df <- data.frame(lfc = lfc[,x], avg=avg)
          model <- lm(lfc ~ avg, df, weights=w)
        })
    } else if (method=="wtlm") {
      models <- lapply(colnames(counts), \(x){
        df <- data.frame(lfc = lfc[,x], avg=avg)
        lqm <- quantile(df$lfc, trimM, na.rm = TRUE)
        uqm <- quantile(df$lfc, 1-trimM, na.rm = TRUE)
        lqa <- quantile(df$avg, trimA, na.rm = TRUE)
        uqa <- quantile(df$avg, 1-trimA, na.rm = TRUE)
        trimRes <- df %>%
          filter(
            lfc >= lqm & lfc <= uqm,
            avg >= lqa & avg <= uqa
          )
        if (method=="wtlm") {
          lm(lfc ~ avg, trimRes, weights=sqrt(2^(trimRes$avg)))
        } else {
          lm(lfc ~ avg, trimRes)
        }
      })
        
        
    }
    
    if (method!="loess") {
        newLFC <- sapply(seq_len(length(models)), 
            \(i) predict(models[[i]], data.frame(avg=avg)))
        res <- (2^(-newLFC))*counts
        colnames(res) <- colnames(counts)
    } else {
        return(res)
    }
}




#' getCounts
#' 
#' @description
#' makes matrix of fragment counts in peaks using a list of bam or bed files
#' @param files a list of filenames for bam or bed files with aligned reads
#' @param ranges a GRanges that contains the ranges of peaks or motifs. If 
#' rowType == "motifs", ranges must contain a mcols called `motif` that contains
#' the name of motif for each row
#' @param rowType using peaks or motifs in rows
#' @param mode if the assays return original counts or weighted counts
#' @param paired paired end data?
#' @param resize a logical to determine if peaks should be resized
#' @param width if resize is `TRUE`, the new width of each peak
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#'  object

getCounts <- function (files,
    atacFrag,
    ranges, # motif matches or peaks, set a parameter to define
    genome,
    species,
    rowType = c("peaks", "motifs"),
    mode = c("total", "weight"),
    paired = TRUE,
    resize = TRUE, 
    width = 300,
    nWidthBins = 30,
    nGCBins = 10,
    cuts = c(0,120,300,500),
    minFrag = 30,
    maxFrag = 3000,
    smooth = c("none", "smooth.2d"),
    aRange = 0,
    peakWeight = c("none", "loess", "lm", "wlm", "tlm", "wtlm"),
    moderating = FALSE,
    singleCell = FALSE,
    ...) {
    
    rowType <- match.arg(rowType, choices=c("peaks", "motifs"))  
    mode <- match.arg(mode, choices=c("total", "weight"))  
    smooth <- match.arg(smooth, c("none", "smooth.2d"))
    peakWeight <- match.arg(peakWeight, 
      choices=c("none", "loess", "lm", "wlm", "tlm", "wtlm"))  
    # get fragments ranges
    if (is.null(atacFrag)) atacFrag <- .importFragments(files)
    
    # sanity check
    .sanityCheck(atacFrag, ranges, type = rowType)
    
    # standard chromosomes
    ranges <- .standardChromosomes(ranges, species = species)
    atacFrag <- lapply(atacFrag, function(dt) {
        gr <- dtToGr(dt)
        gr <- .standardChromosomes(gr, species = species)
        as.data.table(gr)
    })

    
    
    # filter too short or too long fragments
    atacFrag <- .filterFrags(atacFrag, min = minFrag, max = maxFrag)
    
    # match seqLevels
    res <- .matchSeqlevels(atacFrag, ranges)
    atacFrag <- res$atacFrag
    ranges <- res$ranges

    
    if (rowType=='peaks' & resize) {
      ranges <- .resizeRanges(peakRanges = ranges, width = width)
    }
    
    if (rowType=='peaks') {
        asy <- .getOverlapCounts(peakRanges = ranges, 
            atacFrag = atacFrag,
            mode = mode,
            cuts = cuts,
            genome = genome,
            smooth = smooth,
            aRange = aRange,
            species = species,
            nWidthBins = nWidthBins,
            nGCBins = nGCBins,
            peakWeight = peakWeight,
            moderating = moderating,
            singleCell = singleCell)
    } else if (rowType=="motifs"){
        if (mode=="weight") {
          atacFrag <- .weightFragments(atacFrag, genome)
            
        }
        lst <- lapply(names(atacFrag), function(x) {
            frag <- atacFrag[[x]]
            frag$sample <- x
            frag
        })
        atacFrags <- rbindlist(lst)
        asy <- .getInsertionCounts(atacFrags, 
          motifRanges = ranges, 
          mode = mode)
    }
    if (singleCell) {
        cd <- rbindlist(lapply(names(atacFrag), \(x) 
            data.table(sample_id=atacFrag[[x]]$sample[1],
                group_id=atacFrag[[x]]$motif[1],
                barcode=x)))
        cd <- cd[match(cd$barcode, colnames(asy[[1]]))]
        SummarizedExperiment(assays = asy, 
            rowRanges = ranges,
            colData = DataFrame(cd))
    } else {
        SummarizedExperiment(assays = asy, rowRanges = ranges)  
    }

}

