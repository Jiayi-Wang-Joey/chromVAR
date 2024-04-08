suppressPackageStartupMessages({
    source("~/chromVAR/R/getCounts.R")
    source("/mnt/plger/esonder/R/tfbs_pred_comp/package/package/R/TfBindingMultiAssay.R")
    library(motifmatchr)
})

#' @author Emanuel Sonder
.getInsertionProfiles <- function(fragDt, 
    motifRanges,
    margin=100,
    aggFun=sum,
    #minWidth=30,
    #maxWidth=2000,
    chunk=TRUE){
    
    # prep motif data
    motifData <- as.data.table(motifRanges)
    chrLevels <- unique(motifData$seqnames)
    motifLevels <- unique(motifData$motif_id)
    
    # convert to factors (memory usage)
    motifData[,seqnames:=as.integer(factor(seqnames, 
        levels=chrLevels, ordered=TRUE))]
    motifData[,motif_id:=factor(motif_id, 
        levels=motifLevels, ordered=TRUE)]
    
    # determine margins
    motifData[,start_margin:=start-margin]
    motifData[,end_margin:=end+margin]
    
    # determine motif center
    motifData[,motif_center:=floor((end_margin-start_margin)/2)+start_margin]
    
    # convert to factors (memory usage)
    fragDt[,seqnames:=as.integer(factor(seqnames, levels=chrLevels, ordered=TRUE))]
    
    nSamples <- length(unique(fragDt$sample))
    
    setorder(motifData, seqnames)
    setorder(fragDt, seqnames)
    motifData[,motif_match_id:=1:nrow(motifData)]
    motifData <- split(motifData, by="seqnames")
    fragDt <- split(fragDt, by="seqnames")
    
    atacInserts <- mapply(function(md,af){
        
        # convert to granges for faster overlapping
        motifMarginRanges <- dtToGr(md, startCol="start_margin", endCol="end_margin")
        atacStartRanges <- dtToGr(af, startCol="start", endCol="start")
        atacEndRanges <- dtToGr(af, startCol="end", endCol="end")
        
        startHits <- findOverlaps(atacStartRanges, 
            motifMarginRanges, type="within") # check if type within faster or slower
        endHits <- findOverlaps(atacEndRanges, motifMarginRanges, type="within") 
        
        # get overlapping insertion sites
        atacStartInserts <- af[queryHits(startHits), c("sample", "start"), with=FALSE]
        atacEndInserts <-af[queryHits(endHits), c("sample", "end"), with=FALSE]
        setnames(atacStartInserts, "start", "insert")
        setnames(atacEndInserts, "end", "insert")
        
        ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
            rbindlist(list(
                md[subjectHits(startHits), c("motif_center", "start", 
                    "end", "motif_id", "motif_match_id")],
                md[subjectHits(endHits), 
                    c("motif_center", "start", "end", 
                        "motif_id", "motif_match_id")])))
        
        # count insertions around motif
        ai[,rel_pos:=abs(insert-motif_center)]
        ai[,type:=fifelse(insert>=start & insert<=end, 1,0)]
        aiMotif <- subset(ai, type==1)
        aiMargin <- subset(ai, type==0)
        
        aiMotif <- aiMotif[,.(pos_count_sample_match=.N), 
            by=.(motif_match_id, rel_pos, sample, motif_id)]
        aiMotif <- aiMotif[,.(pos_count_sample_match=median(pos_count_sample_match)), 
            by=.(motif_match_id, sample, motif_id)]
        aiMotif$rel_pos <- 0
        
        aiMargin <- aiMargin[,.(pos_count_sample_match=.N), 
            by=.(motif_match_id, rel_pos, sample, motif_id)]
        aiMargin[,rel_pos:=rel_pos-min(rel_pos)+1, by=motif_id]
        rbind(aiMargin, aiMotif)
    }, 
        motifData, 
        fragDt, 
        SIMPLIFY=FALSE)
    
    # combine insertion counts across chromosomes
    atacInserts <- rbindlist(atacInserts)
    atacProfiles <- atacInserts[,.(pos_count_sample=sum(pos_count_sample_match)), 
        by=.(motif_id, rel_pos, sample)]
    atacProfiles <- atacProfiles[,pos_count_global:=sum(pos_count_sample), 
        by=.(motif_id, rel_pos)]
    setorder(atacProfiles, motif_id, rel_pos)
    atacProfiles[,w:=smooth(pos_count_global/sum(pos_count_global),
        twiceit=TRUE), by=motif_id]
    atacProfiles <- atacProfiles[,.(w=first(w)), by=.(rel_pos, motif_id)]
    atacInserts <- merge(atacInserts, 
        atacProfiles[,c("rel_pos", "motif_id", "w"), with=FALSE], 
        by.x=c("motif_id", "rel_pos"),
        by.y=c("motif_id", "rel_pos"))
    atacInserts[,score:=w*pos_count_sample_match]
    # discuss
    atacInserts[,pos_count_sample_match_norm:=
            (pos_count_sample_match/sum(pos_count_sample_match))*1e4, by=sample]
    # calculate per sample motif match scores
    matchScores <- atacInserts[,.(score=sum(score)), 
        by=.(motif_match_id, sample)]
    
    # get back original coordinates
    motifData <- rbindlist(motifData)
    matchScores <- cbind(motifData[matchScores$motif_match_id,
        setdiff(colnames(motifData), "score"), with=FALSE], 
        matchScores)
    matchScores[,seqnames:=chrLevels[seqnames]]
    #matchScores[,motif_name:=motifLevels[motif_id]]
    
    return(matchScores)
}

#' Function to get, from the counts-per-match and for a motif of choice, 
#' a matrix of the counts from the best motif match per peak (0 if none) 
#' @param matchScore a data table containing motif ranges, motif_id, sample 
#' and match score (output by `.getInsertProfile`)
#' @param peakRange a GRange object containing the peak ranges
#' @return a list of counts-per-peaks for motif X

#' @TODO: slow, speed up
.getPeakMatchScore <- function (matchScore,
    peakRange,
    seqCol="seqnames",
    ...
    ) {
    if (!is.null(peakRange$score)) peakRange$score <- NULL
    
    peaks <- as.data.table(peakRange)
    peaks$peakID <- seq_len(nrow(peaks))
    allRes <- lapply(split(matchScore, matchScore$motif_id), \(motifDt){
        sampleRes <- sapply(split(motifDt, motifDt$sample), \(dt) {
            motifRange <- dtToGr(dt, seqCol = seqCol)
            hits <- findOverlaps(motifRange, peakRange, type="within")
            overlaps <- cbind(dt[queryHits(hits),],
                peaks[subjectHits(hits), c("peakID")])
            tmp <- overlaps[, .(scores = max(score, na.rm = TRUE)), by = peakID]
            res <- data.table(peakID = seq_len(nrow(peaks)))
            res <- merge(res, tmp, all=TRUE)
            res[is.na(res)] <- 0
            res$scores
        })
    })
}

#' Function to get the background peak sets and strongest motif match scores in 
#' each peak 
#' @param se a SummarizedExperiment object that containing original fragment count
#' on peak level
#' @param matchScore a data table containing motif ranges, motif_id, sample 
#' and match score (output by `.getInsertProfile`)
#' @param genome BSgenome object
#' @param niterations number of background peaks to sample
#' @return a list containing background peak indices and a matrix containing the 
#' strongest motif match score for peaks per sample
.getBackgroundPeaks <- function (se, 
    matchScore, 
    genome,
    niterations=10) {
    se <- addGCBias(se,genome=genome)
    bg <- getBackgroundPeaks(object=se, niterations=niterations)
    if ("strand" %in% colnames(matchScore)) matchScore$strand <- NULL
    maxScorePeaks <- .genomicRangesMapping(refRanges=rowRanges(se), 
        aggregationFun=max,
        assayTable=matchScore, 
        byCols="sample", 
        seqNamesCol="seqnames", 
        scoreCol="score")
    
    return(list(backgroundPeaks=bg, maxScorePeaks=maxScorePeaks))
}

#' get the motif activity score
#' @param backgroundPeaks a data.table containing background peak indices
#' @param peakMatchScore a list of matrix containing match score on peak level 
#' for each motif
#' @param maxScorePeaks a matrix containing strongest motif match score for 
#' peaks per sample
#' @return a list of matrix containing motif activity score on peak level per sample

.getMotifActivityScore <- function (backgroundPeaks, 
    peakMatchScore, 
    maxScorePeaks) {
    # compute the deviations for each motif
    activityScore <- sapply(peakMatchScore, \(score) {
        # find which peaks containing that motif
        idx <- which(rowSums(score)!=0) 
        # compute background deviations
        dev_bg <- t(sapply(seq_len(ncol(backgroundPeaks)), \(x) {
            bg_i <- backgroundPeaks[idx,x]
            maxPeak <- maxScorePeaks[bg_i,]
            # deviations
            (colSums(maxPeak)-mean(colSums(maxPeak)))/mean(colSums(maxPeak))
        }))
        
        y <- colSums(score[idx,])
        dev_motif <- y-mean(y)/mean(y)
        
        z <- (dev_motif-mean(dev_bg))/sd(dev_bg) #or colMeans and colSds?

        
    })
    return(t(activityScore))
}


computeMotifActivityScore <- function (se,
    atacFrag, 
    peakRange,
    motif,
    species,
    minFrag = 30,
    maxFrag = 3000,
    genome,
    niterations=100,
    ...
    ) {
    # standard chromosomes
    peakRange <- .standardChromosomes(peakRange, species = species)
    atacFrag <- lapply(atacFrag, function(dt) {
        gr <- dtToGr(dt)
        gr <- .standardChromosomes(gr, species = species)
        as.data.table(gr)
    })
    
    
    
    # filter too short or too long fragments
    atacFrag <- .filterFrags(atacFrag, min = minFrag, max = maxFrag)
    
    # match seqLevels
    res <- .matchSeqlevels(atacFrag, peakRange)
    atacFrag <- res$atacFrag
    peakRange <- res$ranges
    
    # get motifRanges
    motif_ix <- matchMotifs(motif, 
        peakRange, 
        genome=genome, 
        out = "positions") 
    
    motifData <- lapply(names(motif_ix), function(x) {
        gr <- motif_ix[[x]]
        mcols(gr)$motif_id <- x
        gr
    })
    motifRanges <- plyranges::bind_ranges(motifData)
    fragDt <- rbindlist(atacFrag, idcol="sample")
    matchScore <- .getInsertionProfiles(fragDt,
        motifRanges=motifRanges)
    peakMatchScore <- .getPeakMatchScore(matchScore, peakRange)
    background <- .getBackgroundPeaks(se, matchScore, genome, niterations)
    activityScore <- .getMotifActivityScore(
        background$backgroundPeaks,
        peakMatchScore, 
        background$maxScorePeaks)
}

    
    
    