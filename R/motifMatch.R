suppressPackageStartupMessages({
    source("~/chromVAR/R/getCounts.R")
    library(motifmatchr)
})

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
        atacProfiles[,w:=w/sum(w), by=motif_id]
        
        #atacProfiles[,motif_name:=motifLevels[motif_id]]
    }
    else
    {
        # uniform weighting
        atacProfiles <- atacInserts[,.(w=1/.N), by=.(motif_id, rel_pos_m)]
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
                md[subjectHits(startHits), c("motif_center", "start", 
                    "end", "motif_id", "motif_match_id")],
                md[subjectHits(endHits), c("motif_center", "start", "end", "motif_id", 
                    "motif_match_id")])))
        
        # count insertions around motif
        ai[,rel_pos:=insert-motif_center]
        ai <- ai[,.(pos_count=.N), by=.(motif_match_id, motif_id, sample, rel_pos)]
        
        ai <- merge(ai, atacProfiles[,c("rel_pos", "motif_id", "w"), with=FALSE], 
            by.x=c("motif_id","rel_pos"), 
            by.y=c("motif_id","rel_pos"), all.x=TRUE)
        ai[,score:=w*pos_count]
        as <- ai[,.(score=sum(score), 
            tot_count=sum(pos_count)), by=.(motif_id, sample)]
        gc()
        # get ai[,rel_pos_m:=fifelse(insert<start, insert-start, 0)]
        #     ai[,rel_pos_m:=fifelse(insert>end, insert-end, rel_pos_m)]
        as
    }, 
        motifData, 
        atacFrag, 
        SIMPLIFY=FALSE)
    
    motifScores <- rbindlist(motifScores)
    motifScores <- motifScores[,.(score=sum(score),
        tot_inserts=sum(tot_count)), 
        by=.(motif_id, sample)]
    
    if(libNorm)
    {
        motifScores[,score:=score/tot_inserts]
    }
    
    #motifScores[,nmotif_name:=motifLevels[motif_id]]
    return(list(ms=motifScores, ap=atacProfiles))
}



.getInserts <- function(atacFrag, 
    coords,
    margin=10,
    aggFun=sum,
    #minWidth=30,
    #maxWidth=2000,
    getType=TRUE,
    chunk=TRUE){
    
    cDt <- as.data.table(coords)
    cDt$id <- seq_len(nrow(cDt))
    #setnames(cDt, "seqnames", "chr")
    cdl <- unique(cDt$seqnames)
    idl <- unique(cDt$id)
    
    # convert to factors (memory usage)
    cDt[,seqnames:=as.integer(factor(seqnames, levels=cdl, ordered=TRUE))]
    cDt[,id:=as.integer(factor(id, levels=idl, ordered=TRUE))]
    
    # determine margins
    cDt[,start_margin:=start-margin]
    cDt[,end_margin:=end+margin]
    
    # determine motif center
    cDt[,center:=floor((end_margin-start_margin)/2)+start_margin]
    
    # convert to factors (memory usage)
    fDt <- copy(atacFrag) 
    fDt <- subset(fDt, seqnames %in% cdl) # frags on different chromosome don't need to be considered
    fDt[,seqnames:=as.integer(factor(seqnames, levels=cdl, ordered=TRUE))]
    
    if(chunk){
        setorder(cDt, seqnames)
        setorder(fDt, seqnames)
        cDts <- split(cDt, by="seqnames")
        fDts <- split(fDt, by="seqnames")
        
        ins <- mapply(function(cDt, fDt){
            
            # convert to granges for faster overlapping
            # coordinate margin ranges
            cr <- dtToGr(cDt, startCol="start_margin", endCol="end_margin", seqCol="seqnames")
            # start ranges
            sr <- dtToGr(fDt, startCol="start", endCol="start", seqCol="seqnames")
            # end ranges
            er <- dtToGr(fDt, startCol="end", endCol="end", seqCol="seqnames")
            
            sHits <- findOverlaps(sr, cr, type="within") 
            eHits <- findOverlaps(er, cr, type="within") 
            
            # get overlapping insertion sites
            sIns <- fDt[queryHits(sHits), c("sample", "start"), with=FALSE]
            eIns <- fDt[queryHits(eHits), c("sample", "end"), with=FALSE]
            setnames(sIns, "start", "insert")
            setnames(eIns, "end", "insert")
            
            ai <- cbind(rbindlist(list(sIns, eIns)),
                rbindlist(list(cDt[subjectHits(sHits), c("center", "id")],
                    cDt[subjectHits(eHits), c("center", "id")])))
            
            # count insertions around motif
            ai[,rel_pos:=abs(insert-center)]}, 
            cDts, 
            fDts, 
            SIMPLIFY=FALSE)
        
        # combine insertion counts across chromosomes
        names(ins) <- cdl
        ins <- rbindlist(ins, idcol="seqnames")
    }
    else
    {
        # convert to granges for faster overlapping
        cr <- dtToGr(cDt, startCol="start_margin", endCol="end_margin", seqCol="seqnames")
        sr <- dtToGr(fDt, startCol="start", endCol="start", seqCol="seqnames")
        er <- dtToGr(fDt, startCol="end", endCol="end", seqCol="seqnames")
        
        sHits <- findOverlaps(sr, cr, type="within") # check if type within faster or slower
        eHits <- findOverlaps(er, cr, type="within") 
        
        # get overlapping insertion sites
        sIns <- fDt[queryHits(sHits), c("sample", "start"), with=FALSE]
        eIns <- fDt[queryHits(eHits), c("sample", "end"), with=FALSE]
        setnames(sIns, "start", "insert")
        setnames(eIns, "end", "insert")
        
        ins <- cbind(rbindlist(list(sIns, eIns)),
            rbindlist(list(cDt[subjectHits(sHits), c("center", "id", "seqnames")],
                cDt[subjectHits(eHits), c("center", "id", "seqnames")])))
        
        # count insertions around motif
        ins[,rel_pos:=abs(insert-center)]
        
        # convert chr naming back 
        ins[,seqnames:=factor(seqnames, levels=1:length(cdl), labels=cdl)]
    }
    
    if(getType)
    {
        # width of side
        mw <- median(width(coords))
        
        ins[,start:=center-floor(mw/2)]
        ins[,end:=center+floor(mw/2)]
        ins[,type:=fifelse(insert>=start & insert<=end, "within", "margin")]
    }
    
    return(ins)
}


.getInsertionBgProfiles <- function(atacFrag, 
    coords,
    #minWidth=30,
    #maxWidth=2000,
    chunk=TRUE){
    
    # prep motif data
    coordData <- as.data.table(coords)
    setnames(coordData, "seqnames", "chr")
    chrLevels <- unique(coordData$chr)
    motifLevels <- unique(coordData$id)
    
    # convert to factors (memory usage)
    
    coordData[,chr:=as.integer(factor(chr, 
        levels=chrLevels, ordered=TRUE))]
    coordData[,id:=as.integer(factor(id, levels=motifLevels, ordered=TRUE))]
    
    # determine motif center
    coordData[,center:=floor((end-start)/2)+start]
    
    # convert to factors (memory usage)
    if("seqnames" %in% colnames(atacFrag)){
        setnames(atacFrag, "seqnames", "chr")
    }
    atacFrag <- copy(atacFrag) #TODO: take out that copy 
    atacFrag[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
    
    # if("seqnames" %in% colnames(atacFrag)){
    #   setnames(atacFrag, "seqnames", "chr")
    # }
    
    nSamples <- length(unique(atacFrag$sample))
    
    setorder(coordData, chr)
    setorder(atacFrag, chr)
    coordData[,coord_id:=1:nrow(coordData)]
    coordData <- split(coordData, by="chr")
    atacFrag <- split(atacFrag, by="chr")
    
    atacInserts <- mapply(function(md,af){
        
        #md[,motif_match_id:=1:nrow(md)]
        
        # convert to granges for faster overlapping
        coordRanges <- dtToGr(md, startCol="start", endCol="end", seqCol="chr")
        atacStartRanges <- dtToGr(af, startCol="start", endCol="start", seqCol="chr")
        atacEndRanges <- dtToGr(af, startCol="end", endCol="end", seqCol="chr")
        
        startHits <- findOverlaps(atacStartRanges, 
            coordRanges, type="within") # check if type within faster or slower
        endHits <- findOverlaps(atacEndRanges, coordRanges, type="within") 
        
        # get overlapping insertion sites
        atacStartInserts <- af[queryHits(startHits), c("sample", "start"), with=FALSE]
        atacEndInserts <-af[queryHits(endHits), c("sample", "end"), with=FALSE]
        setnames(atacStartInserts, "start", "insert")
        setnames(atacEndInserts, "end", "insert")
        
        ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
            rbindlist(list(
                md[subjectHits(startHits), c("center", "start", 
                    "end", "id", "coord_id")],
                md[subjectHits(endHits), c("center", "start", "end", "id", 
                    "coord_id")])))
        
        # count insertions around motif
        ai[,rel_pos:=abs(insert-center)]
        
        ai <- ai[,.(pos_count_sample=.N),  by=.(coord_id, rel_pos, sample, id)]
        ai
    }, 
        coordData, 
        atacFrag, 
        SIMPLIFY=FALSE)
    
    # combine insertion counts across chromosomes
    atacInserts <- rbindlist(atacInserts)
    atacProfiles <- atacInserts[,.(pos_count_sample=sum(pos_count_sample)), by=.(id, rel_pos, sample)]
    atacProfiles <- atacProfiles[,pos_count_global:=sum(pos_count_sample), by=.(id, rel_pos)]
    setorder(atacProfiles, id, rel_pos)
    atacProfiles[,w:=smooth(pos_count_global/sum(pos_count_global),
        twiceit=TRUE), by=id]
    atacProfiles <- atacProfiles[,.(w=first(w)), by=.(rel_pos, id)]
    atacInserts <- merge(atacInserts, 
        atacProfiles[,c("rel_pos", "id", "w"), with=FALSE], 
        by.x=c("id", "rel_pos"),
        by.y=c("id", "rel_pos"))
    #atacInserts[,pos_count_sample_norm:=(pos_count_sample/sum(pos_count_sample))*1e4, by=sample]
    atacInserts[,score:=w*pos_count_sample]
    
    
    # calculate per sample motif match scores
    scores <- atacInserts[,.(score=sum(score)), by=.(coord_id, sample)]
    #scores[,score:=sum(score),by=id]
    
    # get back original coordinates
    coordData <- rbindlist(coordData)
    scores <- cbind(coordData[scores$coord_id,setdiff(colnames(coordData), "score"), with=FALSE], 
        scores)
    scores[,chr:=chrLevels[seqnames]]
    #scores[,motif_name:=motifLevels[motif_id]]
    
    return(scores)
}


#' Mapping & aggregating modalities with genomic coordinates to reference 
#' coordinates.
#'
#' Internal convencience function for mapping different modality scores with 
#' cell type and further labels such as tfs to reference coordinates. The resulting
#' table will have dimension ref coords x byCols (or ref coord x byCols[1] x byCols[2]).
#' ByCols can be for instance cell type labels and/or transcription factor names.
#'
#'@name .genomicRangesMapping
#'@param refRanges GRanges object with reference coordinates
#'@param assayTable modality table to be mapped to the reference coordinates. 
#'Needs to containing genomic coordinates (see args: seqNamesCol, startCol, endCol), and byCols.
#'@param byCols will be the columns /depths of the resulting matrix with dimension 
#' ref coords x byCols (or ref coord x byCols[1] x byCols[2]). 
#' ByCols can be for instance cell type labels and/or transcription factor names.
#' @param seqNamesCol name of the column in motif, atac and chIP-seq data.tables containing
#' the sequence information.
#' @param startCol name of the column in motif, atac and chIP-seq data.tables containing
#' the start coordinate.
#' @param endCol name of the column in motif, atac and chIP-seq data.tables containing
#' the end coordinate.
#' @param scoreCol name of the score column (e.g. motif matching scores, atac fragment counts)
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate
#' if multiple rows of the assayTable overlap a reference coordinate.
#' @param BPPARAM BiocParallel argument either SerialParam() or MulticoreParam(workers=n)
#' @export
.genomicRangesMapping <- function(refRanges, 
    assayTable,
    byCols=c("tf_uniprot_id",
        "cell_type_id"),
    seqNamesCol="chr",
    startCol="start",
    endCol="end",
    scoreCol=NULL,
    calledInCol=NULL,
    aggregationFun=NULL,
    minoverlap=1,
    isReplicated=FALSE,
    BPPARAM=SerialParam()){
    # TODO: - add warning for integer overflows - data.table size
    assayTable <- copy(assayTable)  
    
    # attribute generic names to dimensionalities
    if(length(byCols)==2)
    {
        setnames(assayTable, byCols, c("col_depth", "col_width"))
        byCols <- c("col_depth", "col_width")
        multiTf <- TRUE
    }
    else
    {
        setnames(assayTable, byCols, c("col_width"))
        byCols <- c("col_width")
        multiTf <- FALSE
    }
    
    # get dimensions of tables
    nRefs <- length(refRanges)
    nColsWidth <- length(unique(assayTable$col_width))
    
    # convert to integer for speed-up
    levels <- unique(assayTable$col_width)
    assayTable[,col_width:=as.integer(factor(assayTable$col_width, levels=levels))]
    
    # convert to GRanges for faster overlap finding
    assayRanges <- makeGRangesFromDataFrame(as.data.frame(assayTable),
        keep.extra.columns=TRUE,
        seqnames.field=seqNamesCol,
        start.field=startCol,
        end.field=endCol,
        ignore.strand=TRUE)
    
    # find overlaps with ref. coordinates
    overlapTable <- as.data.table(findOverlaps(assayRanges, refRanges,
        type="within",
        minoverlap=minoverlap,
        ignore.strand=TRUE))
    rm(refRanges, assayRanges)
    
    
    # retrieve tf and cell type ids
    overlapTable <- cbind(overlapTable$subjectHits,
        assayTable[overlapTable$queryHits,
            c(byCols, scoreCol, calledInCol), with=FALSE])
    
    if(multiTf)
    {
        setkey(overlapTable, V1, col_width)
        if(!is.null(calledInCol)) setnames(overlapTable, calledInCol, "calledInCol")
        if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
        overlapTable <- split(overlapTable, by=c("col_depth"))
        
        #TODO:  BPPARAM= serialParam or multiCoreParam() as function arguments
        #TODO:  bplapply(fls[1:3], FUN, BPPARAM = MulticoreParam(), param = param)
        # overlap with ref. coordinates
        overlapTable <- BiocParallel::bplapply(overlapTable, function(table){
            
            if(is.null(aggregationFun))
            {
                # get number of max callers
                #table[, n_max:=length(unique(calledInCol)), 
                #        by=c("col_width")]
                
                if(isReplicated)
                {
                    table[,rep:=tstrsplit(calledInCol, split="-", keep=2)]
                    table[, n_max:=length(unique(calledInCol)), 
                        by=c("col_width", "rep")]
                    
                    table <- table[,.(value=.chIPlabelRule(calledInCol,
                        data.table::first(n_max))),
                        by=c("V1", "col_width", "rep")]
                    table <- table[,.(value=.chIPReplicateLabelRule(value)),
                        by=c("V1", "col_width")]
                }
                else{
                    table[,n_max:=length(unique(calledInCol)), 
                        by=c("col_width")]
                    
                    table <- table[,.(value=.chIPlabelRule(calledInCol,
                        data.table::first(n_max))),
                        by=c("V1", "col_width")]
                }
            }
            else
            {
                table <- table[,.(value=aggregationFun(scoreCol)),
                    by=c("V1", "col_width")]
            }
            
            # one would need to construct a second table here for the neg labels
            
            # convert to sparse matrix
            table <- sparseMatrix(i=table$V1, 
                j=as.integer(table$col_width), # what if non-numeric cell-id
                dims=c(nRefs, nColsWidth),
                x=table$value)
            colnames(table) <- levels
            return(table)},
            BPPARAM=BPPARAM)
        
    }
    else
    {
        # setkeys for speed-up
        overlapTable[,V1:=as.integer(V1)]
        setkey(overlapTable, col_width, V1)
        
        # overlap with ref. coordinates
        if(is.null(aggregationFun)) error("Aggregation function needs to be defined")
        
        # why an lapply 05.04.24
        #overlapTable <- overlapTable[,.("scoreCol"=lapply(.SD, aggregationFun, na.rm=TRUE)),
        #                              by=c("col_width", "V1"),
        #                             .SDcols=scoreCol]
        setnames(overlapTable, scoreCol, "scoreCol")
        overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
            by=c("col_width", "V1")]
        
        # convert to sparse matrix
        overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1, 
            j=overlapTable$col_width,
            dims=c(nRefs, nColsWidth),
            x=overlapTable$scoreCol)
        
        colnames(overlapTable) <- levels
    }
    
    return(overlapTable)
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
    if ("strand" %in% colnames(matchScore)) matchScore$strand <- NULL
    res <- .genomicRangesMapping(peakRange, matchScore, 
        byCols=c("motif_id", "sample"), seqNamesCol = seqCol,
        scoreCol="score", aggregationFun = max)
    lapply(res, \(x) as.matrix(x))

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
# .getBackgroundPeaks <- function (se, 
#     matchScore, 
#     genome,
#     niterations=10) {
#     se <- addGCBias(se,genome=genome)
#     bg <- getBackgroundPeaks(object=se, niterations=niterations)
#     if ("strand" %in% colnames(matchScore)) matchScore$strand <- NULL
#     maxScoreSample <- .genomicRangesMapping(refRanges=rowRanges(se), 
#         aggregationFun=max,
#         assayTable=matchScore, 
#         byCols=c("sample","motif_id"), 
#         seqNamesCol="seqnames", 
#         scoreCol="score")
#     
#     maxIds <- lapply(maxScoreSample, \(x) {
#         dt <- as(x,"TsparseMatrix")
#         #rownames(dt) <- seq_len(nrow(dt))
#         idDt <- data.table(peakID=dt@i+1, 
#             motif_id=colnames(x)[dt@j+1], score=dt@x)
#     })
#     maxId <- rbindlist(maxIds, idcol="sample")
#     maxId[,mean_score:=sum(score)/length(meanScorePeak),by=.(peakID, motif_id)]
#     maxId[,max_score:=max(mean_score),by=peakID]
#     maxDt <- maxId[maxId$max_score==maxId$mean_score,]
#     resDt <- data.table(peakID=seq_len(nrow(se)))
#     resDt <- merge(resDt, maxDt, by = "peakID", all.x=TRUE)
#     maxScorePeaks <- data.table::dcast(resDt, peakID~sample, fill=0, 
#          fun.aggregate = mean, value.var = "score", drop=FALSE)
#     #rownames(maxScorePeaks) <- maxScorePeaks$peakID
#     maxScorePeaks$`NA` <- NULL
#     maxScorePeaks$peakID <- NULL
#     
#     
#     return(list(backgroundPeaks=bg, maxScorePeaks=maxScorePeaks))
# }


.getBackgroundPeaks <- function(se, genome, niterations=100) {
    se <- addGCBias(se,genome=genome)
    bg <- getBackgroundPeaks(object=se, niterations=niterations)
    #bgProfile <- .getInsertionBgProfiles(fragDt, peakRange)
}

#' get the motif activity score
#' @param backgroundPeaks a data.table containing background peak indices
#' @param peakMatchScore a list of matrix containing match score on peak level 
#' for each motif
#' @param maxScorePeaks a matrix containing strongest motif match score for 
#' peaks per sample
#' @return a list of matrix containing motif activity score on peak level per sample

# .getMotifActivityScore <- function (backgroundPeaks, 
#     peakMatchScore, 
#     maxScorePeaks) {
#     # compute the deviations for each motif
#     activityScore <- lapply(peakMatchScore, \(score) {
#         # find which peaks containing that motif
#         idx <- which(rowSums(score)!=0) 
#         # compute background deviations
#         colSum_bg <- sapply(seq_len(ncol(backgroundPeaks)), \(x) {
#             bg_i <- backgroundPeaks[idx,x]
#             maxPeak <- maxScorePeaks[bg_i,]
#             # deviations
#             sMaxPeak <- colSums(maxScorePeaks[bg_i,])
#             #(sMaxPeak-mean(sMaxPeak))/mean(sMaxPeak)
#         })
#         dev_bg <- (colSum_bg-colMeans(colSum_bg))/colMeans(colSum_bg)
#         y <- colSums(score[idx,])
#         dev_motif <- (y-mean(y))/mean(y)
#         
#         z <- (dev_motif-rowMeans(dev_bg))/rowSds(dev_bg) 
#         list(z=z, sumScore_motif=y, dev_bg=dev_bg, sumScore_bg=colSum_bg)
#     })
#     return(activityScore)
# }
.getMotifActivityScore <- function (backgroundPeaks, 
    peakMatchScore, 
    bgScore) {
    
    # compute the deviations for each motif
    activityScore <- lapply(peakMatchScore, function(score) {
        # find which peaks containing that motif
        idx <- which(rowSums(score)!=0) 
        # compute background deviations
        bg_scores <- sapply(seq_len(ncol(backgroundPeaks)), function(x) {
            bg_i <- backgroundPeaks[idx,x]
            peak_set <- bgScore[bg_i,]
            # deviations
            bg_score <- colSums(peak_set)
        })
        
        motif_score <- colSums(score[idx,])
        bg_score <- rowMeans(bg_scores)
        bg_motif_score <- motif_score - rowMeans(bg_scores)
        z <- bg_motif_score/rowSds(bg_scores)
        
        dev_motif <- (motif_score-mean(motif_score))/mean(motif_score)
        dev_bg <- (bg_scores-colMeans(bg_scores))/colMeans(bg_scores)
        dev_z <- (dev_motif - rowMeans(dev_bg))/rowSds(dev_bg)
        list(motif_score=motif_score,bg_motif_score=bg_motif_score,z=z,
            dev_motif, dev_bg=dev_bg,dev_z=dev_z, bg_score=bg_score)

    })
    return(activityScore)
}


computeMotifActivityScore <- function (se,
    atacFrag, 
    #peakRange,
    motif,
    species,
    minFrag = 30,
    maxFrag = 3000,
    genome,
    niterations=100,
    ...
    ) {
    # standard chromosomes
    #rowData(se)$peakID <- seq_len(nrow(se))
    
    # already resized peaks
    se <- filterPeaks(se, non_overlapping = TRUE)
    peakRange <- rowRanges(se)

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
    motif_ix <- motif_ix[!vapply(motif_ix, \(.) length(.)==0, logical(1))]
    motifData <- lapply(names(motif_ix), function(x) {
        gr <- motif_ix[[x]]
        mcols(gr)$motif_id <- x
        gr
    })
    motifRanges <- plyranges::bind_ranges(motifData)
    
    # calculate activity score
    fragDt <- rbindlist(atacFrag, idcol="sample")
    matchScore <- .getInsertionProfiles(fragDt, motifRanges=motifRanges)
    #matchScore <- insertProfile$matchScore
    #globalProfile <- insertProfile$atacProfiles
    peakMatchScore <- .getPeakMatchScore(matchScore, peakRange)
    backgroundPeaks <- .getBackgroundPeaks(se, genome=genome, 
        niterations=niterations)
    mcols(peakRange)$id <- 1
    bgProfile <- .getInsertionBgProfiles(fragDt, peakRange)
    bgScore <- as.matrix(sparseMatrix(i=bgProfile$coord_id, 
        j=as.integer(factor(bgProfile$sample)), x=as.numeric(bgProfile$score),
        dims=dim(se), dimnames=list(rownames(se),colnames(se))))
    
    
    #peakInserts <- .getInserts(fragDt, peakRange)
    #bgDt <- subset(peakInserts, type=="within")[,.(sum_inserts=.N), by=.(id,sample)]
    #bgDt[,norm_inserts:=sum_inserts/sum(sum_inserts)*1e4, by=sample]
    #bgDt[,score:=norm_inserts/width(peakRange)[1]]
    #bgScore <- as.matrix(sparseMatrix(i=bgDt$id,j=as.integer(as.factor(bgDt$sample)),
    #    x=bgDt$score,
    #    dims=dim(se),dimnames=list(rownames(se),colnames(se))))
    
    # deviations, z-score, colSums of motifs and (mean deviations of) backgrounds
    activityScore <- .getMotifActivityScore(
        backgroundPeaks,
        peakMatchScore, 
        bgScore)
    
    motif_score <- t(sapply(activityScore, \(x) x$motif_score))
    bg_score <- t(sapply(activityScore, \(x) x$bg_score))
    return(list(motif_score=motif_score, bg_score=bg_score))
    # group_id <- c(0,0,1,1)
    # lms <- lapply(seq_len(nrow(motif_score)), function(i) {
    #     lm(motif_score[i,] ~ bg_score[i,] + group_id)
    # })
    # names(lms) <- rownames(motif_score)
    # coef <- sapply(lms, function(x) x$coefficients[3])
    
}
# s <- proc.time()
# z <- computeMotifActivityScore(se,
#     atacFrag,
#     motif,
#     species="Homo sapiens",
#     genome=BSgenome.Hsapiens.UCSC.hg38)
# e <- proc.time()
# e-s

