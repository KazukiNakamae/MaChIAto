GetAbsolutePath <- function(relative.path) {
  # Get an absolute path from a given relative.path
  #
  # Args:
  #   relative.path: the charators of relative path
  #
  # Returns:
  #   the charators of absolute path
  return(normalizePath(relative.path))
}

ReverseComplementSeq <- function(seq){
  return(as.character(reverseComplement(DNAString(seq))))
}

FactorToNumeric <- function(factor){
  as.numeric(as.character(factor))
}

SafeVariantCounts <- function(crispr.set){
  if(!is.null(crispr.set)){
    return(variantCounts(crispr.set))
  }else{
    return(NULL)
  }
}

#'@title Make Wt Crispr Set
#'@description Make Wt Crispr Set from reference sequence and mapped bam file 
#'@param untreated.seq untreated sequence flanking out-out indicators on MaChIAto
#'@param index.name name of index
#'@param wt.seq untreated sequence in reference index
#'@param protospacer.seq protospacer sequence
#'@param bam.fn file name (including path) of mapped bam file
#'@param names data names
#'@author Kazuki Nakamae
#'@return The list including mutation infomation classified by mutant type
#'@rdname MakeWtCrisprSet
MakeWtCrisprSet <- function(untreated.seq, index.name, wt.seq, protospacer.seq, bam.fn, names){

  wt.crispr.set <- NULL
  bf <- open(BamFile(bam.fn))
  bam.len <- length(scanBam(bf)[[1]]$qname)
  close(bf)
  if(bam.len == 0){
    # bam file data dosen't exist.
    return(wt.crispr.set)
  }

  tryCatch({
    options(warn=-1)
    untreated.wt.pos <- as.numeric(regexpr(untreated.seq, wt.seq))
    protospacer.untreated.pos <- as.numeric(regexpr(protospacer.seq, untreated.seq))
    wt.guide <- GRanges(seqnames = index.name, strand = c("+"),
                    ranges = IRanges(start = c(untreated.wt.pos)
                      , end = c(untreated.wt.pos + nchar(untreated.seq) - 1 )))
    wt.reference = untreated.seq #guide sequence
    wt.crispr.set <- readsToTarget(bam.fn
                                , wt.guide
                                , reference = wt.reference
                                , target.loc = protospacer.untreated.pos + (nchar(protospacer.seq) - 3) -1
                                , names = names
                                , upstream.snv = 35
                                , downstream.snv = 35
                                , chimeras = "ignore"
                                )
    return(wt.crispr.set)
  }, 
    error = function(e) {
      message(paste("Alignment Process is stopped due to insufficient data. Export NULL in ", names, sep=""))
  },
    silent = TRUE
  )
  options(warn=0)
  return(wt.crispr.set)
}

#'@title Make Knock-in Crispr Set
#'@description Make Knock-in Crispr Set from reference sequence and mapped bam file 
#'@param knockin.seq knock-in sequence flanking out-out indicators on MaChIAto
#'@param index.name name of index
#'@param ki.seq knock-in sequence in reference index
#'@param donor.seq donor sequence
#'@param bam.fn file name (including path) of mapped bam file
#'@param names data names
#'@author Kazuki Nakamae
#'@return The list including mutation infomation classified by mutant type
#'@rdname MakeKiCrisprSet
MakeKiCrisprSet <- function(
  outleft.indicator.seq, outright.indicator.seq
  , inleft.indicator.seq, inright.indicator.seq
  , left.homology.seq, right.homology.seq
  , index.name, precise.knockin.seq, donor.seq, bam.fn, names, ki.mode){

  ki.crispr.set <- NULL
  bf <- open(BamFile(bam.fn))
  bam.len <- length(scanBam(bf)[[1]]$qname)
  close(bf)
  if(bam.len == 0){
    # bam file data dosen't exist.
    return(ki.crispr.set)
  }

  tryCatch({
    options(warn=-1)
    if(ki.mode == "left.precise"){
      # make crispr set
      ref.start = as.numeric(regexpr(outleft.indicator.seq, precise.knockin.seq))
      ref.end = as.numeric(regexpr(inright.indicator.seq, precise.knockin.seq) + nchar(inright.indicator.seq))
      ki.reference = substr(precise.knockin.seq, ref.start, ref.end)
      junction.loc = as.numeric(regexpr(inleft.indicator.seq, ki.reference)) - 1
      # remain.donor.len = nchar(ki.reference) - as.numeric(regexpr(donor.seq, ki.reference))
      ki.guide <- GRanges(seqnames = index.name, strand = c("+")
        ,ranges = IRanges(start = c(ref.start)
        , end = c(ref.end))
      )
      ki.crispr.set <- readsToTarget(bam.fn
                                  , ki.guide
                                  , reference = ki.reference
                                  , target.loc = junction.loc
                                  , names = names
                                  , match.label = "perfect_left_junction"
                                  , upstream.snv = nchar(left.homology.seq)
                                  # , downstream.snv = remain.donor.len
                                  , downstream.snv = 0
                                  , chimeras = "ignore")
      return(ki.crispr.set)
    }else if(ki.mode == "left.imprecise"){
      # make nhej donor
      nhej.donor.seq <- paste(left.extra.seq, left.homology.seq, donor.seq, right.homology.seq, right.extra.seq, sep="")
      precise.knockin.seq <- sub(donor.seq, nhej.donor.seq, precise.knockin.seq)
      donor.seq <- nhej.donor.seq
      # make crispr set
      ref.start = as.numeric(regexpr(outleft.indicator.seq, precise.knockin.seq))
      ref.end = as.numeric(regexpr(inright.indicator.seq, precise.knockin.seq) + nchar(inright.indicator.seq))
      ki.reference = substr(precise.knockin.seq, ref.start, ref.end)
      junction.loc = as.numeric(regexpr(inleft.indicator.seq, ki.reference)) - 1
      remain.donor.len = nchar(ki.reference) - as.numeric(regexpr(donor.seq, ki.reference))
      ki.guide <- GRanges(seqnames = index.name, strand = c("+")
        ,ranges = IRanges(start = c(ref.start)
        , end = c(ref.end))
      )
      ki.crispr.set <- readsToTarget(bam.fn
                                  , ki.guide
                                  , reference = ki.reference
                                  , target.loc = junction.loc
                                  , names = names
                                  , match.label = "duplicated_homology_left_junction"
                                  , upstream.snv = nchar(left.homology.seq) * 2 + nchar(left.extra.seq)
                                  , downstream.snv = remain.donor.len
                                  , chimeras = "ignore")
      return(ki.crispr.set)
    }else if(ki.mode == "rc.left.imprecise"){
      # make reverse complement nhej donor
      nhej.donor.seq <- paste(left.extra.seq, left.homology.seq, donor.seq, right.homology.seq, right.extra.seq, sep="")
      nhej.rc.donor.seq <- ReverseComplementSeq(nhej.donor.seq)
      precise.knockin.seq <- sub(donor.seq, nhej.rc.donor.seq, precise.knockin.seq)
      donor.seq <- nhej.rc.donor.seq
      # make reverse complement indicator
      inleft.rc.indicator.seq <- ReverseComplementSeq(inleft.indicator.seq)
      inright.rc.indicator.seq <- ReverseComplementSeq(inright.indicator.seq)
      # make crispr set
      ref.start = as.numeric(regexpr(outleft.indicator.seq, precise.knockin.seq))
      ref.end = as.numeric(regexpr(inleft.rc.indicator.seq, precise.knockin.seq) + nchar(inleft.rc.indicator.seq))
      ki.reference = substr(precise.knockin.seq, ref.start, ref.end)
      junction.loc = as.numeric(regexpr(inright.rc.indicator.seq, ki.reference)) - 1
      remain.donor.len = nchar(ki.reference) - as.numeric(regexpr(donor.seq, ki.reference))
      ki.guide <- GRanges(seqnames = index.name, strand = c("+")
        ,ranges = IRanges(start = c(ref.start)
        , end = c(ref.end))
      )
      ki.crispr.set <- readsToTarget(bam.fn
                                  , ki.guide
                                  , reference = ki.reference
                                  , target.loc = junction.loc
                                  , names = names
                                  , match.label = "reverse_knockin_left_junction"
                                  , upstream.snv = nchar(left.homology.seq) + nchar(right.homology.seq) + nchar(right.extra.seq)
                                  , downstream.snv = remain.donor.len
                                  , chimeras = "ignore")
      # browser()
      return(ki.crispr.set)
    }else if(ki.mode == "right.precise"){
      # make crispr set
      ref.start = as.numeric(regexpr(inleft.indicator.seq, precise.knockin.seq)) + nchar(inleft.indicator.seq)
      ref.end = as.numeric(regexpr(outright.indicator.seq, precise.knockin.seq))
      ki.reference = substr(precise.knockin.seq, ref.start, ref.end)
      junction.loc = as.numeric(regexpr(inright.indicator.seq, ki.reference)) + nchar(inright.indicator.seq) - 1
      # remain.donor.len = junction.loc
      ki.guide <- GRanges(seqnames = index.name, strand = c("+")
        ,ranges = IRanges(start = c(ref.start)
        , end = c(ref.end))
      )
      ki.crispr.set <- readsToTarget(bam.fn
                                  , ki.guide
                                  , reference = ki.reference
                                  , target.loc = junction.loc
                                  , names = names
                                  , match.label = "perfect_right_junction"
                                  # , upstream.snv = remain.donor.len
                                  , upstream.snv = 0
                                  , downstream.snv = nchar(right.homology.seq)
                                  , chimeras = "ignore")
      return(ki.crispr.set)
    }else if(ki.mode == "right.imprecise"){
      # make nhej donor
      nhej.donor.seq <- paste(left.extra.seq, left.homology.seq, donor.seq, right.homology.seq, right.extra.seq, sep="")
      precise.knockin.seq <- sub(donor.seq, nhej.donor.seq, precise.knockin.seq)
      donor.seq <- nhej.donor.seq
      # make crispr set
      ref.start = as.numeric(regexpr(inleft.indicator.seq, precise.knockin.seq)) + nchar(inleft.indicator.seq)
      ref.end = as.numeric(regexpr(outright.indicator.seq, precise.knockin.seq))
      ki.reference = substr(precise.knockin.seq, ref.start, ref.end)
      junction.loc = as.numeric(regexpr(inright.indicator.seq, ki.reference)) + nchar(inright.indicator.seq) - 1
      remain.donor.len = junction.loc
      ki.guide <- GRanges(seqnames = index.name, strand = c("+")
        ,ranges = IRanges(start = c(ref.start)
        , end = c(ref.end))
      )
      ki.crispr.set <- readsToTarget(bam.fn
                                  , ki.guide
                                  , reference = ki.reference
                                  , target.loc = junction.loc
                                  , names = names
                                  , match.label = "duplicated_homology_right_junction"
                                  , upstream.snv = remain.donor.len
                                  , downstream.snv = nchar(right.homology.seq) * 2 + nchar(right.extra.seq)
                                  , chimeras = "ignore")
      return(ki.crispr.set)
    }else if(ki.mode == "rc.right.imprecise"){
      # make reverse complement nhej donor
      nhej.donor.seq <- paste(left.extra.seq, left.homology.seq, donor.seq, right.homology.seq, right.extra.seq, sep="")
      nhej.rc.donor.seq <- ReverseComplementSeq(nhej.donor.seq)
      precise.knockin.seq <- sub(donor.seq, nhej.rc.donor.seq, precise.knockin.seq)
      donor.seq <- nhej.rc.donor.seq
      # make reverse complement indicator
      inleft.rc.indicator.seq <- ReverseComplementSeq(inleft.indicator.seq)
      inright.rc.indicator.seq <- ReverseComplementSeq(inright.indicator.seq)
      # make crispr set
      ref.start = as.numeric(regexpr(inright.rc.indicator.seq, precise.knockin.seq)) + nchar(inright.rc.indicator.seq)
      ref.end = as.numeric(regexpr(outright.indicator.seq, precise.knockin.seq))
      ki.reference = substr(precise.knockin.seq, ref.start, ref.end)
      junction.loc = as.numeric(regexpr(inleft.rc.indicator.seq, ki.reference)) + nchar(inleft.rc.indicator.seq) - 1
      remain.donor.len = junction.loc
      ki.guide <- GRanges(seqnames = index.name, strand = c("+")
        ,ranges = IRanges(start = c(ref.start)
        , end = c(ref.end))
      )
      ki.crispr.set <- readsToTarget(bam.fn
                                  , ki.guide
                                  , reference = ki.reference
                                  , target.loc = junction.loc
                                  , names = names
                                  , match.label = "reverse_knockin_right_junction"
                                  , upstream.snv = remain.donor.len
                                  , downstream.snv = nchar(left.homology.seq) + nchar(right.homology.seq) + nchar(left.extra.seq)
                                  , chimeras = "ignore")
      return(ki.crispr.set)
    }else{
      stop("Undefined ki.mode.\n")
    }
  }, 
    error = function(e) {
      message(paste("Alignment Process is stopped due to insufficient data. Export NULL in ", names, sep=""))
  },
    silent = TRUE
  )
  options(warn=0)
  return(ki.crispr.set)
}

#'@title Make mutant list 
#'@description Make mutant list from mutant label eg "-1:2D"
#'@param mut.label mutant label vector
#'@author Kazuki Nakamae
#'@return The list including mutation infomation classified by mutant type
#'@rdname MakeMutList
MakeMutList <- function(mut.label){

  # Convert SNV:xB, yB -> SNV:xB, SNV:yB
  for(ind in 1:length(mut.label)){
    which.muttype <- !is.na(charmatch(c("D", "I", "SNV", "M"), mut.label[ind]))
    has.muttype <- sum(which.muttype)
    if((has.muttype < 1) && ind > 1){
      if(pre.muttype.label %in% c("I", "D", "M")){
        mut.label[ind] <- paste(mut.label[ind], pre.muttype.label, sep="")
      }else{
        mut.label[ind] <- paste(pre.muttype.label, ":", mut.label[ind], sep="")
      }
    }else{
      pre.muttype.label <- c("D", "I", "SNV", "M")[which.muttype]

      # if it is special label(eg."no variant"), it will be ignored. 
      if(length(pre.muttype.label) < 1){
        break
      }

    }
  }

  options(warn=-1)
  del.group.ind <- grep("D", mut.label)
  del.list <- list()
  if(length(del.group.ind) > 0){
    for(mut.type in mut.label[del.group.ind]){
      del.type <- gsub("D$", "", mut.type)
      del.vec <- as.numeric(strsplit(del.type, ":")[[1]])
      names(del.vec) <- c("pos", "size")
      temp.list.name <- names(del.list)
      del.list <- c(del.list, list(del.vec))
      names(del.list) <- c(temp.list.name, mut.type)
    }
  }
  ins.group.ind <- grep("I", mut.label)
  ins.list <- list()
  if(length(ins.group.ind) > 0){
    for(mut.type in mut.label[ins.group.ind]){
      ins.type <- gsub("I$", "", mut.type)
      ins.vec <- as.numeric(strsplit(ins.type, ":")[[1]])
      names(ins.vec) <- c("pos", "size")
      temp.list.name <- names(ins.list)
      ins.list <- c(ins.list, list(ins.vec))
      names(ins.list) <- c(temp.list.name, mut.type)
    }
  }
  snv.group.ind <- grep("SNV", mut.label)
  snv.list <- list()
  if(length(snv.group.ind) > 0){
    for(mut.type in mut.label[snv.group.ind]){
      snv.type <- gsub("^SNV:", "", mut.type)
      snv.vec <- c(as.numeric(gsub("[[:alpha:]]", "", snv.type)), gsub("[[:digit:]]", "", snv.type))
      names(snv.vec) <- c("pos", "base")
      temp.list.name <- names(snv.list)
      snv.list <- c(snv.list, list(snv.vec))
      names(snv.list) <- c(temp.list.name, mut.type)
    }
  }
  large.group.ind <- grep("M", mut.label)
  large.list <- list()
  if(length(large.group.ind) > 0){
    for(mut.type in mut.label[large.group.ind]){
      large.type <- gsub("M$", "", mut.type)
      large.vec <- c(as.numeric(large.type))
      names(large.vec) <- c("pos")
      temp.list.name <- names(large.list)
      large.list <- c(large.list, list(large.vec))
      names(large.list) <- c(temp.list.name, mut.type)
    }
  }
  mut.list <- list(deletion = del.list, insertion = ins.list, substitution = snv.list, large.indel = large.list)
  options(warn=0)
  return(mut.list)
}

ExtractConsensusSeqs <- function(ref, consensus.allele, start.pos, pos.vec, indel.size.vec, focused.indel.ind, type = c("left.ki", "right.ki", "mut")){
  type <- match.arg(type)
  if(class(consensus.allele) == "DNAStringSet"){
    warning("only the first element will be used.")
    consensus.allele <- consensus.allele[[1]]
  }
  # extract indel sequence
  focused.pos <- pos.vec[focused.indel.ind]
  focused.indel.size <- indel.size.vec[focused.indel.ind]
  through.indel.ind.vec <- 1:focused.indel.ind
  through.indel.size.sum <- sum(indel.size.vec[through.indel.ind.vec])
  ref.length <- length(ref)
  if(focused.indel.size > 0){ # insertion
    five.seq.end <- focused.pos + through.indel.size.sum - focused.indel.size - start.pos
    insert.seq.start <- focused.pos + through.indel.size.sum - focused.indel.size - start.pos + 1
    insert.seq.end <- focused.pos + through.indel.size.sum - start.pos
    three.seq.start <- focused.pos + through.indel.size.sum - start.pos + 1

    adjust.base.size <- 0
    deletion.size <- 0
  }else{ # deletion
    five.seq.end <- focused.pos + through.indel.size.sum - indel.size.vec[focused.indel.ind] - start.pos
    insert.seq.start <- focused.pos + through.indel.size.sum - indel.size.vec[focused.indel.ind] - start.pos + 1
    insert.seq.end <- focused.pos + through.indel.size.sum - indel.size.vec[focused.indel.ind] - start.pos
    three.seq.start <- focused.pos + through.indel.size.sum - indel.size.vec[focused.indel.ind] - start.pos + 1

    adjust.base.size <- through.indel.size.sum + abs(focused.indel.size)
    deletion.size <- abs(focused.indel.size)
  }
  if(type == "right.ki" & focused.pos > 0){
    five.seq.end <- five.seq.end - 1
    insert.seq.start <- insert.seq.start - 1
    insert.seq.end <- insert.seq.end - 1
    three.seq.start <- three.seq.start - 1
  }else if(type == "mut" & focused.pos > 0){
    five.seq.end <- five.seq.end - 1
    insert.seq.start <- insert.seq.start - 1
    insert.seq.end <- insert.seq.end - 1
    three.seq.start <- three.seq.start - 1
  }
  # browser()
  # detect microhomology
  three.conserved.mh <- ""
  five.conserved.mh <- ""
  if(focused.indel.size <= 0){ # deletion
    split.consensus.seq.list <- list(
      ref.seq = ref,
      consensus.seq = consensus.allele,
      five.seq = subseq(consensus.allele
        , 1
        , five.seq.end
      ),
      insert.seq = DNAString(""),
      deletion.seq = subseq(ref
        , insert.seq.start - adjust.base.size
        , insert.seq.start - adjust.base.size + deletion.size - 1
      ),
      three.seq = subseq(consensus.allele
        , three.seq.start
        , nchar(consensus.allele)
      )
    )

    three.len.lim <- min(nchar(split.consensus.seq.list$deletion.seq), length(split.consensus.seq.list$three.seq))
    three.conserved.seq <- compareStrings(
      subseq(split.consensus.seq.list$deletion.seq, 1, three.len.lim)
      , subseq(split.consensus.seq.list$three.seq, 1, three.len.lim)
    )
    three.conserved.mh <- str_extract(three.conserved.seq, "^[A-Za-z]*")
    five.len.lim <- min(nchar(split.consensus.seq.list$deletion.seq), length(split.consensus.seq.list$five.seq))
    five.conserved.seq <- compareStrings(
      subseq(split.consensus.seq.list$deletion.seq, -five.len.lim, -1)
      , subseq(split.consensus.seq.list$five.seq, -five.len.lim, -1)
    )
    five.conserved.mh <- str_extract(five.conserved.seq, "[A-Za-z]*$")

    if(nchar(three.conserved.mh) > nchar(five.conserved.mh)){
      mh.type.label <- "three.prime.conserved.mmej"
      mh.seq.vec <- c(three.conserved.mh)
      trimmed.seq.vec <- str_replace(as.character(split.consensus.seq.list$deletion.seq), paste0("^", three.conserved.mh), "")
    }else if(nchar(three.conserved.mh) < nchar(five.conserved.mh)){
      mh.type.label <- "five.prime.conserved.mmej"
      mh.seq.vec <- c(five.conserved.mh)
      trimmed.seq.vec <- str_replace(as.character(split.consensus.seq.list$deletion.seq), paste0(five.conserved.mh, "$"), "")
    }else if(nchar(three.conserved.mh) > 0 & nchar(three.conserved.mh) == nchar(five.conserved.mh)){
      mh.type.label <- "both.conserved.conserved"
      mh.seq.vec <- c(three.conserved.mh, five.conserved.mh)
      trimmed.seq.vec <- c(
        str_replace(as.character(split.consensus.seq.list$deletion.seq), paste0("^", three.conserved.mh), "")
        , trimmed.seq.vec <- str_replace(as.character(split.consensus.seq.list$deletion.seq), paste0(five.conserved.mh, "$"), "")
      )
    }else{
      mh.type.label <- "nhej"
      mh.seq.vec <- c("")
      trimmed.seq.vec <- c(as.character(split.consensus.seq.list$deletion.seq))
    }
  }else{
    split.consensus.seq.list <- list(
      ref.seq = ref,
      consensus.seq = consensus.allele,
      five.seq = subseq(consensus.allele
        , 1
        , five.seq.end
      ),
      insert.seq = subseq(consensus.allele
        , insert.seq.start
        , insert.seq.end
      ),
      deletion.seq = DNAString(""),
      three.seq = subseq(consensus.allele
        , three.seq.start
        , nchar(consensus.allele)
      )
    )

    mh.type.label <- "no_deletion"
    mh.seq.vec <- c("")
    trimmed.seq.vec <- c("")
  }
  split.consensus.seq.list <- c(
    split.consensus.seq.list,
    list(
      ej.type.label = mh.type.label,
      mh.seq.vec = mh.seq.vec,
      trimmed.seq.vec = trimmed.seq.vec,
      mononuc.freq.table = oligonucleotideFrequency(DNAStringSet(mh.seq.vec[1]), width = 1),
      dinuc.freq.table = oligonucleotideFrequency(DNAStringSet(mh.seq.vec[1]), width = 2),
      mononuc.trimmed.freq.table = oligonucleotideFrequency(DNAStringSet(trimmed.seq.vec[1]), width = 1),
      dinuc.trimmed.freq.table = oligonucleotideFrequency(DNAStringSet(trimmed.seq.vec[1]), width = 2)
    )
  )
  return(split.consensus.seq.list)
}

#'@title Make total mutant list
#'@description Make mutant list from crispr.set (crispRvariant package)
#'@param crispr.set output of TargetToSeq
#'@author Kazuki Nakamae
#'@return The total mutant list including mutant name ("-1:2D"), mutant count (20) and mutant list (deletion : c(-1, 2)) in variants.count.data, lastly include indel sequence data.
# You can access this infomation using variants.count.total.list[index]
# index:length % 4 == 1 : mut.name
# index:length % 4 == 2 : mut.cnt
# index:length % 4 == 3 : mut.list
# index:length % 4 == 0 : indel.seq
#'@rdname MakeTotalMutList
MakeTotalMutList <- function(crispr.set, mut.type){
  variants.count.data <- SafeVariantCounts(crispr.set)
  if(is.null(variants.count.data)){
    return(list())
  }

  if(nrow(variants.count.data) <= 0){
    return(list())
  }else{
    variants.count.total.list <- list()
    consensus.allele.set <- crispr.set$consensusAlleles()
    for(var.ind in 1:nrow(variants.count.data)){
      #if(var.ind > 1 & nrow(variants.count.data) == 1){
      #  break
      #}
      raw.mut.label <- rownames(variants.count.data)[var.ind]
      mut.label <- strsplit(raw.mut.label, ",")[[1]]
      # consensus seq data (this is insert or deletion only)
      is.indel.label <- !(is.na(str_match(mut.label, "D$")[1,1]) & is.na(str_match(mut.label, "I$")[1,1]))
      if(is.indel.label){
        # make pos.vec
        pos.vec <- as.vector(as.numeric(str_match(mut.label, "-*[0-9]+")))
        # make focused.indel.ind. We condider nearest mutation as a main editing outcome.
        focused.indel.ind <- which(abs(pos.vec) == min(abs(pos.vec)))[1]
        # make indel.size.vec
        del.length.label.vec <- str_match(mut.label, "[0-9]+D$")[,1]
        del.length.label.vec[is.na(del.length.label.vec)] <- 0
        deletion.size.vec <- -as.vector(as.numeric(str_match(del.length.label.vec, "[0-9]+")))
        ins.length.label.vec <- str_match(mut.label, "[0-9]+I$")[,1]
        ins.length.label.vec[is.na(ins.length.label.vec)] <- 0
        insertion.size.vec <- as.vector(as.numeric(str_match(ins.length.label.vec, "[0-9]+")))
        indel.size.vec <- deletion.size.vec + insertion.size.vec

        # consensus sequence is analysied in min positon < x < max positon & indel size is less than 30bp.
        # otherwise, out of range may be occur.
        min.lim.pos <- crispr.set$genome_to_target[as.character(crispr.set$target@ranges@start)] # alignment start point
        # max.lim.pos <- crispr.set$genome_to_target[as.character(crispr.set$target@ranges@start + crispr.set$target@ranges@width - 1)]
        max.lim.pos <- crispr.set$genome_to_target[as.character(crispr.set$target@ranges@start + crispr.set$target@ranges@width - 2)] # alignment end point - 2

        # extract cnsensus sequence and µH
        # the range of microhomology is less and equal to 30bp due to the limitation of data.
        # if((pos.vec[focused.indel.ind] > min.lim.pos) & (max.lim.pos > pos.vec[focused.indel.ind] + abs(indel.size.vec[focused.indel.ind])) & (abs(indel.size.vec[focused.indel.ind]) < 30)){
        # if((pos.vec[focused.indel.ind] > min.lim.pos) & (max.lim.pos > pos.vec[focused.indel.ind] + abs(indel.size.vec[focused.indel.ind]))){
        if(all(pos.vec > min.lim.pos) & all(max.lim.pos > pos.vec + abs(sum(indel.size.vec)))){
          tryCatch({
            concensus.seq.list <- ExtractConsensusSeqs(
              ref = crispr.set$ref
              , consensus.allele = consensus.allele.set[[raw.mut.label]]
              , start.pos = crispr.set$genome_to_target[as.character(crispr.set$target@ranges@start)]
              , pos.vec = pos.vec
              , indel.size.vec = indel.size.vec
              , focused.indel.ind = focused.indel.ind
              , type = mut.type
            )
          }, error = function(e) {
            # browser()
            stop("This is unexpected error. Please get touch in the developer <kazukinakamae@gmail.com>.")
          })
        }else{
          concensus.seq.list <- list(
            ref.seq = DNAString(),
            consensus.seq = DNAString(),
            five.seq = DNAString(),
            insert.seq = DNAString(),
            deletion.seq = DNAString(),
            three.seq = DNAString(),
            ej.type.label = "unidentified",
            mh.seq.vec = c(),
            trimmed.seq.vec = c(),
            mononuc.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 1),
            dinuc.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 2),
            mononuc.trimmed.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 1),
            dinuc.trimmed.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 2)
          )
        }
      }else{
        concensus.seq.list <- list(
          ref.seq = DNAString(),
          consensus.seq = DNAString(),
          five.seq = DNAString(),
          insert.seq = DNAString(),
          deletion.seq = DNAString(),
          three.seq = DNAString(),
          ej.type.label = "no_deletion",
          mh.seq.vec = c(),
          trimmed.seq.vec = c(),
          mononuc.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 1),
          dinuc.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 2),
          mononuc.trimmed.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 1),
          dinuc.trimmed.freq.table = oligonucleotideFrequency(DNAStringSet(""), width = 2)
        )
      }
      # add list
      variants.count.total.list <- c(
        variants.count.total.list,
        list(
          mut.name = raw.mut.label
          , mut.cnt = variants.count.data[var.ind]
          , mut.list = MakeMutList(mut.label)
          , indel.seq = concensus.seq.list
        )
      )
    }
  }
  return(variants.count.total.list)
}

#'@title Make mutant table for plot
#'@description Make mutant table including position and count for plot
#'@param variants.count.total.list The total mutant list
#'@param mut.type Type of mutant
#'@param range.vec Range of position
#'@author Kazuki Nakamae
#'@rdname MakeMutPositionLabelTable
MakeMutPositionLabelTable <- function(variants.count.total.list, mut.type = c("deletion", "insertion", "substitution", "large", "no.overlap.deletion"), range.vec = c(-35:-1,1:35)){
  
  mut.type <- match.arg(mut.type)
  pos.vec <- numeric(0)

  if(length(variants.count.total.list) == 0){
    return(table(factor(, levels = range.vec)))
  }

  for(mut.ind in 1:(length(variants.count.total.list)/4)){

    if(mut.type == "no.overlap.deletion"){

      for(mut in getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], "deletion")){
        if((getElement(mut, "pos") < 0) & (getElement(mut, "pos") + getElement(mut, "size") < 1)){ # The deletion covers with cutting site (position:0) ?
          pos.vec <- c(pos.vec, rep(getElement(mut, "pos"), variants.count.total.list[mut.ind * 4 - 2][[1]]))
        }else if(getElement(mut, "pos") > 0){
          pos.vec <- c(pos.vec, rep(getElement(mut, "pos"), variants.count.total.list[mut.ind * 4 - 2][[1]]))
        }
      }

    }else{

      for(mut in getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], mut.type)){
        pos.vec <- c(pos.vec, rep(getElement(mut, "pos"), variants.count.total.list[mut.ind * 4 - 2][[1]]))
      }

    }
  }

  return(table(factor(pos.vec, levels = range.vec)))
}

MakeSubstitutionNumberReadTable <- function(variants.count.total.list, focus.range = -100:100, range.vec = 0:100, special.label = "", is.special.included = TRUE){
  num.vec <- numeric(0)
  if(length(variants.count.total.list) == 0){
    return(table(factor(, levels = range.vec)))
  }
  for(mut.ind in 1:(length(variants.count.total.list)/4)){
    # get number of special label
    if(variants.count.total.list[mut.ind * 4 - 3][[1]] %in% special.label){
      if(is.special.included){
        num.vec <- c(num.vec, (rep(0, variants.count.total.list[mut.ind * 4 - 2][[1]])))
      }
      next
    }

    # count number of sbstitution among range
    num.cnt <- 0
    for(mut in getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], "substitution")){
      if(getElement(mut, "pos") %in% focus.range){
        num.cnt <- num.cnt + 1
      }
    }
    if(num.cnt > 0){
      num.vec <- c(num.vec, rep(num.cnt, variants.count.total.list[mut.ind * 4 - 2][[1]]))
    }else if(
      !any(
        length(getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], "deletion")) > 0
        , length(getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], "insertion")) > 0
        , length(getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], "large.indel")) > 0
      )
    ){
      num.vec <- c(num.vec, rep(0, variants.count.total.list[mut.ind * 4 - 2][[1]]))
    }
  }
  return(table(factor(num.vec, levels = range.vec)))
}

MakeMutSizeReadTable <- function(variants.count.total.list, focus.range = -100:100, range.vec = -100:100, special.label = "", is.special.included = TRUE){
  size.vec <- numeric(0)
  if(length(variants.count.total.list) == 0){
    return(table(factor(, levels = range.vec)))
  }

  for(mut.ind in 1:(length(variants.count.total.list)/4)){

    # get number of special label
    if(variants.count.total.list[mut.ind * 4 - 3][[1]] %in% special.label){
      if(is.special.included){
        size.vec <- c(size.vec, (rep(0, variants.count.total.list[mut.ind * 4 - 2][[1]])))
      }
      next
    }

    # get number of deletion
    del.size <- 0
    for(del in getElement(variants.count.total.list[mut.ind * 4 -1][[1]], "deletion")){
      if(getElement(del, "pos") %in% focus.range){
        del.size <- del.size + getElement(del, "size")
      }
    }

    # get number of insertion
    ins.size <- 0
    for(ins in getElement(variants.count.total.list[mut.ind * 4 - 1][[1]], "insertion")){
      if(getElement(ins, "pos") %in% focus.range){
        ins.size <- ins.size + getElement(ins, "size")
      }
    }

    indel.size <- ins.size - del.size
    size.vec <- c(size.vec, rep(indel.size, variants.count.total.list[mut.ind * 4 - 2][[1]]))
  }
  return(table(factor(size.vec, levels = range.vec)))
}

MakeKISizeReadTable <- function(variants.count.total.list, mut.type = c("deletion", "insertion"), range.vec = 1:100, special.label = "",
  is.adjust.mode = FALSE, adjustment.applied.pos.range = numeric(0), adjustment.add.size = 0){
  mut.type <- match.arg(mut.type)
  size.vec <- numeric(0)
  if(length(variants.count.total.list) == 0){
    return(table(factor(, levels = range.vec)))
  }
  for(mut.ind in 1:(length(variants.count.total.list)/3)){

    if(!is.adjust.mode){


      # get number of deletion
      del.size <- 0
      for(del in getElement(variants.count.total.list[mut.ind * 3][[1]], "deletion")){
        del.size <- del.size + getElement(del, "size")
      }

      # get number of insertion
      ins.size <- 0
      for(del in getElement(variants.count.total.list[mut.ind * 3][[1]], "insertion")){
        ins.size <- ins.size + getElement(del, "size")
      }

      indel.size <- ins.size - del.size
      size.vec <- c(size.vec, rep(indel.size, variants.count.total.list[mut.ind * 3 - 1][[1]]))
      for(mut in getElement(variants.count.total.list[mut.ind * 3][[1]], "deletion")){
        size.vec <- c(size.vec, rep(getElement(mut, "size"), variants.count.total.list[mut.ind * 3 - 1][[1]]))
      }

    }else{

      if(mut.type == "insertion"){

        # get number of special label
        if(variants.count.total.list[mut.ind * 3 - 2][[1]] %in% special.label){
          # special label == 0 insertion (relatively)
          size.vec <- c(size.vec, (rep(0, variants.count.total.list[mut.ind * 3 - 1][[1]]) + adjustment.add.size))
          next
        }

      }


      # get number
      for(mut in getElement(variants.count.total.list[mut.ind * 3][[1]], mut.type)){

        if(getElement(mut, "pos") %in% adjustment.applied.pos.range){
          # within adjustment.applied.pos.range
          if(mut.type == "insertion"){
            size.vec <- c(size.vec, (rep(getElement(mut, "size"), variants.count.total.list[mut.ind * 3 - 1][[1]]) + adjustment.add.size))
          }else if(mut.type == "deletion"){
            size <- getElement(mut, "size") - adjustment.add.size
            if(!(size %in% range.vec)){
              range.vec <- c(size, range.vec)[order(c(size, range.vec))]
            }
            size.vec <- c(size.vec, (rep(size, variants.count.total.list[mut.ind * 3 - 1][[1]]) ))
            # negative value == insertion
          }
        }else{
          # not within adjustment.applied.pos.range
          size.vec <- c(size.vec, rep(getElement(mut, "size"), variants.count.total.list[mut.ind * 3 - 1][[1]]))
        }

      }

    }

  }
  return(table(factor(size.vec, levels = range.vec)))
}

SaveBarPlot <- function(position.count.table, color.code, file.name, xlab, ylab){
  
  if(sum(position.count.table$Frequency) < 1){
    return()
  }

  options(warn = -1)
  position.count.p <- ggplot(data=position.count.table
    , aes(x=FactorToNumeric(Position), y=Frequency, fill = c("substitution"))
  ) +
  geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = color.code) +
  coord_cartesian(ylim = c(0, max(position.count.table$Frequency) + 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)
      , axis.text.y = element_text(size = 7)
      , axis.title.x = element_text(size = 7)
      , axis.title.y = element_text(size = 7)
    ) +
  xlab(xlab) +
  ylab(ylab)
  # plot(position.count.p)
  ggsave(file = file.name, plot = position.count.p, dpi = 300, width = 0.1 * length(position.count.table$Frequency), height = 3, limitsize = FALSE)
  options(warn = 0)
  message(paste("Save bar plot in ", file.name, sep=""))
}

SavePositionBarPlot <- function(variants.count.total.list, out.dir, mut.type = c("deletion", "insertion", "substitution", "large", "no.overlap.deletion"), range.vec = c(-35:-1,1:35)){
  
  mut.type <- match.arg(mut.type)
  if(mut.type == "deletion"){
    color.code = "#8C4303"
    file.name = file.path(out.dir, "Distribution_of_deletion_labels_on_position.png")
    ylab = "Number of deletion label"
  }else if(mut.type == "no.overlap.deletion"){
    color.code = "#592202"
    file.name = file.path(out.dir, "Distribution_of_no_overlap_deletion_labels_on_position.png")
    ylab = "Number of deletion(no overlap with cutsite) label"
  }else if(mut.type == "insertion"){
    color.code = "#F2B705"
    file.name = file.path(out.dir, "Distribution_of_insertion_labels_on_position.png")
    ylab = "Number of insertion label"
  }else if(mut.type == "substitution"){
    color.code = "#F27D72"
    file.name = file.path(out.dir, "Distribution_of_substitution_labels_on_position.png")
    ylab = "Number of substitution(SNV) label"
  }else if(mut.type == "large"){
    color.code = "#BFB639"
    file.name = file.path(out.dir, "Distribution_of_large_indels_labels_on_position.png")
    ylab = "Number of large indels label"
  }else{
    stop("Unexpected Error")
  }

  mut.table <- data.frame(MakeMutPositionLabelTable(variants.count.total.list
    , mut.type = mut.type
    , range.vec = range.vec)
    , stringsAsFactors=FALSE)
  colnames(mut.table) <- c("Position", "Frequency")
  position.mut.table <- data.frame(Position=mut.table$Position, Frequency=mut.table$Frequency, stringsAsFactors=FALSE)
  if(!is.null(position.mut.table)){
    SaveBarPlot(position.mut.table
      , color.code = color.code
      , file.name = file.name
      , xlab = "Relative position[bp] (Zero is cut site.)"
      , ylab = ylab
    )
  }

}

SaveVariantsData <- function(crispr.set, out.dir, range.vec = -35:35, mut.type){
  is.perfect.process = FALSE 
  tryCatch({
    if(is.null(crispr.set)){
      is.perfect.process = TRUE
      return(is.perfect.process)
    }
    # Save R object
    saveRDS(crispr.set, file = file.path(out.dir, "crispr.set.rds"))
    # variants
    write.csv(as.data.frame(variantCounts(crispr.set)), file.path(out.dir, "variant.counts.csv"))
    # mutation efficiency
    write.csv(as.data.frame(mutationEfficiency(crispr.set)), file.path(out.dir, "mutantation.efficiency.csv"))
    # show alignment
    pdf(file.path(out.dir, "aligned.variants.pdf"), width = nchar(crispr.set$ref)/5, height = length(rownames(variantCounts(crispr.set))) / 3 + 6)
    crispr.set.p <- plotVariants(crispr.set, row.ht.ratio = c(1, 6), col.wdth.ratio = c(8.5, 1)
      , plotAlignments.args = list(top.n = 200, min.insertion.freq = 0, max.insertion.size = 00, highlight.pam = FALSE, line.weight = 1
        ,legend.cols = 3, axis.text.size = 15, tile.height = 0.55, highlight.guide = FALSE, ins.size = 5, plot.text.size = 4)
      , plotFreqHeatmap.args = list(top.n = 200, legend.key.height = ggplot2::unit(1.5, "lines")
        , plot.text.size = 5)
    )
    dev.off()
    message(paste("Save alignment plot in ", file.path(out.dir, "aligned.variants.pdf"), sep=""))
    # save concensus sequences
    write.csv(as.data.frame(consensusSeqs(crispr.set)), file.path(out.dir, "consensus.seqs.csv"))
    # make bar plot of variants
    variants.count.data <- variantCounts(crispr.set)
    colPal1 <- colorRamp(colors(1000), space = "rgb")
    pdf(file.path(out.dir, "variants.barplot.pdf"), width = 10, height = length(rownames(variantCounts(crispr.set))) / 12 + 4)
    variants.count.data.p <- barplotAlleleFreqs(variants.count.data, classify = FALSE, include.table = FALSE
      , bar.colours = rgb(colPal1(sapply(rownames(variants.count.data)
        , function(x){
            color.vec <- sum(strtoi(charToRaw(x), 16L)/2000)
            color.vec <- replace(color.vec, color.vec > 1, 1)
            return(color.vec)
          }
        ))/255
      )
    )
    plot(variants.count.data.p)
    dev.off()
    message(paste("Save variants barplot in ", file.path(out.dir, "variants.barplot.pdf"), sep=""))

    # browser()
    variants.count.total.list <- MakeTotalMutList(crispr.set, mut.type)
    saveRDS(variants.count.total.list, file = file.path(out.dir, "variants.count.total.list.rds"))
    # each mutation type of position plot
    SavePositionBarPlot(variants.count.total.list, out.dir, mut.type = "deletion", range.vec = range.vec)
    SavePositionBarPlot(variants.count.total.list, out.dir, mut.type = "insertion", range.vec = range.vec)
    SavePositionBarPlot(variants.count.total.list, out.dir, mut.type = "substitution", range.vec = range.vec)
    SavePositionBarPlot(variants.count.total.list, out.dir, mut.type = "no.overlap.deletion", range.vec = range.vec)
    SavePositionBarPlot(variants.count.total.list, out.dir, mut.type = "large", range.vec = range.vec)

    is.perfect.process = TRUE
  },error = function(e) {
      message("ERROR!")
      message(e)           
      traceback()
    },
    silent = TRUE
  )
  return(is.perfect.process)
}

SortDataFrame <- function(merge.table, by.vector){
  sortlist <- order(by.vector)
  merge.table <- merge.table[sortlist,]
  return(merge.table)
}

MakeMargeSumPositionTable <- function(merge.table){
  return(data.frame(
    Position=merge.table[, names(merge.table) %in% c("Position")]
    , Frequency=apply(merge.table[, !names(merge.table) %in% c("Position")], MARGIN=1, sum)
      , stringsAsFactors=FALSE
    )
  )
}

MakeMargeSumSizeTable <- function(merge.table){
  return(data.frame(
    Size=merge.table[, names(merge.table) %in% c("Size")]
    , Frequency=apply(merge.table[, !names(merge.table) %in% c("Size")], MARGIN=1, sum)
    , stringsAsFactors=FALSE
    )
  )
}

MakeMargeSumNumberTable <- function(merge.table){
  return(data.frame(
    Number=merge.table[, names(merge.table) %in% c("Number")]
    , Frequency=apply(merge.table[, !names(merge.table) %in% c("Number")], MARGIN=1, sum)
    , stringsAsFactors=FALSE
    )
  )
}

# split variants.count.total.list by enc joining type.
SplitTotalMutListByEjType <- function(variants.count.total.list){
  if(length(variants.count.total.list) == 0){
    return(
      list(mmej.variants.count.total.list = list()
        , nhej.variants.count.total.list = list()
        , unidentified.deletion.variants.count.total.list = list()
        , no.variants.count.total.list = list()
      )
    )
  }
  mmej.variants.count.total.list = list()
  nhej.variants.count.total.list = list()
  unidentified.deletion.variants.count.total.list = list()
  no.variants.count.total.list = list()
  for(mut.ind in 1:(length(variants.count.total.list) / 4)){
    ej.type.label <- variants.count.total.list[mut.ind * 4]$indel.seq$ej.type.label
    if(ej.type.label %in% c("three.prime.conserved.mmej", "five.prime.conserved.mmej", "both.conserved.conserved")){
      mmej.variants.count.total.list <- c(
        mmej.variants.count.total.list
        , variants.count.total.list[mut.ind * 4 - 3]
        , variants.count.total.list[mut.ind * 4 - 2]
        , variants.count.total.list[mut.ind * 4 - 1]
        , variants.count.total.list[mut.ind * 4]
      )
    }else if(ej.type.label %in% "nhej"){
      nhej.variants.count.total.list <- c(
        nhej.variants.count.total.list
        , variants.count.total.list[mut.ind * 4 - 3]
        , variants.count.total.list[mut.ind * 4 - 2]
        , variants.count.total.list[mut.ind * 4 - 1]
        , variants.count.total.list[mut.ind * 4]
      )
    }else if(ej.type.label %in% "unidentified"){
      unidentified.deletion.variants.count.total.list <- c(
        unidentified.deletion.variants.count.total.list
        , variants.count.total.list[mut.ind * 4 - 3]
        , variants.count.total.list[mut.ind * 4 - 2]
        , variants.count.total.list[mut.ind * 4 - 1]
        , variants.count.total.list[mut.ind * 4]
      )
    }else{
      no.variants.count.total.list <- c(
        no.variants.count.total.list
        , variants.count.total.list[mut.ind * 4 - 3]
        , variants.count.total.list[mut.ind * 4 - 2]
        , variants.count.total.list[mut.ind * 4 - 1]
        , variants.count.total.list[mut.ind * 4]
      )
    }
  }
  return(
    list(
      mmej.variants.count.total.list = mmej.variants.count.total.list
      , nhej.variants.count.total.list = nhej.variants.count.total.list
      , unidentified.deletion.variants.count.total.list = unidentified.deletion.variants.count.total.list
      , no.variants.count.total.list = no.variants.count.total.list
    )
  )
}

MakeIsSnReadTable <- function(variants.count.total.list, focus.range = -35:35, substitution.range.vec = 0:100, indels.range.vec = -100:100, special.label = "", is.special.included = TRUE){ # Is: Indels size, Sn: Substitution number
  # split variants.count.total.list
  split.variants.count.total.list <- SplitTotalMutListByEjType(variants.count.total.list)
  mmej.variants.count.total.list <- split.variants.count.total.list$mmej.variants.count.total.list
  nhej.variants.count.total.list <- split.variants.count.total.list$nhej.variants.count.total.list
  unidentified.deletion.variants.count.total.list <- split.variants.count.total.list$unidentified.deletion.variants.count.total.list
  insert.dominant.variants.count.total.list <- split.variants.count.total.list$no.variants.count.total.list

  # indels/size
  indels.size.table <- data.frame(MakeMutSizeReadTable(variants.count.total.list, focus.range = focus.range, range.vec = indels.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(indels.size.table) <- c("Size", "Frequency")
  mmej.indels.size.table <- data.frame(MakeMutSizeReadTable(mmej.variants.count.total.list, focus.range = focus.range, range.vec = indels.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(mmej.indels.size.table) <- c("Size", "Frequency")
  nhej.indels.size.table <- data.frame(MakeMutSizeReadTable(nhej.variants.count.total.list, focus.range = focus.range, range.vec = indels.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(nhej.indels.size.table) <- c("Size", "Frequency")
  unidentified.deletion.indels.size.table <- data.frame(MakeMutSizeReadTable(unidentified.deletion.variants.count.total.list, focus.range = focus.range, range.vec = indels.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(unidentified.deletion.indels.size.table) <- c("Size", "Frequency")
  insert.dominant.indels.size.table <- data.frame(MakeMutSizeReadTable(insert.dominant.variants.count.total.list, focus.range = focus.range, range.vec = indels.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(insert.dominant.indels.size.table) <- c("Size", "Frequency")
  
  # substitution/number
  substitution.number.table <- data.frame(MakeSubstitutionNumberReadTable(variants.count.total.list, focus.range = focus.range, range.vec = substitution.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(substitution.number.table) <- c("Number", "Frequency")
  mmej.substitution.number.table <- data.frame(MakeSubstitutionNumberReadTable(mmej.variants.count.total.list, focus.range = focus.range, range.vec = substitution.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(mmej.substitution.number.table) <- c("Number", "Frequency")
  nhej.substitution.number.table <- data.frame(MakeSubstitutionNumberReadTable(nhej.variants.count.total.list, focus.range = focus.range, range.vec = substitution.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(nhej.substitution.number.table) <- c("Number", "Frequency")
  unidentified.deletion.substitution.number.table <- data.frame(MakeSubstitutionNumberReadTable(unidentified.deletion.variants.count.total.list, focus.range = focus.range, range.vec = substitution.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(unidentified.deletion.substitution.number.table) <- c("Number", "Frequency")
  insert.dominant.substitution.number.table <- data.frame(MakeSubstitutionNumberReadTable(insert.dominant.variants.count.total.list, focus.range = focus.range, range.vec = substitution.range.vec, special.label = special.label, is.special.included = is.special.included), stringsAsFactors=FALSE)
  colnames(insert.dominant.substitution.number.table) <- c("Number", "Frequency")

  is.sn.list <- list(
    total.is.sn.list = list(indels.size.table = indels.size.table, substitution.number.table = substitution.number.table)
    , mmej.is.sn.list = list(indels.size.table = mmej.indels.size.table, substitution.number.table = mmej.substitution.number.table)
    , nhej.is.sn.list = list(indels.size.table = nhej.indels.size.table, substitution.number.table = nhej.substitution.number.table)
    , unidentified.deletion.is.sn.list = list(indels.size.table = unidentified.deletion.indels.size.table, substitution.number.table = unidentified.deletion.substitution.number.table)
    , insert.dominant.is.sn.list = list(indels.size.table = insert.dominant.indels.size.table, substitution.number.table = insert.dominant.substitution.number.table)
  )

  return(is.sn.list)
}

SaveIndelSizeReadBarPlot <- function(indels.size.table, file.name){
  options(warn = -1)
  indels.size.p <- ggplot(data=indels.size.table
    , aes(x=Size, y=Frequency, fill = c("indels"))) +
    geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#D93240") +
    coord_cartesian(ylim = c(0, max(indels.size.table$Frequency) + 100)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    , axis.text.y = element_text(size = 6)
    , axis.title.x = element_text(size = 12)
    , axis.title.y = element_text(size = 12)) +
  xlab("Indels size[bp]") +
  ylab("Number of indels reads")
  # plot(indels.size.p)
  ggsave(file = file.name, plot = indels.size.p, dpi = 300, width = 16, height = 4)
  options(warn = 0)
}

SaveSubstitutionNumberReadBarPlot <- function(substitution.number.table, file.name){
  options(warn = -1)
  substitution.number.p <- ggplot(data=substitution.number.table
    , aes(x=Number, y=Frequency, fill = c("substitution"))) +
    geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#D93240") +
    coord_cartesian(ylim = c(0, max(substitution.number.table$Frequency) + 100)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    , axis.text.y = element_text(size = 6)
    , axis.title.x = element_text(size = 12)
    , axis.title.y = element_text(size = 12)) +
  xlab("Substitution number") +
  ylab("Number of substitution reads")
  # plot(substitution.number.p)
  ggsave(file = file.name, plot = substitution.number.p, dpi = 300, width = 12, height = 6)
  options(warn = 0)
}

MakeAggregatedIsSnTableList <- function(is.sn.list){
  # mutation pie chart
  # classify by indels size
  indels.size.labeled.table <- data.frame(is.sn.list$indels.size.table,
    size.label=apply(is.sn.list$indels.size.table, MARGIN=1, function(indels.size.row){
      size <- as.numeric(indels.size.row[names(indels.size.row) %in% c("Size")])
      if(size == -1){
        return("D1")
      }else if(size == -2){
        return("D2")
      }else if(size < -2){
        return("D>2")
      }else if(size == 1){
        return("I1")
      }else if(size > 1){
        return("I>1")
      }else if(size == 0){
        return("NS")
      }else{
        warning("Unexpected output")
      }
    })
  )
  aggregated.indels.size.table <- aggregate(indels.size.labeled.table$Frequency
    , by=list(size.label=indels.size.labeled.table$size.label)
    , FUN=sum
  )
  # Remove NS(no variants/substitution) label because it will be included in substitution tabels. Avoid double counting.
  colnames(aggregated.indels.size.table) <- c("label", "Frequency")

  # classify by substitutions number
  substitution.number.labeled.table <- data.frame(is.sn.list$substitution.number.table,
    number.label=apply(is.sn.list$substitution.number.table, MARGIN=1, function(substitution.number.row){
      number <- as.numeric(substitution.number.row[names(substitution.number.row) %in% c("Number")])
      if(number == 0){
        return("N")
      }else if(number == 1){
        return("S1")
      }else if(number > 1){
        return("S>1")
      }else{
        warning("Unexpected output")
      }
    })
  )
  aggregated.substitution.number.table <- aggregate(substitution.number.labeled.table$Frequency
    , by=list(number.label=substitution.number.labeled.table$number.label)
    , FUN=sum
  )
  colnames(aggregated.substitution.number.table) <- c("label", "Frequency")
  return(
    list(aggregated.indels.size.table = aggregated.indels.size.table
    ,aggregated.substitution.number.table = aggregated.substitution.number.table)
  )
} 

MakeAggregatedKiIsSnTableList <- function(is.sn.list, direction = c("correct", "reverse"), nhej.insert.length){
  direction <- match.arg(direction)
  if(any(direction %in% c("reverse"))){
    indel.label.func <- function(indels.size.row){
      size <- as.numeric(indels.size.row[names(indels.size.row) %in% c("Size")])
      if(size == nhej.insert.length){
        return("Rd")
      }else if(size > nhej.insert.length){
        return(">Rd")
      }else if(size < nhej.insert.length){
        return("<Rd")
      }else{
        warning("Unexpected output")
      }
    }
    substitution.label.func <- function(substitution.number.row){
      number <- as.numeric(substitution.number.row[names(substitution.number.row) %in% c("Number")])
      if(number == 0){
        return("Rp_S0")
      }else if(number == 1){
        return("Rp_S1")
      }else if(number > 1){
        return("Rp_S>1")
      }else{
        warning("Unexpected output")
      }
    }
  }else{
    indel.label.func <- function(indels.size.row){
      size <- as.numeric(indels.size.row[names(indels.size.row) %in% c("Size")])
      if(size == nhej.insert.length){
        return("Id")
      }else if(size > nhej.insert.length){
        return(">Id")
      }else if((size < nhej.insert.length) & (size > 0)){
        return("Id-Ip")
      }else if(size < 0){
        return("<Ip")
      }else if(size == 0){
        return("Ip")
      }else{
        warning("Unexpected output")
      }
    }
    substitution.label.func <- function(substitution.number.row){
      number <- as.numeric(substitution.number.row[names(substitution.number.row) %in% c("Number")])
      if(number == 0){
        return("P_S0")
      }else if(number == 1){
        return("P_S1")
      }else if(number > 1){
        return("P_S>1")
      }else{
        warning("Unexpected output")
      }
    }
  }
  # mutation pie chart
  # classify by indels size
  indels.size.labeled.table <- data.frame(is.sn.list$indels.size.table,
    size.label = apply(is.sn.list$indels.size.table, MARGIN=1, indel.label.func)
  )
  aggregated.indels.size.table <- aggregate(indels.size.labeled.table$Frequency
    , by=list(size.label=indels.size.labeled.table$size.label)
    , FUN=sum
  )
  # Remove NS(no variants/substitution) label because it will be included in substitution tabels. Avoid double counting.
  colnames(aggregated.indels.size.table) <- c("label", "Frequency")

  # classify by substitutions number
  substitution.number.labeled.table <- data.frame(is.sn.list$substitution.number.table,
    number.label = apply(is.sn.list$substitution.number.table, MARGIN=1, substitution.label.func)
  )
  aggregated.substitution.number.table <- aggregate(substitution.number.labeled.table$Frequency
    , by=list(number.label=substitution.number.labeled.table$number.label)
    , FUN=sum
  )
  colnames(aggregated.substitution.number.table) <- c("label", "Frequency")
  return(
    list(aggregated.indels.size.table = aggregated.indels.size.table
    ,aggregated.substitution.number.table = aggregated.substitution.number.table)
  )
}

SavePieChart <- function(class.count.table, color.vec, elements.name, file.name){

  if(sum(class.count.table$Frequency) < 1){
    return()
  }

    class.count.plot.table <- class.count.table %>%
    arrange(desc(Class)) %>%
    mutate(Lab1.ypos = cumsum(Percentage) - 0.5*Percentage)

    class.count.plot.table <- class.count.plot.table %>%
    arrange(desc(Class)) %>%
    mutate(Lab2.ypos = Lab1.ypos)
    
  # Adjusting label position : this is difficult to make it automatically. We can make plot case by case.
  # browser(expr = basename(file.name) == "aggregated.left.ki.only.i.ri.size.pie.chart.png")
  #for(ind in 2:length(class.count.plot.table$Lab1.ypos)){
  #  if((class.count.plot.table$Lab1.ypos[ind] - class.count.plot.table$Lab1.ypos[ind-1]) < 2){
  #    dec <- (class.count.plot.table$Lab1.ypos[ind-1] - 2)
  #    inc <- (class.count.plot.table$Lab1.ypos[ind-1] + 2)
  #    if(dec <= 0){
  #      #class.count.plot.table$Lab1.ypos[ind-1] <- class.count.plot.table$Lab1.ypos[ind-1]
  #      class.count.plot.table$Lab2.ypos[ind-1] <- class.count.plot.table$Lab2.ypos[ind-1] + 2
  #      class.count.plot.table$Lab2.ypos[ind] <- class.count.plot.table$Lab2.ypos[ind-1] + 3.5
  #    }else if(inc >= 100){
  #      class.count.plot.table$Lab2.ypos[ind] <- 100 - 1.5
  #      class.count.plot.table$Lab2.ypos[ind-1] <- class.count.plot.table$Lab2.ypos[ind] - 3
  #    }else{
  #      class.count.plot.table$Lab2.ypos[ind-1] <- class.count.plot.table$Lab2.ypos[ind-1] + 2
  #      class.count.plot.table$Lab2.ypos[ind] <- class.count.plot.table$Lab2.ypos[ind-1] + 3.5
  #    }
  #  }
  #}
  options(warn = -1)
  class.count.pie.chart <- ggplot(class.count.plot.table, aes(x="", y=Percentage, fill=Class))+
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=color.vec) +
    geom_text(aes(x = 1.30, y=Lab1.ypos, label = paste(round(Percentage, digits = 1), "%", sep = "")), color = "black",  size=5) + 
    geom_text(aes(x = 1.82, y=Lab2.ypos, label = paste(Class, "\n(", Frequency, "", elements.name, ")", sep = "")), color = "black",  size=5) + 
    coord_polar("y", start = pi * 3/2, direction = -1) +
    theme_void()
  # plot(class.count.pie.chart)
  ggsave(file = file.name, plot = class.count.pie.chart, dpi = 300, width = 10, height = 10)
  options(warn = 0)
  message(paste("Save pie chart in ", file.name, sep=""))
}

SaveIndelPlot <- function(is.sn.list, outdir, ki.type = c("mutation")){
  ki.type <- match.arg(ki.type)
  # make indel and substitution plot
  aggregated.is.sn.table.list <- MakeAggregatedIsSnTableList(is.sn.list)
  aggregated.is.sn.table.list$aggregated.indels.size.table <- aggregated.is.sn.table.list$aggregated.indels.size.table[which(aggregated.is.sn.table.list$aggregated.indels.size.table$label != "NS"), ]
  # only substitution number chart
  if(!is.null(aggregated.is.sn.table.list$aggregated.substitution.number.table) & sum(aggregated.is.sn.table.list$aggregated.substitution.number.table$Frequency) > 0){
    
    aggregated.indels.size.class.count.table <- data.frame(
      Class=aggregated.is.sn.table.list$aggregated.substitution.number.table$label
      , Frequency=aggregated.is.sn.table.list$aggregated.substitution.number.table$Frequency
      , Percentage=aggregated.is.sn.table.list$aggregated.substitution.number.table$Frequency / sum(aggregated.is.sn.table.list$aggregated.substitution.number.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲa]Rate_of_", ki.type, "_substitution_number.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲa]Rate_of_", ki.type, "_substitution_number.rds", sep="")))
    SavePieChart(aggregated.indels.size.class.count.table
      , c("#F0F0F2", "#D9BACB", "#D99ABC")
      , "Reads"
      , file.path(outdir, paste("[ⅲa]Rate_of_", ki.type, "_substitution_number_piechart.png", sep="")))
  }
  # only indels size pie chart
  if(!is.null(aggregated.is.sn.table.list$aggregated.indels.size.table) & sum(aggregated.is.sn.table.list$aggregated.indels.size.table$Frequency) > 0){
    
    aggregated.indels.size.class.count.table <- data.frame(
      Class=aggregated.is.sn.table.list$aggregated.indels.size.table$label
      , Frequency=aggregated.is.sn.table.list$aggregated.indels.size.table$Frequency
      , Percentage=aggregated.is.sn.table.list$aggregated.indels.size.table$Frequency / sum(aggregated.is.sn.table.list$aggregated.indels.size.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲb]Rate_of_", ki.type, "_indel_size.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲb]Rate_of_", ki.type, "_indel_size.rds", sep="")))
    SavePieChart(aggregated.indels.size.class.count.table
      , c("#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304")
      , "Reads"
      , file.path(outdir, paste("[ⅲb]Rate_of_", ki.type, "_indel_size_piechart.png", sep="")))
  }

  # substitution number + indels size pie chart
  # make merge aggregated table
  aggregated.mutation.merge.size.number.table <- rbind(aggregated.is.sn.table.list$aggregated.substitution.number.table, aggregated.is.sn.table.list$aggregated.indels.size.table)
  if(!is.null(aggregated.mutation.merge.size.number.table) & sum(aggregated.mutation.merge.size.number.table$Frequency) > 0){
    
    aggregated.indels.size.class.count.table <- data.frame(
      Class=aggregated.mutation.merge.size.number.table$label
      , Frequency=aggregated.mutation.merge.size.number.table$Frequency
      , Percentage=aggregated.mutation.merge.size.number.table$Frequency / sum(aggregated.mutation.merge.size.number.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲc]Rate_of_", ki.type, "_indel_size_and_substitution_number.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲc]Rate_of_", ki.type, "_indel_size_and_substitution_number.rds", sep="")))
    SavePieChart(aggregated.indels.size.class.count.table
      , c("#F0F0F2", "#D9BACB", "#D99ABC", "#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304")
      , "Reads"
      , file.path(outdir, paste("[ⅲc]Rate_of_", ki.type, "_indel_size_and_substitution_number_piechart.png", sep="")))
  }

  # substitution number + indels size pie chart without Unmodified
  # remove N row
  aggregated.mutation.merge.size.number.without.n.table <- aggregated.mutation.merge.size.number.table[which(aggregated.mutation.merge.size.number.table$label != "N"), ]
  if(!is.null(aggregated.mutation.merge.size.number.without.n.table) & sum(aggregated.mutation.merge.size.number.without.n.table$Frequency) > 0){
    
    aggregated.indels.size.class.count.table <- data.frame(
      Class=aggregated.mutation.merge.size.number.without.n.table$label
      , Frequency=aggregated.mutation.merge.size.number.without.n.table$Frequency
      , Percentage=aggregated.mutation.merge.size.number.without.n.table$Frequency / sum(aggregated.mutation.merge.size.number.without.n.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲd]Rate_of_", ki.type, "_indel_size_and_substitution_number(without Unmodified).txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅲd]Rate_of_", ki.type, "_indel_size_and_substitution_number(without Unmodified).rds", sep="")))
    SavePieChart(aggregated.indels.size.class.count.table
      , c("#D9BACB", "#D99ABC", "#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304")
      , "Reads"
      , file.path(outdir, paste("[ⅲd]Rate_of_", ki.type, "_indel_size_and_substitution_number(without Unmodified).png", sep=""))
    )
  }
}

# Bind precise knock-in IsSn list considering precise knock-in as correct size (indels size:0, sunstitution number:0)
BindKiIsSnList <- function(precise.is.sn.list, imprecise.is.sn.list, rev.imprecise.is.sn.list, nhej.insert.length, rc.nhej.insert.length){
  # Standard size is changed from duplicated junction to perfect junction.
  # browser()
  imprecise.is.sn.list$indels.size.table$Size <- factor(-100:100 + nhej.insert.length)
  correct.direction.size.range <- intersect(precise.is.sn.list$indels.size.table$Size, imprecise.is.sn.list$indels.size.table$Size)

  ### correct direction knock-in
  ki.correct.is.sn.list <- precise.is.sn.list
  ki.correct.is.sn.list$indels.size.table <- ki.correct.is.sn.list$indels.size.table[ki.correct.is.sn.list$indels.size.table$Size %in% correct.direction.size.range,]

  ki.correct.is.sn.list$indels.size.table$Frequency <- precise.is.sn.list$indels.size.table[precise.is.sn.list$indels.size.table$Size %in% correct.direction.size.range,]$Frequency +
    imprecise.is.sn.list$indels.size.table[imprecise.is.sn.list$indels.size.table$Size %in% correct.direction.size.range,]$Frequency
  ki.correct.is.sn.list$substitution.number.table$Frequency <- precise.is.sn.list$substitution.number.table$Frequency
  ### reverse direction knock-in
  # Standard size is changed from reverse junction to perfect junction.
  rev.imprecise.is.sn.list$indels.size.table$Size <- factor(-100:100 + rc.nhej.insert.length)
  ki.reverse.is.sn.list <- rev.imprecise.is.sn.list

  return(list(ki.correct.is.sn.list = ki.correct.is.sn.list, ki.reverse.is.sn.list = ki.reverse.is.sn.list))
}

SaveKiIndelPlot <- function(correct.is.sn.list, reverse.is.sn.list, outdir, ki.type = c("left_ki", "right_ki"), nhej.insert.length, rc.nhej.insert.length){
  ki.type <- match.arg(ki.type)
  # make indel and substitution plot
  aggregated.correct.is.sn.table.list <- MakeAggregatedKiIsSnTableList(correct.is.sn.list, "correct", nhej.insert.length)
  aggregated.reverse.is.sn.table.list <- MakeAggregatedKiIsSnTableList(reverse.is.sn.list, "reverse", rc.nhej.insert.length)
  
  # redefine Ip as imprecise knock-in has same length of precise knock-in
  correct.ip.ind <- which(aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label == "Ip")
  correct.ip.freq <- aggregated.correct.is.sn.table.list$aggregated.indels.size.table$Frequency[correct.ip.ind] - sum(aggregated.correct.is.sn.table.list$aggregated.substitution.number.table$Frequency)
  aggregated.correct.is.sn.table.list$aggregated.indels.size.table <- aggregated.correct.is.sn.table.list$aggregated.indels.size.table[-c(correct.ip.ind), ]
  aggregated.correct.is.sn.table.list$aggregated.indels.size.table <- rbind(aggregated.correct.is.sn.table.list$aggregated.indels.size.table, data.frame(label = factor("Ip"), Frequency = correct.ip.freq))
  # only substitution number chart
  if(!is.null(aggregated.correct.is.sn.table.list$aggregated.substitution.number.table)){
    
    aggregated.substitution.number.class.count.table <- data.frame(
      Class=aggregated.correct.is.sn.table.list$aggregated.substitution.number.table$label
      , Frequency=aggregated.correct.is.sn.table.list$aggregated.substitution.number.table$Frequency
      , Percentage=aggregated.correct.is.sn.table.list$aggregated.substitution.number.table$Frequency / sum(aggregated.correct.is.sn.table.list$aggregated.substitution.number.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.substitution.number.class.count.table, file = file.path(outdir, paste("[ⅱa]", ki.type, "_Rate_of_knock-in_substitution_number.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.substitution.number.class.count.table, file = file.path(outdir, paste("[ⅱa]", ki.type, "_Rate_of_knock-in_substitution_number.rds", sep="")))
    SavePieChart(aggregated.substitution.number.class.count.table
      , c("#F29F05", "#7F5124", "#BF7936")
      , "Reads"
      , file.path(outdir, paste("[ⅱa]Rate_of_", ki.type, "_knock-in_substitution_number_piechart.png", sep="")))
  }
  # only indels size pie chart (without precise knock-in)
  aggregated.mutation.merge.size.table <- rbind(aggregated.correct.is.sn.table.list$aggregated.indels.size.table, aggregated.reverse.is.sn.table.list$aggregated.indels.size.table)
  if(!is.null(aggregated.mutation.merge.size.table) & sum(aggregated.mutation.merge.size.table$Frequency) > 0){
  
    aggregated.indels.size.class.count.table <- data.frame(
      Class=aggregated.mutation.merge.size.table$label
      , Frequency=aggregated.mutation.merge.size.table$Frequency
      , Percentage=aggregated.mutation.merge.size.table$Frequency / sum(aggregated.mutation.merge.size.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅱb]", ki.type, "_Rate_of_knock-in_indel_size.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.indels.size.class.count.table, file = file.path(outdir, paste("[ⅱb]", ki.type, "_Rate_of_knock-in_indel_size.rds", sep="")))
    SavePieChart(aggregated.indels.size.class.count.table
      , c("#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
      , "Reads"
      , file.path(outdir, paste("[ⅱb]Rate_of_", ki.type, "_knock-in_indel_size_piechart.png", sep="")))
  }

  # substitution number + indels size pie chart
  # make merge aggregated table
  aggregated.mutation.merge.size.number.table <- rbind(aggregated.correct.is.sn.table.list$aggregated.substitution.number.table, aggregated.mutation.merge.size.table)
  if(!is.null(aggregated.mutation.merge.size.number.table) & sum(aggregated.mutation.merge.size.number.table$Frequency) > 0){
    
    aggregated.merge.size.number.class.count.table <- data.frame(
      Class=aggregated.mutation.merge.size.number.table$label
      , Frequency=aggregated.mutation.merge.size.number.table$Frequency
      , Percentage=aggregated.mutation.merge.size.number.table$Frequency / sum(aggregated.mutation.merge.size.number.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.merge.size.number.class.count.table, file = file.path(outdir, paste("[ⅱc]", ki.type, "_Rate_of_knock-in_substitution_and_indels.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.merge.size.number.class.count.table, file = file.path(outdir, paste("[ⅱc]", ki.type, "_Rate_of_knock-in_substitution_and_indels.rds", sep="")))
    SavePieChart(aggregated.merge.size.number.class.count.table
      , c("#F29F05", "#7F5124", "#BF7936", "#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
      , "Reads"
      , file.path(outdir, paste("[ⅱc]Rate_of_", ki.type, "_knock-in_substitution_and_indels_piechart.png", sep="")))
  }

  # substitution number + indels size pie chart without precise knock-in
  # remove P_S0 row
  aggregated.mutation.merge.size.number.without.p_s0.table <- aggregated.mutation.merge.size.number.table[which(aggregated.mutation.merge.size.number.table$label != "P_S0"), ]
  if(!is.null(aggregated.mutation.merge.size.number.without.p_s0.table) & sum(aggregated.mutation.merge.size.number.without.p_s0.table$Frequency) > 0){
    
    aggregated.merge.size.number.without.p_s0.class.count.table <- data.frame(
      Class=aggregated.mutation.merge.size.number.without.p_s0.table$label
      , Frequency=aggregated.mutation.merge.size.number.without.p_s0.table$Frequency
      , Percentage=aggregated.mutation.merge.size.number.without.p_s0.table$Frequency / sum(aggregated.mutation.merge.size.number.without.p_s0.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.merge.size.number.without.p_s0.class.count.table, file = file.path(outdir, paste("[ⅱd]", ki.type, "_Rate_of_knock-in_substitution_and_indels(without_precise_ki).txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.merge.size.number.without.p_s0.class.count.table, file = file.path(outdir, paste("[ⅱd]", ki.type, "_Rate_of_knock-in_substitution_and_indels(without_precise_ki).rds", sep="")))
    SavePieChart(aggregated.merge.size.number.without.p_s0.class.count.table
      , c("#7F5124", "#BF7936", "#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
      , "Reads"
      , file.path(outdir, paste("[ⅱd]Rate_of_", ki.type, "_knock-in_substitution_and_indels_piechart(without_precise_ki).png", sep=""))
    )
  }
}

RunKiTotalIndelAnalysis <- function(precise.is.sn.list, imprecise.is.sn.list, rev.imprecise.is.sn.list, nhej.insert.length, rc.nhej.insert.length, save.dir, ki.type = c("left_ki", "right_ki")){
  ki.type <- match.arg(ki.type)
  ki.is.sn.list <- BindKiIsSnList(precise.is.sn.list
    , imprecise.is.sn.list
    , rev.imprecise.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )

  # make indel size and substitution number pie chart
  if(!is.null(ki.is.sn.list$ki.correct.is.sn.list) | !is.null(ki.is.sn.list$ki.reverse.is.sn.list)){

    # make indel size bar plot
    if(!is.null(ki.is.sn.list$ki.correct.is.sn.list)){
      saveRDS(ki.is.sn.list$ki.correct.is.sn.list$indels.size.table , file = file.path(save.dir, paste0("[ⅰa]Distribution_of_indel_size_on_correct_", ki.type, "_junction.table.rds")))
      SaveIndelSizeReadBarPlot(ki.is.sn.list$ki.correct.is.sn.list$indels.size.table, file.path(save.dir, paste0("[ⅰa]Distribution_of_indel_size_on_correct_", ki.type, "_junction.png")))
    }
    # make substitution number bar plot
    if(!is.null(ki.is.sn.list$ki.correct.is.sn.list$substitution.number.table)){
      saveRDS(ki.is.sn.list$ki.correct.is.sn.list$substitution.number.table , file = file.path(save.dir, paste0("[ⅰb]Distribution_of_substitution_number_on_correct_", ki.type, "_junction.table.rds")))
      SaveSubstitutionNumberReadBarPlot(ki.is.sn.list$ki.correct.is.sn.list$substitution.number.table, file.path(save.dir, paste0("[ⅰb]Distribution_of_substitution_number_on_correct_", ki.type, "_junction.png")))
    }

    SaveKiIndelPlot(ki.is.sn.list$ki.correct.is.sn.list, ki.is.sn.list$ki.reverse.is.sn.list, outdir = save.dir, ki.type = ki.type, nhej.insert.length = nhej.insert.length, rc.nhej.insert.length = rc.nhej.insert.length)
  }
}

# indels size pie chart by EJ type
SaveInDelByTypePieChart <- function(aggregated.indels.size.table, ej.type, col.vec, type = c("mut", "left_ki", "right_ki"), outdir){
  type <- match.arg(type)
  if(type == "mut"){
    id.number <- "[ⅳ"
  }else if(type %in% c("left_ki", "right_ki")){
    id.number <- "[ⅲ"
  }
  if(!is.null(aggregated.indels.size.table) & sum(aggregated.indels.size.table$Frequency) > 0){
    aggregated.indels.size.class.count.table <- data.frame(
      Class=aggregated.indels.size.table$label
      , Frequency=aggregated.indels.size.table$Frequency
      , Percentage=aggregated.indels.size.table$Frequency / sum(aggregated.indels.size.table$Frequency) * 100
    )
    options(warn=-1) # "appending column names to file" is no problem.
    write.table(aggregated.indels.size.class.count.table, file = file.path(outdir, paste(id.number, "]Rate_of_", ej.type, "_", type, "_indel_size.txt", sep=""))
      , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
    options(warn=0)
    saveRDS(aggregated.indels.size.class.count.table, file = file.path(outdir, paste(id.number, "]Rate_of_", ej.type, "_", type, "_indel_size.rds", sep="")))
    SavePieChart(aggregated.indels.size.class.count.table
      , col.vec
      , "Reads"
      , file.path(outdir, paste(id.number, "]Rate_of_", ej.type, "_", type, "_indel_size_piechart.png", sep="")))
  }
}

# indels size pie chart by indel size
SaveInDelBySizePieChart <- function(
  mmej.aggregated.indels.size.table
  , nhej.aggregated.indels.size.table
  , unidentified.deletion.aggregated.indels.size.table
  , insert.dominant.aggregated.indels.size.table
  , size.label, type = c("mut", "left_ki", "right_ki"), outdir){
    type <- match.arg(type)
    if(type == "mut"){
      id.number <- "[ⅴ"
    }else if(type %in% c("left_ki", "right_ki")){
      id.number <- "[ⅳ"
    }
    if(!is.null(mmej.aggregated.indels.size.table)
      & !is.null(nhej.aggregated.indels.size.table)
      & !is.null(unidentified.deletion.aggregated.indels.size.table)
      & !is.null(insert.dominant.aggregated.indels.size.table)
    ){
    mmej.label <- grep(paste0("_", size.label, "$"), mmej.aggregated.indels.size.table$label, value=TRUE)
    nhej.label <- grep(paste0("_", size.label, "$"), nhej.aggregated.indels.size.table$label, value=TRUE)
    unidentified.deletion.label <- grep(paste0("_", size.label, "$"), unidentified.deletion.aggregated.indels.size.table$label, value=TRUE)
    insert.dominant.label <- grep(paste0("_", size.label, "$"), insert.dominant.aggregated.indels.size.table$label, value=TRUE)
    
    ej.type.by.size.table <- rbind(
      mmej.aggregated.indels.size.table[mmej.aggregated.indels.size.table$label %in% mmej.label,]
      , nhej.aggregated.indels.size.table[nhej.aggregated.indels.size.table$label %in% nhej.label,]
      , unidentified.deletion.aggregated.indels.size.table[unidentified.deletion.aggregated.indels.size.table$label %in% unidentified.deletion.label,]
      , insert.dominant.aggregated.indels.size.table[insert.dominant.aggregated.indels.size.table$label %in% insert.dominant.label,]
    )
    if(!is.null(ej.type.by.size.table)){
      options(warn=-1) # "appending column names to file" is no problem.
      write.table(ej.type.by.size.table, file = file.path(outdir, paste(id.number, "]Rate_of_", size.label, "_", type, "_endjoining_type.txt", sep=""))
        , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
      options(warn=0)
      saveRDS(ej.type.by.size.table, file = file.path(outdir, paste(id.number, "]Rate_of_", size.label, "_", type, "_endjoining_type.rds", sep="")))
      if(sum(ej.type.by.size.table$Frequency) > 0){
        ej.type.by.size.class.count.table <- data.frame(
          Class=ej.type.by.size.table$label
          , Frequency=ej.type.by.size.table$Frequency
          , Percentage=ej.type.by.size.table$Frequency / sum(ej.type.by.size.table$Frequency) * 100
        )
        options(warn=-1) # "appending column names to file" is no problem.
        write.table(ej.type.by.size.table, file = file.path(outdir, paste(id.number, "]Rate_of_", size.label, "_", type, "_endjoining_type.txt", sep=""))
          , quote=FALSE, col.names=TRUE, row.names=FALSE,append=TRUE, sep = ",")
        options(warn=0)
        saveRDS(ej.type.by.size.table, file = file.path(outdir, paste(id.number, "]Rate_of_", size.label, "_", type, "_endjoining_type.rds", sep="")))
        SavePieChart(ej.type.by.size.class.count.table
          , c("#ACF2F2", "#F2A0DC", "#D288F2", "#F2C9F0")
          , "Reads"
          , file.path(outdir, paste(id.number, "]Rate_of_", size.label, "_", type, "_endjoining_type_piechart.png", sep="")))
      }
    }
  }
}

SaveMutEjTypePlot = function(is.sn.list, outdir){
  # make indel and substitution plot
  mmej.aggregated.is.sn.table.list <- MakeAggregatedIsSnTableList(is.sn.list$mmej.is.sn.list)
  mmej.aggregated.is.sn.table.list$aggregated.indels.size.table$label <- paste0("MMEJ_", mmej.aggregated.is.sn.table.list$aggregated.indels.size.table$label)
  nhej.aggregated.is.sn.table.list <- MakeAggregatedIsSnTableList(is.sn.list$nhej.is.sn.list)
  nhej.aggregated.is.sn.table.list$aggregated.indels.size.table$label <- paste0("NHEJ_", nhej.aggregated.is.sn.table.list$aggregated.indels.size.table$label)
  unidentified.deletion.aggregated.is.sn.table.list <- MakeAggregatedIsSnTableList(is.sn.list$unidentified.deletion.is.sn.list)
  unidentified.deletion.aggregated.is.sn.table.list$aggregated.indels.size.table$label <- paste0("Unidentified_Deletion_", unidentified.deletion.aggregated.is.sn.table.list$aggregated.indels.size.table$label)
  insert.dominant.aggregated.is.sn.table.list <- MakeAggregatedIsSnTableList(is.sn.list$insert.dominant.is.sn.list)
  insert.dominant.aggregated.is.sn.table.list$aggregated.indels.size.table$label <- paste0("Insert_Dominant_", insert.dominant.aggregated.is.sn.table.list$aggregated.indels.size.table$label)
  # indels size pie chart by EJ type
  size.label.col.vec <- c("#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304", "#D99ABC")
  SaveInDelByTypePieChart(mmej.aggregated.is.sn.table.list$aggregated.indels.size.table, "MMEJ", size.label.col.vec, "mut", outdir)
  SaveInDelByTypePieChart(nhej.aggregated.is.sn.table.list$aggregated.indels.size.table, "NHEJ", size.label.col.vec, "mut", outdir)
  SaveInDelByTypePieChart(unidentified.deletion.aggregated.is.sn.table.list$aggregated.indels.size.table, "Unidentified_Deletion", size.label.col.vec, "mut", outdir)
  SaveInDelByTypePieChart(insert.dominant.aggregated.is.sn.table.list$aggregated.indels.size.table, "Insert_Dominant", size.label.col.vec, "mut", outdir)
  # indels size pie chart by indel size
  for(size.label in c("D1", "D2", "D>2", "I1", "I>1", "NS")){
    SaveInDelBySizePieChart(
      mmej.aggregated.is.sn.table.list$aggregated.indels.size.table
      , nhej.aggregated.is.sn.table.list$aggregated.indels.size.table
      , unidentified.deletion.aggregated.is.sn.table.list$aggregated.indels.size.table
      , insert.dominant.aggregated.is.sn.table.list$aggregated.indels.size.table
      , size.label, "mut", outdir
    )
  }
}

SaveImKiEjTypePlot = function(null.precise.is.sn.list.set, imprecise.is.sn.list.set, rev.imprecise.is.sn.list.set
  , nhej.insert.length, rc.nhej.insert.length, outdir, ki.type = c("left_ki", "right_ki")){
  ki.type <- match.arg(ki.type)
  # bind imprecise knock-in with reverse imprecise knock-in
  # (precise knock-in is not target of this analysis. but, add null precise knock-in data because BindKiIsSnList require precise knock-in data also.)
  mmej.ki.is.sn.list <- BindKiIsSnList(null.precise.is.sn.list.set$mmej.is.sn.list
    , imprecise.is.sn.list.set$mmej.is.sn.list
    , rev.imprecise.is.sn.list.set$mmej.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )
  nhej.ki.is.sn.list <- BindKiIsSnList(null.precise.is.sn.list.set$nhej.is.sn.list
    , imprecise.is.sn.list.set$nhej.is.sn.list
    , rev.imprecise.is.sn.list.set$nhej.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )
  unidentified.deletion.ki.is.sn.list <- BindKiIsSnList(null.precise.is.sn.list.set$unidentified.deletion.is.sn.list
    , imprecise.is.sn.list.set$unidentified.deletion.is.sn.list
    , rev.imprecise.is.sn.list.set$unidentified.deletion.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )
  insert.dominant.ki.is.sn.list <- BindKiIsSnList(null.precise.is.sn.list.set$insert.dominant.is.sn.list
    , imprecise.is.sn.list.set$insert.dominant.is.sn.list
    , rev.imprecise.is.sn.list.set$insert.dominant.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )

  # make classification by ki length (substitution will be ignored.)
  # In RunKiTotalIndelAnalysis(), we conduct Ip - (P_s0 + P_s1 + P_s>1), but in this case, we dont need it because precise knock-in is not included.
  mmej.aggregated.correct.is.sn.table.list <- MakeAggregatedKiIsSnTableList(mmej.ki.is.sn.list$ki.correct.is.sn.list, "correct", nhej.insert.length)
  mmej.aggregated.reverse.is.sn.table.list <- MakeAggregatedKiIsSnTableList(mmej.ki.is.sn.list$ki.reverse.is.sn.list, "reverse", rc.nhej.insert.length)
  nhej.aggregated.correct.is.sn.table.list <- MakeAggregatedKiIsSnTableList(nhej.ki.is.sn.list$ki.correct.is.sn.list, "correct", nhej.insert.length)
  nhej.aggregated.reverse.is.sn.table.list <- MakeAggregatedKiIsSnTableList(nhej.ki.is.sn.list$ki.reverse.is.sn.list, "reverse", rc.nhej.insert.length)
  unidentified.deletion.aggregated.correct.is.sn.table.list <- MakeAggregatedKiIsSnTableList(unidentified.deletion.ki.is.sn.list$ki.correct.is.sn.list, "correct", nhej.insert.length)
  unidentified.deletion.aggregated.reverse.is.sn.table.list <- MakeAggregatedKiIsSnTableList(unidentified.deletion.ki.is.sn.list$ki.reverse.is.sn.list, "reverse", rc.nhej.insert.length)
  insert.dominant.aggregated.correct.is.sn.table.list <- MakeAggregatedKiIsSnTableList(insert.dominant.ki.is.sn.list$ki.correct.is.sn.list, "correct", nhej.insert.length)
  insert.dominant.aggregated.reverse.is.sn.table.list <- MakeAggregatedKiIsSnTableList(insert.dominant.ki.is.sn.list$ki.reverse.is.sn.list, "reverse", rc.nhej.insert.length)
  # rename
  mmej.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label <- paste0("MMEJ_", mmej.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label)
  nhej.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label <- paste0("NHEJ_", nhej.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label)
  unidentified.deletion.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label <- paste0("Unidentified_Deletion_", unidentified.deletion.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label)
  insert.dominant.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label <- paste0("Insert_Dominant_", insert.dominant.aggregated.correct.is.sn.table.list$aggregated.indels.size.table$label)
  mmej.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label <- paste0("MMEJ_", mmej.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label)
  nhej.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label <- paste0("NHEJ_", nhej.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label)
  unidentified.deletion.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label <- paste0("Unidentified_Deletion_", unidentified.deletion.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label)
  insert.dominant.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label <- paste0("Insert_Dominant_", insert.dominant.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table$label)
  # merge correct direction table with reverse direction table
  mmej.aggregated.im.ki.merge.size.table <- rbind(mmej.aggregated.correct.is.sn.table.list$aggregated.indels.size.table, mmej.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table)
  nhej.aggregated.im.ki.merge.size.table <- rbind(nhej.aggregated.correct.is.sn.table.list$aggregated.indels.size.table, nhej.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table)
  unidentified.deletion.aggregated.im.ki.merge.size.table <- rbind(unidentified.deletion.aggregated.correct.is.sn.table.list$aggregated.indels.size.table, unidentified.deletion.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table)
  insert.dominant.aggregated.im.ki.merge.size.table <- rbind(insert.dominant.aggregated.correct.is.sn.table.list$aggregated.indels.size.table, insert.dominant.aggregated.reverse.is.sn.table.list$aggregated.indels.size.table)
  # indels size pie chart by EJ type
  size.label.col.vec <- c("#C0D904", "#277A8C", "#3E5902", "#60A6A6", "#618C03", "#95D904", "#F2EC9B", "#B3D9C0") # color order  
  SaveInDelByTypePieChart(mmej.aggregated.im.ki.merge.size.table, "MMEJ", size.label.col.vec, ki.type, outdir)
  SaveInDelByTypePieChart(nhej.aggregated.im.ki.merge.size.table, "NHEJ", size.label.col.vec, ki.type, outdir)
  SaveInDelByTypePieChart(unidentified.deletion.aggregated.im.ki.merge.size.table, "Unidentified_Deletion", size.label.col.vec, ki.type, outdir)
  SaveInDelByTypePieChart(insert.dominant.aggregated.im.ki.merge.size.table, "Insert_Dominant", size.label.col.vec, ki.type, outdir)
  # indels size pie chart by indel size
  for(size.label in c("<Ip", ">Id", "Id", "Id-Ip", "Ip", "<Rd", ">Rd", "Rd")){
    SaveInDelBySizePieChart(
      mmej.aggregated.im.ki.merge.size.table
      , nhej.aggregated.im.ki.merge.size.table
      , unidentified.deletion.aggregated.im.ki.merge.size.table
      , insert.dominant.aggregated.im.ki.merge.size.table
      , size.label, ki.type, outdir
    )
  }
}

SaveImKiIndelSizePlot = function(null.precise.is.sn.list.set, imprecise.is.sn.list.set, rev.imprecise.is.sn.list.set
  , nhej.insert.length, rc.nhej.insert.length, outdir, ki.type = c("left_ki", "right_ki")){
  ki.type <- match.arg(ki.type)
  if(ki.type == "left_ki"){
    junction.name <- "left_junction"
  }else{
    junction.name <- "right_junction"
  }
  # bind imprecise knock-in with reverse imprecise knock-in
  # (precise knock-in is not target of this analysis. but, add null precise knock-in data because BindKiIsSnList require precise knock-in data also.)
  mmej.ki.is.sn.list <- BindKiIsSnList(null.precise.is.sn.list.set$mmej.is.sn.list
    , imprecise.is.sn.list.set$mmej.is.sn.list
    , rev.imprecise.is.sn.list.set$mmej.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )
  nhej.ki.is.sn.list <- BindKiIsSnList(null.precise.is.sn.list.set$nhej.is.sn.list
    , imprecise.is.sn.list.set$nhej.is.sn.list
    , rev.imprecise.is.sn.list.set$nhej.is.sn.list
    , nhej.insert.length
    , rc.nhej.insert.length
  )
  if(!is.null(mmej.ki.is.sn.list$ki.correct.is.sn.list$indels.size.table)){
    saveRDS(mmej.ki.is.sn.list$ki.correct.is.sn.list$indels.size.table , file = file.path(outdir, paste0("[iiib]Distribution_of_MMEJ_deletion_size_on_correct_", junction.name, ".table.rds")))
    SaveIndelSizeReadBarPlot(mmej.ki.is.sn.list$ki.correct.is.sn.list$indels.size.table, file.path(outdir, paste0("[iiib]Distribution_of_MMEJ_deletion_size_on_correct_", junction.name, ".png")))
  }
  if(!is.null(nhej.ki.is.sn.list$ki.correct.is.sn.list$indels.size.table)){
    saveRDS(nhej.ki.is.sn.list$ki.correct.is.sn.list$indels.size.table , file = file.path(outdir, paste0("[iiic]Distribution_of_NHEJ_deletion_size_on_correct_", junction.name, ".table.rds")))
    SaveIndelSizeReadBarPlot(nhej.ki.is.sn.list$ki.correct.is.sn.list$indels.size.table, file.path(outdir, paste0("[iiic]Distribution_of_NHEJ_deletion_size_on_correct_", junction.name, ".png")))
  }
  if(!is.null(mmej.ki.is.sn.list$ki.reverse.is.sn.list$indels.size.table)){
    saveRDS(mmej.ki.is.sn.list$ki.reverse.is.sn.list$indels.size.table , file = file.path(outdir, paste0("[iiid]Distribution_of_MMEJ_deletion_size_on_reverse_", junction.name, ".table.rds")))
    SaveIndelSizeReadBarPlot(mmej.ki.is.sn.list$ki.reverse.is.sn.list$indels.size.table, file.path(outdir, paste0("[iiid]Distribution_of_MMEJ_deletion_size_on_reverse_", junction.name, ".png")))
  }
  if(!is.null(nhej.ki.is.sn.list$ki.reverse.is.sn.list$indels.size.table)){
    saveRDS(nhej.ki.is.sn.list$ki.reverse.is.sn.list$indels.size.table , file = file.path(outdir, paste0("[iiie]Distribution_of_NHEJ_deletion_size_on_reverse_", junction.name, ".table.rds")))
    SaveIndelSizeReadBarPlot(nhej.ki.is.sn.list$ki.reverse.is.sn.list$indels.size.table, file.path(outdir, paste0("[iiie]Distribution_of_NHEJ_deletion_size_on_reverse_", junction.name, ".png")))
  }

}


# Make micromology length - count table
MakeMicrohomologyLengthReadTable <- function(variants.count.total.list, collect.type = c("mmej", "nhej"), range.vec = 0:100){
  collect.type <- match.arg(collect.type)
  if(collect.type == "mmej"){
    id.type.vec <- c("three.prime.conserved.mmej", "five.prime.conserved.mmej", "both.prime.conserved.mmej")
  }else if(collect.type == "nhej"){
    id.type.vec <- c("nhej")
  }
  num.vec <- numeric(0)
  if(length(variants.count.total.list) == 0){
    return(table(factor(, levels = range.vec)))
  }
  for(mut.ind in 1:(length(variants.count.total.list)/4)){
    indel.seq.info.list <- variants.count.total.list[mut.ind * 4][[1]]
    if(indel.seq.info.list$ej.type.label %in% id.type.vec){
      num.vec <- c(num.vec
        , (rep(
            nchar(indel.seq.info.list$mh.seq.vec[1])
            , variants.count.total.list[mut.ind * 4 - 2][[1]]
          ))
      )
    }
  }
  return(table(factor(num.vec, levels = range.vec)))
}

# Make trimmed length - count table
MakeTrimmedLengthReadTable <- function(variants.count.total.list, collect.type = c("mmej", "nhej"), range.vec = 0:100){
  collect.type <- match.arg(collect.type)
  if(collect.type == "mmej"){
    id.type.vec <- c("three.prime.conserved.mmej", "five.prime.conserved.mmej", "both.prime.conserved.mmej")
  }else if(collect.type == "nhej"){
    id.type.vec <- c("nhej")
  }
  num.vec <- numeric(0)
  if(length(variants.count.total.list) == 0){
    return(table(factor(, levels = range.vec)))
  }
  for(mut.ind in 1:(length(variants.count.total.list)/4)){
    indel.seq.info.list <- variants.count.total.list[mut.ind * 4][[1]]
    if(indel.seq.info.list$ej.type.label %in% id.type.vec){
      num.vec <- c(num.vec
        , (rep(
            nchar(indel.seq.info.list$trimmed.seq.vec[1])
            , variants.count.total.list[mut.ind * 4 - 2][[1]]
          ))
      )
    }
  }
  return(table(factor(num.vec, levels = range.vec)))
}

MakeMlcTlcReadTable <- function(variants.count.total.list, microhomology.range.vec = 0:100, trimmedseq.range.vec = 0:100){ # Mlc: micromology length - read count, Tlc: trimmed sequence length - read count
  # split variants.count.total.list
  split.variants.count.total.list <- SplitTotalMutListByEjType(variants.count.total.list)
  mmej.variants.count.total.list <- split.variants.count.total.list$mmej.variants.count.total.list
  nhej.variants.count.total.list <- split.variants.count.total.list$nhej.variants.count.total.list
  # micromology length - read count
  mmej.mh.size.table <- data.frame(MakeMicrohomologyLengthReadTable(mmej.variants.count.total.list, collect.type = "mmej"), stringsAsFactors=FALSE)
  colnames(mmej.mh.size.table) <- c("Length", "Frequency")
  nhej.mh.size.table <- data.frame(MakeMicrohomologyLengthReadTable(nhej.variants.count.total.list, collect.type = "nhej"), stringsAsFactors=FALSE)
  colnames(nhej.mh.size.table) <- c("Length", "Frequency")
  
  # trimmed sequence length - read count
  mmej.trim.size.table <- data.frame(MakeTrimmedLengthReadTable(mmej.variants.count.total.list, collect.type = "mmej"), stringsAsFactors=FALSE)
  colnames(mmej.trim.size.table) <- c("Length", "Frequency")
  nhej.trim.size.table <- data.frame(MakeTrimmedLengthReadTable(nhej.variants.count.total.list, collect.type = "nhej"), stringsAsFactors=FALSE)
  colnames(nhej.trim.size.table) <- c("Length", "Frequency")

  mlc.tlc.list <- list(
    mmej.mlc.tlc.list = list(mh.size.table = mmej.mh.size.table, trim.size.table = mmej.trim.size.table)
    , nhej.mlc.tlc.list = list(mh.size.table = nhej.mh.size.table, trim.size.table = nhej.trim.size.table)
  )

  return(mlc.tlc.list)
}

SaveSeqLengthReadBarPlot <- function(length.count.table, file.name, seq.type = c("Microhomology", "Trimmed_seq")){
  seq.type <- match.arg(seq.type)
  options(warn = -1)
  length.count.p <- ggplot(data=length.count.table
    , aes(x=Length, y=Frequency, fill = c(seq.type))) +
    geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#D93240") +
    coord_cartesian(ylim = c(0, max(length.count.table$Frequency) + 100)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    , axis.text.y = element_text(size = 6)
    , axis.title.x = element_text(size = 12)
    , axis.title.y = element_text(size = 12)) +
  xlab(paste0(seq.type, " size[bp]")) +
  ylab("Number of reads")
  # plot(length.count.p)
  ggsave(file = file.name, plot = length.count.p, dpi = 300, width = 16, height = 4)
  options(warn = 0)
}

# Make micromology length - count plot
SaveMlTlScatterPlot <- function(variants.count.total.list, save.dir, type = c("mut", "left_ki", "right_ki")){
  type <- match.arg(type)
  if(type == "mut"){
    id.number <- "[ⅸ"
  }else if(type %in% c("left_ki", "right_ki")){
    id.number <- "[ⅶ"
  }
  id.type.vec <- c("three.prime.conserved.mmej", "five.prime.conserved.mmej", "both.prime.conserved.mmej")
  plot.df <- data.frame(mut.label = character(0), mh.len = numeric(0), trim.len = numeric(0))
  if(length(variants.count.total.list) == 0){
    return()
  }
  for(mut.ind in 1:(length(variants.count.total.list)/4)){
    indel.seq.info.list <- variants.count.total.list[mut.ind * 4][[1]]
    if(indel.seq.info.list$ej.type.label %in% id.type.vec){
      plot.df <- rbind(plot.df
        , data.frame(
          mut.label = variants.count.total.list[mut.ind * 4 - 3][[1]]
          , mh.len = nchar(indel.seq.info.list$mh.seq.vec[1])
          , trim.len = nchar(indel.seq.info.list$trimmed.seq.vec[1])
        )
      )
    }
  }
  saveRDS(plot.df , file = file.path(save.dir, paste0(id.number, "]Distribution_of_", type, "_MMEJ_Microhomology_length_Trimmed_Seq_length.table.rds")))
  if(nrow(plot.df) > 0){
    options(warn=-1)
    ml.tl.scatter.p <- ggplot(data = plot.df, aes(x = mh.len, y = trim.len, color = "#BF9E75")) +
    geom_point() +
    xlim(0, 30) +
    ylim(0, 50) +
    labs(x = "Length of Microhomology [bp]", y = "Length of Trimmed Sequence [bp]") +
    theme(legend.position = "none") +
    geom_smooth(color = "black", size = 0.8, linetype = 2) # loess model
    ggsave(file = file.path(save.dir, paste0(id.number, "]Distribution_of_", type, "_MMEJ_Microhomology_length_Trimmed_Seq_length.png")), plot = ml.tl.scatter.p, dpi = 300, width = 16, height = 4)
    options(warn=0)
  }
}

# Make microhomology sequence-frequency-count plot
SaveMsFreqPlot <- function(variants.count.total.list, save.dir, type = c("mut", "left_ki", "right_ki")){
  type <- match.arg(type)
  if(type == "mut"){
    id.number <- "[ⅹ"
  }else if(type %in% c("left_ki", "right_ki")){
    id.number <- "[ⅷ"
  }
  id.type.vec <- c("three.prime.conserved.mmej", "five.prime.conserved.mmej", "both.prime.conserved.mmej")
  freq.df <- data.frame(
    mut.label = character(0), read.cnt = numeric(0), mh.len = numeric(0), trimmed.len = numeric(0)
    , A = numeric(0), C = numeric(0), G = numeric(0), T = numeric(0)
    , AA = numeric(0), AC = numeric(0), AG = numeric(0), AT = numeric(0)
    , CA = numeric(0), CC = numeric(0), CG = numeric(0), CT = numeric(0)
    , GA = numeric(0), GC = numeric(0), GG = numeric(0), GT = numeric(0)
    , TA = numeric(0), TC = numeric(0), TG = numeric(0), TT = numeric(0)
  )
  if(length(variants.count.total.list) == 0){
    return()
  }
  for(mut.ind in 1:(length(variants.count.total.list)/4)){
    indel.seq.info.list <- variants.count.total.list[mut.ind * 4][[1]]
    if(indel.seq.info.list$ej.type.label %in% id.type.vec){
      # print(cbind(indel.seq.info.list$mononuc.freq.table, indel.seq.info.list$dinuc.freq.table))
      freq.df <- rbind(freq.df
        , data.frame(
          mut.label = variants.count.total.list[mut.ind * 4 - 3][[1]]
          , read.cnt = variants.count.total.list[mut.ind * 4 - 2][[1]]
          , mh.len = nchar(indel.seq.info.list$mh.seq.vec[1])
          , trimmed.len = nchar(indel.seq.info.list$trimmed.seq.vec[1])
          , cbind(indel.seq.info.list$mononuc.freq.table, indel.seq.info.list$dinuc.freq.table)
        )
      )
    }
  }

  if(nrow(freq.df) > 0){
    # aggregated frequency of microhomology nucleotide - read.cnt
    total.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:ncol(freq.df)], levels = colnames(freq.df)[5:ncol(freq.df)])
      ,Frequency = colSums(t(apply(freq.df, MARGIN = 1, function(row){as.numeric(row[5:length(row)]) * as.numeric(row[2])})))
    )
    options(warn=-1)
    total.sum.by.nucleotide.p <- ggplot(data=total.sum.by.nucleotide.table
      , aes(x=Nucleotide, y=Frequency, fill = c("microhomology"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8
        , fill = c("#00931F", "#001298", "#000000", "#FE0000"
          ,"#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
          ,"#009840", "#104277", "#25217D", "#005E6A"
          ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
          ,"#F28322", "#4557BF", "#400101", "#D93223")
      ) +
    coord_cartesian(ylim = c(0, max(total.sum.by.nucleotide.table$Frequency[1:4]) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Nucleotide in microhomology") +
    ylab("Total frequency in reads")
    # plot(substitution.number.p)
    saveRDS(total.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "a]Distribution_of_", type, "_MMEJ_Total_Microhomology_Sequence_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "a]Distribution_of_", type, "_MMEJ_Total_Microhomology_Sequence_Frequency.png")), plot = total.sum.by.nucleotide.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # aggregated frequency of mono microhomology nucleotide per length - read.cnt
    mono.mh.freq.df <- t(apply(freq.df, MARGIN = 1, function(row){
        if(as.numeric(row[3]) == 1){
          as.numeric(row[5:8]) * as.numeric(row[2])
        }else{
          numeric(4)
        }
    }))
    mono.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:8], levels = colnames(freq.df)[5:8])
      ,Frequency = colSums(mono.mh.freq.df)
    )
    options(warn=-1)
    mono.mh.sum.by.nucleotide.p <- ggplot(data=mono.mh.sum.by.nucleotide.table
      , aes(x=Nucleotide, y=Frequency, fill = c("microhomology"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8
        , fill = c("#00931F", "#001298", "#000000", "#FE0000")
      ) +
    coord_cartesian(ylim = c(0, max(mono.mh.sum.by.nucleotide.table$Frequency[1:4]) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Nucleotide in mono-microhomology") +
    ylab("Number of reads")
    # plot(substitution.number.p)
    saveRDS(mono.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "b]Distribution_of_", type, "_MMEJ_Mono_Microhomology_Sequence_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "b]Distribution_of_", type, "_MMEJ_Mono_Microhomology_Sequence_Frequency.png")), plot = mono.mh.sum.by.nucleotide.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # aggregated frequency of di microhomology nucleotide per length - read.cnt
    di.mh.freq.df <- t(apply(freq.df, MARGIN = 1, function(row){
        if(as.numeric(row[3]) == 2){
          as.numeric(row[9:length(row)]) * as.numeric(row[2])
        }else{
          numeric(length(row) - 8)
        }
    }))
    di.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[9:ncol(freq.df)], levels = colnames(freq.df)[9:ncol(freq.df)])
      ,Frequency = colSums(di.mh.freq.df)
    )
    options(warn=-1)
    di.mh.sum.by.nucleotide.p <- ggplot(data=di.mh.sum.by.nucleotide.table
      , aes(x=Nucleotide, y=Frequency, fill = c("microhomology"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8
        , fill = c("#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
          ,"#009840", "#104277", "#25217D", "#005E6A"
          ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
          ,"#F28322", "#4557BF", "#400101", "#D93223")
      ) +
    coord_cartesian(ylim = c(0, max(di.mh.sum.by.nucleotide.table$Frequency[1:16]) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Nucleotide in di-microhomology") +
    ylab("Number of reads")
    # plot(substitution.number.p)
    saveRDS(di.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "c]Distribution_of_", type, "_MMEJ_Di_Microhomology_Sequence_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "c]Distribution_of_", type, "_MMEJ_Di_Microhomology_Sequence_Frequency.png")), plot = di.mh.sum.by.nucleotide.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # Make trimmed length-microhomology sequence-frequency-count 
    trim.len.mh.freq.list <- list()
    for(trimmed.len in 1:max(freq.df$trimmed.len)){
      trim.len.mh.freq.list <- c(trim.len.mh.freq.list, list(t(apply(freq.df, MARGIN = 1, function(row){
          if(as.numeric(row[4]) == trimmed.len){
            as.numeric(row[5:ncol(freq.df)]) * as.numeric(row[2])
          }else{
            numeric(ncol(freq.df) - 4)
          }
      })))
      )
    }
    trim.len.mh.freq.sum.table <- data.frame(matrix(unlist(lapply(trim.len.mh.freq.list, colSums)), nrow=length(lapply(trim.len.mh.freq.list, colSums)), byrow=T))
    colnames(trim.len.mh.freq.sum.table) <- colnames(freq.df[5:ncol(freq.df)])
    # Mono frequency
    mono.trim.len.mh.freq.sum.table <- trim.len.mh.freq.sum.table[, 1:4]
    mono.trim.len.mh.freq.sum.plot.table <- cbind(melt(mono.trim.len.mh.freq.sum.table), trim.len = factor(rep(1:nrow(mono.trim.len.mh.freq.sum.table), ncol(mono.trim.len.mh.freq.sum.table))))
    options(warn=-1)
    mono.trim.len.mh.freq.sum.plot.table %>% group_by(variable) %>% ggplot(aes(x = trim.len, y = value, fill = variable)) +
    geom_bar(stat="identity", position = position_dodge(0.8)) +
    scale_fill_manual(values = c("#00931F", "#001298", "#000000", "#FE0000")) +
    coord_cartesian(ylim = c(0, max(mono.trim.len.mh.freq.sum.plot.table$value) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Length of Intervening Sequence") +
    ylab("Total Frequency in Reads") -> mono.trim.len.mh.freq.sum.plot.p # inside frequency * number of read
    # plot(substitution.number.p)
    saveRDS(mono.trim.len.mh.freq.sum.plot.table , file = file.path(save.dir, paste0(id.number, "d]Distribution_of_", type, "_MMEJ_Intervening_Length_Microhomology_Sequence_Mono_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "d]Distribution_of_", type, "_MMEJ_Intervening_Length_Microhomology_Sequence_Mono_Frequency.png")), plot = mono.trim.len.mh.freq.sum.plot.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # Di frequency
    di.trim.len.mh.freq.sum.table <- trim.len.mh.freq.sum.table[, 5:20]
    di.trim.len.mh.freq.sum.plot.table <- cbind(melt(di.trim.len.mh.freq.sum.table), trim.len = factor(rep(1:nrow(di.trim.len.mh.freq.sum.table), ncol(di.trim.len.mh.freq.sum.table))))
    options(warn=-1)
    di.trim.len.mh.freq.sum.plot.table %>% group_by(variable) %>% ggplot(aes(x = trim.len, y = value, fill = variable)) +
    geom_bar(stat="identity", position = position_dodge(0.8)) +
    scale_fill_manual(values = c("#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
      ,"#009840", "#104277", "#25217D", "#005E6A"
      ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
      ,"#F28322", "#4557BF", "#400101", "#D93223")) +
    coord_cartesian(ylim = c(0, max(di.trim.len.mh.freq.sum.plot.table$value) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Length of Intervening Sequence") +
    ylab("Total Frequency in Reads") -> di.trim.len.mh.freq.sum.plot.p # inside frequency * number of read
    # plot(substitution.number.p)
    saveRDS(di.trim.len.mh.freq.sum.plot.table , file = file.path(save.dir, paste0(id.number, "e]Distribution_of_", type, "_MMEJ_Intervening_Length_Microhomology_Sequence_Di_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "e]Distribution_of_", type, "_MMEJ_Intervening_Length_Microhomology_Sequence_Di_Frequency.png")), plot = di.trim.len.mh.freq.sum.plot.p, dpi = 300, width = 12, height = 6)  
    options(warn=0)
  }else{
    total.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:ncol(freq.df)], levels = colnames(freq.df)[5:ncol(freq.df)])
      ,Frequency = numeric(ncol(freq.df) - 4)
    )
    saveRDS(total.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "a]Distribution_of_", type, "_MMEJ_Total_Microhomology_Sequence_Frequency.table.rds")))
    mono.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:8], levels = colnames(freq.df)[5:8])
      ,Frequency = numeric(8 - 4)
    )
    saveRDS(mono.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "b]Distribution_of_", type, "_MMEJ_Mono_Microhomology_Sequence_Frequency.table.rds")))
    di.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[9:ncol(freq.df)], levels = colnames(freq.df)[9:ncol(freq.df)])
      ,Frequency = numeric(ncol(freq.df) - 8)
    )
    saveRDS(di.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "c]Distribution_of_", type, "_MMEJ_Di_Microhomology_Sequence_Frequency.table.rds")))
    saveRDS(data.frame() , file = file.path(save.dir, paste0(id.number, "d]Distribution_of_", type, "_MMEJ_Intervening_Length_Microhomology_Sequence_Mono_Frequency.table.rds")))
    saveRDS(data.frame() , file = file.path(save.dir, paste0(id.number, "e]Distribution_of_", type, "_MMEJ_Intervening_Length_Microhomology_Sequence_Di_Frequency.table.rds")))
  }
}


# Make microhomology sequence-frequency-count plot
SaveTsFreqPlot <- function(variants.count.total.list, save.dir, type = c("mut", "left_ki", "right_ki")){
  type <- match.arg(type)
  if(type == "mut"){
    id.number <- "[xi"
  }else if(type %in% c("left_ki", "right_ki")){
    id.number <- "[ix"
  }
  id.type.vec <- c("three.prime.conserved.mmej", "five.prime.conserved.mmej", "both.prime.conserved.mmej")
  freq.df <- data.frame(
    mut.label = character(0), read.cnt = numeric(0), mh.len = numeric(0), trimmed.len = numeric(0)
    , A = numeric(0), C = numeric(0), G = numeric(0), T = numeric(0)
    , AA = numeric(0), AC = numeric(0), AG = numeric(0), AT = numeric(0)
    , CA = numeric(0), CC = numeric(0), CG = numeric(0), CT = numeric(0)
    , GA = numeric(0), GC = numeric(0), GG = numeric(0), GT = numeric(0)
    , TA = numeric(0), TC = numeric(0), TG = numeric(0), TT = numeric(0)
  )
  if(length(variants.count.total.list) == 0){
    return()
  }
  for(mut.ind in 1:(length(variants.count.total.list)/4)){
    indel.seq.info.list <- variants.count.total.list[mut.ind * 4][[1]]
    if(indel.seq.info.list$ej.type.label %in% id.type.vec){
      freq.df <- rbind(freq.df
        , data.frame(
          mut.label = variants.count.total.list[mut.ind * 4 - 3][[1]]
          , read.cnt = variants.count.total.list[mut.ind * 4 - 2][[1]]
          , mh.len = nchar(indel.seq.info.list$mh.seq.vec[1])
          , trimmed.len = nchar(indel.seq.info.list$trimmed.seq.vec[1])
          , cbind(indel.seq.info.list$mononuc.trimmed.freq.table, indel.seq.info.list$dinuc.trimmed.freq.table)
        )
      )
    }
  }

  if(nrow(freq.df) > 0){
    # aggregated frequency of microhomology nucleotide - read.cnt
    total.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:ncol(freq.df)], levels = colnames(freq.df)[5:ncol(freq.df)])
      ,Frequency = colSums(t(apply(freq.df, MARGIN = 1, function(row){as.numeric(row[5:length(row)]) * as.numeric(row[2])})))
    )
    options(warn=-1)
    total.sum.by.nucleotide.p <- ggplot(data=total.sum.by.nucleotide.table
      , aes(x=Nucleotide, y=Frequency, fill = c("microhomology"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8
        , fill = c("#00931F", "#001298", "#000000", "#FE0000"
          ,"#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
          ,"#009840", "#104277", "#25217D", "#005E6A"
          ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
          ,"#F28322", "#4557BF", "#400101", "#D93223")
      ) +
    coord_cartesian(ylim = c(0, max(total.sum.by.nucleotide.table$Frequency[1:4]) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Nucleotide in intervening") +
    ylab("Total frequency in reads")
    # plot(substitution.number.p)
    saveRDS(total.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "a]Distribution_of_", type, "_MMEJ_Total_Intervening_Sequence_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "a]Distribution_of_", type, "_MMEJ_Total_Intervening_Sequence_Frequency.png")), plot = total.sum.by.nucleotide.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # aggregated frequency of mono microhomology nucleotide per length - read.cnt
    mono.trimmed.freq.df <- t(apply(freq.df, MARGIN = 1, function(row){
        if(as.numeric(row[3]) == 1){
          as.numeric(row[5:8]) * as.numeric(row[2])
        }else{
          numeric(4)
        }
    }))
    mono.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:8], levels = colnames(freq.df)[5:8])
      ,Frequency = colSums(mono.trimmed.freq.df)
    )
    options(warn=-1)
    mono.mh.sum.by.nucleotide.p <- ggplot(data=mono.mh.sum.by.nucleotide.table
      , aes(x=Nucleotide, y=Frequency, fill = c("intervening"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8
        , fill = c("#00931F", "#001298", "#000000", "#FE0000")
      ) +
    coord_cartesian(ylim = c(0, max(mono.mh.sum.by.nucleotide.table$Frequency[1:4]) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Nucleotide in mono-intervening") +
    ylab("Number of reads")
    # plot(substitution.number.p)
    saveRDS(mono.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "b]Distribution_of_", type, "_MMEJ_Mono_Intervening_Sequence_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "b]Distribution_of_", type, "_MMEJ_Mono_Intervening_Sequence_Frequency.png")), plot = mono.mh.sum.by.nucleotide.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # aggregated frequency of di microhomology nucleotide per length - read.cnt
    di.trimmed.freq.df <- t(apply(freq.df, MARGIN = 1, function(row){
        if(as.numeric(row[3]) == 2){
          as.numeric(row[9:length(row)]) * as.numeric(row[2])
        }else{
          numeric(length(row) - 8)
        }
    }))
    di.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[9:ncol(freq.df)], levels = colnames(freq.df)[9:ncol(freq.df)])
      ,Frequency = colSums(di.trimmed.freq.df)
    )
    options(warn=-1)
    di.mh.sum.by.nucleotide.p <- ggplot(data=di.mh.sum.by.nucleotide.table
      , aes(x=Nucleotide, y=Frequency, fill = c("intervening"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8
        , fill = c("#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
          ,"#009840", "#104277", "#25217D", "#005E6A"
          ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
          ,"#F28322", "#4557BF", "#400101", "#D93223")
      ) +
    coord_cartesian(ylim = c(0, max(di.mh.sum.by.nucleotide.table$Frequency[1:16]) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Nucleotide in di-intervening") +
    ylab("Number of reads")
    # plot(substitution.number.p)
    saveRDS(di.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "c]Distribution_of_", type, "_MMEJ_Di_Intervening_Sequence_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "c]Distribution_of_", type, "_MMEJ_Di_Intervening_Sequence_Frequency.png")), plot = di.mh.sum.by.nucleotide.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # Make trimmed length-microhomology sequence-frequency-count 
    trim.len.mh.freq.list <- list()
    for(trimmed.len in 1:max(freq.df$trimmed.len)){
      trim.len.mh.freq.list <- c(trim.len.mh.freq.list, list(t(apply(freq.df, MARGIN = 1, function(row){
          if(as.numeric(row[4]) == trimmed.len){
            as.numeric(row[5:ncol(freq.df)]) * as.numeric(row[2])
          }else{
            numeric(ncol(freq.df) - 4)
          }
      })))
      )
    }
    trim.len.mh.freq.sum.table <- data.frame(matrix(unlist(lapply(trim.len.mh.freq.list, colSums)), nrow=length(lapply(trim.len.mh.freq.list, colSums)), byrow=T))
    colnames(trim.len.mh.freq.sum.table) <- colnames(freq.df[5:ncol(freq.df)])
    # Mono frequency
    mono.trim.len.mh.freq.sum.table <- trim.len.mh.freq.sum.table[, 1:4]
    mono.trim.len.mh.freq.sum.plot.table <- cbind(melt(mono.trim.len.mh.freq.sum.table), trim.len = factor(rep(1:nrow(mono.trim.len.mh.freq.sum.table), ncol(mono.trim.len.mh.freq.sum.table))))
    options(warn=-1)
    mono.trim.len.mh.freq.sum.plot.table %>% group_by(variable) %>% ggplot(aes(x = trim.len, y = value, fill = variable)) +
    geom_bar(stat="identity", position = position_dodge(0.8)) +
    scale_fill_manual(values = c("#00931F", "#001298", "#000000", "#FE0000")) +
    coord_cartesian(ylim = c(0, max(mono.trim.len.mh.freq.sum.plot.table$value) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Length of Intervening Sequence") +
    ylab("Total Frequency in Reads") -> mono.trim.len.mh.freq.sum.plot.p # inside frequency * number of read
    # plot(substitution.number.p)
    saveRDS(mono.trim.len.mh.freq.sum.plot.table , file = file.path(save.dir, paste0(id.number, "d]Distribution_of_", type, "_MMEJ_Intervening_Length_Intervening_Sequence_Mono_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "d]Distribution_of_", type, "_MMEJ_Intervening_Length_Intervening_Sequence_Mono_Frequency.png")), plot = mono.trim.len.mh.freq.sum.plot.p, dpi = 300, width = 12, height = 6)
    options(warn=0)

    # Di frequency
    di.trim.len.mh.freq.sum.table <- trim.len.mh.freq.sum.table[, 5:20]
    di.trim.len.mh.freq.sum.plot.table <- cbind(melt(di.trim.len.mh.freq.sum.table), trim.len = factor(rep(1:nrow(di.trim.len.mh.freq.sum.table), ncol(di.trim.len.mh.freq.sum.table))))
    options(warn=-1)
    di.trim.len.mh.freq.sum.plot.table %>% group_by(variable) %>% ggplot(aes(x = trim.len, y = value, fill = variable)) +
    geom_bar(stat="identity", position = position_dodge(0.8)) +
    scale_fill_manual(values = c("#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
      ,"#009840", "#104277", "#25217D", "#005E6A"
      ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
      ,"#F28322", "#4557BF", "#400101", "#D93223")) +
    coord_cartesian(ylim = c(0, max(di.trim.len.mh.freq.sum.plot.table$value) + 100)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
      , axis.text.y = element_text(size = 6)
      , axis.title.x = element_text(size = 12)
      , axis.title.y = element_text(size = 12)) +
    xlab("Length of Intervening Sequence") +
    ylab("Total Frequency in Reads") -> di.trim.len.mh.freq.sum.plot.p # inside frequency * number of read
    # plot(substitution.number.p)
    saveRDS(di.trim.len.mh.freq.sum.plot.table , file = file.path(save.dir, paste0(id.number, "e]Distribution_of_", type, "_MMEJ_Intervening_Length_Intervening_Sequence_Di_Frequency.table.rds")))
    ggsave(file = file.path(save.dir, paste0(id.number, "e]Distribution_of_", type, "_MMEJ_Intervening_Length_Intervening_Sequence_Di_Frequency.png")), plot = di.trim.len.mh.freq.sum.plot.p, dpi = 300, width = 12, height = 6)  
    options(warn=0)
  }else{
    total.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:ncol(freq.df)], levels = colnames(freq.df)[5:ncol(freq.df)])
      ,Frequency = numeric(ncol(freq.df) - 4)
    )
    saveRDS(total.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "a]Distribution_of_", type, "_MMEJ_Total_Intervening_Sequence_Frequency.table.rds")))
    mono.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[5:8], levels = colnames(freq.df)[5:8])
      ,Frequency = numeric(8 - 4)
    )
    saveRDS(mono.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "b]Distribution_of_", type, "_MMEJ_Mono_Intervening_Sequence_Frequency.table.rds")))
    di.mh.sum.by.nucleotide.table <- data.frame(
      Nucleotide = factor(colnames(freq.df)[9:ncol(freq.df)], levels = colnames(freq.df)[9:ncol(freq.df)])
      ,Frequency = numeric(ncol(freq.df) - 8)
    )
    saveRDS(di.mh.sum.by.nucleotide.table , file = file.path(save.dir, paste0(id.number, "c]Distribution_of_", type, "_MMEJ_Di_Intervening_Sequence_Frequency.table.rds")))
    saveRDS(data.frame() , file = file.path(save.dir, paste0(id.number, "d]Distribution_of_", type, "_MMEJ_Intervening_Length_Intervening_Sequence_Mono_Frequency.table.rds")))
    saveRDS(data.frame() , file = file.path(save.dir, paste0(id.number, "e]Distribution_of_", type, "_MMEJ_Intervening_Length_Intervening_Sequence_Di_Frequency.table.rds")))
  }
}


CalcKldivergence <- function(x.crispr.set, y.crispr.set){
  if(is.null(x.crispr.set) & is.null(y.crispr.set)){
    return(0)
  }else if(is.null(x.crispr.set) | is.null(y.crispr.set)){
    return(1)
  }
  # prepare vector
  variants.names.vec <- unique(c(rownames(x.crispr.set$cigar_freqs), rownames(y.crispr.set$cigar_freqs)))
  # 0.5 is small pseudocount for avoinding division by zero. This method was inspired by [Felicity Allen, et.al, 2017]
  x.variants.vec <- rep(0.5, length(variants.names.vec))
  names(x.variants.vec) <- variants.names.vec
  x.variants.vec[rownames(x.crispr.set$cigar_freqs)] <- x.variants.vec[rownames(x.crispr.set$cigar_freqs)] + x.crispr.set$cigar_freqs[,1]
  x.variants.vec <- x.variants.vec / sum(x.variants.vec)
  y.variants.vec <- rep(0.5, length(variants.names.vec))
  names(y.variants.vec) <- variants.names.vec
  y.variants.vec[rownames(y.crispr.set$cigar_freqs)] <- y.variants.vec[rownames(y.crispr.set$cigar_freqs)] + y.crispr.set$cigar_freqs[,1]
  y.variants.vec <- y.variants.vec / sum(y.variants.vec)
  return(sum(x.variants.vec * log(x.variants.vec / y.variants.vec)))
}

#################################################################################################################################################################

