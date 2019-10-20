
library(R6)
library(Biostrings)
library(tcltk2)

# Calculate Alignment score between two sequences.
#
# Note:
#   This program was refered to following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   I am grateful to them!
CalcAlign <- R6Class(
  # name
  "CalcAlign",
  # public menber
  public = list(
    
    seq1.DNAStringSet = NA,  # DNAStringSet objects of sequences type 1
    seq2.DNAStringSet = NA,  # DNAStringSet objects of sequences type 2
    len1.num = 0,  # length of alignment sequences type 1
    len2.num = 0,  # length of alignment sequences type 2
    ori1.char = "",  # orientation of alignment sequences type 1
    ori2.char = "",  # orientation of alignment sequences type 2
    seq1.vec = character(0),  # vector of sequences type 1
    seq2.vec = character(0),  # vector of sequences type 2
    is.sense.sense = TRUE,  # alignment between sense sequences. 
                            # If is.sense.sense is FALSE, alignment between sense-antisense sequences is performed.
    align.score.table = numeric(9),
    summary  = NULL,  # summary of all propaties in sequences

    # constructor
    initialize = function(seq1.vec,
                          len1.num = 5,
                          ori1.char = "5' to 3'",
                          seq2.vec,
                          len2.num = 5,
                          ori2.char = "5' to 3'",
                          is.sense.sense = TRUE) {
        # Set sequences as the validation set
        #
        # Args:
        #   seq1.vec: vector of sequences type 1 is to be aligned.
        #   len1.num: length of sequences type 1 is to be aligned.
        #   ori1.char: vector of sequences type 1 is to be aligned.
        #   seq2.vec: vector of sequences type 2 is to be aligned.
        #   len2.num: length of sequences type 2 is to be aligned.
        #   ori2.char: vector of sequences type 2 is to be aligned.
        #   is.sense.sense: alignment between sense sequences. 
        #                   If is.sense.sense is FALSE, alignment between sense-antisense sequences is performed.
        #
        # Returns:
        #   On success, returns TRUE.
        #   On failure, returns FALSE.
      
      tryCatch({
        # Check sequence data
        private$checkSeqVec(seq1.vec)
        private$checkMode(len1.num, "len1.num", "numeric")
        private$checkMode(ori1.char, "ori1.char", "character")
        private$checkSeqVec(seq2.vec)
        private$checkMode(len2.num, "len2.num", "numeric")
        private$checkMode(ori2.char, "ori2.char", "character")
        private$checkMode(is.sense.sense, "is.sense.sense", "logical")
      
        # Set infomation
        self$seq1.DNAStringSet <- DNAStringSet(seq1.vec)
        self$len1.num <- len1.num
        self$ori1.char <- ori1.char
        self$seq2.DNAStringSet <- DNAStringSet(seq2.vec)
        self$len2.num <- len2.num
        self$ori2.char <- ori2.char
        self$is.sense.sense <- is.sense.sense
        
        
        message("Set sequences")
        
        self$seq1.vec <- as.character(self$seq1.DNAStringSet)
        self$seq2.vec <- as.character(self$seq2.DNAStringSet)
        
        # Calculate the optimal global pairwise alignment scores
        # between the seed region and scaffold using the Needleman-Wunsch algorithm
        
        message("Calculate the optimal global pairwise alignment scores")
        
        names(self$align.score.table) <- c("sequence1.char", "aligned sub-sequence1.char",
                                           "sequence2.char", "aligned sub-sequence2.char",
                                           "align.score.global.num","align.score.local.num",
                                           "align.score.overlap.num", "align.score.globallocal.num",
                                           "align.score.localglobal.num") #set name
        pb <- txtProgressBar(min = 1, max = length(self$seq1.vec), style = 3)
        for (i in 1:length(self$seq1.vec)) {
          res <- self$AlignSeqs(seq1.char = self$seq1.vec[i], len1.num = self$len1.num, ori1.char = self$ori1.char,
                                seq2.char = self$seq2.vec[i], len2.num = self$len2.num, ori2.char = self$ori2.char,
                                self$is.sense.sense)
          self$align.score.table <- rbind(self$align.score.table, as.vector(res))
          setTxtProgressBar(pb, i)
        }
        # remove numeric(9) in the first line
        self$align.score.table <- self$align.score.table[-1, ]
        
        message("Record summary")
        self$summary <- as.data.frame(self$align.score.table)
        
        message("Complete !")

      }
      , error = function(e) {
        message(e)
        return(FALSE)
      }
      , silent = TRUE
      )
      return(TRUE)
    },
    # public methods
    ExtractSubSeq = function(seqinfo.list){
      # Extract sub-sequence from seqinfo.list.
      #
      # Args:
      #   seqinfo.list: sub-sequence information
      #     $seq.char:  host-sequence
      #     $ori.char:   derection ;"5' to 3'" | "3' to 5'"
      #     $len.num    length of sub-sequence
      #
      # Returns:
      #   On success, returns characters of sub-sequence.
      #   On failure, returns ERROR MESSAGE.
      #
      # Error handling
      if(!is.list(seqinfo.list) || !is.character(seqinfo.list$ori.char)){
        if(!is.list(seqinfo.list)) stop("seqinfo.list is expected to be list.") 
        else stop("seqinfo.list$ori.char is expected to be characters.")
      }
      if (seqinfo.list$ori.char == "5' to 3'"){
        return(as.character(DNAString(x = seqinfo.list$seq.char,
                                      start = 1,
                                      nchar = seqinfo.list$len.num)
        ))
      }
      else if(seqinfo.list$ori.char == "3' to 5'"){
        return(as.character(DNAString(x = seqinfo.list$seq.char,
                                      start = -seqinfo.list$len.num,
                                      nchar = seqinfo.list$len.num)
        ))
      }
      else stop(paste(seqinfo.list$ori.char, "is not undefined characters. Please use <5' to 3'> or <3' to 5'>", sep = ' '))
    },
    AlignSeqs = function(seq1.char, len1.num, ori1.char, seq2.char, len2.num, ori2.char, is.sense.sense = TRUE){
      # Align between sequences according to given sequence infomation
      #
      # Args:
      #   seq1.char: vector of sequences type 1 is to be aligned.
      #   len1.num: length of sequences type 1 is to be aligned.
      #   ori1.char: vector of sequences type 1 is to be aligned.
      #   seq2.char: vector of sequences type 2 is to be aligned.
      #   len2.num: length of sequences type 2 is to be aligned.
      #   ori2.char: vector of sequences type 2 is to be aligned.
      #   is.sense.sense; alignment between sense sequences. 
      #                   If is.sense.sense is FALSE, alignment between sense-antisense sequences is performed.
      #
      # Returns:
      #   On success, returns data.frame of result.
      #   On failure, returns NULL.
      #
      tryCatch({
        # Error handling
        private$checkMode(seq1.char, "seq1.char", "character")
        private$checkMode(len1.num, "len1.num", "numeric")
        private$checkMode(ori1.char, "ori1.char", "character")
        private$checkMode(seq2.char, "seq1.char", "character")
        private$checkMode(len2.num, "len1.num", "numeric")
        private$checkMode(ori2.char, "ori1.char", "character")
        
        # Create sub-seqeuence infomation to be aligned
        seq1.info.list <- list(
          seq.char = seq1.char,
          ori.char = ori1.char,
          len.num = len1.num
        )
        seq2.info.list <- list(
          seq.char = seq2.char,
          ori.char = ori2.char,
          len.num = len2.num
        )
        
        # Alignment
        subseq1.char <- self$ExtractSubSeq(seq1.info.list)
        subseq2.char <- self$ExtractSubSeq(seq2.info.list)
        if(is.sense.sense){
          subseq1.DNAString <- DNAString(subseq1.char)
          subseq2.DNAString <- DNAString(subseq2.char)
        }else{
          subseq1.DNAString <- reverseComplement(DNAString(subseq1.char))
          subseq2.DNAString <- DNAString(subseq2.char)
          subseq1.char <- as.character(subseq1.DNAString)
        }
        
        
        return(data.frame(sequence1 = seq1.char, subseq1.char, sequence2 = seq2.char, subseq2.char,
                          align.score.global.num = pairwiseAlignment(subseq1.DNAString, subseq2.DNAString, type="global", scoreOnly = TRUE),
                          align.score.local.num = pairwiseAlignment(subseq1.DNAString, subseq2.DNAString, type="local", scoreOnly = TRUE),
                          align.score.overlap.num = pairwiseAlignment(subseq1.DNAString, subseq2.DNAString, type="overlap", scoreOnly = TRUE),
                          align.score.globallocal.num = pairwiseAlignment(subseq1.DNAString, subseq2.DNAString, type="global-local", scoreOnly = TRUE),
                          align.score.localglobal.num = pairwiseAlignment(subseq1.DNAString, subseq2.DNAString, type="local-global", scoreOnly = TRUE),
                          stringsAsFactors = FALSE))
      }
      , error = function(e) {
        message(e)
        return(NULL)
      }
      , silent = TRUE
      )
    }
  ),
  
  private = list(
    checkSeqVec = function(seq.vec) {
      # Check whether seq.vec is vector and this elements is character.
      #
      # Args:
      #   seq.vec: value is to be checked.
      #
      # Returns:
      #   On success, returns 0.
      #   On failure, returns ERROR MESSAGE.
      #
      # Error handling
      if (!is.vector(seq.vec)) {
        #do type identifying for a user input
        if(is.list(seq.vec)){
          Robj = "list"
        } else if (is.matrix(seq.vec)) {
          Robj = "matrix"
        } else if (is.array(seq.vec)) {
          Robj = "array"
        } else if (is.factor(seq.vec)) {
          Robj = "factor"
        } else if (is.data.frame(seq.vec)) {
          Robj = "data.frame"
        } else {
          Robj = "unknown type"
        }
        stop("Argument seq.vec must be a vector: your input was ", Robj, ".")
      }
      for (seq in seq.vec) { #do type checking for elements of user input
        if (!is.character(seq)) {
          stop("Element of argument seq.vec must be a string: element of your input was ", seq, ".")
        }
      }
      return(0)
    },
    checkMode = function(input, argname, mode) {
      # Check type of a user input
      #
      # Args:
      #   input: value is to be checked.
      #   argname: name of argument is to be checked.
      #   mode: type is expected.
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.

      input.mode <- mode(input)
      if (input.mode != mode){
        stop("Argument ", argname, "must be a ", mode, " : your input was ", input.mode, ".\n")
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
  )
)
