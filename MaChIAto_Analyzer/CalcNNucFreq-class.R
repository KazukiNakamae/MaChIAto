
library(R6)
library(ape)

# Calculate N-nucleotide composition.
#
# Note:
#   This program was refered to following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   I am grateful for them!
CalcNNucFreq <- R6Class(
  # name
  "CalcNNucFreq",
  # public menber
  public = list(

    seq  = NA,  # DNAStringSet objects of sequences
    seqfreq.list  = NA,  # list of frequency and position-dependent frequency
    summary  = NULL,  # summary of all propaties in sequences

    exception.vec = numeric(0), # vector of exception position for calculation.

    # constructor
    initialize = function(seq.vec, N, exception.vec = numeric(0)) {
      # Set sequences as the validation set
      #
      # Args:
      #   seq.vec: vector of sequences is to be calculated.
      #   N               :   length of splited nucleotides.
      #   exception.vec   :   insert the position (5' -> 3') expected to be broken.
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.

      tryCatch({
        # Error handling
        # seq.vec
        if (private$getObjType(seq.vec) != "vector") {
          stop("Argument seq.vec must be a vector.")
        }
        for (seq in seq.vec) { #do type checking for elements of user input
          if (!is.character(seq)) {
            stop("Element of argument seq.vec must be a string: element of your input was ", seq, ".")
          }
        }
        # N
        private$checkMode(N, "N", "numeric")
        if(!(N %% 1 == 0) | N < 1){
          stop("Argument N must be natural number.", N, ".")
        }
        # exception.vec
        if (private$getObjType(exception.vec) != "vector") {
          stop("Argument exception must be a vector.")
        }
        for (exception in exception.vec) { #do type checking for elements of user input
          if (!is.numeric(exception)) {
            stop("Element of argument exception.vec must be a numeric: element of your input was ", exception, ".")
          }
        }
        
        
        # Set infomation
        message("Set sequences")
        self$seq <- DNAStringSet(seq.vec) # DNAStringSet objects
        self$exception.vec <- exception.vec

        # Calculate propaties
        message("Calculate nucleotide frequency")
        self$seqfreq.list <- self$CalcSeqFreq(self$seq, exception.vec = self$exception.vec, N = N)
        
        message("Record summary")
        self$summary <- c(list(whole.sequence = as.vector(self$seq)), self$seqfreq.list)
        
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
    CalcSeqFreq = function(seqset, exception.vec = numeric(0), N = 1) {
      # Calculate nucleotide frequency in sequences
      #
      # Args:
      #   seqset: input DNA String set. NOTE : length of sequences must be same.
      #   exception.vec : vector of exception position for calculation.
      #   N: length of splited nucleotides.
      #
      # Returns:
      #   On success, returns list of list of posiosion-dependent and non-positive-dependent nucleotide frequency.
      #   On failure, returns numeric(0).
      
      tryCatch({
        # Error handling
        # seq.table.list
        private$checkMode(seqset, "seqset", "S4")
        private$checkClass(seqset, "seqset", "DNAStringSet")
        # exception.vec
        if (private$getObjType(exception.vec) != "vector") {
          stop("Argument seq.table.list must be a vector.")
        }
        for (exception in exception.vec) {
          private$checkMode(exception, "Element of exception.vec", "numeric")
        }
        # N
        private$checkMode(N, "N", "numeric")
        if(!(N %% 1 == 0) | N < 1){
          stop("Argument N must be natural number. : ", N, ".")
        }

        # calculate exception start position
        exstartpos.vec <- numeric(0)
        for(exception in exception.vec){
          for(exstartpos in (exception - (N - 1)):exception){
            # memolize
            exstartpos.vec = c(exstartpos.vec, exstartpos)
          }
        }
        
        # Calculate N-nucleotide frequency
        startpos <- 1
        endpos <- startpos + N - 1
        pdfreq <- numeric(0)
        freq <- oligonucleotideFrequency(seqset, width = N) * 0
        repeat{
          
          # check parameters
          if(!(startpos %in% exstartpos.vec)){
            # count frequency
            atfreq <- oligonucleotideFrequency(subseq(seqset, start = startpos, end = endpos), width = N)
            # Calculate N-nucleotide frequency
            freq <- freq + atfreq
            # Calculate position-dependent N-nucleotide frequency
            cname <- colnames(atfreq)
            colnames(atfreq) <- paste(startpos, cname, sep = ".")
            pdfreq <- cbind(pdfreq, atfreq)
          }
          
          # update parameters
          startpos <- startpos + 1
          endpos <- startpos + N - 1
          if(endpos > width(seqset)[1]) break
          
        }
        
        res.list = list(freq, pdfreq)
        names(res.list) = c("frequency", "PDfrequency")
        return(res.list)
      }
      , error = function(e) {
        message(e)
        return(numeric(0))
      }
      , silent = TRUE
      )
    }
  ),

  private = list(
    #Table of 12 dinucleotide propaties was refered to CHEN, Wei, et al. iNuc-PhysChem: a sequence-based predictor for identifying nucleosomes via physicochemical properties. PloS one, 2012, 7.10: e47843.
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
    },
    checkClass = function(input, argname, mode) {
      # Check Class of a user input
      #
      # Args:
      #   input: value is to be checked.
      #   argname: name of argument is to be checked.
      #   mode: type is expected.
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.
      
      input.mode <- class(input)
      if (input.mode != mode){
        stop("Argument ", argname, "must be a instance of", mode, " class : your input was ", input.mode, ".\n")
        return(FALSE)
      } else {
        return(TRUE)
      }
    },
    getObjType = function(input) {
      # Get object type of a user input
      #
      # Args:
      #   input: object is to be checked.
      #
      # Returns:
      #   On success, returns name of object type.
      #   On failure, returns null character.
      
      tryCatch({
          if(is.list(input)){
            return("list")
          } else if (is.vector(input)) {
            return("vector")
          } else if (is.matrix(input)) {
            return("matrix")
          } else if (is.array(input)) {
            return("array")
          } else if (is.factor(input)) {
            return("factor")
          } else if (is.data.frame(input)) {
            return("data.frame")
          } else {
            stop("Argument input is not standard R Object", input, ".")
          }
        }, error = function(e) {
          message(e)
          return("")
      }, silent = TRUE
    )}
  )
)
