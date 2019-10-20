
library(R6)
library(stringr)

# Calculate microhomology score and out-of-frame score written by Wei Chen et al.
#
# Note:
#   This program was refered to following reports :
#   BAE, Sangsu, et al. Microhomology-based choice of Cas9 nuclease target sites. Nature methods, 2014, 11.7: 705-706.
#   I am grateful for them!
CalcMHScore <- R6Class(
  # name
  "CalcMHScore",
  # public menber
  public = list(

    seq  = NA,  # DNAStringSet objects of sequences
    seq.vec  = character(0),  # vector of sequences
    seqmh.vec  = numeric(0),  # vector of microhomology score in sequences
    seqframeshift.vec  = numeric(0),  # vector of out-of-frame score in sequences
    summary  = NULL,  # summary of all propaties in sequences

    length_weight = 0.0, #weight for calculating score
    cutpos.vec = numeric(0), #insert the position (5' -> 3') expected to be broken.
    script.dir = "", # directory which has python script

    # constructor
    initialize = function(seq.vec, cutpos.vec, length_weight = 20.0, script.dir = ".") {
      # Set sequences as the validation set
      #
      # Args:
      #   seq.vec: vector of sequences is to be calculated.
      #   length_weight   :   weight for calculating score
      #   cutpos.vec          :   insert the position (5' -> 3') expected to be broken.
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.

      tryCatch({
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
        # Set infomation
        self$seq <- DNAStringSet(seq.vec) # DNAStringSet objects
        
        message("Set sequences")
        
        self$seq.vec <- as.character(self$seq)
        
        # Calculate propaties
        
        message("Calculate propaties")
        
        self$length_weight <- length_weight
        self$cutpos.vec <- cutpos.vec

        # Calculate microhomology score and out-of-frame score
        for (ind in 1:length(self$seq.vec)) {
          res = self$ProfileMH(ind)
          self$seqmh.vec <- c(self$seqmh.vec, res["Microhomology score"])
          self$seqframeshift.vec <- c(self$seqframeshift.vec, res["Out-of-frame score"])
        }
        
        message("Record summary")
        
        self$summary <- data.frame(
          self$seq.vec,
          self$seqmh.vec,
          self$seqframeshift.vec
        )
        colnames(self$summary) <- c("sequence" ,"microhomology.score", "out-of-frame.score")
        
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
    ProfileMH = function(ind) {
      # Calculate microhomology score and out-of-frame score in the given sequence
      #
      # Args:
      #   ind: index of sequence set
      #
      # Returns:
      #   On success, returns data.frame of result.
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(self$seq.vec[ind], "seq", "character")

        res <- system(paste("python", file.path(script.dir, "calcMHscore.py"), self$seq.vec[ind], self$length_weight, self$cutpos.vec), intern = TRUE)

        mhscore <- as.numeric(str_extract(res[8], pattern = "[0-9]+.[0-9]+"))
        ofscore <- as.numeric(str_extract(res[9], pattern = "[0-9]+.[0-9]+"))

        scores <- c(mhscore, ofscore)
        names(scores) <- c("Microhomology score", "Out-of-frame score")

        return(scores)
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
    }
  )
)
