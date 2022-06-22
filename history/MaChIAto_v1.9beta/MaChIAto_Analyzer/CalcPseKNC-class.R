
library(R6)
library(rDNAse)

# Calculate Pseudo k-tuple nucleotide composition (PseKNC) of reports written by Guo, Shou-Hui, et al.
#
# Note:
#   This program was refered to following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   Guo, Shou-Hui, et al. "iNuc-PseKNC: a sequence-based predictor for predicting nucleosome positioning in genomes with pseudo k-tuple nucleotide composition." Bioinformatics 30.11 (2014): 1522-1529.
#   I am grateful to them!
CalcPseKNC <- R6Class(
  # name
  "CalcPseKNC",
  # public menber
  public = list(

    seq  = NA,  # DNAStringSet objects of sequences
    seq.vec  = character(0),  # vector of sequences
    ref.dir = "", # ref.dir: directory path contains reference tables. 
    #Table of 6 phisical Structural propaties was refered to Guo, Shou-Hui, et al. "iNuc-PseKNC: a sequence-based predictor for predicting nucleosome positioning in genomes with pseudo k-tuple nucleotide composition." Bioinformatics 30.11 (2014): 1522-1529.
    prop.table = NULL,
    seqPseKNC.table = numeric(17),  # vector of sequences
    summary  = NULL,  # summary of all propaties in sequences

    # constructor
    initialize = function(seq.vec, lambda = 1, w = 0.05, k = 3, ref.dir) {
      # Set sequences as the validation set
      #
      # Args:
      #   seq.vec: vector of sequences is to be calculated.
      #   lambda: It is one of the parameters for the Pseudo K-tupler Composition.
      #           an integer larger than or equal to 0 and less than or equal to L-2
      #           (L means the length of the shortest sequence in the dataset).
      #           It represents the highest counted rank (or tier) of the correlation along a DNA sequence.
      #           Its default value is 3.
      #   w: It is one of the parameters for the Pseudo K-tupler Composition.
      #      the weight factor ranged from 0 to 1. Its default value is 0.05.
      #   k: It is one of the parameters for the Pseudo K-tupler Composition.
      #      an integer larger than 0 represents the k-tuple. Its default value is 3.
      #   ref.dir: directory path contains reference tables.
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
        self$ref.dir <- ref.dir

        # Load table
        self$prop.table <- read.csv(file.path(self$ref.dir, "Six.Physical.Structural.Prop.csv"), stringsAsFactors = FALSE)
        
        message("Set sequences")
        
        self$seq.vec <- as.character(self$seq)
        
        # Calculate Pseudo k-tuple nucleotide composition using DNA local structural properties
        # which were divided into local translational (rise, slide and shift) and angular (twist, roll and tilt)
        # to generate the PseKNC feature vector
        
        message("Calculate Pseudo k-tuple nucleotide composition")

        names(self$seqPseKNC.table) <- c("sequence", "1_AA", "2_AC",
                                      "3_AG", "4_AT", "5_CA", "6_CC",
                                      "7_CG", "8_CT", "9_GA", "10_GC",
                                      "11_GG", "12_GT", "13_TA", "14_TC",
                                      "15_TG", "16_TT"
        )#set name

        for (seq in seq.vec) {
          res = self$CalcProp(seq, lambda, w, k)
          self$seqPseKNC.table <- rbind(self$seqPseKNC.table, as.vector(res))
        }
        
        # remove numeric(17) in the first line
        self$seqPseKNC.table <- self$seqPseKNC.table[-1, ]
        
        message("Record summary")
        self$summary <- as.data.frame(self$seqPseKNC.table)
        
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
    CalcProp = function(seq, lambda, w, k) {
      # Calculate 16 propaties of dinucleotide in the given sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #   lambda: It is one of the parameters for the Pseudo K-tupler Composition.
      #           an integer larger than or equal to 0 and less than or equal to L-2
      #           (L means the length of the shortest sequence in the dataset).
      #           It represents the highest counted rank (or tier) of the correlation along a DNA sequence.
      #           Its default value is 3.
      #   w: It is one of the parameters for the Pseudo K-tupler Composition.
      #      the weight factor ranged from 0 to 1. Its default value is 0.05.
      #   k: It is one of the parameters for the Pseudo K-tupler Composition.
      #      an integer larger than 0 represents the k-tuple. Its default value is 3.
      #
      # Returns:
      #   On success, returns data.frame of result.
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")
        private$checkMode(lambda, "lambda", "numeric")
        private$checkMode(w, "w", "numeric")
        private$checkMode(k, "k", "numeric")

        #collect propaties for parts of sequence
        customprops <- t(self$prop.table[ ,-1])
        colnames(customprops) <- self$prop.table[ ,"Dinucleotide"]
        PseKNC.table <- extrPseKNC(seq, normalize = TRUE, customprops = customprops, lambda = lambda, w = w, k = k)[-17]

        return(data.frame(sequence = seq, t(PseKNC.table), stringsAsFactors = FALSE))
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
    merges = function(dfs, ...) {
      # merge dataf.rames
      #
      # Args:
      #   dfs: merged data.frames
      #   ...: options of merge()
      #
      # Returns:
      #   On success, returns merged data.frame.
      #   On failure, returns null data.frame.
      
      tryCatch({
          # Error handling
          if (!is.list(dfs)) {
            #do type identifying for a user input
            if(is.data.frame(dfs)){
              Robj = "data.frame"
            } else if (is.matrix(dfs)) {
              Robj = "matrix"
            } else if (is.array(dfs)) {
              Robj = "array"
            } else if (is.factor(dfs)) {
              Robj = "factor"
            } else if (is.vector(dfs)) {
              Robj = "vector"
            } else {
              Robj = "unknown type"
            }
            stop("Argument dfs must be a list: your input was ", Robj, ".")
          }
          
          base <- dfs[1]
          lapply(dfs[-1], function(i) base <<- merge(base, i, ...))
          return(base)
          
        }
        , error = function(e) {
          message(e)
          return(data.frame())
        }
        , silent = TRUE
      )
    }
  )
)
