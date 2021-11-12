
library(R6)

# Calculate physiochemical nucleotide propaties of reports written by Wei Chen et al.
#
# Note:
#   This program was refered to following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   CHEN, Wei, et al. iNuc-PhysChem: a sequence-based predictor for identifying nucleosomes via physicochemical properties. PloS one, 2012, 7.10: e47843.
#   I am grateful to them!
CalcPhyChem <- R6Class(
  # name
  "CalcPhyChem",
  # public menber
  public = list(

    seq  = NA,  # DNAStringSet objects of sequences
    seq.vec  = character(0),  # vector of sequences
    ref.dir = "", # ref.dir: directory path contains reference tables. 
    #Table of 12 dinucleotide propaties was refered to CHEN, Wei, et al. iNuc-PhysChem: a sequence-based predictor for identifying nucleosomes via physicochemical properties. PloS one, 2012, 7.10: e47843.
    prop.table = NULL,
    seqave.table  = numeric(13),  # dataframe of average of the 12 physiochemical nucleotide propaties in sequences
    seqmax.table  = numeric(13),  # dataframe of maximum of the 12 physiochemical nucleotide propaties in sequences
    seqmin.table  = numeric(13),  # dataframe of minimum of the 12 physiochemical nucleotide propaties in sequences
    summary  = NULL,  # summary of all propaties in sequences

    # constructor
    initialize = function(seq.vec, ref.dir) {
      # Set sequences as the validation set
      #
      # Args:
      #   seq.vec: vector of sequences is to be calculated.
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
        self$prop.table <- read.csv(file.path(self$ref.dir, "PhyChemTable.csv"), stringsAsFactors = FALSE)
        
        message("Set sequences")
        
        self$seq.vec <- as.character(self$seq)
        
        # Calculate 12 physiochemical nucleotide propaties
        
        message("Calculate 12 physiochemical nucleotide propaties")
        
        names(self$seqave.table) <- c("sequence", "mean.P(1)", "mean.P(2)",
                                      "mean.P(3)", "mean.P(4)", "mean.P(5)", "mean.P(6)",
                                      "mean.P(7)", "mean.P(8)", "mean.P(9)", "mean.P(10)",
                                      "mean.P(11)", "mean.P(12)"
        )
        names(self$seqmax.table) <- c("sequence", "max.P(1)", "max.P(2)",
                                      "max.P(3)", "max.P(4)", "max.P(5)", "max.P(6)",
                                      "max.P(7)", "max.P(8)", "max.P(9)", "max.P(10)",
                                      "max.P(11)", "max.P(12)"
        )
        names(self$seqmin.table) <- c("sequence", "min.P(1)", "min.P(2)",
                                      "min.P(3)", "min.P(4)", "min.P(5)", "min.P(6)",
                                      "min.P(7)", "min.P(8)", "min.P(9)", "min.P(10)",
                                      "min.P(11)", "min.P(12)"
        )#set name
        for (seq in seq.vec) {
          res = self$CalcProp(seq)
          self$seqave.table <- rbind(self$seqave.table, as.vector(res$Average))
          self$seqmax.table <- rbind(self$seqmax.table, as.vector(res$Maximum))
          self$seqmin.table <- rbind(self$seqmin.table, as.vector(res$Minimum))
        }
        # remove numeric(13) in the first line
        self$seqave.table <- self$seqave.table[-1, ]
        self$seqmax.table <- self$seqmax.table[-1, ]
        self$seqmin.table <- self$seqmin.table[-1, ]
        
        message("Record summary")
        self$summary <- private$merges(
          list(
            as.data.frame(self$seqave.table)
            , as.data.frame(self$seqmax.table)
            , as.data.frame(self$seqmin.table)
          ), by="sequence", sort = FALSE
        )
        
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
    ToNnuc = function(seq, N = 1) {
      # Convert sequence to vector of N-nucleotides
      #
      # Args:
      #   seq: sequence is to be converted.
      #   N: length of splited nucleotides.
      #
      # Returns:
      #   On success, returns vector of N-nucleotides.
      #   On failure, returns character vector : character( <length of "seq"> ).

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")
        private$checkMode(N, "N", "numeric")


        # Convert sequence to vector of dinucleotides.
        nnuc.vec <- character(0)
        nnuclen <- nchar(seq) - N + 1
        for (nnuc_ind in 1:nnuclen){
          nnuc.vec = c(nnuc.vec, substring(seq, nnuc_ind, nnuc_ind + N -1))
        }
        return(nnuc.vec)
      }
      , error = function(e) {
        message(e)
        return(character(length(seq)))
      }
      , silent = TRUE
      )
    },
    CalcProp = function(seq) {
      # Calculate average of 12 propaties of dinucleotide in the given sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns data.frame of result.
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        dinuc.vec <- self$ToNnuc(seq, N = 2)
        resformat = numeric(13)
        names(resformat) <- c("Seq", "P(1)", "P(2)",
                              "P(3)", "P(4)", "P(5)", "P(6)",
                              "P(7)", "P(8)", "P(9)", "P(10)",
                              "P(11)", "P(12)"
        ) #set name
        averes <- maxres <- minres <- resformat


        averes["Seq"] <- maxres["Seq"] <- minres["Seq"] <- seq

        #collect propaties for parts of sequence
        parts.table <- numeric(13)
        for (dinuc in dinuc.vec){
          parts.vec = self$prop.table[self$prop.table$Dinucleotide == dinuc, ]
          parts.vec[1] <- as.character(parts.vec[1])
          parts.vec[2:13] <- as.numeric(parts.vec[2:13])
          parts.table <- rbind(parts.table, parts.vec)
        }
        parts.table <- parts.table[-1, ] # remove numeric(13) in the first line

        #calcutation
        averes[2:13] = apply(parts.table[, 2:13], 2, mean)
        maxres[2:13] = apply(parts.table[, 2:13], 2, max)
        minres[2:13] = apply(parts.table[, 2:13], 2, min)

        return(data.frame(Average = averes, Maximum = maxres, Minimum = minres))
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
