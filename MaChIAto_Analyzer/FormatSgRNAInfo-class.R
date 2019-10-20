
library(R6)
library(Biostrings)

# Format sgRNA sequence information to extract it easily.
#
FormatSgRNAInfo <- R6Class(
  # name
  "FormatSgRNAInfo",
  # public menber
  public = list(
    
    fourtybase.region.surrounding.targetseq.vec = character(0), # 40bp region surrounding the target site.
    protospacer.start.pos.num = 0, # The position where protospacer involved in 40bp region surrounding the target site starts.
    protospacer.end.pos.num = 0, # The position where protospacer involved in 40bp region surrounding the target site ends.
    scaffold.seq.vec = character(0), # The each sgRNA scaffold sequence.
    need.annotation.logical = FALSE, # a logical, if TRUE, the the value of annotation.vec will be validated.
    annotation.vec = character(0), # annotation names is to be given to each 40bp region surrounding the target site.
    
    protospacer.seq.vec = character(0), # protospacer
    fourtybase.region.surrounding.targetseq.df = NULL, # table which contains 0bp region surrounding the target sequences and annotations
    protospacer.seq.df = NULL, # table which contains protospacer sequences and annotations
    scaffold.seq.df = NULL, # table which contains scaffold sequences and annotations
    
    # constructor
    initialize = function(fourtybase.region.surrounding.targetseq.vec = character(0),
                          protospacer.start.pos.num = 12,
                          protospacer.end.pos.num = 31,
                          scaffold.seq.vec = character(0),
                          need.annotation.logical = FALSE,
                          annotation.vec = character(0)
                          ) {
      # Set sequences
      #
      # Args:
      #   fourtybase.region.surrounding.targetseq.vec: 40bp region surrounding the target site.
      #   protospacer.start.pos.num: The position where protospacer involved in 40bp region surrounding the target site starts.
      #   protospacer.end.pos.num: The position where protospacer involved in 40bp region surrounding the target site ends.
      #   scaffold.seq.vec: The each sgRNA scaffold sequence.
      #                   It should be length(scaffold.seq.vec) == length(fourtybase.region.surrounding.targetseq.vec)
      #   need.annotation.logical: a logical, if TRUE, the the value of annotation.vec will be validated.
      #   annotation.vec: annotation names is to be given to each 40bp region surrounding the target site.
      #                   It should be length(annotation.vec) == length(fourtybase.region.surrounding.targetseq.vec)
      #                   If need.annotation.logical == FALSE, the value of annotation.vec will be invalidated.
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.
      
      tryCatch({
        # Check data
        private$checkSeqVec(fourtybase.region.surrounding.targetseq.vec)
        private$checkMode(protospacer.start.pos.num, "protospacer.start.pos.num", "numeric")
        private$checkMode(protospacer.end.pos.num, "protospacer.end.pos.num", "numeric")
        private$checkSeqVec(scaffold.seq.vec)
        stopifnot(length(fourtybase.region.surrounding.targetseq.vec) == length(scaffold.seq.vec))
        private$checkMode(need.annotation.logical, "need.annotation.logical", "logical")
        if(need.annotation.logical){
          private$checkSeqVec(annotation.vec)
          stopifnot(length(fourtybase.region.surrounding.targetseq.vec) == length(annotation.vec))
        }
        
        # Set information
        self$fourtybase.region.surrounding.targetseq.vec = fourtybase.region.surrounding.targetseq.vec
        self$protospacer.start.pos.num = protospacer.start.pos.num
        self$protospacer.end.pos.num = protospacer.end.pos.num
        self$scaffold.seq.vec = scaffold.seq.vec
        self$need.annotation.logical = need.annotation.logical
        if(self$need.annotation.logical){
          self$annotation.vec = annotation.vec
        }else{
          annotation.vec = character(length(self$fourtybase.region.surrounding.targetseq.vec))
        }
        
        # Make protospacer seq set
        self$protospacer.seq.vec = as.character(
          subseq(
            DNAStringSet(self$fourtybase.region.surrounding.targetseq.vec),
            start = self$protospacer.start.pos.num,
            end = self$protospacer.end.pos.num
          )
        )
        
        # Compose tables which contains sequences and annotations
        self$fourtybase.region.surrounding.targetseq.df <- data.frame(whole.seq = self$fourtybase.region.surrounding.targetseq.vec,
                                                                annotation = self$annotation.vec, stringsAsFactors = FALSE)
        self$protospacer.seq.df <- data.frame(protospacer.seq = self$protospacer.seq.vec,
                                        annotation = self$annotation.vec, stringsAsFactors = FALSE)
        self$scaffold.seq.df <- data.frame(scaffold.seq = self$scaffold.seq.vec,
                                     annotation = self$annotation.vec, stringsAsFactors = FALSE)
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
    RetrieveSeqset = function(type.char, annotation.char = "all"){
      # Retrieve desired sequence set
      #
      # Args:
      #   type.char: sequence type {"whole"||"proto"||"scaffold"}
      #   annotation.char : annotation name, if "all", return all sequences which the detaset has.
      #
      # Returns:
      #   On success, returns vector of desired sequences.
      #   On failure, returns ERROR MESSAGE.
      #
      # Error handling
      stopifnot(type.char == "whole" || type.char == "proto" || type.char == "scaffold")
      stopifnot(mode(annotation.char) == "character")
      
      if(annotation.char == "all"){
        if(type.char == "whole"){
          return(self$fourtybase.region.surrounding.targetseq.df$whole.seq)
        }else if(type.char == "proto"){
          return(self$protospacer.seq.df$protospacer.seq)
        }else if(type.char == "scaffold"){
          return(self$scaffold.seq.df$scaffold.seq)
        }else{
          stop("Unexpected Erorr")
        }
      }else{
        if(type.char == "whole"){
          return(self$fourtybase.region.surrounding.targetseq.df$whole.seq[self$fourtybase.region.surrounding.targetseq.df$annotation == annotation.char])
        }else if(type.char == "proto"){
          return(self$protospacer.seq.df$protospacer.seq[self$protospacer.seq.df$annotation == annotation.char])
        }else if(type.char == "scaffold"){
          return(self$scaffold.seq.df$scaffold.seq[self$scaffold.seq.df$annotation == annotation.char])
        }else{
          stop("Unexpected Erorr")
        }
      }
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
