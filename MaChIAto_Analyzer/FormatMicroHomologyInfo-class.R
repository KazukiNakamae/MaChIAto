
library(R6)
library(hash)
library(Biostrings)
# Format sgRNA sequence information to extract it easily.
#
FormatMicroHomologyInfo <- R6Class(
  # name
  "FormatMicroHomologyInfo",
  # public menber
  public = list(
    
    whole.region.microhomology.seq.vec = character(0), # 40bp region surrounding the target site.
    focusing.region.start.pos.num = 0, # The position where focusing.region involved in 40bp region surrounding the target site starts.
    focusing.region.end.pos.num = 0, # The position where focusing.region involved in 40bp region surrounding the target site ends.
    is.left.microhomology.logical = FALSE,
    scaffold.seq.vec = character(0), # The each sgRNA scaffold sequence.
    need.annotation.logical = FALSE, # a logical, if TRUE, the the value of annotation.vec will be validated.
    annotation.vec = character(0), # annotation names is to be given to each 40bp region surrounding the target site.
    
    focusing.region.seq.vec = character(0), # focusing.region
    fourtybase.region.surrounding.targetseq.df = NULL, # table which contains 0bp region surrounding the target sequences and annotations
    focusing.region.seq.df = NULL, # table which contains focusing.region sequences and annotations
    scaffold.seq.df = NULL, # table which contains scaffold sequences and annotations
    
    # constructor
    initialize = function(whole.region.microhomology.seq.vec = character(0),
                          focusing.region.start.pos.num = 1,
                          focusing.region.end.pos.num = 10,
                          is.left.microhomology.logical = TRUE,
                          scaffold.seq.vec = character(0),
                          need.annotation.logical = FALSE,
                          annotation.vec = character(0)
                          ) {
      # Set sequences
      #
      # Args:
      #   whole.region.microhomology.seq.vec: 40bp region surrounding the target site.
      #   focusing.region.start.pos.num: The position where focusing.region involved in 40bp region surrounding the target site starts.
      #   focusing.region.end.pos.num: The position where focusing.region involved in 40bp region surrounding the target site ends.
      #   is.left.microhomology.logical: 
      #   scaffold.seq.vec: The each sgRNA scaffold sequence.
      #                   It should be length(scaffold.seq.vec) == length(whole.region.microhomology.seq.vec)
      #   need.annotation.logical: a logical, if TRUE, the the value of annotation.vec will be validated.
      #   annotation.vec: annotation names is to be given to each 40bp region surrounding the target site.
      #                   It should be length(annotation.vec) == length(whole.region.microhomology.seq.vec)
      #                   If need.annotation.logical == FALSE, the value of annotation.vec will be invalidated.
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.
      
      tryCatch({
        # Check data
        private$checkSeqVec(whole.region.microhomology.seq.vec)
        private$checkMode(focusing.region.start.pos.num, "focusing.region.start.pos.num", "numeric")
        private$checkMode(focusing.region.end.pos.num, "focusing.region.end.pos.num", "numeric")
        private$checkSeqVec(scaffold.seq.vec)
        private$checkMode(is.left.microhomology.logical, "is.left.microhomology.logical", "logical")
        stopifnot(length(whole.region.microhomology.seq.vec) == length(scaffold.seq.vec))
        private$checkMode(need.annotation.logical, "need.annotation.logical", "logical")
        if(need.annotation.logical){
          private$checkSeqVec(annotation.vec)
          stopifnot(length(whole.region.microhomology.seq.vec) == length(annotation.vec))
        }
        # Set information
        self$whole.region.microhomology.seq.vec = whole.region.microhomology.seq.vec
        self$focusing.region.start.pos.num = focusing.region.start.pos.num
        self$focusing.region.end.pos.num = focusing.region.end.pos.num
        self$scaffold.seq.vec = scaffold.seq.vec
        self$is.left.microhomology.logical = is.left.microhomology.logical
        self$need.annotation.logical = need.annotation.logical
        if(self$need.annotation.logical){
            self$annotation.vec = annotation.vec
        }else{
            annotation.vec = character(length(self$whole.region.microhomology.seq.vec))
        }
        
        # Make focusing.region seq set
        self$focusing.region.seq.vec = as.character(
          subseq(
            DNAStringSet(self$whole.region.microhomology.seq.vec),
            start = self$focusing.region.start.pos.num,
            end = self$focusing.region.end.pos.num
          )
        )
        # Compose tables which contains sequences and annotations
        self$fourtybase.region.surrounding.targetseq.df <- data.frame(whole.seq = self$whole.region.microhomology.seq.vec,
                                                                annotation = self$annotation.vec, stringsAsFactors = FALSE)
        self$focusing.region.seq.df <- data.frame(focusing.region.seq = self$focusing.region.seq.vec,
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
      #   type.char: sequence type {"whole"||"focus"||"scaffold"}
      #   annotation.char : annotation name, if "all", return all sequences which the detaset has.
      #
      # Returns:
      #   On success, returns vector of desired sequences.
      #   On failure, returns ERROR MESSAGE.
      #
      # Error handling
      stopifnot(type.char == "whole" || type.char == "focus" || type.char == "scaffold")
      stopifnot(mode(annotation.char) == "character")
      
      if(annotation.char == "all"){
        if(type.char == "whole"){
          return(self$fourtybase.region.surrounding.targetseq.df$whole.seq)
        }else if(type.char == "focus"){
          return(self$focusing.region.seq.df$focusing.region.seq)
        }else if(type.char == "scaffold"){
          return(self$scaffold.seq.df$scaffold.seq)
        }else{
          stop("Unexpected Erorr")
        }
      }else{
        if(type.char == "whole"){
          return(self$fourtybase.region.surrounding.targetseq.df$whole.seq[self$fourtybase.region.surrounding.targetseq.df$annotation == annotation.char])
        }else if(type.char == "focus"){
          return(self$focusing.region.seq.df$focusing.region.seq[self$focusing.region.seq.df$annotation == annotation.char])
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
