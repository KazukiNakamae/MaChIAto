merges <- function(dfs, ...) {
  # merge data.frames
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

calclogodd <- function(pos.mat, neg.mat) {
  # calculate log odd ratio for positive matrix
  #
  # Args:
  #   pos.mat: positive matrix
  #   neg.mat: negative matrix
  #
  # Returns:
  #   On success, returns log odd ratio vector.
  #   On failure, returns null vector.
  
  tryCatch({
    # Error handling
    if (!is.matrix(pos.mat)) {
      #do type identifying for a user input
      if(is.data.frame(pos.mat)){
        Robj = "data.frame"
      } else if (is.list(pos.mat)) {
        Robj = "list"
      } else if (is.array(pos.mat)) {
        Robj = "array"
      } else if (is.factor(pos.mat)) {
        Robj = "factor"
      } else if (is.vector(pos.mat)) {
        Robj = "vector"
      } else {
        Robj = "unknown type"
      }
      stop("Argument neg.mat must be a matrix: your input was ", Robj, ".")
    }
    if (!is.matrix(neg.mat)) {
      #do type identifying for a user input
      if(is.data.frame(neg.mat)){
        Robj = "data.frame"
      } else if (is.list(neg.mat)) {
        Robj = "list"
      } else if (is.array(neg.mat)) {
        Robj = "array"
      } else if (is.factor(neg.mat)) {
        Robj = "factor"
      } else if (is.vector(neg.mat)) {
        Robj = "vector"
      } else {
        Robj = "unknown type"
      }
      stop("Argument neg.mat must be a matrix: your input was ", Robj, ".")
    }
    
    logodd.vec <- numeric(0)
    pos.vec <- apply(pos.mat, 2, sum)
    neg.vec <- apply(neg.mat, 2, sum)
    sum.vec <- pos.vec + neg.vec
    for (property in names(pos.vec)){
      pos.prop <- pos.vec[property] / sum.vec[property]
      neg.prop <- neg.vec[property] / sum.vec[property]
      logodd <- c(log((pos.prop / (1 - pos.prop)) / (neg.prop / (1 - neg.prop))))
      names(logodd) <- property
      logodd.vec <- c(logodd.vec, logodd)
    }
    
    return(logodd.vec)
    
  }
  , error = function(e) {
    message(e)
    return(numeric(0))
  }
  , silent = TRUE
  )
}