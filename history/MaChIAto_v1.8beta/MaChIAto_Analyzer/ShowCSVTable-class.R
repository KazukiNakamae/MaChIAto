
library("R6")

ShowCSVTable <- R6Class(
  # name
  "ShowCSVTable", 
  
  # public menber
  public = list(
    fnList  = NA,  # list of file names
    initialize = function(fnList) {
      self$fnList = fnList
    },
    # get table as dataframe by entering index
    GetDataAsTable = function(fnInd) {
      table <- read.csv(self$fnList[fnInd], header = T, stringsAsFactors = F)
      return(table)
    }
  ),
  
  private = list(
  )
)