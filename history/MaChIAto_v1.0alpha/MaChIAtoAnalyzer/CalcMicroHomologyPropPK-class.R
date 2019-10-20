
# load library
library(R6)
library(Biostrings)
library(feather)

# Calculate sgRNA Propaties using the same procedure as KUAN, Pei Fen, et al (2017).
#
# Note:
#   This program was inspired by following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   I appreciate their great work.

CalcMicroHomologyPropPK <- R6Class(
  # name
  "CalcMicroHomologyPropPK",
  # public menber
  public = list(
    
    whole.region.microhomology.seq.vec = character(0), # 40bp region surrounding the target site.
    focusing.region.start.pos.num = 0, # The position where focusing.region involved in 40bp region surrounding the target site starts.
    focusing.region.end.pos.num = 0, # The position where focusing.region involved in 40bp region surrounding the target site ends.
    scaffold.seq.vec = character(0), # The each sgRNA scaffold sequence.
    is.left.microhomology.logical = FALSE,
    annotation.vec = character(0), # annotation names is to be given to each 40bp region surrounding the target site.
    align.save.file.char = "", # a connection or the name of the file where the alignment score list is saved to or read from.
    script.dir = "", # script.dir: directory path contains scripts.
    ref.dir = "", # ref.dir: directory path contains reference tables. 
    
    # instance
    microhomology.info.FormatMicroHomologyInfo = NULL,
    mononuc.high.CalcNNucFreq = NULL,
    mononuc.low.CalcNNucFreq = NULL,
    dinuc.high.CalcNNucFreq = NULL,
    dinuc.low.CalcNNucFreq = NULL,
    future.thermo.CalcThermo = NULL, 
    future.packer.CalcPacker = NULL, 
    future.phychem.CalcPhyChem = NULL, 
    future.pseknc.CalcPseKNC.1.05.2 = NULL, 
    future.align.CalcAlign = NULL,
    
    # dataset of sequences
    future.thermo.df = NULL, # thermodynamics and secondary structure properties (Thermo)
    future.packer.df = NULL, # DNA secondary structures based on dinucleotide and tetra nucleotide properties (Packer)
    future.phychem.df = NULL, # physiochemical propertiess (PhyChem)
    future.pseknc.df = NULL, # pseudo k-tuple nucleotide composition  (PseKNC)
    future.align.list = NULL, # optimal pairwise alignment  (Align)
    
    # future
    future.freq.tstatic.vec = NULL, # single and dinucleotide frequencies (Freq)
    future.pdmono.odds.ratio.vec = NULL, # position-dependent mono-nucleotide composition (PD Mono)
    future.pddinuc.odds.ratio.vec = NULL, # position-dependent dinucleotide composition (PD Dinuc)
    future.thermo.tstatic.vec = NULL, # thermodynamics and secondary structure properties (Thermo)
    future.packer.tstatic.vec = NULL, # DNA secondary structures based on dinucleotide and tetra nucleotide properties (Packer)
    future.phychem.tstatic.vec = NULL, # physiochemical propertiess (PhyChem)
    future.pseknc.tstatic.vec = NULL, # pseudo k-tuple nucleotide composition  (PseKNC)
    future.align.tstatic.vec = NULL, # optimal pairwise alignment  (Align)
    
    # colnames
    colnames.freq.vec = NULL,
    colnames.thermo.vec = NULL,
    colnames.packer.vec = NULL,
    colnames.phychem.vec = NULL,
    colnames.pseknc.vec = NULL,
    colnames.align.vec = NULL,
    
    # param
    align.range.vec = numeric(0),
    
    # properties is calculated
    all.continuous.properties.df = NULL,
    all.binary.properties.df = NULL,
    all.properties.df = NULL,
    
    # properties sorted by annotation
    high.continuous.properties.df = NULL,
    low.continuous.properties.df = NULL,
    high.properties.df = NULL,
    low.properties.df = NULL,
    
    log.2.odds.ratio.vec = numeric(0),
    t.static.vec = numeric(0),
    
    # constructor
    initialize = function(whole.region.microhomology.seq.vec = character(0),
                          focusing.region.start.pos.num = 1,
                          focusing.region.end.pos.num = 10,
                          is.left.microhomology.logical = TRUE,
                          scaffold.seq.vec = character(0),
                          need.annotation.logical = FALSE,
                          annotation.vec = character(0),
                          align.range.vec = 5:20,
                          align.save.file.char = "",
                          script.dir = "",
                          ref.dir = "") {
      # Set sequences as the validation set
      #
      # Args:
      #   whole.region.microhomology.seq.vec: 40bp region surrounding the target site.
      #   focusing.region.start.pos.num: The position where focusing.region involved in 40bp region surrounding the target site starts.
      #   focusing.region.end.pos.num: The position where focusing.region involved in 40bp region surrounding the target site ends.
      #   scaffold.seq.vec: The each sgRNA scaffold sequence.
      #                   It should be length(scaffold.seq.vec) == length(whole.region.microhomology.seq.vec)
      #   need.annotation.logical: a logical, if TRUE, the the value of annotation.vec will be validated.
      #   annotation.vec: annotation names is to be given to each 40bp region surrounding the target site.
      #                   It should be length(annotation.vec) == length(whole.region.microhomology.seq.vec)
      #                   If need.annotation.logical == FALSE, the value of annotation.vec will be invalidated.
      #   align.range.vec : The range of alighment
      #   align.save.file.char: a connection or the name of the file where the alignment score list is saved to or read from.
      #   ref.dir: directory path contains reference tables. 
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.
      
      tryCatch({
        message("Check data")
        # Check data
        private$checkMode(align.save.file.char, "align.save.file.char", "character")
        private$checkMode(ref.dir, "ref.dir", "character")
        private$checkMode(script.dir, "script.dir", "character")
        stopifnot(is.vector(align.range.vec))
        private$checkMode(focusing.region.start.pos.num, "focusing.region.start.pos.num", "numeric")
        private$checkMode(focusing.region.end.pos.num, "focusing.region.end.pos.num", "numeric")


        # Set data
        self$align.save.file.char <- align.save.file.char
        self$align.range.vec <- align.range.vec
        self$focusing.region.start.pos.num <- focusing.region.start.pos.num
        self$focusing.region.end.pos.num <- focusing.region.end.pos.num
        self$scaffold.seq.vec <- scaffold.seq.vec
        self$script.dir <- script.dir
        self$ref.dir <- ref.dir
        # Load class
        source(file.path(self$script.dir, "FormatMicroHomologyInfo-class.R"))
        source(file.path(self$script.dir, "CalcNNucFreq-class.R"))
        source(file.path(self$script.dir, "CalcThermo-class.R"))
        source(file.path(self$script.dir, "CalcPacker-class.R"))
        source(file.path(self$script.dir, "CalcPhyChem-class.R"))
        source(file.path(self$script.dir, "CalcPseKNC-class.R"))
        source(file.path(self$script.dir, "CalcAlign-class.R"))

        message("Set Microhomology information")
        # create instance of FormatMicroHomologyInfo
        self$microhomology.info.FormatMicroHomologyInfo <- FormatMicroHomologyInfo$new(
          whole.region.microhomology.seq.vec,
          focusing.region.start.pos.num,
          focusing.region.end.pos.num,
          is.left.microhomology.logical,
          scaffold.seq.vec,
          need.annotation.logical,
          annotation.vec)

        
        # check data size
        message("Sequence size")
        for(annotation.char in unique(annotation.vec)){
          message(annotation.char, " : ", length(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("whole", annotation.char)), " sequences")
        }
        message("Total : ", length(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("whole", "all")), " sequences")

        
        # Calculate futures
        private$CalcFutureNN() # Freq, PD-Mono, PD-Dinuc
        private$CalcFutureThermo() # Thermo
        private$CalcFuturePacker() # Packer
        private$CalcFuturePhyChem() # PhyChem
        private$CalcFuturePseKNC() # PseKNC
        private$CalcFutureAlign() # Align

        private$MakeBinaryPropertiesTable() # Make Binary Properties Table

        private$CalcTstatic() # Calc t-static
        
        private$MakePropertiesTable() # Make Properties Table
        
        

        ###################################################################
        
        
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
    ShowFutureAssociation = function() {
      # Show list of t-value or odds ratio
      #
      # Returns:
      #   returns list.
      return(
        list(
          Freq = self$future.freq.tstatic.vec,
          PD.Mono = self$future.pdmono.odds.ratio.vec,
          PD.Dinuc = self$future.pddinuc.odds.ratio.vec,
          Thermo = self$future.thermo.tstatic.vec,
          Packer = self$future.packer.tstatic.vec,
          PhyChem = self$future.phychem.tstatic.vec,
          PseKNC = self$future.pseknc.tstatic.vec,
          Align = self$future.align.tstatic.vec
        )
      )
    },
    ShowDataValue = function(datatype) {
      # Show list of evaluated value per sequence
      #
      # Args:
      #   datatype: data type {"all"|"binary"|"continuous"}
      #
      # Returns:
      #   On success, returns data.frame.
      #   On failure, returns Error.
      
      res.df <- rbind(
        cbind(self$high.properties.df,
              label = rep(1, length(self$high.properties.df[,1])),
              stringsAsFactors = FALSE
        ),
        cbind(self$low.properties.df
              , label = rep(0, length(self$low.properties.df[,1])),
              stringsAsFactors = FALSE
        ),
        stringsAsFactors = FALSE
      )
      if(datatype == "all"){
        return(res.df)
      }else if(datatype == "binary"){
        return(res.df[, c(colnames(self$all.binary.properties.df), "label")])
      }else if(datatype == "continuous"){
        return(res.df[, c(colnames(self$all.continuous.properties.df), "label")])
      }else{
        stop("legal argument are only all or binary or continuous. Please check the input.")
      }
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
        stop("Argument ", argname, " must be a ", mode, " : your input was ", input.mode, ".\n")
        return(FALSE)
      } else {
        return(TRUE)
      }
    },
    calcLogOdds = function(pos.mat, neg.mat) {
      # calculate log odds ratio
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
        stopifnot(is.matrix(pos.mat))
        stopifnot(is.matrix(neg.mat))
        
        # Calculate Lod odds ratios
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
    },
    CalcFutureNN = function() {
      message("--------------------------------")
      message("Calculate futures of nucleotide composition")
      message("--------------------------------")
      # Calc nucleotide composition
      message("Mononucleotide")
      message("high sequence set")
      self$mononuc.high.CalcNNucFreq <- CalcNNucFreq$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("whole", "high"),
                                                              N = 1)
      message("low sequence set")
      self$mononuc.low.CalcNNucFreq <- CalcNNucFreq$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("whole", "low"),
                                                                N = 1)
      message("Dinucleotide")
      message("high sequence set")
      self$dinuc.high.CalcNNucFreq <- CalcNNucFreq$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("whole", "high"),
                                                            N = 2)
      message("low sequence set")
      self$dinuc.low.CalcNNucFreq <- CalcNNucFreq$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("whole", "low"),
                                                              N = 2)
      message("--------------------------------")
      message("--------------------------------")
      message("Calculate log odds ratios of futures that are binary variables,")
      message("--------------------------------")
      
      # get futures(Freq, PD Mono, PD Dinuc)
      # message("single and dinucleotide frequencies (Freq)")
      # self$future.freq.odds.ratio.vec <-  c(private$calcLogOdds(self$mononuc.high.CalcNNucFreq$summary$frequency,
      #                                                 self$mononuc.low.CalcNNucFreq$summary$frequency),
      #                            private$calcLogOdds(self$dinuc.high.CalcNNucFreq$summary$frequency,
      #                                                 self$dinuc.low.CalcNNucFreq$summary$frequency)
      # )
      message("position-dependent mono-nucleotide composition (PD Mono)")
      self$future.pdmono.odds.ratio.vec <-  private$calcLogOdds(self$mononuc.high.CalcNNucFreq$summary$PDfrequency,
                                                                self$mononuc.low.CalcNNucFreq$summary$PDfrequency)
      message("position-dependent dinucleotide composition (PD Dinuc)")
      self$future.pddinuc.odds.ratio.vec <-  private$calcLogOdds(self$dinuc.high.CalcNNucFreq$summary$PDfrequency,
                                                                 self$dinuc.low.CalcNNucFreq$summary$PDfrequency)
    },
    CalcFutureThermo = function() {
      message("--------------------------------")
      message("Calculate futures of thermodynamics and secondary structure properties (Thermo)")
      message("--------------------------------")
      self$future.thermo.CalcThermo  <- CalcThermo$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"), 0.2, self$script.dir, self$ref.dir)
      self$future.thermo.df <- self$future.thermo.CalcThermo$summary
      message("--------------------------------")
    },
    CalcFuturePacker = function() {
      message("--------------------------------")
      message("Calculate futures of DNA secondary structures based on dinucleotide and tetra nucleotide properties (Packer)")
      message("--------------------------------")
      self$future.packer.CalcPacker <- CalcPacker$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"), self$ref.dir)
      self$future.packer.df <- self$future.packer.CalcPacker$summary
      message("--------------------------------")
    },
    CalcFuturePhyChem = function() {
      message("--------------------------------")
      message("Calculate futures of physiochemical propertiess (PhyChem)")
      message("--------------------------------")
      self$future.phychem.CalcPhyChem <- CalcPhyChem$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"), self$ref.dir)
      self$future.phychem.df <- self$future.phychem.CalcPhyChem$summary
      message("--------------------------------")
    },
    CalcFuturePseKNC = function() {
      message("--------------------------------")
      message("Calculate futures of pseudo k-tuple nucleotide composition  (PseKNC)")
      message("--------------------------------")
      self$future.pseknc.CalcPseKNC.1.05.2 <- CalcPseKNC$new(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"),
                                                             lambda = 1, w = 0.5, k = 2, self$ref.dir)
      self$future.pseknc.df <- self$future.pseknc.CalcPseKNC.1.05.2$summary
      message("--------------------------------")
    },
    CalcFutureAlign = function() {
      message("--------------------------------")
      message("Calculate futures of optimal pairwise alignment  (Align)")
      message("--------------------------------")
      if(!file.exists(self$align.save.file.char)){
        self$future.align.list <- list(self$align.range.vec)
        message("Perform optimal pairwise alignments between the scaffold sequence and ...")
        for(k in self$align.range.vec){
          message(k,"bp region")
          message("-----------------------------")

          #if(self$microhomology.info.FormatMicroHomologyInfo$is.left.microhomology.logical){
          #    self$future.align.CalcAlign <- CalcAlign$new(
          #        self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"),
          #        len1.num = k, ori1.char = "3' to 5'",
          #        self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("scaffold", "all"),
          #        len2.num = k, ori2.char = "5' to 3'",
          #        is.sense.sense = TRUE)
          #    self$future.align.list[which(self$align.range.vec %in% k)] <- list(self$future.align.CalcAlign$summary)
          #}else{
          #    self$future.align.CalcAlign <- CalcAlign$new(
          #        self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"),
          #        len1.num = k, ori1.char = "5' to 3'",
          #        self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("scaffold", "all"),
          #        len2.num = k, ori2.char = "5' to 3'",
          #        is.sense.sense = TRUE)
          #    self$future.align.list[which(self$align.range.vec %in% k)] <- list(self$future.align.CalcAlign$summary)
          #}

          if(self$microhomology.info.FormatMicroHomologyInfo$is.left.microhomology.logical){
              self$future.align.CalcAlign <- CalcAlign$new(
                  self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"),
                  len1.num = nchar(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all")[1]),
                  ori1.char = "3' to 5'",
                  self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("scaffold", "all"),
                  len2.num = k, ori2.char = "5' to 3'",
                  is.sense.sense = TRUE)
              self$future.align.list[which(self$align.range.vec %in% k)] <- list(self$future.align.CalcAlign$summary)
          }else{
              self$future.align.CalcAlign <- CalcAlign$new(
                  self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all"),
                  len1.num = nchar(self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "all")[1]),
                  ori1.char = "5' to 3'",
                  self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("scaffold", "all"),
                  len2.num = k, ori2.char = "3' to 5'",
                  is.sense.sense = FALSE)
              self$future.align.list[which(self$align.range.vec %in% k)] <- list(self$future.align.CalcAlign$summary)
          }
            
          message("-----------------------------")
          message("Done.")
        }
        names(self$future.align.list) <- paste(as.character(self$align.range.vec), "bp", sep = "")
        message("Save alignment score list in ", self$align.save.file.char)
        saveRDS(self$future.align.list, self$align.save.file.char)
      }else{
        message("Load alignment score list from the file : ", self$align.save.file.char)
        self$future.align.list <- readRDS(self$align.save.file.char)
        message("Done.")
      }
    },
    MakeBinaryPropertiesTable = function() {
      self$all.binary.properties.df <- rbind(data.frame(
        surrounding.sequence = self$mononuc.high.CalcNNucFreq$summary$whole.sequence,
        sequence = as.character(
          subseq(
            DNAStringSet(self$mononuc.high.CalcNNucFreq$summary$whole.sequence),
            start = self$focusing.region.start.pos.num,
            end = self$focusing.region.end.pos.num
          )
        ),
        #self$mononuc.high.CalcNNucFreq$summary$frequency,
        #self$dinuc.high.CalcNNucFreq$summary$frequency,
        self$mononuc.high.CalcNNucFreq$summary$PDfrequency,
        self$dinuc.high.CalcNNucFreq$summary$PDfrequency,
        stringsAsFactors = FALSE
      ),
      data.frame(
        surrounding.sequence = self$mononuc.low.CalcNNucFreq$summary$whole.sequence,
        sequence = as.character(
          subseq(
            DNAStringSet(self$mononuc.low.CalcNNucFreq$summary$whole.sequence),
            start = self$focusing.region.start.pos.num,
            end = self$focusing.region.end.pos.num
          )
        ),
        #self$mononuc.low.CalcNNucFreq$summary$frequency,
        #self$dinuc.low.CalcNNucFreq$summary$frequency,
        self$mononuc.low.CalcNNucFreq$summary$PDfrequency,
        self$dinuc.low.CalcNNucFreq$summary$PDfrequency,
        stringsAsFactors = FALSE
      ))
    },
    MakePropertiesTable = function(){
      self$all.properties.df <- merge(self$all.binary.properties.df,
                                      self$all.continuous.properties.df,
                                      by = "sequence", all=T
      )
      
      self$high.properties.df <- data.frame(sequence = self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "high"),
                                                 stringsAsFactors = FALSE)
      self$low.properties.df <- data.frame(sequence = self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "low"),
                                                   stringsAsFactors = FALSE)
      self$high.properties.df <- (merge(self$high.properties.df, self$all.properties.df, by="sequence", all=F))
      self$low.properties.df <- (merge(self$low.properties.df, self$all.properties.df, by="sequence", all=F))
    },
    MakeContinuousPropertiesTable = function() {
      # Collection of continuous dataset 
      # freq
      self$colnames.freq.vec <- c(colnames(self$mononuc.high.CalcNNucFreq$summary$frequency),
                                  colnames(self$dinuc.high.CalcNNucFreq$summary$frequency))
      self$all.continuous.properties.df <- data.frame(
        sequence = c(
          as.character(
            subseq(
              DNAStringSet(self$mononuc.high.CalcNNucFreq$summary$whole.sequence),
              start = self$focusing.region.start.pos.num,
              end = self$focusing.region.end.pos.num
            )),
          as.character(
            subseq(
              DNAStringSet(self$mononuc.low.CalcNNucFreq$summary$whole.sequence),
              start = self$focusing.region.start.pos.num,
              end = self$focusing.region.end.pos.num
            ))
        )
        ,stringsAsFactors = FALSE)
      self$all.continuous.properties.df <- cbind(self$all.continuous.properties.df, rbind(self$mononuc.high.CalcNNucFreq$summary$frequency,
                                                                                          self$mononuc.low.CalcNNucFreq$summary$frequency)
      )
      self$all.continuous.properties.df <- cbind(self$all.continuous.properties.df, rbind(self$dinuc.high.CalcNNucFreq$summary$frequency,
                                                                                          self$dinuc.low.CalcNNucFreq$summary$frequency)
      )
      
      # thermo
      self$colnames.thermo.vec <- colnames(self$future.thermo.df)[-1]
      self$all.continuous.properties.df <- (merge(self$all.continuous.properties.df, self$future.thermo.df, by="sequence", all=T))
      # Change a part of colnames
      temp.df <- self$future.packer.df
      # Packer
      colnames(temp.df)[2:7] <- c("min.tetranuc.E" ,
                                  "max.tetranuc.E",
                                  "mean.tetranuc.E",
                                  "min.tetranuc.flex",
                                  "max.tetranuc.flex",
                                  "mean.tetranuc.flex")
      self$colnames.packer.vec <- colnames(temp.df)[-1]
      self$all.continuous.properties.df <- (merge(self$all.continuous.properties.df, temp.df, by="sequence", all=T))
      
      # PhyChem
      temp.df <- self$future.phychem.df
      colnames(temp.df) <- c("sequence", paste(
        c(rep("mean.", 12), rep("max.", 12), rep("min.", 12)),
        c("A.philicity", "base.stacking", "B.DNA twist", "bendability", "DNA.bending.stiffness","DNA.denaturation",
          "duplex.disrupt.energy", "duplex.free.energy", "propeller.twist", "protein.deformation", "protein-DNA.twist","Z.DNA"),
        sep = ""))
      self$colnames.phychem.vec <- colnames(temp.df)[-1]
      self$all.continuous.properties.df <- (merge(self$all.continuous.properties.df, temp.df, by="sequence", all=T))
      
      # PseKNC
      temp.df <- self$future.pseknc.df
      colnames(temp.df) <- c("sequence", paste("PseKNC.", as.character(1:16), ".",
                                               str_extract(colnames(self$future.pseknc.df)[-1], pattern="[ATCG]{2}$"),
                                               sep = "")
      )
      self$colnames.pseknc.vec <- colnames(temp.df)[-1]
      self$all.continuous.properties.df <- (merge(self$all.continuous.properties.df, temp.df, by="sequence", all=T))
      
      # Align
      temp.df <- NULL
      temp.df <- self$future.align.list[["20bp"]]$sequence1 # Start alignment by 20 bp
      for(name in paste(as.character(self$align.range.vec), "bp", sep = "")){
        temp.df <- data.frame(temp.df,
                              self$future.align.list[[name]]$align.score.global.num,
                              stringsAsFactors = FALSE)
      }
      colnames(temp.df) <- c("sequence", paste("Align.", as.character(self$align.range.vec), sep = ""))
      self$colnames.align.vec <- colnames(temp.df)[-1]
      self$all.continuous.properties.df <- (merge(self$all.continuous.properties.df, temp.df, by="sequence", all=T))
      
      # Split it by annotation
      self$high.continuous.properties.df <- data.frame(sequence = self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "high"),
                                                            stringsAsFactors = FALSE)
      self$low.continuous.properties.df <- data.frame(sequence = self$microhomology.info.FormatMicroHomologyInfo$RetrieveSeqset("focus", "low"),
                                                              stringsAsFactors = FALSE)
      self$high.continuous.properties.df <- (merge(self$high.continuous.properties.df, self$all.continuous.properties.df, by="sequence", all=F))
      self$low.continuous.properties.df <- (merge(self$low.continuous.properties.df, self$all.continuous.properties.df, by="sequence", all=F))
    },
    CalcTstatic = function() {
      private$MakeContinuousPropertiesTable() # Make continuous properties table
      
      message("--------------------------------")
      message("Calculate t-statics of futures that are continuous variables")
      message("--------------------------------")
      
      self$t.static.vec <- numeric(0)
      for(colname in colnames(self$high.continuous.properties.df)[-1]){
        temp.t.static <- NA
        tryCatch({
          temp.ttest.res <- t.test(as.numeric(as.character(self$high.continuous.properties.df[[colname]])),
                                   as.numeric(as.character(self$low.continuous.properties.df[[colname]])),
                                   alternative = "two.sided")
          temp.t.static <- temp.ttest.res$statistic
        }
        , error = function(e) {
          temp.t.static <- NA
        }
        , silent = TRUE
        )
        self$t.static.vec <- c(self$t.static.vec, temp.t.static)
      }
      names(self$t.static.vec) <- colnames(self$high.continuous.properties.df)[-1]
      
      message("single and dinucleotide frequencies (Freq)")
      self$future.freq.tstatic.vec <- self$t.static.vec[self$colnames.freq.vec]
      message("thermodynamics and secondary structure properties (Thermo)")
      self$future.thermo.tstatic.vec <- self$t.static.vec[self$colnames.thermo.vec]
      message("DNA secondary structures based on dinucleotide and tetra nucleotide properties (Packer)")
      self$future.packer.tstatic.vec <- self$t.static.vec[self$colnames.packer.vec]
      message("physiochemical propertiess (PhyChem)")
      self$future.phychem.tstatic.vec <- self$t.static.vec[self$colnames.phychem.vec]
      message("pseudo k-tuple nucleotide composition  (PseKNC)")
      self$future.pseknc.tstatic.vec <- self$t.static.vec[self$colnames.pseknc.vec]
      message("optimal pairwise alignment  (Align)")
      self$future.align.tstatic.vec <- self$t.static.vec[self$colnames.align.vec]
      
      message("--------------------------------")
    }
  )
)
