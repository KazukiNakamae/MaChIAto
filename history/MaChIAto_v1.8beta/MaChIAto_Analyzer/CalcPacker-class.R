
library(R6)

# Calculate DNA secondary structures based on dinucleotide and tetra nucleotide propaties of reports written by Martin J Packer et al.
#
# Note:
#   This program was refered to following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A. Sequence-dependent DNA structure: dinucleotide conformational maps. Journal of molecular biology, 2000, 295.1: 71-83.
#   PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A. Sequence-dependent DNA structure: tetranucleotide conformational maps. Journal of molecular biology, 2000, 295.1: 85-103.
#   I am grateful for them!
CalcPacker <- R6Class(
  # name
  "CalcPacker",
  # public menber
  public = list(

    seq  = NA,  # DNAStringSet objects of sequences
    seq.vec  = character(0),  # vector of sequences
    ref.dir = "", # ref.dir: directory path contains reference tables. 
    #Table of tetranucleotide minimum energy and flexibility was refered to PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A. Sequence-dependent DNA structure: tetranucleotide conformational maps. Journal of molecular biology, 2000, 295.1: 85-103.
    tetraenergy = NULL,
    #The flexibility of "TCAG" and "CTGA" are not listed in Table 3 of the paper. I supposed the flexibility of them is 0.
    tetraflex = NULL,
    seqtetraEmin.vec  = character(0),  # vector of the lowest Minimum energy (kJ mol 1) tetranucleotide in sequences
    seqtetraEmax.vec  = character(0),  # vector of the highest Minimum energy (kJ mol 1) tetranucleotide in sequences
    seqtetraEave.vec  = character(0),  # vector of the average Minimum energy (kJ mol 1) tetranucleotide in sequences
    seqtetraFmin.vec  = character(0),  # vector of the lowest Flexibility (kJ mol 1 AÊ 2) for tetranucleotide in sequences
    seqtetraFmax.vec  = character(0),  # vector of the highest Flexibility (kJ mol 1 AÊ 2) for tetranucleotide in sequences
    seqtetraFave.vec  = character(0),  # vector of the average Flexibility (kJ mol 1 AÊ 2) for tetranucleotide in sequences
    seqdiRollmin.vec  = character(0),  # vector of the lowest Rotation (deg.) for dinucleotide in sequences
    seqdiRollmax.vec  = character(0),  # vector of the highest Rotation (deg.) for dinucleotide in sequences
    seqdiRollave.vec  = character(0),  # vector of the average Rotation (deg.) for dinucleotide in sequences
    seqdiTwistmax.vec  = character(0),  # vector of the lowest Twist (deg.) for dinucleotide in sequences
    seqdiTwistmin.vec  = character(0),  # vector of the highest Twist (deg.) for dinucleotide in sequences
    seqdiTwistave.vec  = character(0),  # vector of the average Twist (deg.) for dinucleotide in sequences
    seqdiSlidemin.vec  = character(0),  # vector of the lowest Slide (kJ / mol A^2) for dinucleotide in sequences
    seqdiSlidemax.vec  = character(0),  # vector of the highest Slide (kJ / mol A^2) for dinucleotide in sequences
    seqdiSlideave.vec  = character(0),  # vector of the average Slide (kJ / mol A^2) for dinucleotide in sequences
    seqdiShiftmin.vec  = character(0),  # vector of the lowest Shift (kJ / mol A^2) for dinucleotide in sequences
    seqdiShiftmax.vec  = character(0),  # vector of the highest Shift (kJ / mol A^2) for dinucleotide in sequences
    seqdiShiftave.vec  = character(0),  # vector of the average Shift (kJ / mol A^2) for dinucleotide in sequences
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
        self$tetraenergy <- read.csv(file.path(self$ref.dir, "packer_alltetranuc_energy.csv"), stringsAsFactors = FALSE)
        self$tetraflex <- read.csv(file.path(self$ref.dir, "packer_alltetranuc_flex.csv"), stringsAsFactors = FALSE)
        
        message("Set sequences")
        
        self$seq.vec <- as.character(self$seq)

        message("Calculate propaties")
        
        message("1.Tetranucleotide Energy")
        
        # Tetranucleotide Energy
        for (seq in seq.vec) {
          Eprofile = self$ProfileTetraEnegy(seq)
          self$seqtetraEmin.vec <- c(self$seqtetraEmin.vec, Eprofile["minE"])
          self$seqtetraEmax.vec <- c(self$seqtetraEmax.vec, Eprofile["maxE"])
          self$seqtetraEave.vec <- c(self$seqtetraEave.vec, Eprofile["aveE"])
        }

        message("2.Tetranucleotide Flexibility")
        
        # Tetranucleotide Flexibility
        for (seq in seq.vec) {
          Fprofile = self$ProfileTetraFlex(seq)
          self$seqtetraFmin.vec <- c(self$seqtetraFmin.vec, Fprofile["minF"])
          self$seqtetraFmax.vec <- c(self$seqtetraFmax.vec, Fprofile["maxF"])
          self$seqtetraFave.vec <- c(self$seqtetraFave.vec, Fprofile["aveF"])
        }

        message("3.Dinucleotide Roll")
        
        # dinucleotide Roll
        for (seq in seq.vec) {
          Rollprofile = self$ProfileDinucRoll(seq)
          self$seqdiRollmin.vec <- c(self$seqdiRollmin.vec, Rollprofile["minRoll"])
          self$seqdiRollmax.vec <- c(self$seqdiRollmax.vec, Rollprofile["maxRoll"])
          self$seqdiRollave.vec <- c(self$seqdiRollave.vec, Rollprofile["aveRoll"])
        }

        message("4.Dinucleotide Twist")
        
        # dinucleotide Twist
        for (seq in seq.vec) {
          Twistprofile = self$ProfileDinucTwist(seq)
          self$seqdiTwistmin.vec <- c(self$seqdiTwistmin.vec, Twistprofile["minTwist"])
          self$seqdiTwistmax.vec <- c(self$seqdiTwistmax.vec, Twistprofile["maxTwist"])
          self$seqdiTwistave.vec <- c(self$seqdiTwistave.vec, Twistprofile["aveTwist"])
        }

        message("5.Dinucleotide Slide")
        
        # dinucleotide Slide
        for (seq in seq.vec) {
          Slideprofile = self$ProfileDinucSlide(seq)
          self$seqdiSlidemin.vec <- c(self$seqdiSlidemin.vec, Slideprofile["minSlide"])
          self$seqdiSlidemax.vec <- c(self$seqdiSlidemax.vec, Slideprofile["maxSlide"])
          self$seqdiSlideave.vec <- c(self$seqdiSlideave.vec, Slideprofile["aveSlide"])
        }

        message("6.Dinucleotide Shift")
        
        # dinucleotide Shift
        for (seq in seq.vec) {
          Shiftprofile = self$ProfileDinucShift(seq)
          self$seqdiShiftmin.vec <- c(self$seqdiShiftmin.vec, Shiftprofile["minShift"])
          self$seqdiShiftmax.vec <- c(self$seqdiShiftmax.vec, Shiftprofile["maxShift"])
          self$seqdiShiftave.vec <- c(self$seqdiShiftave.vec, Shiftprofile["aveShift"])
        }

        message("Record summary")
        
        # summary
        self$summary <- data.frame(
          self$seq.vec,
          self$seqtetraEmin.vec,
          self$seqtetraEmax.vec,
          self$seqtetraEave.vec,
          self$seqtetraFmin.vec,
          self$seqtetraFmax.vec,
          self$seqtetraFave.vec,
          self$seqdiRollmin.vec,
          self$seqdiRollmax.vec,
          self$seqdiRollave.vec,
          self$seqdiTwistmin.vec,
          self$seqdiTwistmax.vec,
          self$seqdiTwistave.vec,
          self$seqdiSlidemin.vec,
          self$seqdiSlidemax.vec,
          self$seqdiSlideave.vec,
          self$seqdiShiftmin.vec,
          self$seqdiShiftmax.vec,
          self$seqdiShiftave.vec
        )
        colnames(self$summary) <- c("sequence" ,"min.tetranuc.energy", "max.tetranuc.energy", "mean.tetranuc.energy"
                                    , "min.tetranuc.flexibility", "max.tetranuc.flexibility", "mean.tetranuc.flexibility"
                                    , "min.dinuc.roll", "max.dinuc.roll", "mean.dinuc.roll"
                                    , "min.dinuc.twist", "max.dinuc.twist", "mean.dinuc.twist"
                                    , "min.dinuc.slide", "max.dinuc.slide", "mean.dinuc.slide"
                                    , "min.dinuc.shift", "max.dinuc.shift", "mean.dinuc.shift")
        
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
    ProfileTetraEnegy = function(seq) {
      # Profile tetranucleotide energy of sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns vector of result.([1]:minimum energy/[2]:maximum energy/[3]:average energy)
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        tetranuc.vec <- self$ToNnuc(seq, N = 4)
        minE <- Inf
        maxE <- -Inf
        sumE <- 0
        for (tetranuc in tetranuc.vec){
          E = self$tetraenergy[self$tetraenergy$Tetranucleotide == tetranuc,]$Emin
          if(E < minE) minE <- E
          if(E > maxE) maxE <- E
          sumE = sumE + E
        }
        aveE <- sumE / length(tetranuc.vec)

        res <- c(minE, maxE, aveE)
        names(res) <- c("minE", "maxE", "aveE")

        return(res)
      }
      , error = function(e) {
        message(e)
        return(NULL)
      }
      , silent = TRUE
      )
    },
    ProfileTetraFlex = function(seq) {
      # Profile tetranucleotide flexibility of sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns vector of result.([1]:minimum flexibility/[2]:maximum flexibility/[3]:average flexibility)
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        tetranuc.vec <- self$ToNnuc(seq, N = 4)
        minF <- Inf
        maxF <- -Inf
        sumF <- 0
        for (tetranuc in tetranuc.vec){
          F = self$tetraflex[self$tetraflex$Tetranucleotide == tetranuc,]$Flexibility
          if(F < minF) minF <- F
          if(F > maxF) maxF <- F
          sumF = sumF + F
        }
        aveF <- sumF / length(tetranuc.vec)

        res <- c(minF, maxF, aveF)
        names(res) <- c("minF", "maxF", "aveF")

        return(res)
      }
      , error = function(e) {
        message(e)
        return(NULL)
      }
      , silent = TRUE
      )
    },
    ProfileDinucRoll = function(seq) {
      # Profile dinucleotide roll of sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns vector of result.([1]:minimum roll/[2]:maximum roll/[3]:average roll)
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #Roll
        #NN parameters is cited by following article:
        #PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A.
        #Sequence-dependent DNA structure: dinucleotide conformational maps.
        #Journal of molecular biology, 2000, 295.1: 71-83.
        NNroll <- hash()
        NNroll[c("AA", "TT")] <- 2.3
        NNroll["AT"] <- -8.1
        NNroll["TA"] <- 8.4
        NNroll[c("CA", "TG")] <- 7.4
        NNroll[c("GT", "AC")] <- -2.0
        NNroll[c("CT", "AG")] <- 0.5
        NNroll[c("GA", "TC")] <- 5.0
        NNroll["CG"] <- 6.3
        NNroll["GC"] <- -0.4
        NNroll[c("GG", "CC")] <- 1.4

        dinuc.vec <- self$ToNnuc(seq, N = 2)
        minRoll = Inf
        maxRoll = -Inf
        sumRoll = 0
        for (dinuc in dinuc.vec){
          Roll = NNroll[[dinuc]]
          if(Roll < minRoll) minRoll = Roll
          if(Roll > maxRoll) maxRoll = Roll
          sumRoll = sumRoll + Roll
        }
        aveRoll = sumRoll / length(dinuc.vec)

        res = c(minRoll, maxRoll, aveRoll)
        names(res) = c("minRoll", "maxRoll", "aveRoll")

        return(res)
      }
      , error = function(e) {
        message(e)
        return(NULL)
      }
      , silent = TRUE
      )
    },
    ProfileDinucTwist = function(seq) {
      # Profile dinucleotide twist of sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns vector of result.([1]:minimum twist/[2]:maximum twist/[3]:average twist)
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #Twist
        #NN parameters is cited by following article:
        #PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A.
        #Sequence-dependent DNA structure: dinucleotide conformational maps.
        #Journal of molecular biology, 2000, 295.1: 71-83.
        NNtwist <- hash()
        NNtwist[c("AA", "TT")] <- 37.6
        NNtwist["AT"] <- 39.7
        NNtwist["TA"] <- 34.6
        NNtwist[c("CA", "TG")] <- 32.2
        NNtwist[c("GT", "AC")] <- 35.8
        NNtwist[c("CT", "AG")] <- 35.7
        NNtwist[c("GA", "TC")] <- 38.4
        NNtwist["CG"] <- 33.9
        NNtwist["GC"] <- 37.4
        NNtwist[c("GG", "CC")] <- 35.5

        dinuc.vec <- self$ToNnuc(seq, N = 2)
        minTwist = Inf
        maxTwist = -Inf
        sumTwist = 0
        for (dinuc in dinuc.vec){
          Twist = NNtwist[[dinuc]]
          if(Twist < minTwist) minTwist = Twist
          if(Twist > maxTwist) maxTwist = Twist
          sumTwist = sumTwist + Twist
        }
        aveTwist = sumTwist / length(dinuc.vec)

        res = c(minTwist, maxTwist, aveTwist)
        names(res) = c("minTwist", "maxTwist", "aveTwist")

        return(res)
      }
      , error = function(e) {
        message(e)
        return(NULL)
      }
      , silent = TRUE
      )
    },
    ProfileDinucSlide = function(seq) {
      # Profile dinucleotide slide of sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns vector of result.([1]:minimum slide/[2]:maximum slide/[3]:average slide)
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #Slide
        #NN parameters is cited by following article:
        #PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A.
        #Sequence-dependent DNA structure: dinucleotide conformational maps.
        #Journal of molecular biology, 2000, 295.1: 71-83.
        NNslide <- hash()
        NNslide[c("AA", "TT")] <- 13.72
        NNslide["AT"] <- 11.69
        NNslide["TA"] <- 7.13
        NNslide[c("CA", "TG")] <- 1.35
        NNslide[c("GT", "AC")] <- 9.57
        NNslide[c("CT", "AG")] <- 7.58
        NNslide[c("GA", "TC")] <- 10.28
        NNslide["CG"] <- 4.02
        NNslide["GC"] <- 4.34
        NNslide[c("GG", "CC")] <- 7.36

        dinuc.vec <- self$ToNnuc(seq, N = 2)
        minSlide = Inf
        maxSlide = -Inf
        sumSlide = 0
        for (dinuc in dinuc.vec){
          Slide = NNslide[[dinuc]]
          if(Slide < minSlide) minSlide = Slide
          if(Slide > maxSlide) maxSlide = Slide
          sumSlide = sumSlide + Slide
        }
        aveSlide = sumSlide / length(dinuc.vec)

        res = c(minSlide, maxSlide, aveSlide)
        names(res) = c("minSlide", "maxSlide", "aveSlide")

        return(res)
      }
      , error = function(e) {
        message(e)
        return(NULL)
      }
      , silent = TRUE
      )
    },
    ProfileDinucShift = function(seq) {
      # Profile dinucleotide shift of sequence
      #
      # Args:
      #   seq: sequence is to be profiled.
      #
      # Returns:
      #   On success, returns vector of result.([1]:minimum shift/[2]:maximum shift/[3]:average shift)
      #   On failure, returns NULL.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #Shift
        #NN parameters is cited by following article:
        #PACKER, Martin J.; DAUNCEY, Mark P.; HUNTER, Christopher A.
        #Sequence-dependent DNA structure: dinucleotide conformational maps.
        #Journal of molecular biology, 2000, 295.1: 71-83.
        NNshift <- hash()
        NNshift[c("AA", "TT")] <- 5.35
        NNshift["AT"] <- 1.13
        NNshift["TA"] <- 4.28
        NNshift[c("CA", "TG")] <- 4.61
        NNshift[c("GT", "AC")] <- 9.73
        NNshift[c("CT", "AG")] <- 8.98
        NNshift[c("GA", "TC")] <- 5.44
        NNshift["CG"] <- 12.13
        NNshift["GC"] <- 1.98
        NNshift[c("GG", "CC")] <- 5.51

        dinuc.vec <- self$ToNnuc(seq, N = 2)
        minShift = Inf
        maxShift = -Inf
        sumShift = 0
        for (dinuc in dinuc.vec){
          Shift = NNshift[[dinuc]]
          if(Shift < minShift) minShift = Shift
          if(Shift > maxShift) maxShift = Shift
          sumShift = sumShift + Shift
        }
        aveShift = sumShift / length(dinuc.vec)

        res = c(minShift, maxShift, aveShift)
        names(res) = c("minShift", "maxShift", "aveShift")

        return(res)
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
    }
  )
)
