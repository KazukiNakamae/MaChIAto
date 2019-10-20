
library(R6)
library(hash)
library(stringr)
library(seqinr)

# Calculate Thermodynamics and secendary structure propaties
#
# Note:
#   This program was refered to following reports :
#   KUAN, Pei Fen, et al. A systematic evaluation of nucleotide properties for CRISPR sgRNA design. BMC bioinformatics, 2017, 18.1: 297.
#   WEI, Hairong, et al. A study of the relationships between oligonucleotide properties and hybridization signal intensities from NimbleGen microarray datasets. Nucleic acids research, 2008, 36.9: 2926-2938.
#   SANTALUCIA, John. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proceedings of the National Academy of Sciences, 1998, 95.4: 1460-1465.
#   I am grateful for them!
CalcThermo <- R6Class(
  # name
  "CalcThermo",
  # public menber
  public = list(

    seq  = NA,  # DNAStringSet objects of sequences
    seq.vec  = character(0),  # vector of sequences
    ref.dir = "", # ref.dir: directory path contains reference tables. 
    script.dir = "",
    seqlen.vec  = numeric(0),  # vector of sequences length
    seqCGcont.vec  = numeric(0),  # vector of %CG of sequences
    seqTm.vec  = numeric(0),  # vector of melting temperature of sequences
    seqdS.vec  = numeric(0),  # vector of entoropy change of sequences
    seqdH.vec  = numeric(0),  # vector of entalpy change of sequences
    seqdG.vec  = numeric(0),  # vector of free Gibbs energy change of sequences
    seqMEF.vec  = numeric(0),  # vector of minimum energy folding of sequences
    seqLSL.vec  = numeric(0),  # vector of the longest length of a potential stem-loop of sequence
    seqRepeat.vec  = numeric(0),  # vector of the copy number of a repetitive sequence
    seqlpolyN.vec  = numeric(0),  # vector of the longest polyN of sequence
    summary  = NULL,  # summary of all propaties in sequences

    #[Na+](Molar)
    Naconc = NA,


    # constructor
    initialize = function(seq.vec, Naconc = 0.2, script.dir, ref.dir) {
      # Set sequences as the validation set
      #
      # Args:
      #   seq.vec: vector of sequences is to be calculated.
      #   Naconc: concentration[M] of Na cation (defalt : 0.2)
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
        Naconc.mode = mode(Naconc)
        if (Naconc.mode == "numeric"){
          self$Naconc <- Naconc
        }else {
          stop("Argument Naconc must be a numeric: your input was ", Naconc.mode, ".")
        }

        # Set infomation
        self$seq <- DNAStringSet(seq.vec) # DNAStringSet objects
        self$script.dir <- script.dir
        self$ref.dir <- ref.dir
        
        message("Set sequences")
        
        # nucleotide sequence
        self$seq.vec <- as.character(self$seq)
        
        #propaty
        
        message("Calculate 10 propaties")
        
        message("1...Length")
        
        # length
        self$seqlen.vec <- width(self$seq)

        message("2...%CG")
        
        # CG contents (%CG)
        self$seqCGcont.vec <- as.vector(letterFrequency(self$seq, letters="CG", as.prob = TRUE)) * 100
        
        message("3...Tm")
        
        # melting temrerature (Tm)
        self$seqTm.vec <- 81.5 + 16.6 * (log10(self$Naconc)) + 0.41 * self$seqCGcont.vec - 600 / self$seqlen.vec

        message("4...∆S")
        
        #@TODO:for loop is very low speed.So make calculation using vector
        # entropy change (∆S)
        for(seq in seq.vec){
          self$seqdS.vec <- c(self$seqdS.vec, self$CalcdS(seq))
        }

        message("5...∆H")
        
        # enthalpy change (∆H)
        for(seq in seq.vec){
          self$seqdH.vec <- c(self$seqdH.vec, self$CalcdH(seq))
        }

        message("6...∆G")
        
        # free energy change (∆G)
        for(seq in seq.vec){
          self$seqdG.vec <- c(self$seqdG.vec, self$CalcdG(seq))
        }

        message("7...MEF")
        
        # minimum energy folding (MEF)
        for(seq in seq.vec){
          self$seqMEF.vec <- c(self$seqMEF.vec, self$CalcMEF(seq))
        }

        message("8...LSL")
        
        # length of a potential stem-loop (LSL)
        for(seq in seq.vec) {
          self$seqLSL.vec <- c(self$seqLSL.vec, self$CalcLSL(seq))
        }

        message("9...Repeat")
        
        # repetitive sequence (repeat)
        for(seq in seq.vec) {
          self$seqRepeat.vec <- c(self$seqRepeat.vec, self$CalcRepeat(seq))
        }

        message("10...PolyN")
        
        # longest polyN
        for(seq in seq.vec){
          self$seqlpolyN.vec <- c(self$seqlpolyN.vec, self$CalcpolyN(seq))
        }
        
        message("Record summary")
        
        # summary
        self$summary <- data.frame(
          self$seq.vec,
          self$seqlen.vec,
          self$seqCGcont.vec,
          self$seqTm.vec,
          self$seqdS.vec,
          self$seqdH.vec,
          self$seqdG.vec,
          self$seqMEF.vec,
          self$seqLSL.vec,
          self$seqRepeat.vec,
          self$seqlpolyN.vec
        )
        colnames(self$summary) <- c("sequence" ,"Length", "%CG", "Tm", "∆S", "∆H", "∆G", "MEF", "LSL", "Repeat", "polyN")
       
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
    CalcdS = function(seq) {
      # Calculate entoropy change (∆S) of sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, entoropy change (∆S) of sequence.
      #   On failure, returns -Inf.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #∆S(unified NN,1M NaCl)
        #NN parameters is cited by the following article:
        #SANTALUCIA, John.
        #A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
        #Proceedings of the National Academy of Sciences, 1998, 95.4: 1460-1465.
        dSInit <- hash()
        dSInit[c("A", "T")] <- 4.1
        dSInit[c("C", "G")] <- -2.8
        dSNN <- hash()
        dSNN[c("AA", "TT")] <- -22.2
        dSNN["AT"] <- -20.4
        dSNN["TA"] <- -21.3
        dSNN[c("CA", "TG")] <- -22.7
        dSNN[c("GT", "AC")] <- -22.4
        dSNN[c("CT", "AG")] <- -21.0
        dSNN[c("GA", "TC")] <- -22.2
        dSNN["CG"] <- -27.2
        dSNN["GC"] <- -24.4
        dSNN[c("GG", "CC")] <- -19.9

        # Calculate entoropy change (∆S) of sequence based on ∆S(unified NN,1M NaCl)

        mononuc.vec <- self$ToNnuc(seq, N = 1)
        dinuc.vec <- self$ToNnuc(seq, N = 2)
        dinuclen <- length(dinuc.vec)
        dS <- 0
        names(dS) = c("∆S[cal/kmol]")

        #add initiation parameter
        dS <- dS + dSInit[[mononuc.vec[1]]]
        dS <- dS + dSInit[[mononuc.vec[length(mononuc.vec)]]]

        #add NN parameter
        for (nnuc_ind in 1:dinuclen){
          dS <- dS + dSNN[[dinuc.vec[nnuc_ind]]]
        }

        #consider solt dependence
        dS <- dS + 0.368 * nchar(seq) * log(self$Naconc)

        return(dS)
      }
      , error = function(e) {
        message(e)
        return(-Inf)
      }
      , silent = TRUE
      )
    },
    CalcdH = function(seq) {
      # Calculate entalpy change (∆H) of sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, entalpy change (∆S) of sequence.
      #   On failure, returns -Inf.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #∆H(unified NN,1M NaCl)
        #NN parameters is cited by following article:
        #SANTALUCIA, John.
        #A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
        #Proceedings of the National Academy of Sciences, 1998, 95.4: 1460-1465.
        dHInit <- hash()
        dHInit[c("A", "T")] <- 2.3
        dHInit[c("C", "G")] <- 0.1
        dHNN <- hash()
        dHNN[c("AA", "TT")] <- -7.9
        dHNN["AT"] <- -7.2
        dHNN["TA"] <- -7.2
        dHNN[c("CA", "TG")] <- -8.5
        dHNN[c("GT", "AC")] <- -8.4
        dHNN[c("CT", "AG")] <- -7.8
        dHNN[c("GA", "TC")] <- -8.2
        dHNN["CG"] <- -10.6
        dHNN["GC"] <- -9.8
        dHNN[c("GG", "CC")] <- -8.0

        # Calculate entoropy change (∆S) of sequence based on ∆S(unified NN,1M NaCl)

        mononuc.vec <- self$ToNnuc(seq, N = 1)
        dinuc.vec <- self$ToNnuc(seq, N = 2)
        dinuclen <- length(dinuc.vec)
        dH <- 0
        names(dH) = c("∆H[kcal/mol]")

        #add initiation parameter
        dH <- dH + dHInit[[mononuc.vec[1]]]
        dH <- dH + dHInit[[mononuc.vec[length(mononuc.vec)]]]

        #add NN parameter
        for (nnuc_ind in 1:dinuclen){
          dH <- dH + dHNN[[dinuc.vec[nnuc_ind]]]
        }

        return(dH)
      }
      , error = function(e) {
        message(e)
        return(-Inf)
      }
      , silent = TRUE
      )
    },
    CalcdG = function(seq) {
      # Calculate Gibbs free energy change (∆G) of sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, Gibbs free energy change (∆G) of sequence.
      #   On failure, returns -Inf.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        #∆G(unified NN,1M NaCl)
        #NN parameters is cited by following article:
        #SANTALUCIA, John.
        #A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
        #Proceedings of the National Academy of Sciences, 1998, 95.4: 1460-1465.
        dGInit <- hash()
        dGInit[c("A", "T")] <- 1.03
        dGInit[c("C", "G")] <- 0.98
        dGNN <- hash()
        dGNN[c("AA", "TT")] <- -1.00
        dGNN["AT"] <- -0.88
        dGNN["TA"] <- -0.58
        dGNN[c("CA", "TG")] <- -1.45
        dGNN[c("GT", "AC")] <- -1.44
        dGNN[c("CT", "AG")] <- -1.28
        dGNN[c("GA", "TC")] <- -1.30
        dGNN["CG"] <- -2.17
        dGNN["GC"] <- -2.24
        dGNN[c("GG", "CC")] <- -1.84

        # Calculate Gibbs free energy change (∆G) of sequence based on ∆G(unified NN,1M NaCl)

        mononuc.vec <- self$ToNnuc(seq, N = 1)
        dinuc.vec <- self$ToNnuc(seq, N = 2)
        dinuclen <- length(dinuc.vec)
        dG <- 0
        names(dG) = c("∆G[kcal/mol]")

        #add initiation parameter
        dG <- dG + dGInit[[mononuc.vec[1]]]
        dG <- dG + dGInit[[mononuc.vec[length(mononuc.vec)]]]

        #add NN parameter
        for (nnuc_ind in 1:dinuclen){
          dG <- dG + dGNN[[dinuc.vec[nnuc_ind]]]
        }

        #consider solt dependence
        dG <- dG - 0.114 * nchar(seq) * log(self$Naconc)

        return(dG)
      }
      , error = function(e) {
        message(e)
        return(-Inf)
      }
      , silent = TRUE
      )
    },
    CalcMEF = function(seq) {
      # Calculate minimum energy folding (MEF) of sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, minimum energy folding (MEF) of sequence.
      #   On failure, returns 0.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")
        command <- paste('echo', seq , '|', file.path(self$script.dir, "hybrid-ss-min"), '-n', 'DNA', '--stream', seq = ' ')
        MEF <- as.numeric(system(command, intern = T))

        return(MEF)
      }
      , error = function(e) {
        message(e)
        return(0)
      }
      , silent = TRUE
      )
    },
    CalcLSL = function(seq) {
      # Calculate length of a potential stem-loop (LSL) of sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, length of a potential stem-loop (LSL) of sequence.
      #   On failure, returns 0.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        command <- paste('echo'
                         , seq
                         , '|'
                         # , file.path(self$script.dir, "palindrome")
                         , file.path("palindrome") # If EMBOSS is globally installed.
                         , '-filter'
                         , '-minpallen'
                         , '1'
                         , '-maxpallen'
                         , nchar(seq)
                         , '-gaplimit'
                         , nchar(seq)
                         , '-nummismatches'
                         , '0'
                         , seq = ' ')
        res <- system(command, intern = T, ignore.stderr = TRUE)
        len.vec = numeric(0)
        for(palInd in 13:length(res)) {
          position <- unlist(str_extract_all(res[palInd], "[0-9]+")) #positon of palindrome : [1]start / [2]end
          if(length(position) == 2) {
            len.vec <- c(len.vec, abs(as.numeric(position[2]) - (as.numeric(position[1])) + 1))
          } else { #no palindrome
            len.vec <- c(len.vec, 0)
          }
        }
        LSL <- max(len.vec)

        return(LSL)
      }
      , error = function(e) {
        message(e)
        return(0)
      }
      , silent = TRUE
      )
    },
    CalcRepeat = function(seq) {
      # Calculate the most copy number of a repetitive sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, the most copy number of a repetitive sequence.
      #   On failure, returns 0.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        temp.dir <- file.path(self$script.dir, "temp")
        system(paste("mkdir", temp.dir, sep = " "))
        #make .fasta file
        write.fasta(seq
                    , names = "seq"
                    , file.out = file.path(temp.dir, "seq.fa")
                    , open = "w"
                    , nbchar = 60
                    , as.string = FALSE)
        #find copy number of a repetitive sequence
        res = system(paste(file.path(self$script.dir, "trf409")
          , file.path(temp.dir, "seq.fa")
          , "2", "5", "7", "80", "10", "0", "2000", "-ngs", "-h", sep = " "), intern = TRUE)
        system(paste("rm", "-r", temp.dir, sep = " ")) # clean
        mostcn <- 0
        if(length(res) >= 2){ # repeat seq > 1
          cn.vec <- numeric(0)
          for(repind in seq(2, length(res), 2)){ #find the MOST copy number of a repetitive sequence
            res.vec <- unlist(strsplit(res[repind], split = " "))
            cn.vec <- c(cn.vec, as.numeric(res.vec[4]))
          }
          mostcn <- max(cn.vec) # the most copy number in repeat seqs
        }

        return(mostcn)
      }
      , error = function(e) {
        message(e)
        return(0)
      }
      , silent = TRUE
      )
    },
    CalcpolyN = function(seq) {
      # Calculate longest polyN of sequence
      #
      # Args:
      #   seq: sequence is to be calculated.
      #
      # Returns:
      #   On success, longest polyN of sequence.
      #   On failure, returns 0.

      tryCatch({
        # Error handling
        private$checkMode(seq, "seq", "character")

        base.vec = unlist(strsplit(seq, split = ""))
        prebase = ""
        polycnt = 1
        polycnt.vec = numeric(0)
        for (base in base.vec) {
          if (prebase == "") {
            prebase = base
          } else if (prebase == base) {
            polycnt = polycnt + 1
          } else {
            polycnt.vec = c(polycnt.vec, polycnt)
            polycnt = 1
            prebase = base
          }
        }
        polycnt.vec = c(polycnt.vec, polycnt)

        maxpolyN = max(polycnt.vec)

        return(maxpolyN)
      }
      , error = function(e) {
        message(e)
        return(0)
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
