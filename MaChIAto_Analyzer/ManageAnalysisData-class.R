
library(R6)
# Manage analysis data of all features about the response variable
#
ManageAnalysisData <- R6Class(
  # name
  "ManageAnalysisData",
  # public menber
  public = list(
    
    value.set = NULL, # The table that includes the rate of response variable and sequence information {data.frame}.
    response.variable = "", # The esponse variable thta the analysis focuses on {character}.
    target.type = "",
    #   "bmh" mode that calculate values of both homology arms
    #   "lmh" mode that calculate values of Left homology arm
    #   "rmh" mode that calculate values of Right homology arm
    #   "protospacer" mode that calculate ones of target site.
    #   "fourtybase.region.surrounding.targetseq" mode that calculate ones of 40 bp srrounding cut site (+17 bp in protospacer).
    should.convert.log = FALSE, # TRUE: the values of response.variable is converted into the log value {Logical}.
    extra.table.list = NULL,  # list of feture table that user add {data.frame}.
    should.remove.max.min = FALSE,  # TRUE: The feature of max/min XXX value is removed in the analysis {logical}.
    save.dir = "",  # saving directory {character}.
    name = "",  # data name {character}.

    fltr.value.set = NULL, # adjusted value.set

    can.convert.log = FALSE,  # TRUE: the values of response.variable cannot be converted into the log value {logical}.

    full.extra.table = NULL, # full feture table that user add {data.frame}.

    microhomology.seqset = NULL, # target microhomology sequence dataset {data.frame}
    focus.length = 0, # range that it aligns {numeric}
    sgRNA.seqset = NULL, #  sgRNA sequence dataset including scaffold {data.frame}
    label.vec = character(0), #  label list according to response variable {vector}
    is.left.both.microhomology.logical = FALSE, # flag represent left microhomology and both microhomology are target {logical}
    out.align.file = "", # filename into which alignment result is saved
    out.prop.file = "", # filename into which the generated instance of CalcMicroHomologyPropPK is saved
    
    target.seqset = NULL, # target sequence dataset {data.frame}
    is.protospacer.logical = FALSE, # flag represent protospacer is target {logical}

    full.continuous.table = NULL, # table has all continuous variables {data.frame}
    continuous.feature.name = character(0), # 
    response.value.by.seq.table = NULL, # table has information of response variable and sequence {data.frame}

    full.feature.vec.list = list(), # list has feature names per feature group {data.frame}
    
    # constructor
    initialize = function(value.set = NULL,
                          response.variable = "",
                          target.type = "",
                          should.convert.log = FALSE,
                          extra.table.list = NULL,
                          should.remove.max.min = FALSE,
                          save.dir = "./",
                          name = ""
                          ) {
      # Set sequences
      #
      # Args:
      #   value.set: The table that includes the rate of response variable and sequence information {data.frame}.
      #   response.variable: The esponse variable thta the analysis focuses on {character}.
      #   target.type:
      #   "bmh" mode that calculate values of both homology arms
      #   "lmh" mode that calculate values of Left homology arm
      #   "rmh" mode that calculate values of Right homology arm
      #   "protospacer" mode that calculate ones of target site.
      #   should.convert.log: TRUE: the values of response.variable is converted into the log value {logical}.
      #   extra.table.list: list of feture table that user add {list}.
      #   should.remove.max.min: TRUE: The feature of max/min XXX value is removed in the analysis {logical}.
      #   save.dir: saving directory {character}; default : current directory
      #
      # Returns:
      #   On success, returns TRUE.
      #   On failure, returns FALSE.
      tryCatch({
        # Check value.set
        stopifnot(is.data.frame(value.set))

        # Check response.variable
        private$checkMode(response.variable, "response.variable", "character")
        response.variable.list <- colnames(value.set)
        stopifnot(response.variable %in% response.variable.list)
        # Check target.type
        private$checkMode(target.type, "target.type", "character")
        target.type.list <- c("bmh", "rmh", "lmh", "ebmh", "ermh", "elmh", "protospacer")
        stopifnot(target.type %in% target.type.list)

        # Check should.convert.log
        private$checkMode(should.convert.log, "should.convert.log", "logical")
        
        # Check value.set
        stopifnot(is.list(extra.table.list))

        # Check should.remove.max.min
        private$checkMode(should.remove.max.min, "should.remove.max.min", "logical")

        # Check save.dir
        private$checkMode(save.dir, "save.dir", "character")

        # Check name
        private$checkMode(name, "name", "character")

        # Set information
        self$value.set <- value.set
        self$response.variable <- response.variable
        self$target.type <- target.type
        self$extra.table.list <- extra.table.list
        self$should.remove.max.min <- should.remove.max.min
        self$save.dir <- save.dir
        self$name <- name
        self$should.convert.log <- should.convert.log

        # Load class
        source(file.path(script.dir, "CalcMicroHomologyPropPK-class.R"))
        source(file.path(script.dir, "CalcSgRNAPropPK-class.R"))

        # Filtering
        # self$fltr.value.set <- filter(self$value.set, self$response.variable > 0) # If it uses ">= 0", log(0) can be observed.
        self$fltr.value.set <- self$value.set

        # Make label vector
        message("Run clustering using k-means...")
        label.vec <- self$MakeKmeansLabel(self$fltr.value.set[, self$response.variable])

        self$focus.length <- nchar(self$fltr.value.set[1, self$target.type])
        self$sgRNA.seqset <- paste(self$fltr.value.set[, "protospacer"], self$fltr.value.set[, "scaffold"], sep="")
        self$label.vec <- label.vec
        self$out.align.file <- file.path(self$save.dir, paste0(self$name, "_", self$target.type, "_", response.variable,  "_align.rds"))
        self$out.prop.file <- file.path(self$save.dir, paste0(self$name, "_", self$target.type, "_", response.variable, "_data.rds"))

        message("-----------------------------------------------------------------------------------")
        message(paste0("target sequence: ", self$target.type))
        message(paste0("response.variable: ", self$response.variable))
        message(paste0("sample name: ", name))
        message("Calcurate features")
        if(self$target.type %in% c("bmh", "rmh", "lmh", "ebmh", "ermh", "elmh")){

          self$microhomology.seqset <- self$fltr.value.set[, self$target.type]
          if(self$target.type %in% c("bmh", "lmh", "ebmh", "elmh")){
            self$is.left.both.microhomology.logical <- TRUE
          }else{
            self$is.left.both.microhomology.logical <- FALSE
          }
          my.CalcXPropPK <- self$CalcMHProp()

        }else if(self$target.type %in% c("protospacer")){

          # self$target.seqset <- self$fltr.value.set[, self$target.type]
          # self$is.protospacer.logical <- TRUE
          self$target.seqset <- self$fltr.value.set[, "fourtybase.region.surrounding.targetseq"]

          my.CalcXPropPK <- self$CalcTargetProp()

        }else{
          stop("Target type is undefined.", " : ", self$target.type, ".\n")
        }

        message("Merge the extra tables")
        for(list.ind in 1:length(extra.table.list)){
          if(list.ind == 1){
            self$full.extra.table <- extra.table.list[[list.ind]]
          }else{
            self$full.extra.table <- merge(self$full.extra.table, extra.table.list[[list.ind]], by = "group.name", all=TRUE)
          }
        }

        # continuous feature table
        message("Make the table of continuous values")
        self$full.continuous.table <- self$MergeExtraData(my.CalcXPropPK$ShowDataValue("continuous"), self$full.extra.table, self$target.type)
        write.csv(self$full.continuous.table, file(file.path(self$save.dir,"full.continuous.features.table.utf8.csv"), encoding="UTF-8"))
        
        # make name list of continuous feature
        message("Make table that has continuous values")
        self$continuous.feature.name <- colnames(self$full.continuous.table)[which(!(colnames(self$full.continuous.table) %in% c("sequence", self$target.type, "group.name", "label")))]

        # make sequecence - response table
        message("Make response table")
        self$response.value.by.seq.table <- data.frame(
          sequence = as.vector(merge(self$full.continuous.table, self$fltr.value.set, by = self$target.type, all=FALSE)[, self$target.type]),
          response.value = as.vector(merge(self$full.continuous.table, self$fltr.value.set, by = self$target.type, all=FALSE)[, self$response.variable]),
          stringsAsFactors=FALSE
        )

        # In this case, "max" and "min" parameter is removed because a diversity of these parameter is likely to be too low in the samples.#####
        if(self$should.remove.max.min){
          message("Remove variables that represents max/min value")
          small.continuous.feature.name.ind.vec <- !as.logical(
              table(
                factor(
                  c(grep("max", self$continuous.feature.name), grep("min", self$continuous.feature.name))
                , levels=1:length(self$continuous.feature.name)
                )
              )
            )
          self$continuous.feature.name <- self$continuous.feature.name[small.continuous.feature.name.ind.vec]
        }

        if(!self$can.convert.log){
          message(paste0("Make scattter plot between ", self$response.variable, " and each continuous feature"))
        }else{
          message(paste0("Make scattter plot between log(", self$response.variable, " + 1e-100) and each continuous feature"))
          # browser(expr = (self$response.variable == "Accuracy"))
        }
        self$PlotFeature()

        # Recognize base feature group name
        base.feature.vec.list <- list(
            SgRNA.Align = my.CalcXPropPK$colnames.align.vec
            , Freq = my.CalcXPropPK$colnames.freq.vec
            , Packer = my.CalcXPropPK$colnames.packer.vec
            , Phychem = my.CalcXPropPK$colnames.phychem.vec
            , Pseknc = my.CalcXPropPK$colnames.pseknc.vec
            , Thermo = my.CalcXPropPK$colnames.thermo.vec
        )
        
        # Recognize extra feature group name
        self$full.feature.vec.list <- base.feature.vec.list
        for(extra.table.ind in 1:length(extra.table.list)){
          temp.feature.name.vec <- colnames(extra.table.list[[extra.table.ind]])
          temp.feature.name.vec <- temp.feature.name.vec[!(temp.feature.name.vec %in% c("group.name", self$target.type))]
          if(names(extra.table.list)[extra.table.ind] == "Genome_Property" &
            any(colnames(self$full.extra.table) %in% "KO.Editing.Efficiency")){
            temp.feature.name.vec <- c(temp.feature.name.vec, "KO.Editing.Efficiency")
          }
          self$full.feature.vec.list <- c(self$full.feature.vec.list, list(temp.feature.name.vec))
          temp.name.vec <- names(self$full.feature.vec.list)[!names(self$full.feature.vec.list) %in% ""]
          names(self$full.feature.vec.list) <- c(temp.name.vec, names(extra.table.list)[extra.table.ind])
        }
        self$full.feature.vec.list <- c(self$full.feature.vec.list, list(All=unlist(self$full.feature.vec.list, use.names=FALSE)))
        # browser()

        # It makes the visibility of software decrease.
        # message("You can check a name of feature group...")
        # message("name of feature group -> feature name in the group")
        # for(feature.group.ind in 1:length(self$full.feature.vec.list)){
        #   message("----------------------------------------------------------------")
        #   message(paste0(names(self$full.feature.vec.list)[feature.group.ind], " -> ", self$full.feature.vec.list[feature.group.ind]))
        #   message("----------------------------------------------------------------")
        # }

      }
      , error = function(e) {
        message(e)
        browser()
        return(FALSE)
      }
      , silent = TRUE
      )
      return(TRUE)
    },
    # public methods
    MakeKmeansLabel = function(vec){
      # Make high/low cluster using k-means method
      #
      # Args:
      #   vec: vector of response variable
      #
      # Returns:
      #   On success, returns vector of labels including "high", "low" and "out".
      #   On failure, returns ERROR MESSAGE.
      # 
      # Cautions:
      #   Only k=2 is analilable in the recent version.
      tryCatch({
        num.vec <- as.numeric(vec)
        if(sd(num.vec) == 0){
          message("Poor diversity!!!")
          message("Add labels randomly. However, there is no use analying using these labels")
          num.vec.len <- length(num.vec)
          high.n <- round(num.vec.len / 2)
          return(c(rep("high", high.n), rep("low", num.vec.len - high.n)))
        }

        k <- 2
        # if(length(num.vec) > 10){
        #  max.k <- 10
        #}else{
        #  max.k <- length(num.vec) - 1
        #}

        while(TRUE){
          #should.stop.loop <- FALSE
          #if(k > max.k){
          #  message("Clustering is difficult. Set k=2...")
          #  k <- 2
          #  should.stop.loop <- TRUE
          #}else{
          #  message(paste0("Set k=", k, "..."))
          #}

          should.stop.loop <- TRUE # The other script cannot follow the cluster data with k > 2. We turn off the loop algorithm.
          message(paste0("Set k=", k, "..."))
          set.seed(1)
          # Make clusters
          kmeans.data <- kmeans(num.vec, k)
          # Recognize order of center value in the clusters
          center.vec <- kmeans.data$centers[,1]
          max.center.ind <- which.max(center.vec)
          second.center <- max(center.vec[center.vec != max(center.vec)])
          second.center.ind <- which(center.vec == second.center)
          if(k > 2){
            third.center <- max(center.vec[-c(max.center.ind, second.center.ind)])
            third.center.ind <- which(center.vec == third.center)
          }

          # Check whether the between_SS / total_SS is 80% or more.
          message(paste0("between_SS / total_SS = ", (kmeans.data$betweenss/kmeans.data$totss) * 100, "%"))
          # browser()
          if((kmeans.data$betweenss / kmeans.data$totss >= 0.8) | should.stop.loop){
            if(kmeans.data$betweenss / kmeans.data$totss >= 0.8){
              message("Good fit!")
            }else{
              message("The clustering is not good. You can use other dataset when you want to get clear data.")
            }
            # NOTE: this function can make large high group to merge two cluster.
            # This balance can be adjusted in the further updates (but, it is not TODO now).
            if(k == 2){
              high1.ind <- max.center.ind
              low.ind <- second.center.ind
              break
            }else{
              high1.ind <- max.center.ind
              high2.ind <- second.center.ind
              low.ind <- third.center.ind
              break
            }
          }else{
            k <- k + 1
            next
          }
        }

        # Attach labels
        label.vec <- replace(kmeans.data$cluster, kmeans.data$cluster == high1.ind, "high")
        if(k == 2){
          label.vec <- replace(label.vec, kmeans.data$cluster == low.ind, "low")
        }else{
          label.vec <- replace(label.vec, kmeans.data$cluster == high2.ind, "high")
          label.vec <- replace(label.vec, kmeans.data$cluster == low.ind, "low")
        }
        if(k > 3){
          out.ind <- (1:k)[!(1:k %in% c(high1.ind, high2.ind, low.ind))]
          label.vec <- replace(label.vec, kmeans.data$cluster %in% out.ind, "out")
        }
        return(label.vec)
      }, 
      error = function(e) { 
        browser()
        stop()
      }, silent = TRUE)
    },

    CalcMHProp = function(){
      # execute CalcMicroHomologyPropPK$new() with proper parameters
      #
      # Returns:
      #   instance of CalcMicroHomologyPropPK

      ### global variable
      # script.dir
      # ref.dir
      ###

      if(!file.exists(file.path(self$out.prop.file))){
          prop.info.CalcMicroHomologyPropPK <- CalcMicroHomologyPropPK$new(
              whole.region.microhomology.seq.vec = self$microhomology.seqset,
              is.left.both.microhomology.logical = self$is.left.both.microhomology.logical,
              scaffold.seq.vec = self$sgRNA.seqset,
              need.annotation.logical = TRUE,
              annotation.vec = self$label.vec,
              align.range.vec = seq(20, min(nchar(self$sgRNA.seqset)) , 20),
              align.save.file.char = file.path(self$out.align.file),
              script.dir = script.dir,
              ref.dir = ref.dir
          )
          saveRDS(prop.info.CalcMicroHomologyPropPK, file.path(self$out.prop.file))
      }else{
          prop.info.CalcMicroHomologyPropPK <- readRDS(file.path(self$out.prop.file))
      }
      return(prop.info.CalcMicroHomologyPropPK)
    },

    CalcTargetProp = function(){
      # execute CalcSgRNAPropPK$new() with proper parameters
      #
      # Returns:
      #   instance of CalcSgRNAPropPK

      ### global variable
      # script.dir
      # ref.dir
      ###

      if(!file.exists(file.path(self$out.prop.file))){
          # if(self$is.protospacer.logical){
          #   focusing.region.start.pos.num = 1
          #   focusing.region.end.pos.num = self$focus.length
          # }else{
            focusing.region.start.pos.num <- 11
            focusing.region.end.pos.num <- 30
          # }

          prop.info.CalcSgRNAPropPK <- CalcSgRNAPropPK$new(
              fourtybase.region.surrounding.targetseq.vec = self$target.seqset,
              protospacer.start.pos.num = focusing.region.start.pos.num,
              protospacer.end.pos.num = focusing.region.end.pos.num,
              # is.left.both.microhomology.logical = self$is.left.both.microhomology.logical,
              scaffold.seq.vec = self$sgRNA.seqset,
              need.annotation.logical = TRUE,
              annotation.vec = self$label.vec,
              align.range.vec = seq(20, nchar(self$sgRNA.seqset[1]) , 20),
              align.save.file.char = file.path(self$out.align.file),
              script.dir = script.dir,
              ref.dir = ref.dir
          )
          saveRDS(prop.info.CalcSgRNAPropPK, file.path(self$out.prop.file))
      }else{
          prop.info.CalcSgRNAPropPK <- readRDS(file.path(self$out.prop.file))
      }
      return(prop.info.CalcSgRNAPropPK)
    },

    MergeExtraData = function(continuous.set, extra.set, target.type){
      # Merge extra dataset with value dataset
      #
      # Args:
      #   continuous.set: a part of value table: it includes all continuous features {data.frame}
      #   extra.set : extra table that user added {data.frame}
      #   target.type: the focused sequence target
      #   "bmh" mode that calculate values of both homology arms
      #   "lmh" mode that calculate values of Left homology arm
      #   "rmh" mode that calculate values of Right homology arm
      #   "protospacer" mode that calculate ones of target site.
      #
      # Returns:
      #   The merged table.
      temp.df <- data.frame("sequence" = as.vector(extra.set[, target.type]), stringsAsFactors = FALSE)
      for(name in colnames(extra.set)){
        # NOTE: The old version has filtering function that add features developer want to analize into the value dataset.
        #       This was unflexible. This version adds all features in extra table into the value dataset.
        temp.df <- cbind(temp.df, extra.set[,name])
        colnames(temp.df)[length(temp.df)] <- name
      }

      sum.df <- merge(continuous.set[, 1:length(continuous.set) - 1], temp.df, by ="sequence", all=T)
      sum.df <- merge(sum.df, continuous.set[,c("sequence", "label")], by ="sequence", all=T)
      if(any(is.na(sum.df$label))){
        return(sum.df[-which(is.na(sum.df$label)), ])
      }else{
        return(sum.df)
      }
    },        

    PlotFeature = function(){
      # Plot Scatter between reposonse and continuous features and Box plot in three groups (High, Medium, Low)
      #
      # Returns:
      #   x, y value. if this cannot be transformed into log value in the case, return NA
        MakePlotValue = function(scat.data.df, feature.name, response.name, can.normalization){
          response.name <- "response.value"

          x.vec <- as.numeric(scat.data.df[, feature.name])

          # Convert value into log value (Optional)
          if(self$should.convert.log){
            if(min(x.vec[complete.cases(x.vec)]) > 0){ # MEMO: NA/Inf value will be removed. In the time, they are ignored.
              x.vec <- log(x.vec)
              self$can.convert.log <- TRUE
            }else{
              message(paste0(feature.name, " cannot be calculated with log-transformed values."))
              return(NA)
            }
          }

          # normalization
          x.min <- min(x.vec[complete.cases(x.vec)])
          x.max <- max(x.vec[complete.cases(x.vec)])
          y.vec <- as.numeric(scat.data.df[, response.name])
          y.min <- min(y.vec[complete.cases(y.vec)])
          y.max <- max(y.vec[complete.cases(y.vec)])

          if(can.normalization){ # turn on normalization
            if(all(is.na(scale(x.vec)[, 1]))){ # if all values are same or NA
              plot.value.x <- rep(0, length(x.vec))
            }else{
              plot.value.x <- scale(x.vec, center = x.min, scale = (x.max - x.min))[, 1]
            }

            if(all(is.na(scale(y.vec)[, 1]))){ # if all values are same or NA
              plot.value.y <- rep(0, length(y.vec))
            }else{
              plot.value.y <- scale(y.vec, center = y.min, scale = (y.max - y.min))[, 1]
            }
            norm.label <- "_norm"
          }else{
            plot.value.x <- x.vec
            plot.value.y <-  y.vec
            norm.label <- ""
          }

          return(list(plot.value.x = plot.value.x, plot.value.y = plot.value.y, norm.label = norm.label))
        }
        PlotScattter = function(feature.name, order, can.normalization){
          tryCatch({
            scat.data.df <- merge(self$response.value.by.seq.table, self$full.continuous.table, by="sequence", all=FALSE)

            temp.list <- MakePlotValue(scat.data.df, feature.name, response.name, can.normalization)
            if(length(temp.list) == 1){
              if(is.na(temp.list)){ # make double nest to avoid warning
                return()
              }
            }
            plot.value.x <- temp.list$plot.value.x
            plot.value.y <- temp.list$plot.value.y
            norm.label <- temp.list$norm.label

            plot.value.df <- data.frame(
                x = plot.value.x,
                y = plot.value.y,
                stringsAsFactors = FALSE)
            rownames(plot.value.df) <- scat.data.df$group.name
            plot.value.df <- plot.value.df[complete.cases(plot.value.df), ] # remove row contains NAs
            plot.value.df <- plot.value.df[!is.infinite(rowSums(plot.value.df)),] # remove row contains Inf or Inf
            if((length(plot.value.df$x) < 3) | (length(plot.value.df$y) < 3)){# not enough finite observations
              return(0)
            }
            if(sd(plot.value.df$x) != 0 & sd(plot.value.df$y) != 0){
              cor.res <- cor(plot.value.df$x, plot.value.df$y)
            }else{
              cor.res <- NA
            }
            tryCatch({
              if(sd(plot.value.df$x) != 0 & sd(plot.value.df$y) != 0){
                cor.test.res <- cor.test(plot.value.df$x, plot.value.df$y, method="pearson")
              }
            }, 
            error = function(e) {
              cor.test.res <- NA
            },
            silent = TRUE
            )

            # Make plot
            options(warn=-1)
            res.plot <- ggplot(plot.value.df, aes(x=x,y=y)) + 
            geom_point()+
            geom_smooth(method=lm, se=FALSE, formula = y ~ x)+
            labs(x = feature.name, y = self$response.variable, title = cor.res)
            options(warn=0)
            if(is.na(cor.res)){
              cor.label <- "NA"
            }else if(abs(cor.res) > 0.4){
              print(paste(order, ".", feature.name, " vs ",self$response.variable, "=", cor.res, sep = ""))
              cor.label <- paste0("_R=", cor.res)
            }else{
              cor.label <- ""
            }

            # Save
            modified.feature.name <- gsub("%", "rate", feature.name)
            print(paste0(self$response.variable, "_feature", order, "_", modified.feature.name, cor.label))

            ggsave(file.path(self$save.dir, paste(self$response.variable, "_feature", order, "_", modified.feature.name, norm.label, cor.label, ".png", sep = "")), plot = res.plot, device = NULL, path = NULL,
            scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
            dpi = 150, limitsize = TRUE)

            saved.data.table <- plot.value.df
            colnames(saved.data.table) <- c(modified.feature.name, self$response.variable)
            
            write.csv(saved.data.table,
              file.path(self$save.dir,
                paste(self$response.variable, "_feature", order, "_", modified.feature.name, norm.label, cor.label, "_data.csv", sep = ""))
            )
            write.csv(data.frame(Cor = cor.res, pvalue = cor.test.res$p.value),
              file.path(self$save.dir,
                paste(self$response.variable, "_feature", order, "_", modified.feature.name, norm.label, cor.label, "_statics.csv", sep = ""))
            )

          }, 
          error = function(e) {
            browser()
          },
          silent = TRUE)
        }
        PlotBox = function(feature.name, order, can.normalization){
          tryCatch({

            scat.data.df <- merge(self$response.value.by.seq.table, self$full.continuous.table, by="sequence", all=FALSE)

            temp.list <- MakePlotValue(scat.data.df, feature.name, response.name, can.normalization)
            if(length(temp.list) == 1){
              if(is.na(temp.list)){ # make double nest to avoid warning
                return()
              }
            }
            feature.value.vec <- temp.list$plot.value.x
            response.value.vec <- temp.list$plot.value.y
            norm.label <- temp.list$norm.label
            plot.value.df <- data.frame(
                feature = feature.value.vec,
                response = response.value.vec,
                stringsAsFactors = FALSE)
            rownames(plot.value.df) <- scat.data.df$group.name
            plot.value.df <- plot.value.df[complete.cases(plot.value.df), ] # remove row contains NAs
            plot.value.df <- plot.value.df[!is.infinite(rowSums(plot.value.df)),] # remove row contains Inf or Inf
            if((length(plot.value.df$feature) < 3) | (length(plot.value.df$response) < 3)){# not enough finite observations
              return(0)
            }

            # make rank
            response.rank.vec <- rank(-plot.value.df$response)

            p.value.label <- paste("")
            modified.feature.name <- gsub("%", "rate", feature.name)
            if(length(unique(response.rank.vec)) > 2){
              group.n <- round(length(response.rank.vec) / 3)
              high.vec <- response.rank.vec <= group.n
              low.vec <- response.rank.vec > length(response.rank.vec) - group.n
              response.group.vec <- response.rank.vec
              response.group.vec[high.vec] <- "High"
              response.group.vec[low.vec] <- "Low"
              response.group.vec[!(high.vec | low.vec)] <- "Medium"
              feature.resgroup.df <- data.frame(sample = rownames(plot.value.df), value = plot.value.df$feature, group = factor(response.group.vec, levels=c("Low", "Medium", "High")))
            
              # make box plot
              # browser(expr = (feature.name == "TT"))
              box.p <- ggplot(feature.resgroup.df, aes(x=group, y=value)) +
                geom_boxplot(fill = c("#0487D9", "#04BF45", "#F21905"), color = "black", alpha = 0.5, show.legend = FALSE, outlier.shape = NA) + 
                geom_jitter(position=position_jitter(width=0.3, height=0), alpha = 1) +
                ggtitle(paste0(feature.name, " per group of ", self$response.variable)) +
                theme(axis.text.x = element_text(size = 10)
                  , axis.text.y = element_text(size = 10)
                  , axis.title.x = element_text(size = 12)
                  , axis.title.y = element_text(size = 12)
                )

              # t-test
              if(sd(plot.value.df$feature) != 0 & sum(high.vec, na.rm = TRUE) > 1 & sum(low.vec, na.rm = TRUE) > 1){

                t.test.res <- t.test(plot.value.df$feature[high.vec], plot.value.df$feature[low.vec])
                print(paste0(feature.name, " among groups of ", self$response.variable, ": the p-value for the difference between the High group and the Low group is ", t.test.res$p.value))
                if(t.test.res$p.value < 0.05){
                  p.value.label <- paste0("_pvalue=", as.character(t.test.res$p.value))
                }
                # save statistics
                n1 <- length(plot.value.df$feature[high.vec])
                n2 <- length(plot.value.df$feature[low.vec])
                sd1 <- sd(plot.value.df$feature[high.vec])
                sd2 <- sd(plot.value.df$feature[low.vec])
                mean1 <- mean(plot.value.df$feature[high.vec])
                mean2 <- mean(plot.value.df$feature[low.vec])
                cohend <- cohen.d(plot.value.df$feature[high.vec], plot.value.df$feature[low.vec])
                write.csv(data.frame(samplesize1 = n1,
                    samplesize2 = n2,
                    sd1 = sd1,
                    sd2 = sd2,
                    mean1 = mean1,
                    mean2 = mean2,
                    t = t.test.res$statistic,
                    df = t.test.res$parameter,
                    pvalue = t.test.res$p.value,
                    effectsize = cohend$estimate,
                    magnitude = cohend$magnitude
                  ),
                  file.path(file.path(self$save.dir, paste0(self$response.variable, "_feature", order, "_", modified.feature.name, norm.label, p.value.label, "_boxplot_statics.csv")))
                )
    
              }else{

                print(paste0(modified.feature.name, " among groups of ", self$response.variable, ": the p-value for the difference between the High group and the Low group is NA"))

              }

            }else{
              
              message("There is a lot of ties in ", self$response.variable, ". Make one box plot using all values.")
              feature.resgroup.df <- data.frame(sample = rownames(plot.value.df), value = plot.value.df$feature, group = factor(c("All"), levels=c("All")))
              box.p <- ggplot(feature.resgroup.df, aes(x=group, y=value)) +
                geom_boxplot(fill = c("grey"), color = "black", alpha = 0.5, show.legend = FALSE, outlier.shape = NA) + 
                geom_jitter(width = 0.25, size = 1, alpha = 1) +
                ggtitle(paste0(feature.name, " per group of ", self$response.variable)) +
                theme(axis.text.x = element_text(size = 10)
                  , axis.text.y = element_text(size = 10)
                  , axis.title.x = element_text(size = 12)
                  , axis.title.y = element_text(size = 12)
                )

            }

            # save plot
            ggsave(file.path(self$save.dir, paste0(self$response.variable, "_feature", order, "_", modified.feature.name, norm.label, p.value.label, "_boxplot.png")),
                plot = box.p, device = NULL, path = NULL,
                scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                dpi = 150, limitsize = TRUE)
            write.csv(data.frame(sample = feature.resgroup.df$sample, value = feature.resgroup.df$value, group = feature.resgroup.df$group, response = plot.value.df$response),
              paste(file.path(self$save.dir, paste0(self$response.variable, "_feature", order, "_", modified.feature.name, norm.label, p.value.label, "_boxplot_data.csv")))
            )
            
          }, 
          error = function(e) {
            browser()
          },
          silent = TRUE)
        }
        for(order in 1:length(self$continuous.feature.name)){
            feature.name <- self$continuous.feature.name[1:length(self$continuous.feature.name)][order]
            PlotScattter(feature.name, order, can.normalization = FALSE)
            PlotBox(feature.name, order, can.normalization = FALSE)
            PlotScattter(feature.name, order, can.normalization = TRUE)
            PlotBox(feature.name, order, can.normalization = TRUE)
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
