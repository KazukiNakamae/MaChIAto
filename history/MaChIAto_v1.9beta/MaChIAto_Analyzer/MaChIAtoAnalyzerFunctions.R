GetAbsolutePath = function(relative.path) {
  # Get an absolute path from a given relative.path
  #
  # Args:
  #   relative.path: the charators of relative path
  #
  # Returns:
  #   the charators of absolute path
  return(normalizePath(relative.path))
}

PCbiplot = function(PC, x="PC1", y="PC2") {
    # PC being a prcomp object
    data <- data.frame(obsnames=row.names(PC$x), PC$x)
    plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
    plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
    datapc <- transform(datapc,
            v1 = .7 * mult * (get(x)),
            v2 = .7 * mult * (get(y))
            )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
    plot
}

SubtractCtrlDF = function(df){

    SubtractCtrl = function(vec){
        options(warn=-1) # is.na(as.numeric(vec[1])) put warning(). It is, however, no problem.
        if(is.na(as.numeric(vec[1]))){
            return(vec)
        }else{
            vec <- as.numeric(vec)
            return(vec - vec[1])
        }
        options(warn=0)
    }

    return(apply(df, MARGIN=2, SubtractCtrl))

}

# MakeValueSet = function(df, sample.type){
#    
#     for(name in names(df)){
#         if(!exists("sum.df")){
#             sum.df <- df[[name]][df[[name]][,"sample.type"] == sample.type,]
#         }else{
#             sum.df <- rbind(
#                 sum.df,
#                 df[[name]][df[[name]][,"sample.type"] == sample.type,]
#             )
#         }
#     }
#     return(sum.df)
# }

MakeDummyValueSet = function(value.set){
  dummy.value.set <- value.set
  for(colname in colnames(value.set)){
      if(is.numeric(value.set[, colname])){
          dummy.value.set[, colname] <- rep(0, length(value.set[, colname]))
      }else if(colname == "sample.type"){
          dummy.value.set[, colname] <- rep("untreated", length(value.set[, colname]))
      }else if(colname == "sample.name"){
          dummy.value.set[, colname] <- rep("machiato_dummy_sample", length(value.set[, colname]))
      }else{
          dummy.value.set[, colname] <- value.set[, colname]
      }
  }
  return(dummy.value.set)
}

CreateManageAnalysisData = function(
  value.set,
  response.variable,
  target.type,
  should.convert.log,
  extra.table.list,
  save.dir,
  name
  ){
  if(nrow(value.set) < 100){ # too low sample numbers # POINT!!!
      should.remove.max.min <- TRUE
  }else{
      should.remove.max.min <- FALSE
  }
  return(
    ManageAnalysisData$new(
      value.set = value.set, 
      response.variable = response.variable,
      target.type = target.type,
      should.convert.log = should.convert.log,
      extra.table.list = extra.table.list,
      should.remove.max.min = should.remove.max.min,
      save.dir = save.dir,
      name = name
    )
  )
}

FltrByQuantile = function(df, indicator){
  return(df[(as.numeric(df[[indicator]]) > quantile(as.numeric(df[[indicator]]), 0.00)), ])
}

MakeMHScoreTable = function(fltr.extra.table, target.type, script.dir, output.dir){
  bmh.vec <- fltr.extra.table[, "bmh"]
  cor.rmh.vec <- fltr.extra.table[, "rmh"]
  cor.lmh.vec <- fltr.extra.table[, "lmh"]
  rev.rmh.vec <- as.character(reverseComplement(DNAStringSet(fltr.extra.table[, "rmh"])))
  rev.lmh.vec <- as.character(reverseComplement(DNAStringSet(fltr.extra.table[, "lmh"])))
  ebmh.vec <- fltr.extra.table[, "ebmh"]
  cor.ermh.vec <- fltr.extra.table[, "ermh"]
  cor.elmh.vec <- fltr.extra.table[, "elmh"]
  rev.ermh.vec <- as.character(reverseComplement(DNAStringSet(fltr.extra.table[, "ermh"])))
  rev.elmh.vec <- as.character(reverseComplement(DNAStringSet(fltr.extra.table[, "elmh"])))
  message("Microhomology score between the left homology sequence and the right homology sequence")
  cc.mh.score.CalcMHScore <- CalcMHScore$new(fltr.extra.table[, "bmh"]
      , cutpos.vec = nchar(fltr.extra.table[, "lmh"])
      , length_weight = 20
      , script.dir = script.dir
      , save.dir = output.dir
  )
  cc.mh.score.CalcMHScore.summary <- as.data.frame(cc.mh.score.CalcMHScore$summary)
  colnames(cc.mh.score.CalcMHScore.summary) <- c("bmh", "Endogenous.MP.Score", "Endogenous.MPOF.Score")
  message("Microhomology score between the left homology sequence and the right donor homology RC sequence")
  cr.mh.score.CalcMHScore <- CalcMHScore$new(paste0(cor.lmh.vec, rev.ermh.vec)
      , cutpos.vec = nchar(cor.rmh.vec)
      , length_weight = 20
      , script.dir = script.dir
      , save.dir = output.dir
  )
  cr.mh.score.CalcMHScore.summary <- as.data.frame(cr.mh.score.CalcMHScore$summary)
  colnames(cr.mh.score.CalcMHScore.summary) <- c("seq", "LeftRevKI.MP.Score", "LeftRevKI.MPOF.Score")
  message("Microhomology score between the left donor homology RC sequence and the right homology sequence")
  rc.mh.score.CalcMHScore <- CalcMHScore$new(paste0(rev.elmh.vec, cor.rmh.vec)
      , cutpos.vec = nchar(rev.elmh.vec)
      , length_weight = 20
      , script.dir = script.dir
      , save.dir = output.dir
  )
  rc.mh.score.CalcMHScore.summary <- as.data.frame(rc.mh.score.CalcMHScore$summary)
  colnames(rc.mh.score.CalcMHScore.summary) <- c("seq", "RightRevKI.MP.Score", "RightRevKI.MPOF.Score")
  message("Microhomology score between the left homology sequence and then left donor homology sequence")
  bl.mh.score.CalcMHScore <- CalcMHScore$new(paste0(cor.lmh.vec, cor.elmh.vec)
      , cutpos.vec = nchar(cor.lmh.vec)
      , length_weight = 20
      , script.dir = script.dir
      , save.dir = output.dir
  )
  bl.mh.score.CalcMHScore.summary <- as.data.frame(bl.mh.score.CalcMHScore$summary)
  colnames(bl.mh.score.CalcMHScore.summary) <- c("seq", "LeftKI.MP.Score", "LeftKI.MPOF.Score")
  message("Microhomology score between the right donor homology sequence and the right homology sequence")
  br.mh.score.CalcMHScore <- CalcMHScore$new(paste0(cor.ermh.vec, cor.rmh.vec)
      , cutpos.vec = nchar(cor.ermh.vec)
      , length_weight = 20
      , script.dir = script.dir
      , save.dir = output.dir
  )
  br.mh.score.CalcMHScore.summary <- as.data.frame(br.mh.score.CalcMHScore$summary)
  colnames(br.mh.score.CalcMHScore.summary) <- c("seq", "RightKI.MP.Score", "RightKI.MPOF.Score")
  bind.table <- cbind(cc.mh.score.CalcMHScore.summary, cr.mh.score.CalcMHScore.summary, rc.mh.score.CalcMHScore.summary, bl.mh.score.CalcMHScore.summary, br.mh.score.CalcMHScore.summary)
  bind.table <- bind.table[, which(!(colnames(bind.table) %in% "seq"))]
  if(target.type == "bmh"){
    bind.table <- cbind(data.frame(group.name=fltr.extra.table$group.name), bind.table)
  }else{
    pre.bind.table <- cbind(fltr.extra.table[, c(target.type, "group.name")], bind.table)
    bind.table <- subset(pre.bind.table, select = -c(bmh))
  }

  return(list(bind.table))
}

MakeKmeansAnnotationSetDF = function(vec){
  tryCatch({
    num.vec <- as.numeric(vec)
    if(sd(num.vec) == 0){
      return(character(length(num.vec)))
    }
    k <- 2
    while(TRUE){
      set.seed(1)
      kmeans.data <- kmeans(num.vec, k)
      # identify size index
      size.vec <- kmeans.data$size
      max.size.ind <- which.max(size.vec)
      if(length(size.vec[size.vec!=max(size.vec)]) == 0){ # size is equal
        second.size.ind <- which(size.vec == size.vec[max.size.ind])[-max.size.ind]
      }else{ # size is not equal
        second.size <- max(size.vec[size.vec!=max(size.vec)])
        second.size.ind <- which(size.vec == second.size)
      }
      if(size.vec[second.size.ind] / size.vec[max.size.ind] > 0.25 &
        size.vec[second.size.ind] >= 5){
        if(kmeans.data$centers[,1][max.size.ind] > kmeans.data$centers[,1][second.size.ind]){
          high.ind <- max.size.ind
          low.ind <- second.size.ind
          break
        }else if(kmeans.data$centers[,1][max.size.ind] < kmeans.data$centers[,1][second.size.ind]){
          high.ind <- second.size.ind
          low.ind <- max.size.ind
          break
        }else(
          next
        )
      }
      k <- k + 1
      if(k > 10){
        return(rep(NA , length(num.vec)))
      }
    }
    annotation.vec <- replace(kmeans.data$cluster, kmeans.data$cluster==high.ind, "high")
    annotation.vec <- replace(annotation.vec, kmeans.data$cluster==low.ind, "low")
    if(k > 2){
      out.ind <- (1:k)[!(1:k %in% c(low.ind, high.ind))]
      annotation.vec <- replace(annotation.vec, kmeans.data$cluster %in% out.ind, "out")
    }
    return(annotation.vec)
  }, 
  error = function(e) { 
    return(vec)
  }, silent = TRUE)
}

Removeoutliers = function(value.set, response){
  outliers <- boxplot(value.set[, response], plot=FALSE)$out
  if(length(outliers)){
    value.set <- value.set[-which(value.set[, response] %in% outliers),]
  }
  return(value.set)
}

PlotBarPlot = function(vec, type, title){
  temp.df <- data.frame(feature.name = names(vec), value = vec)
  if(type == "mono.grad"){
    temp.df %>% 
      ggplot(aes(reorder(feature.name, value), value)) + 
      geom_col(aes(fill = value)) + 
      scale_fill_gradient(low = "gray", 
                          high = "black") + 
      coord_flip() + 
      labs(x = "Value", y = "Future", title = title) +
      theme(axis.text = element_text(size = 8, colour="black"))  -> p
  }else if(type == "black"){
    temp.df %>% 
      ggplot(aes(reorder(feature.name, value), value)) + 
      geom_col() +
      coord_flip() + 
      labs(x = "Future") +
      theme_classic() +
      theme(axis.text = element_text(size = 8, colour="black")) -> p
  }
  return(p)
}

# Multiple plot function (http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

MakeCorPlot = function(x, y, x.name, y.name, save.dir){
  cor <- cor(x, y)
  df <- data.frame(x=x, y=y, stringsAsFactors = FALSE)
  # Make plot
  res.plot <- ggplot(df, aes(x=x,y=y)) + 
    geom_point()+
    geom_smooth(method=lm, se=FALSE)+
    labs(x = x.name, y = y.name, title = cor)
  ggsave(file.path(save.dir, paste(x.name, "_vs_", y.name, ".png", sep = "")), plot = res.plot, device = NULL, path = NULL,
    scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
    dpi = 300, limitsize = TRUE)
  return(res.plot)
}

MakeCorCompPlot = function(ManageAnalysisData.1, ManageAnalysisData.2, feature.group.ind, save.dir, response.name, sample.name.vec, should.convert.log){
  # generate boxplots regarding correlation value for comparison
  #
  # Args:
  #   ManageAnalysisData.1: data contains feature value {ManageAnalysisData}
  #   ManageAnalysisData.2: data contains feature value to be compared {ManageAnalysisData}
  #   feature.group.ind: index of feature group {number}
  #   save.dir: saving directory {character}
  #   response.name: name of response value {character}
  #   sample.name.vec: names of ManageAnalysisData.1 and ManageAnalysisData.2 {vector}
  #   should.convert.log
  #
  # Returns:
  #   On success, returns TRUE.
  #   On failure, returns NA.

  tryCatch({
    feature.data.df1 <- ManageAnalysisData.1$full.continuous.table
    feature.data.df2 <- ManageAnalysisData.2$full.continuous.table
    response.continuous.df1 <- ManageAnalysisData.1$response.value.by.seq.table
    response.continuous.df2 <- ManageAnalysisData.2$response.value.by.seq.table
    feature.name.vec = ManageAnalysisData.1$full.feature.vec.list[[feature.group.ind]]
    feature.set.name = names(ManageAnalysisData.1$full.feature.vec.list)[feature.group.ind]

    CalcCor = function(feature.data.df, response.continuous.df, feature.name, should.convert.log){
      scat.data.df <- merge(feature.data.df, response.continuous.df, by="sequence", all=F)
      response.name <- colnames(scat.data.df)[length(scat.data.df)]

      x.vec <- as.numeric(scat.data.df[, feature.name])

      # Convert value into log value (Optional)
      if(should.convert.log){
        
        if(min(x.vec[complete.cases(x.vec)]) > 0){ # MEMO: NA/Inf value will be removed. In the time, they are ignored.
          x.vec <- log(x.vec)
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

      if(all(is.na(scale(x.vec))[, 1])){ # if all values are same or NA
        plot.value.x <- rep(0, length(x.vec))
      }else{
        plot.value.x <- scale(x.vec, center = x.min, scale = (x.max - x.min))[, 1]
      }

      if(all(is.na(scale(y.vec))[, 1])){ # if all values are same or NA
        plot.value.y <- rep(0, length(y.vec))
      }else{
        plot.value.y <- scale(y.vec, center = y.min, scale = (y.max - y.min))[, 1]
      }

      # Calc cor
      plot.value.df <- data.frame(
          x = plot.value.x,
          y = plot.value.y,
          stringsAsFactors = FALSE)
      plot.value.df <- plot.value.df[complete.cases(plot.value.df), ] # remove row contains NAs # Point!!!
      plot.value.df <- plot.value.df[!is.infinite(rowSums(plot.value.df)),] # remove row contains Inf or Inf # Point!!!
      if((length(plot.value.df$x) < 3) | (length(plot.value.df$y) < 3)){# not enough finite observations
        return(NA)
      }
      if(sd(plot.value.df$x) != 0 & sd(plot.value.df$y) != 0){
        cor.res <- cor(plot.value.df$x, plot.value.df$y)
      }else{
        cor.res <- NA
      }

      return(cor.res)
    }
    
    # merge cor data
    cor1.vec <- numeric(0)
    cor2.vec <- numeric(0)
    for(feature.name in feature.name.vec){

      temp1 <- c(CalcCor(feature.data.df1, response.continuous.df1, feature.name, should.convert.log))
      if(length(temp1) == 1){
        if(is.na(temp1)){ # make double nest to avoid warning
          next
        }
      }

      temp2 <- c(CalcCor(feature.data.df2, response.continuous.df2, feature.name, should.convert.log))
      if(length(temp2) == 1){
        if(is.na(temp2)){ # make double nest to avoid warning
          next
        }
      }

      names(temp1) <- feature.name
      cor1.vec <- c(cor1.vec, temp1)
      names(temp2) <- feature.name
      cor2.vec <- c(cor2.vec, temp2)

    }
    if( (length(cor1.vec) > 2) & (length(cor2.vec) > 2) ){

      # prepare plot
      options(warn=-1)
      t.test.res <- t.test(x = cor1.vec, y = cor2.vec, paired = TRUE)
      ks.test.res <- ks.test(x = cor1.vec, y = cor2.vec)
      options(warn=0)
      plot.cor.df <- data.frame(
          cor = c(cor1.vec, cor2.vec),
          group = as.factor(c(rep(sample.name.vec[1], length(cor1.vec)), rep(sample.name.vec[2], length(cor2.vec)))),
          feature.name = c(names(cor1.vec), names(cor2.vec)),
          stringsAsFactors = FALSE)
      plot.cor.df <- plot.cor.df[complete.cases(plot.cor.df), ]
      groupall.n <- nrow(plot.cor.df)
      group1.n <- nrow(plot.cor.df[plot.cor.df$group == sample.name.vec[1],])
      group2.n <- groupall.n - group1.n
      is_outlier = function(x) {
        return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
      }
      outliner.vec <- c(
        rep(TRUE, nrow(plot.cor.df[plot.cor.df$group == sample.name.vec[1],])),
        !is_outlier(plot.cor.df[plot.cor.df$group == sample.name.vec[2],]$cor)
      )
      plot.cor.df$feature.name[outliner.vec] <- ""

      p1 <- ggplot(plot.cor.df, aes(x=group, y=cor)) +
        geom_boxplot(fill = c("#D96A6A", "#F2DA91"), color = "black", alpha = 1, show.legend = FALSE, outlier.shape = NA) + 
        geom_jitter(width = 0.25, aes(colour=group), size = 0.5, alpha = 1) +
        ggtitle(paste0(response.name, "_", feature.set.name)) +
        geom_text(aes(label = feature.name), size = 3) +
        scale_color_manual(values = c("#8C3A32", "#F2A88D"))+
        ylim(-1, 1)
      ggsave(file.path(save.dir, paste(response.name, "_", feature.set.name, "_Tp.value=", t.test.res$p.value, "_KSp.value="
          , ks.test.res$p.value, "_n1=", group1.n, "_n2=", group2.n, "_labeled.png", sep = ""))
          , plot = p1, device = NULL, path = NULL,
        scale = 1, width = 3, height = 6, units = c("in", "cm", "mm"),
        dpi = 300, limitsize = TRUE)
      p2 <- ggplot(plot.cor.df, aes(x=group, y=cor)) +
        geom_boxplot(fill = c("#D96A6A", "#F2DA91"), color = "black", alpha = 1, show.legend = FALSE, outlier.shape = NA) + 
        geom_jitter(width = 0.25, aes(colour=group), size = 0.5, alpha = 1) +
        ggtitle(paste0(response.name, "_", feature.set.name)) +
        scale_color_manual(values = c("#8C3A32", "#F2A88D"))+
        ylim(-1, 1)
      ggsave(file.path(save.dir, paste(response.name, "_", feature.set.name, "_Tp.value=", t.test.res$p.value, "_KSp.value="
          , ks.test.res$p.value, "_n1=", group1.n, "_n2=", group2.n, ".png", sep = ""))
          , plot = p2, device = NULL, path = NULL,
        scale = 1, width = 3, height = 6, units = c("in", "cm", "mm"),
        dpi = 150, limitsize = TRUE)
      message(paste0(
          "Saved as ",
          file.path(save.dir,
            paste0(response.name, "_", feature.set.name, "_Tp.value=", t.test.res$p.value, "_KSp.value="
            , ks.test.res$p.value, "_n1=", group1.n, "_n2=", group2.n, ".png")
          )
        ))

    }else{

      message("Insufficient amount of data.")

    }
  }
  , error = function(e) {
    browser()
    message(e)
    message("Could not generate boxplot. The system skiped it")
    # browser()
    return(FALSE)
  }
  , silent = TRUE
  )
  return(TRUE)
}

CompDoubleData = function(ManageAnalysisData.1, ManageAnalysisData.2, value.name.1, value.name.2, save.dir, response.name, should.convert.log){
  # perform comparison between two data set
  #
  # Args:
  #   ManageAnalysisData.1: data contains feature value {ManageAnalysisData}
  #   ManageAnalysisData.2: data contains feature value to be compared {ManageAnalysisData}
  #   value.name.1: name of value.set {character}
  #   value.name.2: name of value.set {character}
  #   save.dir: saving directory {character}
  #   response.name: name of response value {character}
  #
  # Returns:
  #   On success, returns TRUE.
  #   On failure, returns FALSE.
  tryCatch({
      message(paste0("Perform comparison between ", value.name.1, " and ", value.name.2))
      for(feature.group.ind in 1:length(ManageAnalysisData.1$full.feature.vec.list)){
        message("-------------------------------------------------------")
        message(paste0("Feature group: ", names(ManageAnalysisData.1$full.feature.vec.list)[feature.group.ind]))
        # browser()
        MakeCorCompPlot(ManageAnalysisData.1 = ManageAnalysisData.1,
          ManageAnalysisData.2 = ManageAnalysisData.2,
          feature.group.ind = feature.group.ind,
          save.dir = save.dir,
          response.name = response.name,
          sample.name.vec = c(value.name.1, value.name.2),
          should.convert.log = should.convert.log
          )
      }
    }
  , error = function(e) {
    browser()
    message(e)
    return(FALSE)
  }
  , silent = TRUE
  )
  return(TRUE)
}

CopyFile = function(file.1, file.2){
  # copy file.1 as file.2
  #
  # Args:
  #   file.1: copied file {character}
  #   file.2: file to be generated {character}
  system(paste("cp", file.1, file.2, sep=" "))
}

ExtractResopnseName = function(value.set){
  # extract response name list from value.set 
  #
  # Args:
  #   value.set: table contains response values and sequence infromation {data.frame}
  response.name.vec <- colnames(value.set)
  response.name.vec <- response.name.vec[!(response.name.vec %in% c(
    "name", "lmh", "rmh", "bmh",
    "elmh", "ermh", "ebmh",
    "protospacer", "scaffold", "fourtybase.region.surrounding.targetseq",
    "sample.name", "group.name", "sample.type"
    ))]
    return(response.name.vec)
}

#PlotTopRankScattter = function(feature.data.df, response.continuous.df, rank.vec, min.rank, sample.type){
#    PlotScattter = function(feature.data.df, response.continuous.df, top.feature.name, rank){
#      tryCatch({# urgent deal
#            scat.data.df <- merge(feature.data.df, response.continuous.df, by="sequence", all=F)
#            response.name <- colnames(scat.data.df)[length(scat.data.df)]
#
#            if(all(is.na(scale(as.numeric(scat.data.df[, top.feature.name])))[, 1])){ # if all values are same or NA
#              plot.value.x <- rep(0, length(as.numeric(scat.data.df[, top.feature.name])))
#            }else{
#              plot.value.x <- scale(as.numeric(scat.data.df[, top.feature.name]))[, 1]
#            }
#
#            if(all(is.na(scale(as.numeric(scat.data.df[, response.name])))[, 1])){ # if all values are same or NA
#              plot.value.y <- rep(0, length(as.numeric(scat.data.df[, response.name])))
#            }else{
#              plot.value.y <- scale(as.numeric(scat.data.df[, response.name]))[, 1]
#            }
#
#            plot.value.df <- data.frame(
#                x = plot.value.x,
#                y = plot.value.y,
#                stringsAsFactors = FALSE)
#            plot.value.df <- plot.value.df[complete.cases(plot.value.df), ] # remove row contains NAs
#            plot.value.df <- plot.value.df[!is.infinite(rowSums(plot.value.df)),] # remove row contains Inf or Inf
#            if((length(plot.value.df$x) < 3) | (length(plot.value.df$y) < 3)){# not enough finite observations
#              return(0)
#            }
#            cor.res <- cor(plot.value.df$x, plot.value.df$y)
#            tryCatch({
#              cor.test.res <- cor.test(plot.value.df$x, plot.value.df$y, method="pearson")
#            }, 
#            error = function(e) {
#            },
#            silent = TRUE
#            )
#
#            # Make plot
#            res.plot <- ggplot(plot.value.df, aes(x=x,y=y)) + 
#            geom_point()+
#            geom_smooth(method=lm, se=FALSE)+
#            labs(x = top.feature.name, y = response.name, title = cor.res)
#            if(is.na(cor.res)){
#              cor.label <- "NA"
#            }else if(abs(cor.res) > 0.4){
#              print(paste(rank, ".", top.feature.name, " vs ",response.name, "=", cor.res, sep = ""))
#              cor.label <- paste0("R=", cor.res)
#            }else{
#              cor.label <- ""
#            }
#
#            # Save
#            print(paste0(sample.type, "_", response.name, "_feature", rank, "_", top.feature.name, "_", cor.label))
#            ggsave(file.path(save.dir, paste(sample.type, "_", response.name, "_feature", rank, "_", top.feature.name, "_", cor.label, ".png", sep = "")), plot = res.plot, device = NULL, path = NULL,
#            scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
#            dpi = 300, limitsize = TRUE)
#            write.csv(data.frame(Cor=cor.res, pvalue=cor.test.res$p.value),
#            file.path(save.dir, paste(sample.type, "_", response.name, "_feature", rank, "_", top.feature.name, "_", cor.label, "_statics.csv", sep = "")))
#            write.csv(data.frame(plot.value.df$x, plot.value.df$y, col.names=c(top.feature.name, response.name), row.names=rownames(feature.data.df)),
#            file.path(save.dir, paste(sample.type, "_", response.name, "_feature", rank, "_", top.feature.name, "_", cor.label, "_data.csv", sep = "")))
#      }, 
#      error = function(e) {
#      },
#      silent = TRUE)
#    }
#    for(rank in 1:min.rank){
#        top.feature.name <- rank.vec[1:min.rank][rank]
#        PlotScattter(feature.data.df, response.continuous.df, top.feature.name, rank)
#    }
#}

SavePCABiplot = function(value.set , feature.name, file.name){
  pdf(file.name, width = 10, height = 8)
  value.onlynumeric.set <- value.set[, feature.name]
  value.onlynumeric.set <- value.onlynumeric.set[, colSums(value.onlynumeric.set != 0) > 0]
  value.set.rpca = prcomp(x = value.onlynumeric.set, scale=TRUE)
  value.set.rpca.loading <- sweep(value.set.rpca$rotation, MARGIN=2, value.set.rpca$sdev, FUN="*")
  par(las=1)
  plot(x = NULL, type = "n", xlab = "PC1", ylab = "PC2", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty="n")
  axis(side = 1, at = seq(-1,1,0.2), tck = 1.0, lty = "dotted", lwd=0.5, col = "#dddddd", labels = expression(-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0))
  axis(side = 2, at = seq(-1,1,0.2), tck=1.0, lty="dotted", lwd=0.5, col = "#dddddd", labels=expression(-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0))
  for(i in 1:length(value.set.rpca.loading[,1]))
  {
    arrows(0, 0, value.set.rpca.loading[i, 1], value.set.rpca.loading[i, 2], col="#ff8c00", length=0.1)
  }
  pointLabel(x = value.set.rpca.loading[, 1], y = value.set.rpca.loading[, 2], labels = rownames(value.set.rpca.loading), cex = 0.7)
  box(bty = "l")

  importance.df <- as.data.frame(summary(value.set.rpca)$importance)["Proportion of Variance",]
  rownames(importance.df) <- c("Proportion_of_Variance")
  importance.df.melt <- melt(importance.df)
  colnames(importance.df.melt) <- c("PC", "Proportion_of_Variance")
  g <- ggplot(importance.df.melt, aes(x = PC, y = Proportion_of_Variance, fill = PC))
  g <- g + geom_bar(stat = "identity")
  g <- g + scale_y_continuous(labels = percent)
  g <- g + theme(axis.text.x = element_text(size = 10, colour="black", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10, colour="black")
  )
  plot(g)
  dev.off()
}

SaveTSNEplot = function(value.set , feature.name, file.name){
  #pdf(file.name, width = 10, height = 8)
  value.set.unique <- unique(value.set)
  value.set.matrix <- as.matrix(value.set.unique[, feature.name])
  set.seed(1)
  tsne.res <- Rtsne(value.set.matrix, perplexity = 5)
  plot(tsne.res$Y, col=rainbow(length(value.set.unique$group.name)))
  #plot(g)
  #dev.off()
}

CV = function(){
  # Leave-One-Out法で3次多項式モデル・9次多項式モデルそれぞれの
  # 交差検証データへの予測値を算出する
  lm3_vec<-rep(0,16)
  for (i in 1:16){
      tmp<-lm(y~V1+V2+V3,d[-i,c(1:3,10)])
      lm3_vec[i]<-predict(tmp,newdata=d[i,1:3])
  }
  lm9_vec<-rep(0,16)
  for (i in 1:16){
      tmp<-lm(y~.,d[-i,])
      lm9_vec[i]<-predict(tmp,newdata=d[i,-10])
  }
  # 学習データと見比べてみる
  plot(x,y,cex=4,xlim=c(0,8),ylim=c(-80,120))
  par(new=T)
  plot(x,lm3_vec,cex=4,xlim=c(0,8),ylim=c(-80,120),col='red')
  par(new=T)
  plot(x,lm9_vec,cex=4,xlim=c(0,8),ylim=c(-80,120),col='blue')


  #######
  # wleによる選択
  result <- wle.cp(log(Volume) ~ log(Girth) + log(Height), data = trees)  # cp 統計量
  result <- wle.cv(log(Volume) ~ log(Girth) + log(Height), data = trees)  # CV
  result <- wle.aic(log(Volume) ~ log(Girth) + log(Height), data = trees) # AIC
  summary(result)

  ######
  # DAAG
  # library(DAAG)
  cv.lm(df=mydata, fit, m=3) # 3 fold cross-validation

  #######
  # Calculate Relative Importance for Each Predictor
  library(relaimpo)
  calc.relimp(fit,type=c("lmg","last","first","pratt"),
    rela=TRUE)
  # Bootstrap Measures of Relative Importance (1000 samples) 
  boot <- boot.relimp(fit, b = 1000, type = c("lmg", 
    "last", "first", "pratt"), rank = TRUE, 
    diff = TRUE, rela = TRUE)
  booteval.relimp(boot) # print result
  plot(booteval.relimp(boot,sort=TRUE)) # plot result

  ######
  # elastic model
  # library(eNetXplorer)
}