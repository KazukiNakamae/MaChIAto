GetAbsolutePath <- function(relative.path) {
  # Get an absolute path from a given relative.path
  #
  # Args:
  #   relative.path: the charators of relative path
  #
  # Returns:
  #   the charators of absolute path
  return(normalizePath(relative.path))
}

PCbiplot <- function(PC, x="PC1", y="PC2") {
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

SubtractCtrlDF <- function(df){

    SubtractCtrl <- function(vec){
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

MakeValueSet <- function(df, sample.type){
   
    for(name in names(df)){
        if(!exists("sum.df")){
            sum.df <- df[[name]][df[[name]][,"sample.type"] == sample.type,]
        }else{
            sum.df <- rbind(
                sum.df,
                df[[name]][df[[name]][,"sample.type"] == sample.type,]
            )
        }
    }
    return(sum.df)
}

FltrByQuantile <- function(df, indicator){
  return(df[(as.numeric(df[[indicator]]) > quantile(as.numeric(df[[indicator]]), 0.00)), ])
}

MakeMHScoreTable <- function(fltr.extra.table, script.dir){
  bmh.vec <- fltr.extra.table[, "bmh"]
  cor.rmh.vec <- fltr.extra.table[, "rmh"]
  cor.lmh.vec <- fltr.extra.table[, "lmh"]
  rev.rmh.vec <- as.character(reverseComplement(DNAStringSet(fltr.extra.table[, "rmh"])))
  rev.lmh.vec <- as.character(reverseComplement(DNAStringSet(fltr.extra.table[, "lmh"])))
  cc.mh.score.CalcMHScore <- CalcMHScore$new(fltr.extra.table[, "bmh"]
      , cutpos.vec = nchar(fltr.extra.table[, "lmh"])
      , length_weight = 20
      , script.dir = script.dir
  )
  cc.mh.score.CalcMHScore.summary <- as.data.frame(cc.mh.score.CalcMHScore$summary)
  colnames(cc.mh.score.CalcMHScore.summary) <- c("bmh", "Endogenous.MP.Score", "Endogenous.MPOF.Score")
  cr.mh.score.CalcMHScore <- CalcMHScore$new(paste0(cor.lmh.vec, rev.rmh.vec)
      , cutpos.vec = nchar(cor.rmh.vec)
      , length_weight = 20
      , script.dir = script.dir
  )
  cr.mh.score.CalcMHScore.summary <- as.data.frame(cr.mh.score.CalcMHScore$summary)
  colnames(cr.mh.score.CalcMHScore.summary) <- c("seq", "LeftRevKI.MP.Score", "LeftRevKI.MPOF.Score")
  rc.mh.score.CalcMHScore <- CalcMHScore$new(paste0(rev.lmh.vec, cor.rmh.vec)
      , cutpos.vec = nchar(rev.lmh.vec)
      , length_weight = 20
      , script.dir = script.dir
  )
  rc.mh.score.CalcMHScore.summary <- as.data.frame(rc.mh.score.CalcMHScore$summary)
  colnames(rc.mh.score.CalcMHScore.summary) <- c("seq", "RightRevKI.MP.Score", "RightRevKI.MPOF.Score")
  bl.mh.score.CalcMHScore <- CalcMHScore$new(paste0(cor.lmh.vec, cor.lmh.vec)
      , cutpos.vec = nchar(cor.lmh.vec)
      , length_weight = 20
      , script.dir = script.dir
  )
  bl.mh.score.CalcMHScore.summary <- as.data.frame(bl.mh.score.CalcMHScore$summary)
  colnames(bl.mh.score.CalcMHScore.summary) <- c("seq", "LeftKI.MP.Score", "LeftKI.MPOF.Score")
  br.mh.score.CalcMHScore <- CalcMHScore$new(paste0(cor.rmh.vec, cor.rmh.vec)
      , cutpos.vec = nchar(cor.rmh.vec)
      , length_weight = 20
      , script.dir = script.dir
  )
  br.mh.score.CalcMHScore.summary <- as.data.frame(br.mh.score.CalcMHScore$summary)
  colnames(br.mh.score.CalcMHScore.summary) <- c("seq", "RightKI.MP.Score", "RightKI.MPOF.Score")
  bind.table <- cbind(cc.mh.score.CalcMHScore.summary, cr.mh.score.CalcMHScore.summary, rc.mh.score.CalcMHScore.summary, bl.mh.score.CalcMHScore.summary, br.mh.score.CalcMHScore.summary)
  bind.table <- bind.table[, which(!(colnames(bind.table) %in% "seq"))]
  return(bind.table)
}

MakeKmeansAnnotationSetDF <- function(vec){
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

MakeKmeansAnnotationSetDF2 <- function(vec){
  tryCatch({
    num.vec <- as.numeric(vec)
    if(sd(num.vec) == 0){
      return(character(length(num.vec)))
    }
    k <- 3
    while(TRUE){
      if(k > 10){
        return(rep(NA , length(num.vec)))
      }
      set.seed(1)
      kmeans.data <- kmeans(num.vec, k)
      center.vec <- kmeans.data$centers[,1]
      max.center.ind <- which.max(center.vec)
      second.center <- max(center.vec[center.vec!=max(center.vec)])
      second.center.ind <- which(center.vec == second.center)
      third.center <- max(center.vec[-c(max.center.ind, second.center.ind)])
      third.center.ind <- which(center.vec == third.center)

      if(kmeans.data$size[third.center.ind] >= 5){
        high1.ind <- max.center.ind
        high2.ind <- second.center.ind
        low.ind <- third.center.ind
        break
      }else{
        k <- k + 1
        next
      }
    }
    annotation.vec <- replace(kmeans.data$cluster, kmeans.data$cluster==high1.ind, "high")
    annotation.vec <- replace(annotation.vec, kmeans.data$cluster==high2.ind, "high")
    annotation.vec <- replace(annotation.vec, kmeans.data$cluster==low.ind, "low")
    if(k > 3){
      out.ind <- (1:k)[!(1:k %in% c(high1.ind, high2.ind, low.ind))]
      annotation.vec <- replace(annotation.vec, kmeans.data$cluster %in% out.ind, "out")
    }
    return(annotation.vec)
  }, 
  error = function(e) { 
    return(vec)
  }, silent = TRUE)
}

Removeoutliers <- function(value.set, response){
  outliers <- boxplot(value.set[, response], plot=FALSE)$out
  if(length(outliers)){
    value.set <- value.set[-which(value.set[, response] %in% outliers),]
  }
  return(value.set)
}

CalcMHProp <- function(microhomology.seqset, focus.length, scaffoldset, annotationset, is.left.microhomology.logical, out.align.file, out.prop.file){
    if(!file.exists(file.path(output.dir, out.prop.file))){
        if(is.left.microhomology.logical){
            focusing.region.end.pos.num = nchar(microhomology.seqset[1])
            focusing.region.start.pos.num = focusing.region.end.pos.num - (focus.length - 1)
        }else{
            focusing.region.start.pos.num = 1
            focusing.region.end.pos.num = focus.length
        }
        prop.info.CalcMicroHomologyPropPK <- CalcMicroHomologyPropPK$new(
            whole.region.microhomology.seq.vec = microhomology.seqset,
            focusing.region.start.pos.num = focusing.region.start.pos.num,
            focusing.region.end.pos.num = focusing.region.end.pos.num,
            is.left.microhomology.logical = is.left.microhomology.logical,
            scaffold.seq.vec = scaffoldset,
            need.annotation.logical = TRUE,
            annotation.vec = annotationset,
            align.range.vec = seq(20, nchar(scaffoldset[1]) , 20),
            #for debug
            #align.range.vec = seq(20, nchar(scaffoldset[1]) , nchar(scaffoldset[1])),
            align.save.file.char = file.path(output.dir, out.align.file),
            script.dir = script.dir,
            ref.dir = ref.dir
        )
        saveRDS(prop.info.CalcMicroHomologyPropPK, file.path(output.dir, out.prop.file))
    }else{
        prop.info.CalcMicroHomologyPropPK <- readRDS(file.path(output.dir, out.prop.file))
    }
    return(prop.info.CalcMicroHomologyPropPK)
}

MergeExtraData <- function(main.df, ext.df, seq.name, value.names.vec){
    temp.df <- data.frame("sequence" = as.vector(ext.df[,seq.name]), stringsAsFactors = FALSE)
    for(name in colnames(ext.df)){
        if(name %in% value.names.vec){
            temp.df <- cbind(temp.df, ext.df[,name])
            colnames(temp.df)[length(temp.df)] <- name
        }
    }
    sum.df <- merge(main.df[, 1:length(main.df) - 1], temp.df, by ="sequence", all=F)
    sum.df <- merge(sum.df, main.df[,c("sequence", "label")], by ="sequence", all=F)
    return(sum.df)
}

PlotBarPlot <- function(vec, type, title){
  temp.df <- data.frame(future.name = names(vec), value = vec)
  if(type == "mono.grad"){
    temp.df %>% 
      ggplot(aes(reorder(future.name, value), value)) + 
      geom_col(aes(fill = value)) + 
      scale_fill_gradient(low = "gray", 
                          high = "black") + 
      coord_flip() + 
      labs(x = "Value", y = "Future", title = title) +
      theme(axis.text = element_text(size = 8, colour="black"))  -> p
  }else if(type == "black"){
    temp.df %>% 
      ggplot(aes(reorder(future.name, value), value)) + 
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
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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

MakeCorPlot = function(x, y, x.name, y.name, figure.dir){
  cor <- cor(x, y)
  df <- data.frame(x=x, y=y, stringsAsFactors = FALSE)
  # Make plot
  res.plot <- ggplot(df, aes(x=x,y=y)) + 
    geom_point()+
    geom_smooth(method=lm, se=FALSE)+
    labs(x = x.name, y = y.name, title = cor)
  ggsave(file.path(figure.dir, paste(x.name, "_vs_", y.name, ".png", sep = "")), plot = res.plot, device = NULL, path = NULL,
    scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
    dpi = 300, limitsize = TRUE)
  return(res.plot)
}

MakeCorCompPlot = function(future.data.df1, future.data.df2, response.continuous.df1, response.continuous.df2, future.name.vec, output.dir, figure.dir, future.set.name, samplename, group.name.vec){
  CalcCor = function(future.data.df, response.continuous.df, future.name){
    scat.data.df <- merge(future.data.df, response.continuous.df, by="sequence", all=F)
    response.name <- colnames(scat.data.df)[length(scat.data.df)]
    print(future.name)
    # browser()
    plot.value.df <- data.frame(
      x = as.numeric(scat.data.df[, future.name]),
      y = as.numeric(scat.data.df[, response.name]),
      stringsAsFactors = FALSE)
    norm.plot.value.df<- as.data.frame(scale(plot.value.df), stringsAsFactors = FALSE)
    if((length(norm.plot.value.df$x) < 3) | (length(norm.plot.value.df$y) < 3)){# not enough finite observations
      return(0)
    }
    cor <- cor(norm.plot.value.df$x, norm.plot.value.df$y)
    return(cor)
  }
  # merge cor data
  cor1.vec <- numeric(0)
  cor2.vec <- numeric(0)
  for(future.name in future.name.vec){
    temp1 <- c(CalcCor(future.data.df1, response.continuous.df1, future.name))
    names(temp1) <- future.name
    cor1.vec <- c(cor1.vec, temp1)
    temp2 <- c(CalcCor(future.data.df2, response.continuous.df2, future.name))
    names(temp2) <- future.name
    cor2.vec <- c(cor2.vec, temp2)
  }
  # browser()
  # prepare plot
  t.test.res <- t.test(x = cor1.vec, y = cor2.vec, paired = TRUE)
  ks.test.res <- ks.test(x = cor1.vec, y = cor2.vec)
  plot.cor.df <- data.frame(
      cor = c(cor1.vec, cor2.vec),
      group = as.factor(c(rep(group.name.vec[1], length(cor1.vec)), rep(group.name.vec[2], length(cor2.vec)))),
      future.name = c(names(cor1.vec), names(cor2.vec)),
      stringsAsFactors = FALSE)
  plot.cor.df <- plot.cor.df[complete.cases(plot.cor.df), ]
  groupall.n <- nrow(plot.cor.df)
  group1.n <- nrow(plot.cor.df[plot.cor.df$group == group.name.vec[1],])
  group2.n <- groupall.n - group1.n
  is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }
  outliner.vec <- c(
    rep(TRUE, nrow(plot.cor.df[plot.cor.df$group == group.name.vec[1],])),
    !is_outlier(plot.cor.df[plot.cor.df$group == group.name.vec[2],]$cor)
  )
  plot.cor.df$future.name[outliner.vec] <- ""
  # plot
  #p <- ggpaired(plot.cor.df, x = "group", y = "cor", fill = c("#8C3A32", "#F2A88D")
  #  , line.color = "#011826", line.size = 0.4, title = paste0(samplename, "_", future.set.name, "_p.value=", t.test.res$p.value, "_n1=", group1.n, "_n2=", group2.n)
  #  , ylab = "Pearson Correlation")+
  #  stat_compare_means(paired = TRUE)+
  #  geom_text(aes(label = future.name), size = 3)+
  #  ylim(-1, 1)
  p1 <- ggplot(plot.cor.df, aes(x=group, y=cor)) +
    geom_boxplot(fill = c("#D96A6A", "#F2DA91"), color = "black", alpha = 1, show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(width = 0.25, aes(colour=group), size = 0.5, alpha = 1) +
    ggtitle(paste0(samplename, "_", future.set.name)) +
    geom_text(aes(label = future.name), size = 3) +
    scale_color_manual(values = c("#8C3A32", "#F2A88D"))+
    ylim(-1, 1)
  # browser()
  ggsave(file.path(figure.dir, paste(samplename, "_", future.set.name, "_Tp.value=", t.test.res$p.value, "_KSp.value="
      , ks.test.res$p.value, "_n1=", group1.n, "_n2=", group2.n, "_labeled.png", sep = ""))
      , plot = p1, device = NULL, path = NULL,
    scale = 1, width = 3, height = 6, units = c("in", "cm", "mm"),
    dpi = 300, limitsize = TRUE)
  p2 <- ggplot(plot.cor.df, aes(x=group, y=cor)) +
    geom_boxplot(fill = c("#D96A6A", "#F2DA91"), color = "black", alpha = 1, show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(width = 0.25, aes(colour=group), size = 0.5, alpha = 1) +
    ggtitle(paste0(samplename, "_", future.set.name)) +
    scale_color_manual(values = c("#8C3A32", "#F2A88D"))+
    ylim(-1, 1)
  # browser()
  ggsave(file.path(figure.dir, paste(samplename, "_", future.set.name, "_Tp.value=", t.test.res$p.value, "_KSp.value="
      , ks.test.res$p.value, "_n1=", group1.n, "_n2=", group2.n, ".png", sep = ""))
      , plot = p2, device = NULL, path = NULL,
    scale = 1, width = 3, height = 6, units = c("in", "cm", "mm"),
    dpi = 300, limitsize = TRUE)
}

PlotTopRankScattter = function(future.data.df, response.continuous.df, rank.vec, min.rank, sample.type){
    PlotScattter = function(future.data.df, response.continuous.df, top.future.name, rank){
      tryCatch({# urgent deal
            scat.data.df <- merge(future.data.df, response.continuous.df, by="sequence", all=F)
            response.name <- colnames(scat.data.df)[length(scat.data.df)]
            plot.value.df <- data.frame(
                x = as.numeric(scat.data.df[, top.future.name]),
                y = as.numeric(scat.data.df[, response.name]),
                stringsAsFactors = FALSE)
            plot.value.df <- plot.value.df[complete.cases(plot.value.df), ] # remove row contains NAs
            plot.value.df <- plot.value.df[!is.infinite(rowSums(plot.value.df)),] # remove row contains Inf or Inf
            norm.plot.value.df<- as.data.frame(scale(plot.value.df), stringsAsFactors = FALSE)
            if((length(norm.plot.value.df$x) < 3) | (length(norm.plot.value.df$y) < 3)){# not enough finite observations
              return(0)
            }
            cor <- cor(norm.plot.value.df$x, norm.plot.value.df$y)
            # Make plot
            res.plot <- ggplot(norm.plot.value.df, aes(x=x,y=y)) + 
            geom_point()+
            geom_smooth(method=lm, se=FALSE)+
            labs(x = top.future.name, y = response.name, title = cor)
            if(abs(cor) > 0.4){
              print(paste(rank, ".", top.future.name, " vs ",response.name, "=", cor, sep = ""))
              cor.label <- paste0("_cor=", cor)
            }else{
              cor.label <- ""
            }
            # Save
            ggsave(file.path(figure.dir, paste(rank, ".", sample.type, ".", top.future.name, " vs ",response.name, cor.label, ".png", sep = "")), plot = res.plot, device = NULL, path = NULL,
            scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
            dpi = 300, limitsize = TRUE)
      }, 
      error = function(e) {
      },
      silent = TRUE)
    }
    for(rank in 1:min.rank){
        top.future.name <- rank.vec[1:min.rank][rank]
        PlotScattter(future.data.df, response.continuous.df, top.future.name, rank)
    }
}

SavePCABiplot = function(value.set , future.name, file.name){
  pdf(file.name, width = 10, height = 8)
  value.onlynumeric.set <- value.set[, future.name]
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

SaveTSNEplot = function(value.set , future.name, file.name){
  #pdf(file.name, width = 10, height = 8)
  value.set.unique <- unique(value.set)
  value.set.matrix <- as.matrix(value.set.unique[, future.name])
  set.seed(1)
  # browser()
  tsne.res <- Rtsne(value.set.matrix, perplexity = 5)
  plot(tsne.res$Y, col=rainbow(length(value.set.unique$group.name)))
  #plot(g)
  #dev.off()
}

CV <- function(){
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