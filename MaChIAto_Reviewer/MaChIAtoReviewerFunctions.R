
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

SaveTable <- function(table, file.name.without.ext){
  options(warn = -1)
  if(is.list(table) == TRUE){
    write.table(data.frame(matrix(unlist(table), nrow=length(table), byrow=T)) , file = paste0(file.name.without.ext, ".txt"), quote=FALSE, col.names=TRUE, row.names=TRUE,append=TRUE, sep = ",")
  }else{
    write.table(table , file = paste0(file.name.without.ext, ".txt"), quote=FALSE, col.names=TRUE, row.names=TRUE,append=TRUE, sep = ",")
  }
  options(warn = 0)
  saveRDS(table , file = paste0(file.name.without.ext, ".rds"))
}

SavePieChart <- function(class.count.table, color.vec, elements.name, file.name){

  if(sum(class.count.table$Frequency) < 1){
    return()
  }

    class.count.plot.table <- class.count.table %>%
    arrange(desc(Class)) %>%
    mutate(Lab1.ypos = cumsum(Percentage) - 0.5*Percentage)

    class.count.plot.table <- class.count.plot.table %>%
    arrange(desc(Class)) %>%
    mutate(Lab2.ypos = Lab1.ypos)
    
  options(warn = -1)
  class.count.pie.chart <- ggplot(class.count.plot.table, aes(x="", y=Percentage, fill=Class))+
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=color.vec) +
    geom_text(aes(x = 1.30, y=Lab1.ypos, label = paste(round(Percentage, digits = 1), "%", sep = "")), color = "black",  size=5) + 
    geom_text(aes(x = 1.82, y=Lab2.ypos, label = paste(Class, "\n(", Frequency, "", elements.name, ")", sep = "")), color = "black",  size=5) + 
    coord_polar("y", start = pi * 3/2, direction = -1) +
    theme_void()
  # plot(class.count.pie.chart)
  ggsave(file = file.name, plot = class.count.pie.chart, dpi = 350, width = 10, height = 10)
  options(warn = 0)
}

# make summary of rate barplot
MakeSummaryClassBarplot <- function(rate.table.mat, condition.label, col.vec, label.x, save.path){
  melt.crispresso.rate.table.list <- melt(rate.table.mat)
  colnames(melt.crispresso.rate.table.list) <- c("Class", "Sample", "Percentage")

  melt.crispresso.rate.table.list.name.df <- t(sapply(melt.crispresso.rate.table.list$Sample, function(x){
    strsplit(as.character(x), "-")[[1]]
  }))
  colnames(melt.crispresso.rate.table.list.name.df) <- c("Target", "Label")
  crispresso.rate.barplot.df <- cbind(melt.crispresso.rate.table.list, melt.crispresso.rate.table.list.name.df)

  # ordering
  target.eff.rank <- rank(filter(crispresso.rate.barplot.df, Class == "Precise knock-in", Label == condition.label)$Percentage)
  crispresso.rate.barplot.df$Target <-
    factor(crispresso.rate.barplot.df$Target,
      levels = filter(crispresso.rate.barplot.df, Class == "Precise knock-in", Label == condition.label)$Target[order(target.eff.rank)])

  options(warn = -1)
  crispresso.rate.barplot <- ggplot(data = crispresso.rate.barplot.df) +
    geom_bar(aes(y = Percentage, x = Label, fill = Class), stat="identity", width=0.50) +
    facet_grid(~Target, scale='free_x') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='gray50')) +
    scale_fill_manual(values=col.vec) +
    labs(x = label.x)
  ggsave(file = save.path, plot = crispresso.rate.barplot, dpi = 350, width = 40, height = 3, limitsize = FALSE)
  options(warn = 0)
}

ExtractEndChar <- function(char.vec, len){
  substr(char.vec, nchar(char.vec) - (len - 1), nchar(char.vec))
}

ExtractStartChar <- function(char.vec, len){
  strtrim(char.vec, nchar(char.vec) - len)
}

MakeClassRateVec <- function(rate.table.mat, class.name, sample.label){
  condition.rate.vec <- rate.table.mat[class.name,
    ExtractEndChar(colnames(rate.table.mat),
    nchar(sample.label)) == sample.label
  ]
  names(condition.rate.vec) <- ExtractStartChar(names(condition.rate.vec), nchar(sample.label) + 1)
  return(condition.rate.vec)
}

MakeCrisprsetList <- function(res.file.df, res.file.path){
  crisprset.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(file.exists(file.path(row["aligner.res.dir.path"], res.file.path))){
      res <- list(readRDS(file.path(row["aligner.res.dir.path"], res.file.path))$cigar_freqs)
      names(res) <- row["sample.name"]
    }else{
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], res.file.path), " was not found."))
      res <- NULL
    }
    return(res)
  })
  return(Filter(Negate(is.null), crisprset.list)) # return list without NULL element
}

MakeCrisprsetListSimple <- function(res.file.df, res.file.path){
  crisprset.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(file.exists(file.path(row["aligner.res.dir.path"], res.file.path))){
      res <- readRDS(file.path(row["aligner.res.dir.path"], res.file.path))$cigar_freqs
    }else{
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], res.file.path), " was not found."))
      res <- NULL
    }
    return(res)
  })
  return(Filter(Negate(is.null), crisprset.list)) # return list without NULL element
}

MakeSpecificCrisprsetList <- function(res.file.df, res.file.path, condition.label){
  crisprset.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] %in% c(condition.label) & file.exists(file.path(row["aligner.res.dir.path"], res.file.path))){
      res <- list(readRDS(file.path(row["aligner.res.dir.path"], res.file.path))$cigar_freqs)
      names(res) <- row["sample.name"]
    }else if(row["sample.label"] %in% c(condition.label)){
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], res.file.path), " was not found."))
      res <- NULL
    }else{
      res <- NULL
    }
    return(res)
  })
  return(Filter(Negate(is.null), crisprset.list)) # return list without NULL element
}


CalcKLdivergence <- function(x.cigar_freqs, y.cigar_freqs){
  if(is.null(x.cigar_freqs) & is.null(y.cigar_freqs)){
    return(0)
  }else if(is.null(x.cigar_freqs) | is.null(y.cigar_freqs)){
    return(Inf)
  }
  CalcKLValue <- function(x.cigar_freqs, y.cigar_freqs){
    # prepare vector
    variants.names.vec <- unique(c(rownames(x.cigar_freqs), rownames(y.cigar_freqs)))
    # 0.5 is small pseudocount for avoinding division by zero. This method was inspired by [Felicity Allen, et.al, 2017]
    x.variants.vec <- rep(0.5, length(variants.names.vec))
    names(x.variants.vec) <- variants.names.vec
    x.variants.vec[rownames(x.cigar_freqs)] <- x.variants.vec[rownames(x.cigar_freqs)] + x.cigar_freqs[,1]
    x.variants.vec <- x.variants.vec / sum(x.variants.vec)
    y.variants.vec <- rep(0.5, length(variants.names.vec))
    names(y.variants.vec) <- variants.names.vec
    y.variants.vec[rownames(y.cigar_freqs)] <- y.variants.vec[rownames(y.cigar_freqs)] + y.cigar_freqs[,1]
    y.variants.vec <- y.variants.vec / sum(y.variants.vec)
    return(sum(x.variants.vec * log(x.variants.vec / y.variants.vec)))
  }
  MakePositionIndependentVariantsFreq <- function(cigar_freqs){
    rownames(cigar_freqs) <- gsub("-*", "", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub("[0-9]*[ATCG]", "S", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub("M:", ";;;", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub("[a-zA-Z0-9]*:", "", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub(";;;", "M:", rownames(cigar_freqs))
    aggregated.freq.table <- aggregate(reads ~ label, data = data.frame(label = rownames(cigar_freqs),  reads = cigar_freqs[,1]), sum)
    position.independent.cigar_freqs <- matrix(
      aggregated.freq.table$reads,
      ncol = 1
    )
    rownames(position.independent.cigar_freqs) <- aggregated.freq.table$label
    return(position.independent.cigar_freqs)
  }

  position.independent.x.cigar_freqs <- MakePositionIndependentVariantsFreq(x.cigar_freqs)
  position.independent.y.cigar_freqs <- MakePositionIndependentVariantsFreq(y.cigar_freqs)
  return(
    CalcKLValue(position.independent.x.cigar_freqs, position.independent.y.cigar_freqs) +
    CalcKLValue(position.independent.y.cigar_freqs, position.independent.x.cigar_freqs)
  )
}

CalcKLdivergence2 <- function(x.cigar_freqs, y.indel.data){
  if(is.null(x.cigar_freqs) & is.null(y.indel.data)){
    return(0)
  }else if(is.null(x.cigar_freqs) | is.null(y.indel.data)){
    return(Inf)
  }
  CalcKLValue <- function(x.cigar_freqs, y.cigar_freqs){
    # prepare vector
    variants.names.vec <- unique(c(rownames(x.cigar_freqs), rownames(y.cigar_freqs)))
    # 0.5 is small pseudocount for avoinding division by zero. This method was inspired by [Felicity Allen, et.al, 2017]
    x.variants.vec <- rep(0.5, length(variants.names.vec))
    names(x.variants.vec) <- variants.names.vec
    x.variants.vec[rownames(x.cigar_freqs)] <- x.variants.vec[rownames(x.cigar_freqs)] + x.cigar_freqs[,1]
    x.variants.vec <- x.variants.vec / sum(x.variants.vec)
    y.variants.vec <- rep(0.5, length(variants.names.vec))
    names(y.variants.vec) <- variants.names.vec
    y.variants.vec[rownames(y.cigar_freqs)] <- y.variants.vec[rownames(y.cigar_freqs)] + y.cigar_freqs[,1]
    y.variants.vec <- y.variants.vec / sum(y.variants.vec)
    return(sum(x.variants.vec * log(x.variants.vec / y.variants.vec)))
  }
  MakePositionIndependentVariantsFreq <- function(cigar_freqs){
    rownames(cigar_freqs) <- gsub("-*", "", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub("[0-9]*[ATCG]", "S", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub("M:", ";;;", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub("[a-zA-Z0-9]*:", "", rownames(cigar_freqs))
    rownames(cigar_freqs) <- gsub(";;;", "M:", rownames(cigar_freqs))
    aggregated.freq.table <- aggregate(reads ~ label, data = data.frame(label = rownames(cigar_freqs),  reads = cigar_freqs[,1]), sum)
    position.independent.cigar_freqs <- matrix(
      aggregated.freq.table$reads,
      ncol = 1
    )
    rownames(position.independent.cigar_freqs) <- aggregated.freq.table$label
    return(position.independent.cigar_freqs)
  }
  position.independent.x.cigar_freqs <- MakePositionIndependentVariantsFreq(x.cigar_freqs)
  #browser()
  x.indel.data <- matrix(position.independent.x.cigar_freqs[which(rownames(position.independent.x.cigar_freqs) %in% rownames(y.indel.data)),], ncol = 1)
  rownames(x.indel.data) <- rownames(position.independent.x.cigar_freqs)[which(rownames(position.independent.x.cigar_freqs) %in% rownames(y.indel.data))]
  return(
    CalcKLValue(x.indel.data, y.indel.data) +
    CalcKLValue(y.indel.data, x.indel.data)
  )
}

SaveVariantsBarPlot <- function(cnt.simple.list, file.name.without.ext){
  variants.cnt.vec <- numeric(0)
  variants.name.vec <- character(0)
  for(ind in 1:length(cnt.simple.list)){
    variants.cnt.vec <- c(variants.cnt.vec, cnt.simple.list[[ind]])
    variants.name.vec <- c(variants.name.vec, rownames(cnt.simple.list[[ind]]))
  }
  variants.cnt.mat <- matrix(variants.cnt.vec, ncol=1)
  rownames(variants.cnt.mat) <- variants.name.vec
  # make aggregated table
  aggregated.variants.cnt.vec <- sapply(by(variants.cnt.mat, rownames(variants.cnt.mat),colSums), identity)
  aggregated.variants.cnt.mat <- matrix(aggregated.variants.cnt.vec, ncol=1)
  rownames(aggregated.variants.cnt.mat) <- names(aggregated.variants.cnt.vec)
  # make aggregated summary table
  ind.snv.aggregated.variants.cnt.mat <- grep("SNV", rownames(aggregated.variants.cnt.mat))
  ind.novariant.aggregated.variants.cnt.mat <- grep("no variant", rownames(aggregated.variants.cnt.mat))
  ind.indel.aggregated.variants.cnt.mat <- (1:nrow(aggregated.variants.cnt.mat))[-c(ind.snv.aggregated.variants.cnt.mat, ind.novariant.aggregated.variants.cnt.mat)]
  ind.majorindel.aggregated.variants.cnt.mat <- which(aggregated.variants.cnt.mat[ind.indel.aggregated.variants.cnt.mat] > 10000)
  ind.minorindel.aggregated.variants.cnt.mat <- ind.indel.aggregated.variants.cnt.mat[-c(ind.majorindel.aggregated.variants.cnt.mat)]
  variants.read.table <- data.frame(
    LABELS = c(
      "No variants",
      "SNV",
      "Minor InDels",
      names(aggregated.variants.cnt.mat[ind.majorindel.aggregated.variants.cnt.mat,])
    ),
    READS = c(
      sum(aggregated.variants.cnt.mat[ind.novariant.aggregated.variants.cnt.mat]),
      sum(aggregated.variants.cnt.mat[ind.snv.aggregated.variants.cnt.mat]),
      sum(aggregated.variants.cnt.mat[ind.minorindel.aggregated.variants.cnt.mat]),
      aggregated.variants.cnt.mat[ind.majorindel.aggregated.variants.cnt.mat]
    ),
    GROUPS = rep("ALL_CRISPResso", length(c("No variants", "SNV", "Minor InDels", names(aggregated.variants.cnt.mat[ind.majorindel.aggregated.variants.cnt.mat,]))))
  )
  SaveTable(variants.read.table, paste0(file.name.without.ext, ".read"))
  variants.read.plot <- ggplot(variants.read.table, aes(x = GROUPS, y = READS, fill = LABELS)) +
    geom_bar(stat = "identity", colour = "black")
  ggsave(file = paste0(file.name.without.ext, ".read.png"), plot = variants.read.plot, dpi = 300, width = 10, height = 10)

  variants.variants.table <- data.frame(
    LABELS = c(
      "No variants",
      "SNV",
      "Minor InDels",
      "Major InDels"),
    VARIANTS = c(
      length(ind.novariant.aggregated.variants.cnt.mat),
      length(ind.snv.aggregated.variants.cnt.mat),
      length(ind.minorindel.aggregated.variants.cnt.mat),
      length(ind.majorindel.aggregated.variants.cnt.mat)
    ),
    GROUPS = rep("ALL_CRISPResso", length(c("No variants", "SNV", "Minor InDels", "Major InDels")))
  )
  SaveTable(variants.variants.table, paste0(file.name.without.ext, ".variants"))
  variants.variants.plot <- ggplot(variants.variants.table, aes(x = GROUPS, y = VARIANTS, fill = LABELS)) +
    geom_bar(stat = "identity", colour = "black")
  ggsave(file = paste0(file.name.without.ext, ".variants.png"), plot = variants.variants.plot, dpi = 300, width = 10, height = 10)
}

DistinguishVariantsByUnmodified <- function(cigar_freqs){
  match.no.variants.vec <- str_match(rownames(cigar_freqs), "no variant")
  match.snv.vec <- str_match(rownames(cigar_freqs), "SNV")
  has.non.indel <- c(which(match.no.variants.vec == "no variant"), which(match.snv.vec == "SNV"))
  has.non.indel <- has.non.indel[complete.cases(has.non.indel)]
  if(length(has.non.indel) > 0){
    indel.non.indel.cnt.vec <- c(length(cigar_freqs[has.non.indel,]) ,length(cigar_freqs[-has.non.indel,])) # length() : variants type / sum() : read
  }else{
    indel.non.indel.cnt.vec <- c(0 ,length(cigar_freqs[,1]))# length() : variants type / sum() : read
  }
  names(indel.non.indel.cnt.vec) <- c("No variant or SNV", "Indel or Knock-in")
  return(indel.non.indel.cnt.vec)
}

SummaryVariantsByUnmodified <- function(res.file.df, path.name){
  crisprset.list <- MakeCrisprsetList(res.file.df, path.name)
  unmodified.snv.mat <- matrix(, nrow = 2)
  for(ind in 1:length(crisprset.list)){
    temp.sample.name <- names(crisprset.list[[ind]])
    cigar_freqs <- getElement(crisprset.list[[ind]], temp.sample.name)
    unmodified.snv.mat <- cbind(unmodified.snv.mat,
      DistinguishVariantsByUnmodified(cigar_freqs))
    colnames(unmodified.snv.mat)[ind + 1] <- temp.sample.name
  }
  unmodified.snv.mat <- unmodified.snv.mat[,-1]
  rowsums.unmodified.snv.vec <- rowSums(unmodified.snv.mat)
  return(rowsums.unmodified.snv.vec)
}

MakePreciseKnockinLengthVec <- function(res.file.df, wt.file.path, precise.knockin.file.path){
  knockin.length.vec <- apply(res.file.df, MARGIN = 1, function(row){
    if(file.exists(file.path(row["aligner.res.dir.path"], wt.file.path))){
      wt.dna <- readDNAStringSet(file.path(row["aligner.res.dir.path"], wt.file.path))[[1]]
      precise.knockin.dna <- readDNAStringSet(file.path(row["aligner.res.dir.path"], precise.knockin.file.path))[[1]]
      res <- length(precise.knockin.dna) - length(wt.dna)
    }else{
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], wt.file.path), " was not found."))
      res <- NA
    }
    return(res)
  })
  return(Filter(Negate(is.na), knockin.length.vec)) # return list without NULL element
}

#'@title Make mutant list 
#'@description Make mutant list from mutant label eg "-1:2D"
#'@param mut.label mutant label vector
#'@author Kazuki Nakamae
#'@return The list including mutation infomation classified by mutant type
#'@rdname MakeMutList
MakeMutList <- function(mut.label){

  # Convert SNV:xB, yB -> SNV:xB, SNV:yB
  for(ind in 1:length(mut.label)){
    which.muttype <- !is.na(charmatch(c("D", "I", "SNV", "M"), mut.label[ind]))
    has.muttype <- sum(which.muttype)
    if((has.muttype < 1) && ind > 1){
      if(pre.muttype.label %in% c("I", "D", "M")){
        mut.label[ind] <- paste(mut.label[ind], pre.muttype.label, sep="")
      }else{
        mut.label[ind] <- paste(pre.muttype.label, ":", mut.label[ind], sep="")
      }
    }else{
      pre.muttype.label <- c("D", "I", "SNV", "M")[which.muttype]

      # if it is special label(eg."no variant"), it will be ignored. 
      if(length(pre.muttype.label) < 1){
        break
      }

    }
  }

  options(warn=-1)
  del.group.ind <- grep("D", mut.label)
  del.list <- list()
  if(length(del.group.ind) > 0){
    for(mut.type in mut.label[del.group.ind]){
      del.type <- gsub("D$", "", mut.type)
      del.vec <- as.numeric(strsplit(del.type, ":")[[1]])
      names(del.vec) <- c("pos", "size")
      temp.list.name <- names(del.list)
      del.list <- c(del.list, list(del.vec))
      names(del.list) <- c(temp.list.name, mut.type)
    }
  }
  ins.group.ind <- grep("I", mut.label)
  ins.list <- list()
  if(length(ins.group.ind) > 0){
    for(mut.type in mut.label[ins.group.ind]){
      ins.type <- gsub("I$", "", mut.type)
      ins.vec <- as.numeric(strsplit(ins.type, ":")[[1]])
      names(ins.vec) <- c("pos", "size")
      temp.list.name <- names(ins.list)
      ins.list <- c(ins.list, list(ins.vec))
      names(ins.list) <- c(temp.list.name, mut.type)
    }
  }
  snv.group.ind <- grep("SNV", mut.label)
  snv.list <- list()
  if(length(snv.group.ind) > 0){
    for(mut.type in mut.label[snv.group.ind]){
      snv.type <- gsub("^SNV:", "", mut.type)
      snv.vec <- c(as.numeric(gsub("[[:alpha:]]", "", snv.type)), gsub("[[:digit:]]", "", snv.type))
      names(snv.vec) <- c("pos", "base")
      temp.list.name <- names(snv.list)
      snv.list <- c(snv.list, list(snv.vec))
      names(snv.list) <- c(temp.list.name, mut.type)
    }
  }
  large.group.ind <- grep("M", mut.label)
  large.list <- list()
  if(length(large.group.ind) > 0){
    for(mut.type in mut.label[large.group.ind]){
      large.type <- gsub("M$", "", mut.type)
      large.vec <- c(as.numeric(large.type))
      names(large.vec) <- c("pos")
      temp.list.name <- names(large.list)
      large.list <- c(large.list, list(large.vec))
      names(large.list) <- c(temp.list.name, mut.type)
    }
  }
  mut.list <- list(deletion = del.list, insertion = ins.list, substitution = snv.list, large.indel = large.list)
  options(warn=0)
  return(mut.list)
}

DistinguishVariantsByPreciseKnockin <- function(cigar_freqs, precise.knockin.length){

  mut.list.matrix <- sapply(rownames(cigar_freqs), function(raw.mut.label){
    MakeMutList(strsplit(raw.mut.label, ",")[[1]])
    })
  indel.size.vec <- apply(mut.list.matrix, MARGIN = 2, function(col){
    indel.size <- 0
    if(length(col$insertion) > 0){
      indel.size <- indel.size + sum(unlist(col$insertion)[seq(2, length(unlist(col$insertion)), 2)])
    }
    if(length(col$deletion) > 0){
      indel.size <- indel.size - sum(unlist(col$deletion)[seq(2, length(unlist(col$deletion)), 2)])
    }
    return(indel.size)
  })
  precise.knockin.cnt <- length(which(indel.size.vec %in% precise.knockin.length))
  non.precise.knockin.cnt <- length(indel.size.vec) - precise.knockin.cnt
  precise.knockin.non.precise.knockin.cnt.vec <- c(precise.knockin.cnt, non.precise.knockin.cnt)

  names(precise.knockin.non.precise.knockin.cnt.vec) <- c("precise size", "Others")
  return(precise.knockin.non.precise.knockin.cnt.vec)
}

SummaryVariantsByPreciseKnockin <- function(res.file.df, res.file.path, wt.file.path, precise.knockin.file.path){
  crisprset.list <- MakeCrisprsetList(res.file.df, res.file.path)
  knockin.length.vec <- MakePreciseKnockinLengthVec(res.file.df, wt.file.path, precise.knockin.file.path)
  precise.knockin.mat <- matrix(, nrow = 2)
  for(ind in 1:length(crisprset.list)){
    temp.sample.name <- names(crisprset.list[[ind]])
    cigar_freqs <- getElement(crisprset.list[[ind]], temp.sample.name)
    precise.knockin.mat <- cbind(precise.knockin.mat,
      DistinguishVariantsByPreciseKnockin(cigar_freqs, knockin.length.vec[ind]))
    colnames(precise.knockin.mat)[ind + 1] <- temp.sample.name
  }
  precise.knockin.mat <- precise.knockin.mat[,-1]
  rowsums.precise.knockin.vec <- rowSums(precise.knockin.mat)
  return(rowsums.precise.knockin.vec)
}

SafeGetElement <- function(x, name){
  res <- NULL
  tryCatch({
    res <- getElement(x, name)
  }, 
    error = function(e) {
    res <- NULL
  },
    silent = TRUE
  )
  return(res)
}

MakeAlluvialDiagrams <- function(melt.df, color.vec, save.path.name, axis1.title, axis2.title, focus.type, unit = "kinds"){

  # If there is zero element, making plot will fail. To avoid it, we add 0.1 into the element.
  melt.df$value[melt.df$value == 0] <- 0.1 

  options(warn = -1)
  alluvial.plot <- ggplot(as.data.frame(melt.df), aes(y = value, axis1 = Var1, axis2 = Var2)) +
  geom_alluvium(aes(fill = Var2) , width = 1/12) +
  geom_stratum(width = 1/3, fill = color.vec, color = "white") +
  geom_text(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c(axis1.title, axis2.title), expand = c(.05, .05)) +
  scale_y_continuous(breaks=seq(0, sum(melt.df$value), round(sum(melt.df$value)/10))) +
  ggtitle(paste0("Comparison between ", axis1.title, " and ", axis2.title, " in ", focus.type)) +
  labs(x = "", y = paste0("Sum of identified variants in samples [", unit, "]"), fill = axis2.title) +
  theme(
    axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 10)
    , axis.title.y = element_text(size = 12)
    , legend.title = element_text(size = 12)
    , legend.text = element_text(size = 12))
  ggsave(file = save.path.name, plot = alluvial.plot, dpi = 350, width = 8, height = 10)
  options(warn = 0)
}

MakeBumpChart <- function(all.rank.table, group.table, standard.label, label.vec, save.path.name){
  rank.table <- all.rank.table[, label.vec]
  malt.rank.table <- melt(rank.table)
  rank.by.standard.vec <- rownames(rank.table[order(rank.table[, standard.label]), ])
  malt.rank.table$Var1 <- factor(malt.rank.table$Var1, level = rank.by.standard.vec)
  standard.group.vec <- group.table[, standard.label]
  line.color.vec <- c(
    rep("#F21905", length(standard.group.vec[standard.group.vec == "High"])),
    #rep("#F27405", length(standard.group.vec[standard.group.vec == 2])),
    #rep("#D9CB04", length(standard.group.vec[standard.group.vec == 3])),
    rep("#04BF45", length(standard.group.vec[standard.group.vec == "Medium"])),
    rep("#0487D9", length(standard.group.vec[standard.group.vec == "Low"]))
  )
  #line.color.vec <- rev(rainbow(2 * nrow(rank.table)))[1:nrow(rank.table)] # it is used in no group.
  options(warn = -1)
  bump.plot <- ggplot(malt.rank.table, aes(factor(Var2, level = label.vec), value,
      group = Var1, colour = Var1, label = Var1)) +
    geom_line(aes(color = Var1), alpha = 0.3, size = 2) +
    geom_point(aes(color = Var1), alpha = 1, size = 4) +
    geom_point(color = "#FFFFFF", size = 1) +
    coord_cartesian(ylim = c(1, nrow(rank.table))) +
    scale_y_reverse(breaks = 1:nrow(rank.table)) +
    labs(x = "", y = "Rank") +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values = line.color.vec)
  ggsave(file = save.path.name, plot = bump.plot, dpi = 350, width = 10, height = 10)
  options(warn = 0)
}

MakePercentageBoxPlot <- function(value.df, col.vec, ylim.vec, save.path.name){
  melt.df <- melt(value.df)
  options(warn = -1)
  box.plot <- ggplot(melt.df, aes(x=variable, y=value, fill = factor(variable))) + 
    geom_boxplot(width = 0.5, alpha = 0.6, show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(width = 0.25, alpha=0.9) +
    coord_cartesian(ylim = ylim.vec) +
    scale_fill_manual(values = col.vec) +
    theme(legend.position = "none") +
    labs(x = "", y = "Percentage[%]") +
    scale_y_continuous(breaks=seq(0, ylim.vec[2], 5))
  ggsave(file = save.path.name, plot = box.plot, dpi = 350, width = 10, height = 10)
  options(warn = 0)
}

MakeFoldBoxPlot <- function(value.df, col.vec, is.negative.vec, ylim.vec, save.path.name){
  value.df[,is.negative.vec] <- value.df[, is.negative.vec]^(-1)
  melt.df <- melt(value.df)
  options(warn = -1)
  box.plot <- ggplot(melt.df, aes(x=variable, y=log2(value), fill = factor(variable))) + 
    geom_boxplot(width = 0.5, alpha = 0.6, show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(width = 0.25, alpha=0.9) +
    coord_cartesian(ylim = ylim.vec) +
    scale_fill_manual(values = col.vec) +
    theme(legend.position = "none") +
    labs(x = "", y = "Fold [log2]") +
    scale_y_continuous(breaks=seq(ylim.vec[1], ylim.vec[2], 1))
  ggsave(file = save.path.name, plot = box.plot, dpi = 350, width = 10, height = 10)
  options(warn = 0)
}


MakePairedStaticalProfile <- function(vec1, vec2, save.path.name){
  res.t.test <- t.test(vec1, vec2, paired=T)
  sink(save.path.name)
  print("Paired t-test")
  if(!is.null(res.t.test)){
    print(paste0("t-statistic : ", res.t.test$statistic))
    print(paste0("degrees of freedom : ", res.t.test$parameter))
    print(paste0("p-value : ", res.t.test$p.value))
  }else{
    print("t-test was failed...")
  }
  print(paste0("mean of condition1 : ", mean(vec1)))
  print(paste0("SD of condition1 : ", sd(vec1)))
  print(paste0("mean of condition2 : ", mean(vec2)))
  print(paste0("SD of condition2 : ", sd(vec2)))
  print(paste0("pearson correlation : ", cor(vec1, vec2)))
  sink()
}

MakeGroupSampleName <- function(group.table, focus.variable, group.name, condition.label){
  name.vec <- names(group.table[,focus.variable])[group.table[,focus.variable] == group.name]
  name.vec <- paste0(name.vec, "-", condition.label)
  return(name.vec)
}

MakeRateSumTable <- function(res.file.df, table.data.path, focus.sample.name.vec){

  table.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] %in% c(condition1.label, condition2.label) & file.exists(file.path(row["aligner.res.dir.path"], table.data.path))){
      return(readRDS(file.path(row["aligner.res.dir.path"], table.data.path)))
    }else if(row["sample.label"] %in% c(condition1.label, condition2.label)){
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], table.data.path), " was not found."))
      return(NULL)
    }else{
      return(NULL)
    }
  })

  table.name.vec <- res.file.df$sample.name[!sapply(table.list, is.null)]

  table.list <- Filter(Negate(is.null), table.list)

  names(table.list) <- table.name.vec

  table.sum <- table.list[[1]]
  table.sum$Frequency <- numeric(length(table.sum$Frequency))
  for(ind in 1:length(table.list)){
    if((names(table.list)[[ind]] %in% focus.sample.name.vec) & (sum(table.list[[ind]]$Frequency) > 0)){
      table.sum$Frequency <- table.sum$Frequency +
        table.list[[ind]]$Frequency / sum(table.list[[ind]]$Frequency)
    }
  }
  return(table.sum)
}

MakeRateSumList <- function(res.file.df, table.data.path, focus.sample.name.vec){
  
  table.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] %in% c(condition1.label, condition2.label) & file.exists(file.path(row["aligner.res.dir.path"], table.data.path))){
      return(readRDS(file.path(row["aligner.res.dir.path"], table.data.path)))
    }else if(row["sample.label"] %in% c(condition1.label, condition2.label)){
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], table.data.path), " was not found."))
      return(NULL)
    }else{
      return(NULL)
    }
  })

  table.name.vec <- res.file.df$sample.name[!sapply(table.list, is.null)]

  table.list <- Filter(Negate(is.null), table.list)

  names(table.list) <- table.name.vec
  table.sum <- data.frame(variable = character(0), value = numeric(0), trim.len = numeric(0))
  for(ind in 1:length(table.list)){
    if((names(table.list)[[ind]] %in% focus.sample.name.vec) & (sum(table.list[[ind]]$value) > 0)){
      table.list[[ind]]$value <- table.list[[ind]]$value / sum(table.list[[ind]]$value)
      table.sum <- rbind(table.sum, table.list[[ind]])
    }
  }
  return(aggregate(table.sum$value, by=list(Var1 = table.sum$variable, Var2 = table.sum$trim.len), FUN=sum))
}

MakeGroupRateSumTable <- function(res.file.df, table.data.path, group.table, focus.variable, condition.label){
  # TODO:x.label.vec will be automatically set
  group.rate.sum.table <- rbind(
    High = MakeRateSumTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "High", condition.label)
    )$Frequency,
    Medium = MakeRateSumTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "Medium", condition.label)
    )$Frequency,
    Low = MakeRateSumTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "Low", condition.label)
    )$Frequency
  )

  colnames(group.rate.sum.table) <- as.character(MakeRateSumTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "High", condition.label)
    )[,1])
  return(group.rate.sum.table)
}

MakeGroupRateSumList <- function(res.file.df, table.data.path, group.table, focus.variable, condition.label){
  group.rate.sum.list <- list(
    High = MakeRateSumList(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "High", condition.label)
    ),
    Medium = MakeRateSumList(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "Medium", condition.label)
    ),
    Low = MakeRateSumList(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "Low", condition.label)
    )
  )
  return(group.rate.sum.list)
}

MakeAggregatedTable <- function(res.file.df, table.data.path, focus.sample.name.vec){

  table.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] %in% c(condition1.label, condition2.label) & file.exists(file.path(row["aligner.res.dir.path"], table.data.path))){
      return(readRDS(file.path(row["aligner.res.dir.path"], table.data.path)))
    }else if(row["sample.label"] %in% c(condition1.label, condition2.label)){
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], table.data.path), " was not found."))
      return(NULL)
    }else{
      return(NULL)
    }
  })

  table.name.vec <- res.file.df$sample.name[!sapply(table.list, is.null)]

  table.list <- Filter(Negate(is.null), table.list)

  names(table.list) <- table.name.vec

  table.sum <- data.frame(mh.len = numeric(0), trim.len = numeric(0))
  for(ind in 1:length(table.list)){
    if(names(table.list)[[ind]] %in% focus.sample.name.vec){
      table.sum <- rbind(table.sum, table.list[[ind]])
    }
  }
  return(table.sum)
}

MakeGroupAggregatedTable <- function(res.file.df, table.data.path, group.table, focus.variable, condition.label){
  group.list <- list(
    High = MakeAggregatedTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "High", condition.label)
    ),
    Medium = MakeAggregatedTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "Medium", condition.label)
    ),
    Low = MakeAggregatedTable(
      res.file.df,
      table.data.path,
      MakeGroupSampleName(group.table, focus.variable, "Low", condition.label)
    )
  )
  return(group.list)
}


SaveBarGroupPlot <- function(count.merge.table, file.name, xlab, ylab, group.name, width.margin){

  if(sum(count.merge.table["High", ]) < 1 &
    sum(count.merge.table["Low", ]) < 1 &
    sum(count.merge.table["Medium", ]) < 1){
    message("There is no data.")
    return()
  }
  
  melt.df <- melt(count.merge.table)
  count.p <- ggplot(data = melt.df, aes(x=Var2, y=value ,fill=Var1)) +
  geom_bar(aes(fill=Var1),stat="identity",position="identity") +
  scale_fill_manual(values = alpha(c("#F21905", "#04BF45", "#0487D9"), .5)) +
  geom_point(color = "black", size = 1, alpha = 1.0) +
  geom_point(aes(color = Var1), size = 0.8, alpha = 1.0) +
  coord_cartesian(ylim = c(0, max(melt.df$value))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)
      , axis.text.y = element_text(size = 5)
      , axis.title.x = element_text(size = 8)
      , axis.title.y = element_text(size = 8)
      , legend.title = element_text(size = 8)
      , legend.text = element_text(size = 8)
    ) +
  xlab(xlab) +
  ylab(ylab) +
  guides(fill=FALSE) +
  labs(color = group.name)
  options(warn = -1)
  colnames.count.merge.table <- as.numeric(colnames(count.merge.table))
  if(!is.na(colnames.count.merge.table[1])){ # is the x-axis is continuous value ?
    count.p <- count.p + scale_x_continuous(breaks = colnames.count.merge.table)
  }

  ggsave(file = file.name, plot = count.p, dpi = 350, width = width.margin * length(count.merge.table["High", ]), height = 3, limitsize = FALSE)
  options(warn = 0)
  message(paste("Save bar plot in ", file.name, sep=""))
}

SaveMultiBarGroupPlot <- function(count.merge.list, file.name, title, ylab, group.name, width.margin){
  if(length(count.merge.list$High$x) < 1 &
    length(count.merge.list$Low$x) < 1 &
    length(count.merge.list$Medium$x) < 1){
    message("There is no data.")
    return()
  }

  class.count.merge.table <- rbind(
    cbind(count.merge.list$High, class = rep("High", nrow(count.merge.list$High))),
    cbind(count.merge.list$Medium, class = rep("Medium", nrow(count.merge.list$Medium))),
    cbind(count.merge.list$Low, class = rep("Low", nrow(count.merge.list$Low)))
  )
  
  if(length(unique(class.count.merge.table$Var1)) == 4){
    col.vec <- c("#00931F", "#001298", "#000000", "#FE0000")
  }else if(length(unique(class.count.merge.table$Var1)) == 16){
    col.vec <- c("#3B7302", "#97BF04", "#6F7302", "#D9CC1E"
      ,"#009840", "#104277", "#25217D", "#005E6A"
      ,"#594011", "#8C7870", "#0D0D0D", "#401D1A"
      ,"#F28322", "#4557BF", "#400101", "#D93223")
  }
  
  options(warn=-1)
  class.count.merge.table %>% group_by(Var1) %>% ggplot(aes(x = Var2, y = x, fill = Var1)) +
  geom_bar(stat="identity", position = position_dodge(0.8)) +
  facet_grid(class ~ .) +
  scale_fill_manual(values = col.vec) +
  coord_cartesian(ylim = c(0, max(class.count.merge.table$x))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    , axis.text.y = element_text(size = 6)
    , axis.title.x = element_text(size = 12)
    , axis.title.y = element_text(size = 12)) +
  labs(title = paste0(title, " : ", group.name, " Group"), x = "Length of intervening sequence", y = ylab, fill = "Nucleotide") -> plot.p # inside frequency * number of read
  ggsave(file = file.name, plot = plot.p, dpi = 300, width = width.margin * length(class.count.merge.table$Var2), height = 6, limitsize = FALSE)
  options(warn=0)
  message(paste("Save bar plot in ", file.name, sep=""))
}


SaveScatterGroupPlot <- function(count.merge.table, file.name, xlab, group.name){

  if(length(count.merge.table$High$mh.len) < 1 &
    length(count.merge.table$Low$mh.len) < 1 &
    length(count.merge.table$Medium$mh.len) < 1){
    message("There is no data.")
    return()
  }
  
  class.count.merge.table <- rbind(
    cbind(count.merge.table$High, class = rep("High", nrow(count.merge.table$High))),
    cbind(count.merge.table$Medium, class = rep("Medium", nrow(count.merge.table$Medium))),
    cbind(count.merge.table$Low, class = rep("Low", nrow(count.merge.table$Low)))
  )
  aggregated.count.table <- aggregate(class.count.merge.table[,c("mh.len")],
    list(mh.len.class = class.count.merge.table$mh.len, trim.len.class = class.count.merge.table$trim.len, group.class = class.count.merge.table$class),
    length
  )
  
  options(warn=-1)
  count.p <- ggplot(data=aggregated.count.table, aes(mh.len.class, trim.len.class)) +
    geom_smooth(aes(color = group.class), size = 0.8, linetype = 2, se = FALSE, method = 'loess', alpha = 0.5, show.legend = FALSE) + # loess model
    geom_point(aes(size = x, fill = group.class), shape = 21) +
    scale_size_area(max_size = 12, name = "Reads")+
    scale_fill_manual(values = alpha(c("#F21905", "#04BF45", "#0487D9"), .5)) +
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    labs(x = "Length of microhomology [bp]", y = "Length of intervening Sequence [bp]", fill = group.name, title = xlab)

  ggsave(file = file.name, plot = count.p, dpi = 350, width = 12, height = 6, limitsize = FALSE)
  options(warn=0)
  message(paste("Save bar plot in ", file.name, sep=""))
}

MakePositionRateSumTablePlots <- function(table.data.path, type, file.name.without.ext.vec){

  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  ############################ Condition1 Reads

  condition1.position.by.condition1.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition1.position.by.condition1.group.table, file.name.without.ext.vec[1])
  SaveBarGroupPlot(
    condition1.position.by.condition1.group.table,
    paste0(file.name.without.ext.vec[1], ".png"),
    paste0(type, " Start position (Direction : ->)"),
    "Cumulative ratio for condition1 reads [labels]",
    "Condition1",
    0.1
  )

  condition1.position.by.condition2.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition1.position.by.condition2.group.table, file.name.without.ext.vec[2])
  SaveBarGroupPlot(
    condition1.position.by.condition2.group.table,
    paste0(file.name.without.ext.vec[2], ".png"),
    paste0(type, " Start position (Direction : ->)"),
    "Cumulative ratio for condition1 reads [labels]",
    "Condition2",
    0.1
  )

  condition1.position.by.enhancement.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition1.position.by.enhancement.group.table, file.name.without.ext.vec[3])
  SaveBarGroupPlot(
    condition1.position.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[3], ".png"),
    paste0(type, " Start position (Direction : ->)"),
    "Cumulative ratio for condition1 reads [labels]",
    "Condition2/Conditon1 in precise knock-in",
    0.1
  )

  ############################ Condition2 Reads

  condition2.position.by.condition1.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition2.position.by.condition1.group.table, file.name.without.ext.vec[4])
  SaveBarGroupPlot(
    condition2.position.by.condition1.group.table,
    paste0(file.name.without.ext.vec[4], ".png"),
    paste0(type, " Start position (Direction : ->)"),
    "Cumulative ratio for condition2 reads [labels]",
    "Condition1",
    0.1
  )

  condition2.position.by.condition2.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition2.position.by.condition2.group.table, file.name.without.ext.vec[5])
  SaveBarGroupPlot(
    condition2.position.by.condition2.group.table,
    paste0(file.name.without.ext.vec[5], ".png"),
    paste0(type, " Start position (Direction : ->)"),
    "Cumulative ratio for condition2 reads [labels]",
    "Condition2",
    0.1
  )

  condition2.position.by.enhancement.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition2.position.by.enhancement.group.table, file.name.without.ext.vec[6])
  SaveBarGroupPlot(
    condition2.position.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[6], ".png"),
    paste0(type, " Start position (Direction : ->)"),
    "Cumulative ratio for condition2 reads [labels]",
    "Condition2/Conditon1 in precise knock-in",
    0.1
  )

}

SavePeriodogramGroupPlot <- function(count.merge.table, is.negative.only = TRUE, file.name, group.name){

  if(sum(count.merge.table["High", ]) < 1 &
    sum(count.merge.table["Low", ]) < 1 &
    sum(count.merge.table["Medium", ]) < 1){
    message("There is no data.")
    return()
  }
  
  # define analysis range
  non.zero.group.table <- count.merge.table[, !apply(count.merge.table == 0, 2, any)]
  min.ind <- head(colnames(non.zero.group.table), n=1)
  if(is.negative.only){
    max.ind <- -1
  }else{
    max.ind <- tail(colnames(non.zero.group.table), n=1)
    if(is.null(max.ind)){
      message("It is insufficient amount for frequency analysis.")
      return()
    }
  }
  focus.range.vec <- as.character(min.ind:max.ind)[as.character(min.ind:max.ind) %in% colnames(count.merge.table)]
  extract.indel.group.table <- count.merge.table[, focus.range.vec]
  # calculate
  if(ncol(extract.indel.group.table) < 5){
    message("It is insufficient amount for frequency analysis.")
    return()
  }
  high.spec.list <- spec.pgram(extract.indel.group.table["High",], spans = c(3,3), plot = FALSE)
  high.spec.table <- data.frame(freq = high.spec.list$freq, spec = high.spec.list$spec, class = rep("High", length(high.spec.list$freq)))
  medium.spec.list <- spec.pgram(extract.indel.group.table["Medium",], spans = c(3,3), plot = FALSE)
  medium.spec.table <- data.frame(freq = medium.spec.list$freq, spec = medium.spec.list$spec, class = rep("Medium", length(medium.spec.list$freq)))
  low.spec.list <- spec.pgram(extract.indel.group.table["Low",], spans = c(3,3), plot = FALSE)
  low.spec.table <- data.frame(freq = low.spec.list$freq, spec = low.spec.list$spec, class = rep("Low", length(low.spec.list$freq)))
  merge.spec.table <- rbind(high.spec.table, medium.spec.table, low.spec.table)

  options(warn = -1)
  periodogram.p <- ggplot(aes(y = spec, x = freq, colour = class), data = merge.spec.table, stat="identity") +
  geom_line() +
  scale_colour_manual(values = c("#F21905", "#04BF45", "#0487D9")) +
  geom_vline(xintercept = 1/(2:(1/min(merge.spec.table$freq))), linetype = "dashed", color = "gray", size = 0.3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)
      , axis.text.y = element_text(size = 5)
      , axis.title.x = element_text(size = 8)
      , axis.title.y = element_text(size = 8)
      , legend.title = element_text(size = 8)
      , legend.text = element_text(size = 8)
    ) +
  xlab("Frequency [1/Period]") +
  ylab("Spectrum") +
  labs(color = group.name)

  ggsave(file = file.name, plot = periodogram.p, dpi = 350, width = 6, height = 3, limitsize = FALSE)
  options(warn = 0)
  message(paste("Save bar plot in ", file.name, sep=""))
}


MakeRateSumTablePlots <- function(table.data.path, x.axis.title, file.name.without.ext.vec, is.indel.size.analysis){

  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  ############################ Condition1 Reads

  condition1.indel.size.by.condition1.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition1.indel.size.by.condition1.group.table, file.name.without.ext.vec[1])
  SaveBarGroupPlot(
    condition1.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[1], ".png"),
    x.axis.title,
    "Cumulative ratio for condition1 reads [reads]",
    "Condition1",
    0.1
  )
  SavePeriodogramGroupPlot(
    condition1.indel.size.by.condition1.group.table,
    is.indel.size.analysis,
    paste0(file.name.without.ext.vec[1], "_periodogram", ".png"),
    "Condition1"
  )

  condition1.indel.size.by.condition2.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition1.indel.size.by.condition2.group.table, file.name.without.ext.vec[2])
  SaveBarGroupPlot(
    condition1.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[2], ".png"),
    x.axis.title,
    "Cumulative ratio for condition1 reads [reads]",
    "Condition2",
    0.1
  )
  SavePeriodogramGroupPlot(
    condition1.indel.size.by.condition2.group.table,
    is.indel.size.analysis,
    paste0(file.name.without.ext.vec[2], "_periodogram", ".png"),
    "Condition2"
  )

  condition1.indel.size.by.enhancement.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition1.indel.size.by.enhancement.group.table, file.name.without.ext.vec[3])
  SaveBarGroupPlot(
    condition1.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[3], ".png"),
    x.axis.title,
    "Cumulative ratio for condition1 reads [reads]",
    "Condition2/Conditon1 in precise knock-in",
    0.1
  )
  SavePeriodogramGroupPlot(
    condition1.indel.size.by.enhancement.group.table,
    is.indel.size.analysis,
    paste0(file.name.without.ext.vec[3], "_periodogram", ".png"),
    "Condition2/Conditon1 in precise knock-in"
  )

  ############################ Condition2 Reads

  condition2.indel.size.by.condition1.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition2.indel.size.by.condition1.group.table, file.name.without.ext.vec[4])
  SaveBarGroupPlot(
    condition2.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[4], ".png"),
    x.axis.title,
    "Cumulative ratio for condition2 reads [reads]",
    "Condition1",
    0.1
  )
  SavePeriodogramGroupPlot(
    condition2.indel.size.by.condition1.group.table,
    is.indel.size.analysis,
    paste0(file.name.without.ext.vec[4], "_periodogram", ".png"),
    "Condition1"
  )

  condition2.indel.size.by.condition2.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition2.indel.size.by.condition2.group.table, file.name.without.ext.vec[5])
  SaveBarGroupPlot(
    condition2.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[5], ".png"),
    x.axis.title,
    "Cumulative ratio for condition2 reads [reads]",
    "Condition2",
    0.1
  )
  SavePeriodogramGroupPlot(
    condition2.indel.size.by.condition2.group.table,
    is.indel.size.analysis,
    paste0(file.name.without.ext.vec[5], "_periodogram", ".png"),
    "Condition2"
  )

  condition2.indel.size.by.enhancement.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition2.indel.size.by.enhancement.group.table, file.name.without.ext.vec[6])
  SaveBarGroupPlot(
    condition2.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[6], ".png"),
    x.axis.title,
    "Cumulative ratio for condition2 reads [reads]",
    "Condition2/Conditon1 in precise knock-in",
    0.1
  )
  SavePeriodogramGroupPlot(
    condition2.indel.size.by.enhancement.group.table,
    is.indel.size.analysis,
    paste0(file.name.without.ext.vec[6], "_periodogram", ".png"),
    "Condition2/Conditon1 in precise knock-in"
  )

}


MakeMaxSumClassTable <- function(res.file.df, table.data.path, focus.sample.name.vec){

  table.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] %in% c(condition1.label, condition2.label) &
      file.exists(file.path(row["aligner.res.dir.path"], table.data.path))
      ){
      if(sum(readRDS(file.path(row["aligner.res.dir.path"], table.data.path))$Frequency) > 0){ # total is non-zero.
        temp.table <- readRDS(file.path(row["aligner.res.dir.path"], table.data.path))
        if(!any(colnames(temp.table) %in% c("Percentage"))){ # there is no "Percentage" label.
          temp.table <- cbind(temp.table, Percentage = temp.table$Frequency /sum(temp.table$Frequency) * 100)
        }
        if(any(colnames(temp.table) %in% c("label"))){ # there is "label" label.
          colnames(temp.table)[which(colnames(temp.table) %in% "label")] <- "Class" # rename
        }
        return(temp.table)
      }else{
        message(paste0("Ignore ", row["sample.name"], " because each value is zero."))
        return(NULL)
      }
    }else if(row["sample.label"] %in% c(condition1.label, condition2.label)){
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], table.data.path), " was not found."))
      return(NULL)
    }else{
      return(NULL)
    }
  })

  table.name.vec <- res.file.df$sample.name[!sapply(table.list, is.null)]

  table.list <- Filter(Negate(is.null), table.list)

  names(table.list) <- table.name.vec
  
  table.sum <- table.list[[1]]
  table.sum$Frequency <- numeric(length(table.sum$Frequency))
  table.sum$Percentage <- numeric(length(table.sum$Percentage))
  for(ind in 1:length(table.list)){
    if((names(table.list)[[ind]] %in% focus.sample.name.vec) & (sum(table.list[[ind]]$Frequency) > 0)){
      
      max.ind <- which(table.list[[ind]]$Frequency == max(table.list[[ind]]$Frequency))
      table.sum$Frequency[max.ind] <- table.sum$Frequency[max.ind] + 1
    }
  }
  table.sum$Percentage <- table.sum$Frequency / sum(table.sum$Frequency) * 100
  return(table.sum)
}

MakeReadsSumClassTable <- function(res.file.df, table.data.path, focus.sample.name.vec){

  table.list <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] %in% c(condition1.label, condition2.label) &
      file.exists(file.path(row["aligner.res.dir.path"], table.data.path))
      ){
      if(sum(readRDS(file.path(row["aligner.res.dir.path"], table.data.path))$Frequency) > 0){ # total is non-zero.
        temp.table <- readRDS(file.path(row["aligner.res.dir.path"], table.data.path))
        if(!any(colnames(temp.table) %in% c("Percentage"))){ # there is no "Percentage" label.
          temp.table <- cbind(temp.table, Percentage = temp.table$Frequency /sum(temp.table$Frequency) * 100)
        }
        if(any(colnames(temp.table) %in% c("label"))){ # there is "label" label.
          colnames(temp.table)[which(colnames(temp.table) %in% "label")] <- "Class" # rename
        }
        return(temp.table)
      }else{
        message(paste0("Ignore ", row["sample.name"], " because each value is zero."))
        return(NULL)
      }
    }else if(row["sample.label"] %in% c(condition1.label, condition2.label)){
      message(paste0("Ignore ", row["sample.name"], " because ", file.path(row["aligner.res.dir.path"], table.data.path), " was not found."))
      return(NULL)
    }else{
      return(NULL)
    }
  })

  table.name.vec <- res.file.df$sample.name[!sapply(table.list, is.null)]

  table.list <- Filter(Negate(is.null), table.list)

  names(table.list) <- table.name.vec
  
  table.sum <- table.list[[1]]
  table.sum$Frequency <- numeric(length(table.sum$Frequency))
  table.sum$Percentage <- numeric(length(table.sum$Percentage))
  for(ind in 1:length(table.list)){
    if((names(table.list)[[ind]] %in% focus.sample.name.vec) & (sum(table.list[[ind]]$Frequency) > 0)){
      
      table.sum$Frequency <- table.sum$Frequency +
        table.list[[ind]]$Frequency #/ sum(table.list[[ind]]$Frequency)

    }
  }
  table.sum$Percentage <- table.sum$Frequency / sum(table.sum$Frequency) * 100
  return(table.sum)
}


MakeMaxSumInDelTablePlots <- function(table.data.path, file.name.without.ext.vec,
  type = c("indel.size.indel", "ej.indel.size.indel", "indel.size.impreciseki", "ej.indel.size.impreciseki", "id-ip.ej.indel.size.impreciseki")){
  type <- match.arg(type)
  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  if(type == "ej.indel.size.indel"){
    col.vec <- c("#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304", "#D99ABC")
  }else if(type == "indel.size.indel"){
    col.vec <- c("#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304")
  }else if(type == "indel.size.impreciseki"){
    col.vec <- c("#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
  }else if(type == "ej.indel.size.impreciseki"){
    # col.vec <- c("#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
    col.vec <- c("#C0D904", "#277A8C", "#3E5902", "#60A6A6", "#618C03", "#95D904", "#F2EC9B", "#B3D9C0")
  }else if(type == "id-ip.ej.indel.size.impreciseki"){
    col.vec <- c("#ACF2F2", "#F2A0DC", "#D288F2", "#F2C9F0")
  }
  
  ############################ Condition1 Reads

  ### condition1.precise.knockin group

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "High")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[1], "_high_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Medium")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[1], "_medium_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Low")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[1], "_low_group.png")
  )

  ### condition2.precise.knockin group

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "High")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[2], "_high_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[2], "_medium_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Low")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[2], "_low_group.png")
  )

  ### enhancement.1.2.precise.knockin group

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "High")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[3], "_high_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[3], "_medium_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Low")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[3], "_low_group.png")
  )

  ############################ Condition2 Reads

  ### condition1.precise.knockin group

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "High")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[4], "_high_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Medium")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[4], "_medium_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Low")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[4], "_low_group.png")
  )

  ### condition2.precise.knockin group

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "High")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[5], "_high_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[5], "_medium_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Low")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[5], "_low_group.png")
  )

  ### enhancement.1.2.precise.knockin group

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "High")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[6], "_high_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[6], "_medium_group.png")
  )

  SavePieChart(
    MakeMaxSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Low")),
    col.vec, "Genes", paste0(file.name.without.ext.vec[6], "_low_group.png")
  )
}

MakeReadsSumInDelTablePlots <- function(table.data.path, file.name.without.ext.vec,
  type = c("indel.size.indel", "ej.indel.size.indel", "indel.size.impreciseki", "ej.indel.size.impreciseki", "id-ip.ej.indel.size.impreciseki")){
  type <- match.arg(type)

  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  if(type == "ej.indel.size.indel"){
    col.vec <- c("#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304", "#D99ABC")
  }else if(type == "indel.size.indel"){
    col.vec <- c("#F2ECD8", "#BF9E75", "#8C4A32", "#D9CC14", "#8C8304")
  }else if(type == "indel.size.impreciseki"){
    col.vec <- c("#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
  }else if(type == "ej.indel.size.impreciseki"){
    # col.vec <- c("#C0D904", "#3E5902", "#618C03", "#95D904", "#F2EC9B", "#277A8C", "#60A6A6", "#B3D9C0")
    col.vec <- c("#C0D904", "#277A8C", "#3E5902", "#60A6A6", "#618C03", "#95D904", "#F2EC9B", "#B3D9C0")
  }else if(type == "id-ip.ej.indel.size.impreciseki"){
    col.vec <- c("#ACF2F2", "#F2A0DC", "#D288F2", "#F2C9F0")
  }
  
  ############################ Condition1 Reads

  ### condition1.precise.knockin group

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "High")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[1], "_high_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Medium")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[1], "_medium_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Low")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[1], "_low_group.png")
  )

  ### condition2.precise.knockin group

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "High")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[2], "_high_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[2], "_medium_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Low")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[2], "_low_group.png")
  )

  ### enhancement.1.2.precise.knockin group

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "High")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[3], "_high_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[3], "_medium_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition1.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Low")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[3], "_low_group.png")
  )

  ############################ Condition2 Reads

  ### condition1.precise.knockin group

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "High")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[4], "_high_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Medium")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[4], "_medium_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition1.precise.knockin",
        group.name = "Low")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[4], "_low_group.png")
  )

  ### condition2.precise.knockin group

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "High")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[5], "_high_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[5], "_medium_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "condition2.precise.knockin",
        group.name = "Low")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[5], "_low_group.png")
  )

  ### enhancement.1.2.precise.knockin group

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "High")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[6], "_high_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Medium")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[6], "_medium_group.png")
  )

  SavePieChart(
    MakeReadsSumClassTable(
      fltr.res.file.df, table.data.path,
      MakeGroupSampleName(
        group.table = machiato.group.table,
        condition.label = condition2.label,
        focus.variable = "enhancement.1.2.precise.knockin",
        group.name = "Low")),
    col.vec, "Reads", paste0(file.name.without.ext.vec[6], "_low_group.png")
  )
}

MakeScatterPlots <- function(table.data.path, x.axis.title, file.name.without.ext.vec){

  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  ############################ Condition1 Reads

  condition1.indel.size.by.condition1.group.table <- MakeGroupAggregatedTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition1.indel.size.by.condition1.group.table, file.name.without.ext.vec[1])
  SaveScatterGroupPlot(
    condition1.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[1], ".png"),
    x.axis.title,
    "Condition1"
  )

  condition1.indel.size.by.condition2.group.table <- MakeGroupAggregatedTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition1.indel.size.by.condition2.group.table, file.name.without.ext.vec[2])
  SaveScatterGroupPlot(
    condition1.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[2], ".png"),
    x.axis.title,
    "Condition2"
  )

  condition1.indel.size.by.enhancement.group.table <- MakeGroupAggregatedTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition1.indel.size.by.enhancement.group.table, file.name.without.ext.vec[3])
  SaveScatterGroupPlot(
    condition1.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[3], ".png"),
    x.axis.title,
    "Condition2/Conditon1 in precise knock-in"
  )

  ############################ Condition2 Reads

  condition2.indel.size.by.condition1.group.table <- MakeGroupAggregatedTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition2.indel.size.by.condition1.group.table, file.name.without.ext.vec[4])
  SaveScatterGroupPlot(
    condition2.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[4], ".png"),
    x.axis.title,
    "Condition1"
  )

  condition2.indel.size.by.condition2.group.table <- MakeGroupAggregatedTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition2.indel.size.by.condition2.group.table, file.name.without.ext.vec[5])
  SaveScatterGroupPlot(
    condition2.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[5], ".png"),
    x.axis.title,
    "Condition2"
  )

  condition2.indel.size.by.enhancement.group.table <- MakeGroupAggregatedTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition2.indel.size.by.enhancement.group.table, file.name.without.ext.vec[6])
  SaveScatterGroupPlot(
    condition2.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[6], ".png"),
    x.axis.title,
    "Condition2/Conditon1 in precise knock-in"
  )

}


MakeSeqRateSumTablePlots <- function(table.data.path, x.axis.title, y.axis.unit, file.name.without.ext.vec, width.margin.1, width.margin.2){

  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  ############################ Condition1 Reads

  condition1.indel.size.by.condition1.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition1.indel.size.by.condition1.group.table, file.name.without.ext.vec[1])
  SaveBarGroupPlot(
    condition1.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[1], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition1",
    width.margin.1
  )

  condition1.indel.size.by.condition2.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition1.indel.size.by.condition2.group.table, file.name.without.ext.vec[2])
  SaveBarGroupPlot(
    condition1.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[2], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2",
    width.margin.1
  )

  condition1.indel.size.by.enhancement.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition1.indel.size.by.enhancement.group.table, file.name.without.ext.vec[3])
  SaveBarGroupPlot(
    condition1.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[3], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2/Conditon1 in precise knock-in",
    width.margin.2
  )

  ############################ Condition2 Reads

  condition2.indel.size.by.condition1.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition2.indel.size.by.condition1.group.table, file.name.without.ext.vec[4])
  SaveBarGroupPlot(
    condition2.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[4], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition1",
    width.margin.1
  )

  condition2.indel.size.by.condition2.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition2.indel.size.by.condition2.group.table, file.name.without.ext.vec[5])
  SaveBarGroupPlot(
    condition2.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[5], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2",
    width.margin.1
  )

  condition2.indel.size.by.enhancement.group.table <- MakeGroupRateSumTable(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition2.indel.size.by.enhancement.group.table, file.name.without.ext.vec[6])
  SaveBarGroupPlot(
    condition2.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[6], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2/Conditon1 in precise knock-in",
    width.margin.2
  )

}

MakeSeqLengthRateSumTablePlots <- function(table.data.path, x.axis.title, y.axis.unit, file.name.without.ext.vec, width.margin.1, width.margin.2){

  ### NOTE : Global variables is used ###.
  # fltr.res.file.df
  # machiato.group.table
  # condition1.label
  # condition2.label

  ############################ Condition1 Reads

  condition1.indel.size.by.condition1.group.table <- MakeGroupRateSumList(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition1.indel.size.by.condition1.group.table, file.name.without.ext.vec[1])
  SaveMultiBarGroupPlot(
    condition1.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[1], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition1",
    width.margin.1
  )

  condition1.indel.size.by.condition2.group.table <- MakeGroupRateSumList(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition1.indel.size.by.condition2.group.table, file.name.without.ext.vec[2])
  SaveMultiBarGroupPlot(
    condition1.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[2], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2",
    width.margin.1
  )

  condition1.indel.size.by.enhancement.group.table <- MakeGroupRateSumList(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition1.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition1.indel.size.by.enhancement.group.table, file.name.without.ext.vec[3])
  SaveMultiBarGroupPlot(
    condition1.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[3], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2/Conditon1 in precise knock-in",
    width.margin.2
  )

  ############################ Condition2 Reads

  condition2.indel.size.by.condition1.group.table <- MakeGroupRateSumList(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition1.precise.knockin")
  SaveTable(condition2.indel.size.by.condition1.group.table, file.name.without.ext.vec[4])
  SaveMultiBarGroupPlot(
    condition2.indel.size.by.condition1.group.table,
    paste0(file.name.without.ext.vec[4], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition1",
    width.margin.1
  )

  condition2.indel.size.by.condition2.group.table <- MakeGroupRateSumList(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "condition2.precise.knockin")
  SaveTable(condition2.indel.size.by.condition2.group.table, file.name.without.ext.vec[5])
  SaveMultiBarGroupPlot(
    condition2.indel.size.by.condition2.group.table,
    paste0(file.name.without.ext.vec[5], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2",
    width.margin.1
  )

  condition2.indel.size.by.enhancement.group.table <- MakeGroupRateSumList(fltr.res.file.df, table.data.path, machiato.group.table, condition.label = condition2.label, focus.variable = "enhancement.1.2.precise.knockin")
  SaveTable(condition2.indel.size.by.enhancement.group.table, file.name.without.ext.vec[6])
  SaveMultiBarGroupPlot(
    condition2.indel.size.by.enhancement.group.table,
    paste0(file.name.without.ext.vec[6], ".png"),
    x.axis.title,
    paste0("Cumulative ratio for condition1 reads [", y.axis.unit, "]"),
    "Condition2/Conditon1 in precise knock-in",
    width.margin.2
  )

}

MakeHeatMap <- function(predicted.data.list, detected.data.list, y.axis.title, x.axis.title, file.name.excluding.ext, condition.label){
  
  # make matrix for plot
  
  predicted.comp.kl.mat <- matrix(NA, ncol = length(predicted.data.list))
  colnames(predicted.comp.kl.mat) <- names(predicted.data.list)
  rownames(predicted.comp.kl.mat) <- "temp"
  for(ind in 1:length(detected.data.list)){
    predicted.comp.kl.vec <- numeric(0)
    machiato.nhej.crisprset.name <- names(detected.data.list[[ind]])
    #if(substr(machiato.nhej.crisprset.name, nchar(machiato.nhej.crisprset.name), nchar(machiato.nhej.crisprset.name)) %in% c(condition2.label)){
    #  next
    #}
    for(temp.name in names(predicted.data.list)){
      predicted.comp.kl.vec <- c(predicted.comp.kl.vec, CalcKLdivergence2(getElement(detected.data.list[[ind]], machiato.nhej.crisprset.name), getElement(predicted.data.list, temp.name)))
      if(length(predicted.comp.kl.vec) < 2){
        names(predicted.comp.kl.vec) <- c(temp.name)
      }else{
        names(predicted.comp.kl.vec) <- c(names(predicted.comp.kl.vec)[1:length(predicted.comp.kl.vec)-1], temp.name)
      }
    }
    predicted.comp.kl.mat <- rbind(predicted.comp.kl.mat, predicted.comp.kl.vec)
    # just below sentence was written due to misunderstanding relationship between colnames and rownames.
    # actually, i should change colnames because colnames mean predicted score. but i changed rownames.
    # i need a few change to fix it. but it works correctly without dealing with it. So it remains.
    rownames(predicted.comp.kl.mat) <- c(rownames(predicted.comp.kl.mat)[1:nrow(predicted.comp.kl.mat)-1], gsub(paste0("-", condition.label, "$"), "", c(machiato.nhej.crisprset.name)))
  }
  # rownames(predicted.comp.kl.mat) # MaChIAto
  # colnames(predicted.comp.kl.mat) # Predicted Score
  predicted.comp.kl.mat <- predicted.comp.kl.mat[-1, ]
  predicted.comp.kl.mat <- predicted.comp.kl.mat[, which(colnames(predicted.comp.kl.mat) %in% paste0(rownames(predicted.comp.kl.mat), "-", condition.label))]

  reorder.ind.vec <- unlist(sapply(paste0(rownames(predicted.comp.kl.mat), "-", condition.label), function(x){
    return(which(colnames(predicted.comp.kl.mat) %in% x))
  }))
  predicted.comp.kl.mat <- predicted.comp.kl.mat[, reorder.ind.vec]
  
  predicted.comp.kl.p.mat <- predicted.comp.kl.mat # just renamed same matrix
  colnames(predicted.comp.kl.p.mat) <- gsub(paste0("-", condition.label, "$"), "", colnames(predicted.comp.kl.p.mat))
  
  # make heat map
  options(warn=-1)
  indelphi.comp.kl.p <- ggplot(melt(log2(predicted.comp.kl.p.mat[,ncol(predicted.comp.kl.p.mat):1]))
    , aes(x=Var1, y=Var2, fill=value)) +
  geom_raster(aes(fill = value)) +
  labs(title ="Symmetrized KL divergence"
    , x = x.axis.title
    , y = y.axis.title) +
  theme(axis.ticks = element_blank()
    , axis.text.x = element_text(angle = 270, hjust = 0, size=7)
    , axis.text.y = element_text(size=7)
    , plot.title = element_text(size = 10)
    , axis.title.x = element_text(size = 10)
    , axis.title.y = element_text(size = 10)
  ) +
  scale_fill_distiller(palette = "Blues", name = "Symmetrized KL Divergence (log2)", limits = c(min(log2(predicted.comp.kl.p.mat)), max(log2(predicted.comp.kl.p.mat))))

  SaveTable(predicted.comp.kl.p.mat, file.name.excluding.ext)
  ggsave(file = paste0(file.name.excluding.ext, ".png")
    , plot = indelphi.comp.kl.p, dpi = 350, width = 11, height = 10)


  # make boxplot
  match.value.vec <- numeric(0)
  nonmatch.value.vec <- numeric(0)
  for(row.ind in 1:nrow(predicted.comp.kl.mat)){
    for(col.ind in 1:ncol(predicted.comp.kl.mat)){
      if(row.ind == col.ind){
        match.value.vec <- c(match.value.vec, log2(predicted.comp.kl.mat)[row.ind, col.ind])
      }else{
        nonmatch.value.vec <- c(nonmatch.value.vec, log2(predicted.comp.kl.mat)[row.ind, col.ind])
      }
    }
  }

  boxplot.df <- data.frame(variable = c(rep("corresponding.target", length(match.value.vec)), rep("noncorresponding.target", length(nonmatch.value.vec))),
    value = c(match.value.vec, nonmatch.value.vec))
  
  options(warn = -1)
  ylim.vec <- c(min(boxplot.df$value), max(boxplot.df$value))
  box.plot <- ggplot(boxplot.df, aes(x=variable, y=value, fill = factor(variable))) + 
    geom_boxplot(width = 0.5, alpha = 0.8, show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(width = 0.05, alpha=0.4, colour = "#F29F05") +
    coord_cartesian(ylim = ylim.vec) +
    scale_fill_manual(values = c("#0528F2", "#F24405")) +
    theme(axis.ticks = element_blank()
      , axis.text.x = element_text(size=15)
      , axis.text.y = element_text(size=15)
      , axis.title.y = element_text(size = 20)
      , legend.position = "none"
    ) +
    labs(x = "", y = "Symmetrized KL Divergence (log2)") +
    scale_y_continuous(breaks=seq(round(ylim.vec[1] - 1), round(ylim.vec[2] + 1), 1))
  
  ggsave(file = paste0(file.name.excluding.ext, "_boxplot.png"), plot = box.plot, dpi = 350, width = 10, height = 10)
  options(warn = 0)

  # run test
  res.t.test <- t.test(match.value.vec, nonmatch.value.vec, paired = FALSE)
  sink(paste(file.name.excluding.ext, "_ttest.txt"))
  print("t-test")
  if(!is.null(res.t.test)){
    print(paste0("t-statistic : ", res.t.test$statistic))
    print(paste0("degrees of freedom : ", res.t.test$parameter))
    print(paste0("p-value : ", res.t.test$p.value))
  }else{
    print("t-test was failed...")
  }
  print(paste0("number of corresponding.target : ", length(match.value.vec)))
  print(paste0("mean of corresponding.target : ", mean(match.value.vec)))
  print(paste0("SD of corresponding.target : ", sd(match.value.vec)))
  print(paste0("number of noncorresponding.target : ", length(nonmatch.value.vec)))
  print(paste0("mean of noncorresponding.target : ", mean(nonmatch.value.vec)))
  print(paste0("SD of noncorresponding.target : ", sd(nonmatch.value.vec)))
  sink()

  options(warn=0)
}

MakeAlluvialDiagramsCRISPRessoMaChIAto <- function(){
  ### Make Read Total all table
  all.dataframe.reads.list <- apply(fltr.res.file.df, MARGIN = 1, function(row){
    if(file.exists(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"))){
      temp.all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
      part.name.vec <- c(
        "NHEJ",
        "UNMODIFIED",
        "HDR",
        "X.Reads",
        "CRISPResso_reclassification_labels"
      )
      temp.part.df <- temp.all.df[, part.name.vec]
      return(aggregate(X.Reads ~ ., data = temp.part.df, sum))
    }else{
      return(NULL)
    }
  })
  all.dataframe.reads.list <- Filter(Negate(is.null), all.dataframe.reads.list) # return list without NULL element
  total.all.reads.dataframe <- all.dataframe.reads.list[[1]]
  for(ind in 2:length(all.dataframe.reads.list)){
    total.all.reads.dataframe <- rbind(total.all.reads.dataframe, all.dataframe.reads.list[[ind]])
    total.all.reads.dataframe <- aggregate(X.Reads ~ ., data = total.all.reads.dataframe, sum)
  }
  # add CRISPResso_original_classification_labels
  total.all.reads.dataframe <- rbind(
    cbind(filter(total.all.reads.dataframe, UNMODIFIED == "True"), CRISPResso_original_classification_labels = "Unmodified"),
    cbind(filter(total.all.reads.dataframe, NHEJ == "True"), CRISPResso_original_classification_labels = "NHEJ"),
    cbind(filter(total.all.reads.dataframe, UNMODIFIED == "False" & NHEJ == "False" & HDR == "False"), CRISPResso_original_classification_labels = "Mixed HDR-NHEJ"),
    cbind(filter(total.all.reads.dataframe, HDR == "True"), CRISPResso_original_classification_labels = "HDR")
  )

  ############################################

  ### Make Variants Total all table
  all.dataframe.variants.list <- apply(fltr.res.file.df, MARGIN = 1, function(row){
    if(file.exists(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"))){
      temp.all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
      temp.all.df <- cbind(temp.all.df, X.Variants = 1)
      part.name.vec <- c(
        "NHEJ",
        "UNMODIFIED",
        "HDR",
        "X.Variants",
        "CRISPResso_reclassification_labels"
      )
      temp.part.df <- temp.all.df[, part.name.vec]
      return(aggregate(X.Variants ~ ., data = temp.part.df, sum))
    }else{
      return(NULL)
    }
  })
  all.dataframe.variants.list <- Filter(Negate(is.null), all.dataframe.variants.list) # return list without NULL element
  total.all.variants.dataframe <- all.dataframe.variants.list[[1]]
  for(ind in 2:length(all.dataframe.variants.list)){
    total.all.variants.dataframe <- rbind(total.all.variants.dataframe, all.dataframe.variants.list[[ind]])
    total.all.variants.dataframe <- aggregate(X.Variants ~ ., data = total.all.variants.dataframe, sum)
  }
  # add CRISPResso_original_classification_labels
  total.all.variants.dataframe <- rbind(
    cbind(filter(total.all.variants.dataframe, UNMODIFIED == "True"), CRISPResso_original_classification_labels = "Unmodified"),
    cbind(filter(total.all.variants.dataframe, NHEJ == "True"), CRISPResso_original_classification_labels = "NHEJ"),
    cbind(filter(total.all.variants.dataframe, UNMODIFIED == "False" & NHEJ == "False" & HDR == "False"), CRISPResso_original_classification_labels = "Mixed HDR-NHEJ"),
    cbind(filter(total.all.variants.dataframe, HDR == "True"), CRISPResso_original_classification_labels = "HDR")
  )



  ### plot
  melt.template.table <- melt(table(total.all.variants.dataframe[, c("CRISPResso_reclassification_labels", "CRISPResso_original_classification_labels")]))
  colnames(melt.template.table) <- c("Var2", "Var1", "value")

  # Make variants count table
  total.all.variant.dataframe <- aggregate(X.Variants ~ ., data = total.all.variants.dataframe[c("X.Variants", "CRISPResso_reclassification_labels", "CRISPResso_original_classification_labels")], sum)
  total.all.dataframe.variantcnt.melt.table <- melt.template.table
  for(var2 in unique(total.all.dataframe.variantcnt.melt.table$Var2)){
    for(var1 in unique(total.all.dataframe.variantcnt.melt.table$Var1)){
      temp.value <- filter(total.all.variant.dataframe, CRISPResso_original_classification_labels == var1 & CRISPResso_reclassification_labels == var2)$X.Variants
      if(length(temp.value) > 0){ # this is number
        total.all.dataframe.variantcnt.melt.table[total.all.dataframe.variantcnt.melt.table$Var1 == var1 & total.all.dataframe.variantcnt.melt.table$Var2 == var2, ]$value <- temp.value
      }else{
        total.all.dataframe.variantcnt.melt.table[total.all.dataframe.variantcnt.melt.table$Var1 == var1 & total.all.dataframe.variantcnt.melt.table$Var2 == var2, ]$value <- 0
      }
    }
  }
  SaveTable(total.all.dataframe.variantcnt.melt.table, file.path(comparison.analysis.dir, "[extra]Alluvial_Diagrams_of_CRISPResso_MaChIAto_Reads"))
  MakeAlluvialDiagrams(total.all.dataframe.variantcnt.melt.table,
    c("#F20C36", "#A60321", "#548C1C", "#F2E2CE", "gray50", "#F20C36", "#A60321", "#548C1C", "#F2E2CE"),
    file.path(comparison.analysis.dir, "[extra]Alluvial_Diagrams_of_CRISPResso_MaChIAto_Variants.png"),
    "CRISPResso Classification", "MaChIAto Classification", "all allele",
    "kinds")

  # Make reads count table
  total.all.read.dataframe <- aggregate(X.Reads ~ ., data = total.all.reads.dataframe[c("X.Reads", "CRISPResso_reclassification_labels", "CRISPResso_original_classification_labels")], sum)
  total.all.dataframe.readcnt.melt.table <- melt.template.table
  for(var2 in unique(total.all.dataframe.readcnt.melt.table$Var2)){
    for(var1 in unique(total.all.dataframe.readcnt.melt.table$Var1)){
      temp.value <- filter(total.all.read.dataframe, CRISPResso_original_classification_labels == var1 & CRISPResso_reclassification_labels == var2)$X.Reads
      if(length(temp.value) > 0){ # this is number
        total.all.dataframe.readcnt.melt.table[total.all.dataframe.readcnt.melt.table$Var1 == var1 & total.all.dataframe.readcnt.melt.table$Var2 == var2, ]$value <- temp.value
      }else{
        total.all.dataframe.readcnt.melt.table[total.all.dataframe.readcnt.melt.table$Var1 == var1 & total.all.dataframe.readcnt.melt.table$Var2 == var2, ]$value <- 0
      }
    }
  }
  SaveTable(total.all.dataframe.readcnt.melt.table, file.path(comparison.analysis.dir, "[extra]Alluvial_Diagrams_of_CRISPResso_MaChIAto_Reads"))
  MakeAlluvialDiagrams(total.all.dataframe.readcnt.melt.table,
    c("#F20C36", "#A60321", "#548C1C", "#F2E2CE", "gray50", "#F20C36", "#A60321", "#548C1C", "#F2E2CE"),
    file.path(comparison.analysis.dir, "[extra]Alluvial_Diagrams_of_CRISPResso_MaChIAto_Reads.png"),
    "CRISPResso Classification", "MaChIAto Classification", "all allele",
    "reads")
}