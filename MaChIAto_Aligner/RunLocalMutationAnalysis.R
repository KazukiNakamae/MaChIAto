
RunLocalMutationAnalysis <- function(bam.datatype.dir, result.datatype.dir, index.dir, seq.dir){

  ### Load sequence
  wt.seq <- as.character(readDNAStringSet(file.path(index.dir, "wt.fa")))
  protospacer.seq <- as.character(readDNAStringSet(file.path(seq.dir, "sgRNA.fa")))
  donor.seq <- as.character(readDNAStringSet(file.path(seq.dir, "donor.fa")))
  inleft.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "In-Left_indicator_sequence.fa")))
  inright.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "In-Right_indicator_sequence.fa")))
  left.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Leftside_homology_arm.fa")))
  outleft.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Out-Left_indicator_sequence.fa")))
  outright.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Out-Right_indicator_sequence.fa")))
  right.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Rightside_homology_arm.fa")))
  untreated.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Untreated_sequence.fa")))

  # bam dir
  bam.muttype.dir <- file.path(bam.datatype.dir, "mut_type")
  bam.other_type.dir <- file.path(bam.datatype.dir, "other_type")
  # result dir
  alignment.analysis.dir <- file.path(result.datatype.dir, "[Ⅰ]alignment_analysis")
  untreated.dir <- file.path(alignment.analysis.dir, "untreated")
  deletion.dir <- file.path(alignment.analysis.dir, "deletion")
  substitution.dir <- file.path(alignment.analysis.dir, "substitution")
  insertion.dir <- file.path(alignment.analysis.dir, "insertion")
  unexpected_mutations.dir <- file.path(alignment.analysis.dir, "unexpected_mutations")
  mutation.analysis.dir <- file.path(result.datatype.dir, "[Ⅱa]analysis_of_mutation_from_endogenous_locus")

  message("Make directory for saving the analysis data")
  for(dirname in c(bam.muttype.dir, alignment.analysis.dir, untreated.dir, deletion.dir, bam.other_type.dir
    , substitution.dir, insertion.dir, unexpected_mutations.dir
    , mutation.analysis.dir)){
    dir.create(file.path(dirname), showWarnings=FALSE)
  }
  ###################################################################################################################################################

  message("Start local mutation alignment...")
  ### Mutant
  ### untreated
  if(!file.exists(file.path(untreated.dir, "crispr.set.rds"))){
    untreated.crispr.set <- MakeWtCrisprSet(
    untreated.seq
    , "wt"
    , wt.seq
    , protospacer.seq
    , file.path(bam.muttype.dir, "untreated.bam")
    , "untreated")
    if(!is.null(untreated.crispr.set))SaveVariantsData(untreated.crispr.set, untreated.dir, range.vec = -35:35, mut.type = "mut")
  }else{
    untreated.crispr.set <- readRDS(file.path(untreated.dir, "crispr.set.rds"))
  }

  ### deletion
  if(!file.exists(file.path(deletion.dir, "crispr.set.rds"))){
    overall_edited_deletion.crispr.set <- MakeWtCrisprSet(
      untreated.seq
      , "wt"
      , wt.seq
      , protospacer.seq
      , file.path(bam.muttype.dir, "overall_edited_deletion.bam")
      , "overall_edited_deletion")
    if(!is.null(overall_edited_deletion.crispr.set))SaveVariantsData(overall_edited_deletion.crispr.set, deletion.dir, range.vec = -35:35, mut.type = "mut")
  }else{
    overall_edited_deletion.crispr.set <- readRDS(file.path(deletion.dir, "crispr.set.rds"))
  }

  ### substitution
  if(!file.exists(file.path(substitution.dir, "crispr.set.rds"))){
    overall_edited_substitution.crispr.set <- MakeWtCrisprSet(
      untreated.seq
      , "wt"
      , wt.seq
      , protospacer.seq
      , file.path(bam.muttype.dir, "overall_edited_substitution.bam")
      , "overall_edited_substitution")
    if(!is.null(overall_edited_substitution.crispr.set))SaveVariantsData(overall_edited_substitution.crispr.set, substitution.dir, range.vec = -35:35, mut.type = "mut")
  }else{
    overall_edited_substitution.crispr.set <- readRDS(file.path(substitution.dir, "crispr.set.rds"))
  }

  ### insertion
  if(!file.exists(file.path(insertion.dir, "crispr.set.rds"))){
    overall_edited_insertion.crispr.set <- MakeWtCrisprSet(
      untreated.seq
      , "wt"
      , wt.seq
      , protospacer.seq
      , file.path(bam.muttype.dir, "overall_edited_insertion.bam")
      , "overall_edited_insertion")
    if(!is.null(overall_edited_insertion.crispr.set))SaveVariantsData(overall_edited_insertion.crispr.set, insertion.dir, range.vec = -35:35, mut.type = "mut")
  }else{
    overall_edited_insertion.crispr.set <- readRDS(file.path(insertion.dir, "crispr.set.rds"))
  }

  ### unexpected_mutations
  if(!file.exists(file.path(unexpected_mutations.dir, "crispr.set.rds"))){
    unexpected_mutations.crispr.set <- MakeWtCrisprSet(
      untreated.seq
      , "wt"
      , wt.seq
      , protospacer.seq
      , file.path(bam.muttype.dir, "unexpected_mutations.bam")
      , "unexpected_mutations")
    if(!is.null(unexpected_mutations.crispr.set))SaveVariantsData(unexpected_mutations.crispr.set, unexpected_mutations.dir, range.vec = -35:35, mut.type = "mut")
  }else{
    unexpected_mutations.crispr.set <- readRDS(file.path(unexpected_mutations.dir, "crispr.set.rds"))
  }

  ### aggregated mutation alasysis
  message("Run aggregated_mutation_position_analysis...")

  # untreated/substitution
  untreated.variants.count.data <- MakeTotalMutList(untreated.crispr.set, "mut")
  untreated.substitution.table <- data.frame(MakeMutPositionLabelTable(untreated.variants.count.data, mut.type = "substitution", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(untreated.substitution.table) <- c("Position", "untreated.substitution.Freq")

  overall_edited_substitution.variants.count.data <- MakeTotalMutList(overall_edited_substitution.crispr.set, "mut")
  # overall_edited_substitution/substitution
  overall_edited_substitution.substitution.table <- data.frame(MakeMutPositionLabelTable(overall_edited_substitution.variants.count.data, mut.type = "substitution", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_substitution.substitution.table) <- c("Position", "overall_edited_substitution.substitution.Freq")
  # overall_edited_substitution/deletion
  overall_edited_substitution.deletion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_substitution.variants.count.data, mut.type = "deletion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_substitution.deletion.table) <- c("Position", "overall_edited_substitution.deletion.Freq")
  # overall_edited_substitution/no.overlap.deletion
  overall_edited_substitution.no.overlap.deletion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_substitution.variants.count.data, mut.type = "no.overlap.deletion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_substitution.no.overlap.deletion.table) <- c("Position", "overall_edited_substitution.no.overlap.deletion.Freq")
  # overall_edited_substitution/insertion
  overall_edited_substitution.insertion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_substitution.variants.count.data, mut.type = "insertion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_substitution.insertion.table) <- c("Position", "overall_edited_substitution.insertion.Freq")

  overall_edited_deletion.variants.count.data <- MakeTotalMutList(overall_edited_deletion.crispr.set, "mut")
  # overall_edited_deletion/deletion
  overall_edited_deletion.deletion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_deletion.variants.count.data, mut.type = "deletion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_deletion.deletion.table) <- c("Position", "overall_edited_deletion.deletion.Freq")
  # overall_edited_deletion/no.overlap.deletion
  overall_edited_deletion.no.overlap.deletion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_deletion.variants.count.data, mut.type = "no.overlap.deletion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_deletion.no.overlap.deletion.table) <- c("Position", "overall_edited_deletion.no.overlap.deletion.Freq")
  # overall_edited_deletion/insertion
  overall_edited_deletion.insertion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_deletion.variants.count.data, mut.type = "insertion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_deletion.insertion.table) <- c("Position", "overall_edited_deletion.insertion.Freq")

  overall_edited_insertion.variants.count.data <- MakeTotalMutList(overall_edited_insertion.crispr.set, "mut")
  # overall_edited_insertion/deletion
  overall_edited_insertion.deletion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_insertion.variants.count.data, mut.type = "deletion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_insertion.deletion.table) <- c("Position", "overall_edited_insertion.deletion.Freq")
  # overall_edited_insertion/no.overlap.deletion
  overall_edited_insertion.no.overlap.deletion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_insertion.variants.count.data, mut.type = "no.overlap.deletion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_insertion.no.overlap.deletion.table) <- c("Position", "overall_edited_insertion.no.overlap.deletion.Freq")
  # overall_edited_insertion/insertion
  overall_edited_insertion.insertion.table <- data.frame(MakeMutPositionLabelTable(overall_edited_insertion.variants.count.data, mut.type = "insertion", range.vec = c(-35:-1,1:35)), stringsAsFactors=FALSE)
  colnames(overall_edited_insertion.insertion.table) <- c("Position", "overall_edited_insertion.insertion.Freq")


  # insertion (overall_edited_deletion + overall_edited_insertion + overall_edited_substitution)
  # merge
  mutation.insertion.merge.table <- merge(overall_edited_deletion.insertion.table, overall_edited_insertion.insertion.table, by = "Position", all=TRUE)
  mutation.insertion.merge.table <- merge(mutation.insertion.merge.table, overall_edited_substitution.insertion.table, by = "Position", all=TRUE)
  mutation.insertion.merge.table <- SortDataFrame(mutation.insertion.merge.table, mutation.insertion.merge.table$Position)
  # sum merge table
  insertion.mutation.position.table <- MakeMargeSumPositionTable(mutation.insertion.merge.table)
  if(!is.null(insertion.mutation.position.table)){
    mutation.insertion.p <- ggplot(data=insertion.mutation.position.table
      , aes(x=Position, y=Frequency, fill = c("insertion"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#D93240") +
      coord_cartesian(ylim = c(0, max(insertion.mutation.position.table$Frequency) + 100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Relative Position[bp]") +
    ylab("Number of insertion label")
    # plot(mutation.insertion.p)
    saveRDS(insertion.mutation.position.table, file = file.path(mutation.analysis.dir, "[ⅰa]edited_mutation.insertion.position.table.rds"))
    ggsave(file = file.path(mutation.analysis.dir, "[ⅰa]edited_mutation.insertion.position.png"), plot = mutation.insertion.p, dpi = 300, width = 12, height = 6)
  }


  # deletion (overall_edited_deletion + overall_edited_insertion + overall_edited_substitution)
  # merge
  mutation.deletion.merge.table <- merge(overall_edited_deletion.deletion.table, overall_edited_insertion.deletion.table, by = "Position", all=TRUE)
  mutation.deletion.merge.table <- merge(mutation.deletion.merge.table, overall_edited_substitution.deletion.table, by = "Position", all=TRUE)
  mutation.deletion.merge.table <- SortDataFrame(mutation.deletion.merge.table, mutation.deletion.merge.table$Position)
  # sum merge table
  deletion.mutation.position.table <- MakeMargeSumPositionTable(mutation.deletion.merge.table)
  if(!is.null(deletion.mutation.position.table)){
    mutation.deletion.deletion.p <- ggplot(data=deletion.mutation.position.table
      , aes(x=Position, y=Frequency, fill = c("deletion"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#2E7324") +
      coord_cartesian(ylim = c(0, max(deletion.mutation.position.table$Frequency) + 100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Relative Position[bp]") +
    ylab("Number of deletion label")
    # plot(mutation.deletion.deletion.p)
    saveRDS(deletion.mutation.position.table, file = file.path(mutation.analysis.dir, "[ⅰb]edited_mutation.deletion.position.table.rds"))
    ggsave(file = file.path(mutation.analysis.dir, "[ⅰb]edited_mutation.deletion.position.png"), plot = mutation.deletion.deletion.p, dpi = 300, width = 12, height = 6)
  }


  # no.overlap.deletion (overall_edited_deletion + overall_edited_insertion + overall_edited_substitution)
  # merge
  mutation.deletion.nooverlap.merge.table <- merge(overall_edited_deletion.no.overlap.deletion.table, overall_edited_insertion.no.overlap.deletion.table, by = "Position", all=TRUE)
  mutation.deletion.nooverlap.merge.table <- merge(mutation.deletion.nooverlap.merge.table, overall_edited_substitution.no.overlap.deletion.table, by = "Position", all=TRUE)
  mutation.deletion.nooverlap.merge.table <- SortDataFrame(mutation.deletion.nooverlap.merge.table, mutation.deletion.nooverlap.merge.table$Position)
  # sum merge table
  mutation.deletion.nooverlap.position.table <- MakeMargeSumPositionTable(mutation.deletion.nooverlap.merge.table)
  if(!is.null(mutation.deletion.nooverlap.position.table)){
    mutation.deletion.no.overlap.deletion.p <- ggplot(data=mutation.deletion.nooverlap.position.table
      , aes(x=Position, y=Frequency, fill = c("no.overlap.deletion"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#1C4012") +
      coord_cartesian(ylim = c(0, max(mutation.deletion.nooverlap.position.table$Frequency) + 100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Relative Position[bp]") +
    ylab("Number of no overlap deletion label")
    # plot(mutation.deletion.no.overlap.deletion.p)
    saveRDS(mutation.deletion.nooverlap.position.table, file = file.path(mutation.analysis.dir, "[ⅰc]edited_mutation.deletion.no.overlap.position.table.rds"))
    ggsave(file = file.path(mutation.analysis.dir, "[ⅰc]edited_mutation.deletion.no.overlap.position.png"), plot = mutation.deletion.no.overlap.deletion.p, dpi = 300, width = 12, height = 6)
  }

  # untreated.substitution barplot
  untreated.substitution.position.table <- data.frame(Position=untreated.substitution.table$Position, Frequency=untreated.substitution.table$untreated.substitution.Freq, stringsAsFactors=FALSE)
  if(!is.null(untreated.substitution.position.table)){
    untreated.substitution.p <- ggplot(data=untreated.substitution.position.table
      , aes(x=Position, y=Frequency, fill = c("substitution"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#D95F18") +
      coord_cartesian(ylim = c(0, sum(untreated.substitution.position.table$Frequency))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Relative Position[bp]") +
    ylab("Number of substitution label")
    # plot(untreated.substitution.p)
    saveRDS(untreated.substitution.position.table, file = file.path(mutation.analysis.dir, "[ⅰd]untreated.substitution.position.table.rds"))
    ggsave(file = file.path(mutation.analysis.dir, "[ⅰd]untreated.substitution.position.png"), plot = untreated.substitution.p, dpi = 300, width = 12, height = 6)
  }

  # mutation.substitution barplot
  mutation.substitution.position.table <- data.frame(Position=overall_edited_substitution.substitution.table$Position, Frequency=overall_edited_substitution.substitution.table$overall_edited_substitution.substitution.Freq, stringsAsFactors=FALSE)
  if(!is.null(mutation.substitution.position.table)){
    mutation.substitution.p <- ggplot(data=mutation.substitution.position.table
      , aes(x=Position, y=Frequency, fill = c("substitution"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#D95F18") +
      coord_cartesian(ylim = c(0, max(mutation.substitution.position.table$Frequency) + 100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Relative Position[bp]") +
    ylab("Number of substitution label")
    # plot(mutation.substitution.p)
    saveRDS(mutation.substitution.position.table, file = file.path(mutation.analysis.dir, "[ⅰe]mutation.substitution.position.table.rds"))
    ggsave(file = file.path(mutation.analysis.dir, "[ⅰe]mutation.substitution.position.png"), plot = mutation.substitution.p, dpi = 300, width = 12, height = 6)
  }

  # substitution barplot (untreated + mutation.)
  # merge
  both.substitution.merge.table <- merge(untreated.substitution.table, overall_edited_substitution.substitution.table, by = "Position", all=TRUE)
  both.substitution.merge.table <- SortDataFrame(both.substitution.merge.table, both.substitution.merge.table$Position)
  # sum merge table
  both.substitution.position.table <- MakeMargeSumPositionTable(both.substitution.merge.table)
  if(!is.null(both.substitution.position.table)){
    both.substitution.p <- ggplot(data=both.substitution.position.table
      , aes(x=Position, y=Frequency, fill = c("substitution"))) +
      geom_bar(stat="identity", width=0.5, alpha = 0.8, fill = "#591414") +
      coord_cartesian(ylim = c(0, max(both.substitution.position.table$Frequency) + 100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Relative Position[bp]") +
    ylab("Number of substitution label")
    # plot(both.substitution.p)
    saveRDS(both.substitution.position.table, file = file.path(mutation.analysis.dir, "[ⅰf]both.substitution.position.table.rds"))
    ggsave(file = file.path(mutation.analysis.dir, "[ⅰf]both.substitution.position.png"), plot = both.substitution.p, dpi = 300, width = 12, height = 6)
  }

  ##############################################################################################################################################
  message("Run aggregated mutation indels size analysis...")
  
  ### mutation position/size analysis

  # size batplot

  # indels (overall_edited_deletion + overall_edited_insertion + overall_edited_substitution)
  indels.variants.count.data <- c(
    overall_edited_substitution.variants.count.data
    , overall_edited_deletion.variants.count.data
    , overall_edited_insertion.variants.count.data
  )

  # focus.range = -35:35
  mutation.is.sn.list <- MakeIsSnReadTable(
    indels.variants.count.data
    , focus.range = -35:35
    , substitution.range.vec = 0:100
    , indels.range.vec = -100:100
    , special.label = "no variant"
    , is.special.included = TRUE
  )
  # mutation.is.sn.list$total.is.sn.list$indels.size.table
  # mutation.is.sn.list$total.is.sn.list$substitution.number.table

  # make mutation size plot
  if(!is.null(mutation.is.sn.list$total.is.sn.list$indels.size.table)){
    saveRDS(mutation.is.sn.list$total.is.sn.list$indels.size.table , file = file.path(mutation.analysis.dir, "[ⅱa]Distribution_of_indel_size_on_target_site.table.rds"))
    SaveIndelSizeReadBarPlot(mutation.is.sn.list$total.is.sn.list$indels.size.table, file.path(mutation.analysis.dir, "[ⅱa]Distribution_of_indel_size_on_target_site.png"))
  }

  # make substitution number plot
  if(!is.null(mutation.is.sn.list$total.is.sn.list$substitution.number.table)){
    saveRDS(mutation.is.sn.list$total.is.sn.list$substitution.number.table , file = file.path(mutation.analysis.dir, "[ⅱb]Distribution_of_substitution_number_on_target_site.table.rds"))
    SaveSubstitutionNumberReadBarPlot(mutation.is.sn.list$total.is.sn.list$substitution.number.table, file.path(mutation.analysis.dir, "[ⅱb]Distribution_of_substitution_number_on_target_site.png"))
  }

  # Make indel size and substitution pie chart
  if(!is.null(mutation.is.sn.list)){
    SaveIndelPlot(mutation.is.sn.list$total.is.sn.list, outdir = mutation.analysis.dir, ki.type = "mutation")
  }

  # Make ej.type rate pie chart
  SaveMutEjTypePlot(mutation.is.sn.list, mutation.analysis.dir)

  # Make deletion length-count plot
  if(!is.null(mutation.is.sn.list$mmej.is.sn.list$indels.size.table)){
    saveRDS(mutation.is.sn.list$mmej.is.sn.list$indels.size.table , file = file.path(mutation.analysis.dir, "[ⅵa]Distribution_of_MMEJ_indel_size_on_target_site.table.rds"))
    SaveIndelSizeReadBarPlot(mutation.is.sn.list$mmej.is.sn.list$indels.size.table, file.path(mutation.analysis.dir, "[ⅵa]Distribution_of_MMEJ_indel_size_on_target_site.png"))
  }
  if(!is.null(mutation.is.sn.list$nhej.is.sn.list$indels.size.table)){
    saveRDS(mutation.is.sn.list$nhej.is.sn.list$indels.size.table , file = file.path(mutation.analysis.dir, "[ⅵb]Distribution_of_NHEJ_indel_size_on_target_site.table.rds"))
    SaveIndelSizeReadBarPlot(mutation.is.sn.list$nhej.is.sn.list$indels.size.table, file.path(mutation.analysis.dir, "[ⅵb]Distribution_of_NHEJ_indel_size_on_target_site.png"))
  }

  #prepare microhomology/trimmed_seq length table
  mutation.mlc.tlc.list <- MakeMlcTlcReadTable(indels.variants.count.data, microhomology.range.vec = 0:100, trimmedseq.range.vec = 0:100)
  mutation.mmej.mlc.tlc.list <- mutation.mlc.tlc.list$mmej.mlc.tlc.list
  mutation.nhej.mlc.tlc.list <- mutation.mlc.tlc.list$nhej.mlc.tlc.list

  # Make microhomology length-count plot
  if(!is.null(mutation.mmej.mlc.tlc.list$mh.size.table)){
    saveRDS(mutation.mmej.mlc.tlc.list$mh.size.table , file = file.path(mutation.analysis.dir, "[ⅶa]Distribution_of_MMEJ_microhomology_length_on_target_site.table.rds"))
    SaveSeqLengthReadBarPlot(
      mutation.mmej.mlc.tlc.list$mh.size.table
      , file.path(mutation.analysis.dir, "[ⅶa]Distribution_of_MMEJ_microhomology_length_on_target_site.png")
      , seq.type = "Microhomology"
    )
  }
  if(!is.null(mutation.nhej.mlc.tlc.list$mh.size.table)){
    saveRDS(mutation.nhej.mlc.tlc.list$mh.size.table , file = file.path(mutation.analysis.dir, "[ⅶb]Distribution_of_NHEJ_microhomology_length_on_target_site.table.rds"))
    SaveSeqLengthReadBarPlot(
      mutation.nhej.mlc.tlc.list$mh.size.table
      , file.path(mutation.analysis.dir, "[ⅶb]Distribution_of_NHEJ_microhomology_length_on_target_site.png")
      , seq.type = "Microhomology"
    )
  }

  # Make trimmed length-count plot
  if(!is.null(mutation.mmej.mlc.tlc.list$trim.size.table)){
    saveRDS(mutation.mmej.mlc.tlc.list$trim.size.table , file = file.path(mutation.analysis.dir, "[ⅷa]Distribution_of_MMEJ_Trimmed_Seq_length_on_target_site.table.rds"))
    SaveSeqLengthReadBarPlot(
      mutation.mmej.mlc.tlc.list$trim.size.table
      , file.path(mutation.analysis.dir, "[ⅷa]Distribution_of_MMEJ_Trimmed_Seq_length_on_target_site.png")
      , seq.type = "Trimmed_seq"
    )
  }
  if(!is.null(mutation.nhej.mlc.tlc.list$trim.size.table)){
    saveRDS(mutation.nhej.mlc.tlc.list$trim.size.table , file = file.path(mutation.analysis.dir, "[ⅷb]Distribution_of_NHEJ_Trimmed_Seq_length_on_target_site.table.rds"))
    SaveSeqLengthReadBarPlot(
      mutation.nhej.mlc.tlc.list$trim.size.table
      , file.path(mutation.analysis.dir, "[ⅷb]Distribution_of_NHEJ_Trimmed_Seq_length_on_target_site.png")
      , seq.type = "Trimmed_seq"
    )
  }

  # Make trimmed length - microhomology length plot
  SaveMlTlScatterPlot(indels.variants.count.data, mutation.analysis.dir, "mut")

  # Make microhomology sequence-frequency-count plot
  SaveMsFreqPlot(indels.variants.count.data, mutation.analysis.dir, "mut")

  # Make intervening sequence-frequency-count plot
  SaveTsFreqPlot(indels.variants.count.data, mutation.analysis.dir, "mut")

  return(TRUE)
}