
RunLocalKIAnalysis <- function(bam.datatype.dir, result.datatype.dir, index.dir, seq.dir){

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

  # junction length
  left.precise.junction.range <- -nchar(left.homology.seq):-1
  left.precise.leftpart.junction.range <- left.precise.junction.range[1:(length(left.precise.junction.range)/2)]
  left.precise.rightpart.junction.range <- left.precise.junction.range[(length(left.precise.junction.range)/2+1):length(left.precise.junction.range)]
  left.endjoining.junction.range <- -(nchar(left.homology.seq) + nchar(left.homology.seq)):-1
  rc.left.endjoining.junction.range <- -(nchar(left.homology.seq) + nchar(right.homology.seq)):-1
  right.precise.junction.range <- 1:nchar(right.homology.seq)
  right.leftpart.precise.junction.range <- right.precise.junction.range[1:(length(right.precise.junction.range)/2)]
  right.rightpart.precise.junction.range <- right.precise.junction.range[(length(right.precise.junction.range)/2+1):length(right.precise.junction.range)]
  right.endjoining.junction.range <- 1:(nchar(right.homology.seq) + nchar(right.homology.seq))
  rc.right.endjoining.junction.range <- 1:(nchar(right.homology.seq) + nchar(left.homology.seq))

  # bam dir
  bam.ki_type.dir <- file.path(bam.datatype.dir, "ki_type")
  bam.rc_ki_type.dir <- file.path(bam.datatype.dir, "rc_ki_type")
  # result dir
  alignment.analysis.dir <- file.path(result.datatype.dir, "[Ⅰ]alignment_analysis")
  left_precise_ki.dir <- file.path(alignment.analysis.dir, "left_precise_ki")
  left_imprecise_ki.dir <- file.path(alignment.analysis.dir, "left_imprecise_ki")
  left_imprecise_ki_complex.dir <- file.path(alignment.analysis.dir, "left_imprecise_ki_complex")
  left_other_mutation.dir <- file.path(alignment.analysis.dir, "left_other_mutation")
  rc_left_imprecise_ki.dir <- file.path(alignment.analysis.dir, "rc_left_imprecise_ki")
  rc_left_imprecise_ki_complex.dir <- file.path(alignment.analysis.dir, "rc_left_imprecise_ki_complex")
  rc_left_other_mutation.dir <- file.path(alignment.analysis.dir, "rc_left_other_mutation")
  right_precise_ki.dir <- file.path(alignment.analysis.dir, "right_precise_ki")
  right_imprecise_ki.dir <- file.path(alignment.analysis.dir, "right_imprecise_ki")
  right_imprecise_ki_complex.dir <- file.path(alignment.analysis.dir, "right_imprecise_ki_complex")
  right_other_mutation.dir <- file.path(alignment.analysis.dir, "right_other_mutation")
  rc_right_imprecise_ki.dir <- file.path(alignment.analysis.dir, "rc_right_imprecise_ki")
  rc_right_imprecise_ki_complex.dir <- file.path(alignment.analysis.dir, "rc_right_imprecise_ki_complex")
  rc_right_other_mutation.dir <- file.path(alignment.analysis.dir, "rc_right_other_mutation")
  knock.in.analysis.dir <- file.path(result.datatype.dir, "[Ⅱb]analysis_of_knock-in_junction")

  message("Make directory for saving the analysis data")
  for(dirname in c(alignment.analysis.dir, bam.rc_ki_type.dir, left_precise_ki.dir, left_imprecise_ki.dir
    , left_imprecise_ki_complex.dir, bam.ki_type.dir, left_other_mutation.dir, rc_left_imprecise_ki.dir, rc_left_imprecise_ki_complex.dir
    , rc_left_other_mutation.dir, right_precise_ki.dir, right_imprecise_ki.dir, right_imprecise_ki_complex.dir
    , right_other_mutation.dir, rc_right_imprecise_ki.dir, rc_right_imprecise_ki_complex.dir, rc_right_other_mutation.dir
    , knock.in.analysis.dir)){
    dir.create(file.path(dirname), showWarnings=FALSE)
  }
  
  ###################################################################################################################################################

  message("Start local mutation alignment...")

  ### left precise knock-in
  if(!file.exists(file.path(left_precise_ki.dir, "crispr.set.rds"))){
    left_precise_ki.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "mmej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "LEFT_PRECISE_KNOCK_IN.bam")
      , "left_precise_ki"
      , "left.precise")
    if(!is.null(left_precise_ki.crispr.set))SaveVariantsData(left_precise_ki.crispr.set, left_precise_ki.dir, range.vec = left.precise.junction.range, mut.type = "left.ki")
  }else{
    left_precise_ki.crispr.set <- readRDS(file.path(left_precise_ki.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(left_imprecise_ki.dir, "crispr.set.rds"))){
    left_imprecise_ki.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "LEFT_IMPRECISE_KNOCK_IN.bam")
      , "left_imprecise_ki"
      , "left.imprecise")
    if(!is.null(left_imprecise_ki.crispr.set))SaveVariantsData(left_imprecise_ki.crispr.set, left_imprecise_ki.dir, range.vec = left.endjoining.junction.range, mut.type = "left.ki")
  }else{
    left_imprecise_ki.crispr.set <- readRDS(file.path(left_imprecise_ki.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(left_imprecise_ki_complex.dir, "crispr.set.rds"))){
    left_imprecise_ki_complex.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "LEFT_IMPRECISE_KNOCK_IN_COMPLEX.bam")
      , "left_imprecise_ki_complex"
      , "left.imprecise")
    if(!is.null(left_imprecise_ki_complex.crispr.set))SaveVariantsData(left_imprecise_ki_complex.crispr.set, left_imprecise_ki_complex.dir, range.vec = left.endjoining.junction.range, mut.type = "left.ki")
  }else{
    left_imprecise_ki_complex.crispr.set <- readRDS(file.path(left_imprecise_ki_complex.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(left_other_mutation.dir, "crispr.set.rds"))){
    left_other_mutation.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "LEFT_OTHER_MUATIONS.bam")
      , "left_other_mutation"
      , "left.imprecise")
    if(!is.null(left_other_mutation.crispr.set))SaveVariantsData(left_other_mutation.crispr.set, left_other_mutation.dir, range.vec = left.endjoining.junction.range, mut.type = "left.ki")
  }else{
    left_other_mutation.crispr.set <- readRDS(file.path(left_other_mutation.dir, "crispr.set.rds"))
  }
  
  if(!file.exists(file.path(rc_left_imprecise_ki.dir, "crispr.set.rds"))){
    rc_left_imprecise_ki.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_rc_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.rc_ki_type.dir, "rc_LEFT_IMPRECISE_KNOCK_IN.bam")
      , "rc_left_imprecise_ki"
      , "rc.left.imprecise")
    if(!is.null(rc_left_imprecise_ki.crispr.set))SaveVariantsData(rc_left_imprecise_ki.crispr.set, rc_left_imprecise_ki.dir, range.vec = rc.left.endjoining.junction.range, mut.type = "left.ki")
  }else{
    rc_left_imprecise_ki.crispr.set <- readRDS(file.path(rc_left_imprecise_ki.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(rc_left_imprecise_ki_complex.dir, "crispr.set.rds"))){
    rc_left_imprecise_ki_complex.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_rc_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.rc_ki_type.dir, "rc_LEFT_IMPRECISE_KNOCK_IN_COMPLEX.bam")
      , "rc_left_imprecise_ki_complex"
      , "rc.left.imprecise")
    if(!is.null(rc_left_imprecise_ki_complex.crispr.set))SaveVariantsData(rc_left_imprecise_ki_complex.crispr.set, rc_left_imprecise_ki_complex.dir, range.vec = rc.left.endjoining.junction.range, mut.type = "left.ki")
  }else{
    rc_left_imprecise_ki_complex.crispr.set <- readRDS(file.path(rc_left_imprecise_ki_complex.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(rc_left_other_mutation.dir, "crispr.set.rds"))){
    rc_left_other_mutation.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_rc_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.rc_ki_type.dir, "rc_LEFT_OTHER_MUATIONS.bam")
      , "rc_left_other_mutation"
      , "rc.left.imprecise")
    if(!is.null(rc_left_other_mutation.crispr.set))SaveVariantsData(rc_left_other_mutation.crispr.set, rc_left_other_mutation.dir, range.vec = rc.left.endjoining.junction.range, mut.type = "left.ki")
  }else{
    rc_left_other_mutation.crispr.set <- readRDS(file.path(rc_left_other_mutation.dir, "crispr.set.rds"))
  }


  ### right precise knock-in
  if(!file.exists(file.path(right_precise_ki.dir, "crispr.set.rds"))){
    right_precise_ki.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "mmej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "RIGHT_PRECISE_KNOCK_IN.bam")
      , "right_precise_ki"
      , "right.precise")
    if(!is.null(right_precise_ki.crispr.set))SaveVariantsData(right_precise_ki.crispr.set, right_precise_ki.dir, range.vec = right.precise.junction.range, mut.type = "right.ki")
  }else{
    right_precise_ki.crispr.set <- readRDS(file.path(right_precise_ki.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(right_imprecise_ki.dir, "crispr.set.rds"))){
    right_imprecise_ki.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "RIGHT_IMPRECISE_KNOCK_IN.bam")
      , "right_imprecise_ki"
      , "right.imprecise")
    if(!is.null(right_imprecise_ki.crispr.set))SaveVariantsData(right_imprecise_ki.crispr.set, right_imprecise_ki.dir, range.vec = right.endjoining.junction.range, mut.type = "right.ki")
  }else{
    right_imprecise_ki.crispr.set <- readRDS(file.path(right_imprecise_ki.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(right_imprecise_ki_complex.dir, "crispr.set.rds"))){
    right_imprecise_ki_complex.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX.bam")
      , "right_imprecise_ki_complex"
      , "right.imprecise")
    if(!is.null(right_imprecise_ki_complex.crispr.set))SaveVariantsData(right_imprecise_ki_complex.crispr.set, right_imprecise_ki_complex.dir, range.vec = right.endjoining.junction.range, mut.type = "right.ki")
  }else{
    right_imprecise_ki_complex.crispr.set <- readRDS(file.path(right_imprecise_ki_complex.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(right_other_mutation.dir, "crispr.set.rds"))){
    right_other_mutation.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.ki_type.dir, "RIGHT_OTHER_MUATIONS.bam")
      , "right_other_mutation"
      , "right.imprecise")
    if(!is.null(right_other_mutation.crispr.set))SaveVariantsData(right_other_mutation.crispr.set, right_other_mutation.dir, range.vec = right.endjoining.junction.range, mut.type = "right.ki")
  }else{
    right_other_mutation.crispr.set <- readRDS(file.path(right_other_mutation.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(rc_right_imprecise_ki.dir, "crispr.set.rds"))){
    rc_right_imprecise_ki.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_rc_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.rc_ki_type.dir, "rc_RIGHT_IMPRECISE_KNOCK_IN.bam")
      , "rc_right_imprecise_ki"
      , "rc.right.imprecise")
    if(!is.null(rc_right_imprecise_ki.crispr.set))SaveVariantsData(rc_right_imprecise_ki.crispr.set, rc_right_imprecise_ki.dir, range.vec = rc.right.endjoining.junction.range, mut.type = "right.ki")
  }else{
    rc_right_imprecise_ki.crispr.set <- readRDS(file.path(rc_right_imprecise_ki.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(rc_right_imprecise_ki_complex.dir, "crispr.set.rds"))){
    rc_right_imprecise_ki_complex.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_rc_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.rc_ki_type.dir, "rc_RIGHT_IMPRECISE_KNOCK_IN_COMPLEX.bam")
      , "rc_right_imprecise_ki_complex"
      , "rc.right.imprecise")
    if(!is.null(rc_right_imprecise_ki_complex.crispr.set))SaveVariantsData(rc_right_imprecise_ki_complex.crispr.set, rc_right_imprecise_ki_complex.dir, range.vec = rc.right.endjoining.junction.range, mut.type = "right.ki")
  }else{
    rc_right_imprecise_ki_complex.crispr.set <- readRDS(file.path(rc_right_imprecise_ki_complex.dir, "crispr.set.rds"))
  }

  if(!file.exists(file.path(rc_right_other_mutation.dir, "crispr.set.rds"))){
    rc_right_other_mutation.crispr.set <- MakeKiCrisprSet(
      outleft.indicator.seq
      , outright.indicator.seq
      , inleft.indicator.seq
      , inright.indicator.seq
      , left.homology.seq
      , right.homology.seq
      , "nhej_rc_ki"
      , mmej.ki.seq
      , donor.seq
      , file.path(bam.rc_ki_type.dir, "rc_RIGHT_OTHER_MUATIONS.bam")
      , "rc_right_other_mutation"
      , "rc.right.imprecise")
    if(!is.null(rc_right_other_mutation.crispr.set))SaveVariantsData(rc_right_other_mutation.crispr.set, rc_right_other_mutation.dir, range.vec = rc.right.endjoining.junction.range, mut.type = "right.ki")
  }else{
    rc_right_other_mutation.crispr.set <- readRDS(file.path(rc_right_other_mutation.dir, "crispr.set.rds"))
  }

  message("Local alignment is Completed.")

  ##########################################################################################################################################################################
  ### left junction position/size analysis
  message("Run knock-in indels size analysis...")

  ##########################################################################---------------------------------------------########

  ### left knock-in
  left_precise_ki.variants.count.data <- MakeTotalMutList(left_precise_ki.crispr.set, "left.ki")
  # focus.range = left.precise.junction.range
  left_precise_ki.is.sn.list <- MakeIsSnReadTable(
    left_precise_ki.variants.count.data
    , focus.range = left.precise.junction.range
    , substitution.range.vec = 0:100
    , indels.range.vec = -100:100 + nchar(left.homology.seq) + nchar(right.homology.seq)
    , special.label = "perfect_left_junction"
    , is.special.included = TRUE
  )

  left_imprecise_ki.variants.count.data <- MakeTotalMutList(left_imprecise_ki.crispr.set, "left.ki")
  left_imprecise_ki_complex.variants.count.data <- MakeTotalMutList(left_imprecise_ki_complex.crispr.set, "left.ki")
  left_other_mutation_ki.variants.count.data <- MakeTotalMutList(left_other_mutation.crispr.set, "left.ki")
  left_imprecise_integration.variants.count.data <- c(
    left_imprecise_ki.variants.count.data
    , left_imprecise_ki_complex.variants.count.data
    , left_other_mutation_ki.variants.count.data
  )
  # focus.range = left.endjoining.junction.range
  left_imprecise_integration.is.sn.list <- MakeIsSnReadTable(
    left_imprecise_integration.variants.count.data
    , focus.range = left.endjoining.junction.range
    , substitution.range.vec = 0:100
    , indels.range.vec = -100:100
    , special.label = "duplicated_homology_left_junction"
    , is.special.included = TRUE
  )

  rc_left_imprecise_ki.variants.count.data <- MakeTotalMutList(rc_left_imprecise_ki.crispr.set, "left.ki")
  rc_left_imprecise_ki_complex.variants.count.data <- MakeTotalMutList(rc_left_imprecise_ki_complex.crispr.set, "left.ki")
  rc_left_other_mutation_ki.variants.count.data <- MakeTotalMutList(rc_left_other_mutation.crispr.set, "left.ki")
  rc_left_imprecise_integration.variants.count.data <- c(
    rc_left_imprecise_ki.variants.count.data
    , rc_left_imprecise_ki_complex.variants.count.data
    , rc_left_other_mutation_ki.variants.count.data
  )
  # focus.range = rc.left.endjoining.junction.range
  rc_left_imprecise_integration.is.sn.list <- MakeIsSnReadTable(
    rc_left_imprecise_integration.variants.count.data
    , focus.range = rc.left.endjoining.junction.range
    , substitution.range.vec = 0:100
    , indels.range.vec = -100:100
    , special.label = "reverse_knockin_left_junction"
    , is.special.included = TRUE
  )

  # The pattrn : donor insert & LEFT_NO_DETECTED is not anallized because MaChIAto considers no indicator sequence as unclassified samples.
  # And donor insert & LEFT_NO_DETECTED has right in-out indicator according to MaChIAto algorithm. This region is analized.

  ### right knock-in
  right_precise_ki.variants.count.data <- MakeTotalMutList(right_precise_ki.crispr.set, "right.ki")
  # focus.range = right.precise.junction.range
  right_precise_ki.is.sn.list <- MakeIsSnReadTable(
    right_precise_ki.variants.count.data
    , focus.range = right.precise.junction.range
    , substitution.range.vec = 0:100
    #, indels.range.vec = -100:100 + nchar(right.homology.seq) + nchar(right.homology.seq)
    , indels.range.vec = -100:100 + nchar(right.homology.seq) + nchar(left.homology.seq)
    , special.label = "perfect_right_junction"
    , is.special.included = TRUE
  )

  right_imprecise_ki.variants.count.data <- MakeTotalMutList(right_imprecise_ki.crispr.set, "right.ki")
  right_imprecise_ki_complex.variants.count.data <- MakeTotalMutList(right_imprecise_ki_complex.crispr.set, "right.ki")
  right_other_mutation_ki.variants.count.data <- MakeTotalMutList(right_other_mutation.crispr.set, "right.ki")
  right_imprecise_integration.variants.count.data <- c(
    right_imprecise_ki.variants.count.data
    , right_imprecise_ki_complex.variants.count.data
    , right_other_mutation_ki.variants.count.data
  )
  # focus.range = right.endjoining.junction.range
  right_imprecise_integration.is.sn.list <- MakeIsSnReadTable(
    right_imprecise_integration.variants.count.data
    , focus.range = right.endjoining.junction.range
    , substitution.range.vec = 0:100
    , indels.range.vec = -100:100
    , special.label = "duplicated_homology_right_junction"
    , is.special.included = TRUE
  )

  rc_right_imprecise_ki.variants.count.data <- MakeTotalMutList(rc_right_imprecise_ki.crispr.set, "right.ki")
  rc_right_imprecise_ki_complex.variants.count.data <- MakeTotalMutList(rc_right_imprecise_ki_complex.crispr.set, "right.ki")
  rc_right_other_mutation_ki.variants.count.data <- MakeTotalMutList(rc_right_other_mutation.crispr.set, "right.ki")
  rc_right_imprecise_integration.variants.count.data <- c(
    rc_right_imprecise_ki.variants.count.data
    , rc_right_imprecise_ki_complex.variants.count.data
    , rc_right_other_mutation_ki.variants.count.data
  )
  # focus.range = rc.right.endjoining.junction.range
  rc_right_imprecise_integration.is.sn.list <- MakeIsSnReadTable(
    rc_right_imprecise_integration.variants.count.data
    , focus.range = rc.right.endjoining.junction.range
    , substitution.range.vec = 0:100
    , indels.range.vec = -100:100
    , special.label = "reverse_knockin_right_junction"
    , is.special.included = TRUE
  )

  # The pattrn : donor insert & RIGHT_NO_DETECTED is not anallized because MaChIAto considers no indicator sequence as unclassified samples.
  # And donor insert & RIGHT_NO_DETECTED has left in-out indicator according to MaChIAto algorithm. This region is analized.

  ##########################################################################---------------------------------------------########

  ##########################################################################---------------------------------------------########

  ### left knock-in
  RunKiTotalIndelAnalysis(
    precise.is.sn.list = left_precise_ki.is.sn.list$total.is.sn.list
    , imprecise.is.sn.list = left_imprecise_integration.is.sn.list$total.is.sn.list
    , rev.imprecise.is.sn.list = rc_left_imprecise_integration.is.sn.list$total.is.sn.list
    # , nhej.insert.length = nchar(left.homology.seq) + nchar(left.homology.seq)
    # , rc.nhej.insert.length = nchar(left.homology.seq) + nchar(right.homology.seq)
    , nhej.insert.length = nchar(left.homology.seq) + nchar(left.extra.seq)
    , rc.nhej.insert.length = nchar(right.homology.seq) + nchar(right.extra.seq)
    , save.dir = knock.in.analysis.dir
    , ki.type = "left_ki"
  )

  ### right knock-in
  RunKiTotalIndelAnalysis(
    precise.is.sn.list = right_precise_ki.is.sn.list$total.is.sn.list
    , imprecise.is.sn.list = right_imprecise_integration.is.sn.list$total.is.sn.list
    , rev.imprecise.is.sn.list = rc_right_imprecise_integration.is.sn.list$total.is.sn.list
    # , nhej.insert.length = nchar(right.homology.seq) + nchar(right.homology.seq)
    # , rc.nhej.insert.length = nchar(left.homology.seq) + nchar(right.homology.seq)
    , nhej.insert.length = nchar(right.homology.seq) + nchar(right.extra.seq)
    , rc.nhej.insert.length = nchar(left.homology.seq) + nchar(left.extra.seq)
    , save.dir = knock.in.analysis.dir
    , ki.type = "right_ki"
  )

  ##########################################################################---------------------------------------------########
  # MMEJ analysis
  ##########################################################################---------------------------------------------########

  # Make ej.type rate pie chart
  # left junction
  SaveImKiEjTypePlot(
    null.precise.is.sn.list.set = MakeIsSnReadTable(
      list()
      , focus.range = left.precise.junction.range
      , substitution.range.vec = 0:100
      # , indels.range.vec = -100:100 + nchar(left.homology.seq) + nchar(right.homology.seq)
      , indels.range.vec = -100:100 + nchar(left.homology.seq) + nchar(right.homology.seq) + nchar(left.extra.seq) + nchar(right.extra.seq)
      , special.label = "perfect_left_junction"
      , is.special.included = TRUE
    )
    , imprecise.is.sn.list.set = left_imprecise_integration.is.sn.list
    , rev.imprecise.is.sn.list.set = rc_left_imprecise_integration.is.sn.list
    #, nhej.insert.length = nchar(left.homology.seq) + nchar(left.homology.seq)
    #, rc.nhej.insert.length = nchar(left.homology.seq) + nchar(right.homology.seq)
    , nhej.insert.length = nchar(left.homology.seq) + nchar(left.extra.seq)
    , rc.nhej.insert.length = nchar(right.homology.seq) + nchar(right.extra.seq)
    , outdir = knock.in.analysis.dir
    , ki.type = "left_ki")
  # right jucntion
  SaveImKiEjTypePlot(
    null.precise.is.sn.list.set = MakeIsSnReadTable(
      list()
      , focus.range = right.precise.junction.range
      , substitution.range.vec = 0:100
      # , indels.range.vec = -100:100 + nchar(right.homology.seq) + nchar(right.homology.seq)
      , indels.range.vec = -100:100 + nchar(left.homology.seq) + nchar(right.homology.seq) + nchar(left.extra.seq) + nchar(right.extra.seq)
      , special.label = "perfect_right_junction"
      , is.special.included = TRUE
    )
    , imprecise.is.sn.list.set = right_imprecise_integration.is.sn.list
    , rev.imprecise.is.sn.list.set = rc_right_imprecise_integration.is.sn.list
    # , nhej.insert.length = nchar(right.homology.seq) + nchar(right.homology.seq)
    # , rc.nhej.insert.length = nchar(right.homology.seq) + nchar(right.homology.seq)
    , nhej.insert.length = nchar(right.homology.seq) + nchar(right.extra.seq)
    , rc.nhej.insert.length = nchar(left.homology.seq) + nchar(left.extra.seq)
    , outdir = knock.in.analysis.dir
    , ki.type = "right_ki")

  # End joining indel size plot
  # left junction
  SaveImKiIndelSizePlot(
    null.precise.is.sn.list.set = MakeIsSnReadTable(
      list()
      , focus.range = left.precise.junction.range
      , substitution.range.vec = 0:100
      , indels.range.vec = -100:100 + nchar(left.homology.seq) + nchar(right.homology.seq) + nchar(left.extra.seq) + nchar(right.extra.seq)
      , special.label = "perfect_left_junction"
      , is.special.included = TRUE
    )
    , imprecise.is.sn.list.set = left_imprecise_integration.is.sn.list
    , rev.imprecise.is.sn.list.set = rc_left_imprecise_integration.is.sn.list
    , nhej.insert.length = nchar(left.homology.seq) + nchar(left.extra.seq)
    , rc.nhej.insert.length = nchar(right.homology.seq) + nchar(right.extra.seq)
    , outdir = knock.in.analysis.dir
    , ki.type = "left_ki")
  # right junction
  SaveImKiIndelSizePlot(
    null.precise.is.sn.list.set = MakeIsSnReadTable(
      list()
      , focus.range = right.precise.junction.range
      , substitution.range.vec = 0:100
      , indels.range.vec = -100:100 + nchar(left.homology.seq) + nchar(right.homology.seq) + nchar(left.extra.seq) + nchar(right.extra.seq)
      , special.label = "perfect_right_junction"
      , is.special.included = TRUE
    )
    , imprecise.is.sn.list.set = right_imprecise_integration.is.sn.list
    , rev.imprecise.is.sn.list.set = rc_right_imprecise_integration.is.sn.list
    , nhej.insert.length = nchar(right.homology.seq) + nchar(right.extra.seq)
    , rc.nhej.insert.length = nchar(left.homology.seq) + nchar(left.extra.seq)
    , outdir = knock.in.analysis.dir
    , ki.type = "right_ki")

  # Prepare microhomology/trimmed_seq length table
  # left jucntion
  left_imprecise_integration.bothside.variants.count.data <- c(left_imprecise_integration.variants.count.data, rc_left_imprecise_integration.variants.count.data)
  left.imprecise.ki.mlc.tlc.list <- MakeMlcTlcReadTable(left_imprecise_integration.bothside.variants.count.data, microhomology.range.vec = 0:100, trimmedseq.range.vec = 0:100)
  left.imprecise.mmej.mlc.tlc.list <- left.imprecise.ki.mlc.tlc.list$mmej.mlc.tlc.list
  left.imprecise.nhej.mlc.tlc.list <- left.imprecise.ki.mlc.tlc.list$nhej.mlc.tlc.list
  # right jucntion
  right_imprecise_integration.bothside.variants.count.data <- c(right_imprecise_integration.variants.count.data, rc_right_imprecise_integration.variants.count.data)
  right.imprecise.ki.mlc.tlc.list <- MakeMlcTlcReadTable(right_imprecise_integration.bothside.variants.count.data, microhomology.range.vec = 0:100, trimmedseq.range.vec = 0:100)
  right.imprecise.mmej.mlc.tlc.list <- right.imprecise.ki.mlc.tlc.list$mmej.mlc.tlc.list
  right.imprecise.nhej.mlc.tlc.list <- right.imprecise.ki.mlc.tlc.list$nhej.mlc.tlc.list

  # Make trimmed length-count plot
  # left jucntion
  if(!is.null(left.imprecise.mmej.mlc.tlc.list$mh.size.table)){
    saveRDS(left.imprecise.mmej.mlc.tlc.list$mh.size.table , file = file.path(knock.in.analysis.dir, "[ⅴa]Distribution_of_MMEJ_microhomology_length_on_left_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      left.imprecise.mmej.mlc.tlc.list$mh.size.table
      , file.path(knock.in.analysis.dir, "[ⅴa]Distribution_of_MMEJ_microhomology_length_on_left_junction.png")
      , seq.type = "Microhomology"
    )
  }
  if(!is.null(left.imprecise.nhej.mlc.tlc.list$mh.size.table)){
    saveRDS(left.imprecise.nhej.mlc.tlc.list$mh.size.table , file = file.path(knock.in.analysis.dir, "[ⅴb]Distribution_of_NHEJ_microhomology_length_on_left_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      left.imprecise.nhej.mlc.tlc.list$mh.size.table
      , file.path(knock.in.analysis.dir, "[ⅴb]Distribution_of_NHEJ_microhomology_length_on_left_junction.png")
      , seq.type = "Microhomology"
    )
  }
  # right jucntion
  if(!is.null(right.imprecise.mmej.mlc.tlc.list$mh.size.table)){
    saveRDS(right.imprecise.mmej.mlc.tlc.list$mh.size.table , file = file.path(knock.in.analysis.dir, "[ⅴa]Distribution_of_MMEJ_microhomology_length_on_right_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      right.imprecise.mmej.mlc.tlc.list$mh.size.table
      , file.path(knock.in.analysis.dir, "[ⅴa]Distribution_of_MMEJ_microhomology_length_on_right_junction.png")
      , seq.type = "Microhomology"
    )
  }
  if(!is.null(right.imprecise.nhej.mlc.tlc.list$mh.size.table)){
    saveRDS(right.imprecise.nhej.mlc.tlc.list$mh.size.table , file = file.path(knock.in.analysis.dir, "[ⅴb]Distribution_of_NHEJ_microhomology_length_on_right_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      right.imprecise.nhej.mlc.tlc.list$mh.size.table
      , file.path(knock.in.analysis.dir, "[ⅴb]Distribution_of_NHEJ_microhomology_length_on_right_junction.png")
      , seq.type = "Microhomology"
    )
  }

  # Make trimmed length-count plot
  # left jucntion
  if(!is.null(left.imprecise.mmej.mlc.tlc.list$trim.size.table)){
    saveRDS(left.imprecise.mmej.mlc.tlc.list$trim.size.table , file = file.path(knock.in.analysis.dir, "[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_left_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      left.imprecise.mmej.mlc.tlc.list$trim.size.table
      , file.path(knock.in.analysis.dir, "[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_left_junction.png")
      , seq.type = "Trimmed_seq"
    )
  }
  if(!is.null(left.imprecise.nhej.mlc.tlc.list$trim.size.table)){
    saveRDS(left.imprecise.nhej.mlc.tlc.list$trim.size.table , file = file.path(knock.in.analysis.dir, "[ⅵb]Distribution_of_NHEJ_Trimmed_Seq_length_on_left_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      left.imprecise.nhej.mlc.tlc.list$trim.size.table
      , file.path(knock.in.analysis.dir, "[ⅵb]Distribution_of_NHEJ_Trimmed_Seq_length_on_left_junction.png")
      , seq.type = "Trimmed_seq"
    )
  }
  # right jucntion
  if(!is.null(right.imprecise.mmej.mlc.tlc.list$trim.size.table)){
    saveRDS(right.imprecise.mmej.mlc.tlc.list$trim.size.table , file = file.path(knock.in.analysis.dir, "[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_right_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      right.imprecise.mmej.mlc.tlc.list$trim.size.table
      , file.path(knock.in.analysis.dir, "[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_right_junction.png")
      , seq.type = "Trimmed_seq"
    )
  }
  if(!is.null(right.imprecise.nhej.mlc.tlc.list$trim.size.table)){
    saveRDS(right.imprecise.nhej.mlc.tlc.list$trim.size.table , file = file.path(knock.in.analysis.dir, "[ⅵb]Distribution_of_NHEJ_Trimmed_Seq_length_on_right_junction.table.rds"))
    SaveSeqLengthReadBarPlot(
      right.imprecise.nhej.mlc.tlc.list$trim.size.table
      , file.path(knock.in.analysis.dir, "[ⅵb]Distribution_of_NHEJ_Trimmed_Seq_length_on_right_junction.png")
      , seq.type = "Trimmed_seq"
    )
  }

  # Make trimmed length - microhomology length plot
  # left jucntion
  SaveMlTlScatterPlot(left_imprecise_integration.bothside.variants.count.data, knock.in.analysis.dir, "left_ki")
  # right jucntion
  SaveMlTlScatterPlot(right_imprecise_integration.bothside.variants.count.data, knock.in.analysis.dir, "right_ki")

  # Make microhomology sequence-frequency-count plot
  # left jucntion
  SaveMsFreqPlot(left_imprecise_integration.bothside.variants.count.data, knock.in.analysis.dir, "left_ki")
  # right jucntion
  SaveMsFreqPlot(right_imprecise_integration.bothside.variants.count.data, knock.in.analysis.dir, "right_ki")

  # Make intervening sequence-frequency-count plot
  # left jucntion
  SaveTsFreqPlot(left_imprecise_integration.bothside.variants.count.data, knock.in.analysis.dir, "left_ki")
  # right jucntion
  SaveTsFreqPlot(right_imprecise_integration.bothside.variants.count.data, knock.in.analysis.dir, "right_ki")


  return(TRUE)
}

