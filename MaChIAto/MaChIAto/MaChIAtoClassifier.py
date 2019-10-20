# -*- coding: utf8 -*-
"""
MaChIAto - Kazuki Nakamae 2019
Software pipeline for the analysis of outcomes of MMEJ-assisted gene knock-in using CRISPR-Cas9 from deep sequencing data
https://github.com/Kazuki-Nakamae/MaChIAto
"""

from MaChIAtoShared import *
from MaChIAtoFunctions import should_stop, reverse_complement, findall_both, is_in

class MaChIAtoClassifier:
    """read and count sequence from the paticular pandas dataframe.

    This reads Aligned_Sequence from the following pandas dataframe
    ---
    Aligned_Sequence Reference_Sequence NHEJ   UNMODIFIED HDR    n_deleted n_inserted n_mutated #Reads  %Reads
    <object>         <object>           <bool> <bool>     <bool> <float64> <float64>  <int64>   <int64> <float64>
    ...
    ---
    Attributes:
    in_pddf: A pandas dataframe given as input argument.
    out_pddf: A pandas dataframe optimized for MaChIAto processing.

    example of out_pddf
    >out_pddf.iloc[1]
    Aligned_Sequence                    TTACATGATGCAGAAAGTTGATATCCCTCCGCTTCTTACTCTTTTT... : CRISPResso output
    Reference_Sequence                  TTACATGATGCAGAAAGTTGATATCCCTCCGCTTCTTACTCTTTTT... : CRISPResso output
    NHEJ                                                                            False : CRISPResso output
    UNMODIFIED                                                                      False : CRISPResso output
    HDR                                                                              True : CRISPResso output
    n_deleted                                                                           2 : CRISPResso output
    n_inserted                                                                         77 : CRISPResso output
    n_mutated                                                                           0 : CRISPResso output
    #Reads                                                                           8090 : CRISPResso output
    %Reads                                                                        6.52614 : CRISPResso output
    OutOut_Aligned_Sequence             TGAAGAAGCTGTGGCAAAAGCTGATAAGCTGGCTGAAGAGCATTCA... : Partial sequence, which is extracted from Aligned_Sequence, between out indicators
    has_outout_indicators                                                            True : Dose the sequence have both out-out indicators ?
    has_outout_rc_indicators                                                        False : Dose the reverse-complement (rc) sequence have both out-out indicators ?
    is_untreated                                                                    False : Is the sequence untreated sequence ?
    is_perfect_outout_ki                                                             True : Is the sequence perfect knock-in sequence ?
    outout_insert_end_point                                                           176 : The acceptable position under the right outout threshold in accordance with perfect knock-in reference sequence
    outout_insert_start_point                                                         177 : The acceptable position under the left outout threshold in accordance with perfect knock-in reference sequence
    outout_left_match_len                                                             176 : The acceptable length under the right outout threshold
    outout_outleft_align_score                                                          1 : The alignment score when <outout_left_match_len> leftside sequences are aligned.
    outout_outright_align_score                                                         1 : The alignment score when <outout_right_match_len> rightside sequences are aligned.
    outout_right_match_len                                                              1 : The acceptable length under the left outout threshold
    Leftside_InOut_Aligned_Sequence     TGAAGAAGCTGTGGCAAAAGCTGATAAGCTGGCTGAAGAGCATTCA... : Partial sequence, which is extracted from OutOut_Aligned_Sequence, between left inout indicators
    has_leftinout_indicators                                                         True : Dose the sequence have left-out indicator and left-in indicator
    Rightside_InOut_Aligned_Sequence    CGATAAGTAAGAGGGGTCTTTGTCCTCTGTACTGTCTCTCTCCTTG... : Partial sequence, which is extracted from OutOut_Aligned_Sequence, between right inout indicators
    has_rightinout_indicators                                                        True : Dose the sequence have right-out indicator and right-in indicator
    is_nearly_leftinout_ki                                                          True : Is the sequence nearly leftside perfect knock-in sequence ?
    leftinout_insert_end_point                                                         59 : The acceptable position under the left left-inout threshold in accordance with left perfect knock-in reference sequence
    leftinout_insert_start_point                                                       60 : The acceptable position under the right left-inout threshold in accordance with left perfect knock-in reference sequence
    leftinout_left_match_len                                                           59 : The acceptable length under the left left-inout threshold
    leftinout_outleft_align_score                                                       1 : The alignment score when <leftinout_left_match_len> leftside sequences are aligned.
    leftinout_outright_align_score                                                      1 : The alignment score when <leftinout_right_match_len> rightside sequences are aligned.
    leftinout_right_match_len                                                           1 : The acceptable length under the right left-inout threshold
    is_nearly_rightinout_ki                                                         True : Is the sequence nearly rightside perfect knock-in sequence ?
    rightinout_insert_end_point                                                        59 : The acceptable position under the left right-inout threshold in accordance with right perfect knock-in reference sequence
    rightinout_insert_start_point                                                      60 : The acceptable position under the right right-inout threshold in accordance with right perfect knock-in reference sequence
    rightinout_left_match_len                                                          59 : The acceptable length under the left right-inout threshold
    rightinout_outleft_align_score                                                      1 : The alignment score when <rightinout_left_match_len> leftside sequences are aligned.
    rightinout_outright_align_score                                                     1 : The alignment score when <rightinout_right_match_len> rightside sequences are aligned.
    rightinout_right_match_len                                                          1 : The acceptable length under the right right-inout threshold

    """

    def __init__(self, in_pddf = None, seq_info_dict = None):
        """ Inits MaChIAtoClassifier. """
        self.in_pddf = None
        self.out_pddf = None
        self.temp_pddf = None # This is used for optimization.
        self.cnt_dict = {
            "total" : None,
            "has_outout_indicators" : None,
            "has_outout_rc_indicators" : None,
            "has_outout_reversible_indicators" : None,
            "has_outout_any_indicators" : None,
            "has_outout_no_indicators" : None,
            "is_nonspecific_outout_indicators" : None,
            "is_untreated" : None,
            "is_treated" : None,
            "is_perfect_outout_ki" : None,
            "has_leftinout_indicators": None,
            "has_multi_leftinout_indicators": None,
            "is_perfect_ki_has_leftinout_indicators" : None,
            "is_nearly_leftinout_ki" : None,
            "has_rightinout_indicators" : None,
            "has_multi_rightinout_indicators" : None,
            "is_perfect_ki_has_rightinout_indicators" : None,
            "is_nearly_rightinout_ki" : None,
            "no_outout_indicators" : None,
            "both_inout_indicators" : None,
            "leftonly_inout_indicators" : None,
            "rightonly_inout_indicators" : None,
            "both_rc_inout_indicators" : None,
            "leftonly_rc_inout_indicators" : None,
            "rightonly_rc_inout_indicators" : None,
            "no_inout_indicators" : None,
            "no_leftinout_indicators" : None,
            "left_perfect_ki" : None,
            "left_nearly_ki" : None,
            "left_nhej_ki" : None,
            "left_imperfect_ki_complex" : None,
            "other_left_mutations" : None,
            "no_rightinout_indicators" : None,
            "right_perfect_ki" : None,
            "right_nearly_ki" : None,
            "right_nhej_ki" : None,
            "right_imperfect_ki_complex" : None,
            "other_right_mutations" : None,
            "no_leftinout_rc_indicators" : None,
            "left_rc_nhej_ki" : None,
            "left_rc_imperfect_ki_complex" : None,
            "other_left_rc_mutations" : None,
            "no_rightinout_rc_indicators" : None,
            "right_rc_nhej_ki" : None,
            "right_rc_imperfect_ki_complex" : None,
            "other_right_rc_mutations" : None,
            "donor_insertion" : None,
            "untreated" : None,
            "unexpected_mutations" : None,
            "overall_edited_substitution" : None,
            "overall_edited_insertion" : None,
            "overall_edited_deletion" : None,
            "LEFT_PRECISE_KNOCK_IN" : None,
            "LEFT_IMPRECISE_KNOCK_IN" : None,
            "LEFT_IMPRECISE KNOCK_IN_COMPLEX" : None,
            "LEFT_OTHER_MUATIONS" : None,
            "LEFT_NO_DETECTED" : None,
            "RIGHT_PRECISE_KNOCK_IN" : None,
            "RIGHT_IMPRECISE_KNOCK_IN" : None,
            "RIGHT_IMPRECISE KNOCK_IN_COMPLEX" : None,
            "RIGHT_OTHER_MUATIONS" : None,
            "RIGHT_NO_DETECTED" : None,
            "PRECISE_KNOCK_IN" : None,
            "IMPRECISE_KNOCK_IN" : None,
            "INSERT_EDITING" : None,
            "DELETION_EDITING" : None,
            "SUBSTITUTION_EDITING" : None,
            "UNTREATED" : None,
            "OTHERS" : None}
        self.read_cnt_dict = {
            "total" : None,
            "has_outout_indicators" : None,
            "has_outout_rc_indicators" : None,
            "has_outout_reversible_indicators" : None,
            "has_outout_any_indicators" : None,
            "has_outout_no_indicators" : None,
            "is_nonspecific_outout_indicators" : None,
            "is_untreated" : None,
            "is_treated" : None,
            "is_perfect_outout_ki" : None,
            "has_leftinout_indicators": None,
            "has_multi_leftinout_indicators": None,
            "is_perfect_ki_has_leftinout_indicators" : None,
            "is_nearly_leftinout_ki" : None,
            "has_rightinout_indicators" : None,
            "has_multi_rightinout_indicators" : None,
            "is_perfect_ki_has_rightinout_indicators" : None,
            "is_nearly_rightinout_ki" : None,
            "no_outout_indicators" : None,
            "both_inout_indicators" : None,
            "leftonly_inout_indicators" : None,
            "rightonly_inout_indicators" : None,
            "both_rc_inout_indicators" : None,
            "leftonly_rc_inout_indicators" : None,
            "rightonly_rc_inout_indicators" : None,
            "no_inout_indicators" : None,
            "no_leftinout_indicators" : None,
            "left_perfect_ki" : None,
            "left_nearly_ki" : None,
            "left_nhej_ki" : None,
            "left_imperfect_ki_complex" : None,
            "other_left_mutations" : None,
            "no_rightinout_indicators" : None,
            "right_perfect_ki" : None,
            "right_nearly_ki" : None,
            "right_nhej_ki" : None,
            "right_imperfect_ki_complex" : None,
            "other_right_mutations" : None,
            "no_leftinout_rc_indicators" : None,
            "left_rc_nhej_ki" : None,
            "left_rc_imperfect_ki_complex" : None,
            "other_left_rc_mutations" : None,
            "no_rightinout_rc_indicators" : None,
            "right_rc_nhej_ki" : None,
            "right_rc_imperfect_ki_complex" : None,
            "other_right_rc_mutations" : None,
            "donor_insertion" : None,
            "untreated" : None,
            "unexpected_mutations" : None,
            "overall_edited_substitution" : None,
            "overall_edited_insertion" : None,
            "overall_edited_deletion" : None,
            "LEFT_PRECISE_KNOCK_IN" : None,
            "LEFT_IMPRECISE_KNOCK_IN" : None,
            "LEFT_IMPRECISE KNOCK_IN_COMPLEX" : None,
            "LEFT_OTHER_MUATIONS" : None,
            "LEFT_NO_DETECTED" : None,
            "RIGHT_PRECISE_KNOCK_IN" : None,
            "RIGHT_IMPRECISE_KNOCK_IN" : None,
            "RIGHT_IMPRECISE KNOCK_IN_COMPLEX" : None,
            "RIGHT_OTHER_MUATIONS" : None,
            "RIGHT_NO_DETECTED" : None,
            "PRECISE_KNOCK_IN" : None,
            "IMPRECISE_KNOCK_IN" : None,
            "INSERT_EDITING" : None,
            "DELETION_EDITING" : None,
            "SUBSTITUTION_EDITING" : None,
            "UNTREATED" : None,
            "OTHERS" : None
        }

        self.seq_info_dict = seq_info_dict.copy()

        if self.seq_info_dict["optimization_mode"] in "OFF":
            info('INPUT :')
            [info('%s : %s' % (str(key), str(self.seq_info_dict[key]))) for key in self.seq_info_dict]

        ###############################################################
        # Internal Setting
        self.oneside_deletion_length = 30
        ###############################################################

        ###############################################################
        if self.seq_info_dict["optimization_mode"] in "OFF":
            info("Reading data table...(1)")
        ###############################################################
        # input pandas dataframe
        self.in_pddf = in_pddf.copy()
        # make sequences without "-"
        self.out_pddf = in_pddf.replace(r'-+', '', regex = True).copy()

        ###############################################################
        if self.seq_info_dict["optimization_mode"] in "OFF":
            info("Defining specific sequences for identifying alleles...(2)")
        ###############################################################

        # define pattrn for reconizing KO/KI alleles
        self._define_sequence_indentifers()

        if self.seq_info_dict["optimization_mode"] in "OFF":
            info('Out-Left indicator sequence : %s' % self.left_out_indicator_seq)
            info('Leftside junction : %s' % self.left_junction_seq)
            info('Rightside junction : %s' % self.right_junction_seq)
            info('Out-Right indicator sequence : %s' % self.right_out_indicator_seq)
            info('Leftside homology arm : %s' % self.left_homologyarm_seq)
            info('Rightside homology arm : %s' % self.right_homologyarm_seq)
            info('Untreated sequence : %s' % self.untreated_seq)
            info('In-Left indicator sequence : %s' % self.left_in_indicator_seq)
            info('In-Right indicator sequence :%s' % self.right_in_indicator_seq)
            info('Left cutout sequence :%s' % self.left_cutout_seq)
            info('Right cutout sequence :%s' % self.right_cutout_seq)
            info('Left cutout maximum deletion sequence :%s' % self.left_cutout_maxdel_seq)
            info('Right cutout maximum deletion sequence :%s' % self.right_cutout_maxdel_seq)

        ###############################################################
        if self.seq_info_dict["optimization_mode"] in "OFF":
            info("Making reference sequences using the specific sequences...(3)")
        ###############################################################

        self._make_references()

        ###############################################################
        if self.seq_info_dict["optimization_mode"] in "OFF":
            info("Prepairing for recognization...(4)")
        ###############################################################

        # make regex for recognization
        self._make_regex()

        ###############################################################
        if self.seq_info_dict["optimization_mode"] in "OFF":
            info("Classifying alleles...(5)")
        ###############################################################

        self._scan_whole_seq() # When optimizing parameter, this only method is to be used.
        if self.seq_info_dict["optimization_mode"] in "OFF":
            self._scan_left_seq()
            self._scan_right_seq()
            self._scan_left_rc_seq()
            self._scan_right_rc_seq()
            self._scan_mutation_seq()
            self._add_global_classification_labels()

    #--- Constructor END ------------------------------------------------------------------------------------#

    def _define_sequence_indentifers(self):
        """define specific sequences for identifying alleles considering parameters """
        try:
            # identify outside indicators and junctions
            perfect_ki_regex_str = "([atcgATCG]{"\
                + str(self.seq_info_dict["left_out_indicator_len"])\
                + "})([atcgATCG]{"\
                + str(self.seq_info_dict["left_junction_len"])\
                + "})" + self.seq_info_dict["donor_seq"]\
                + "([atcgATCG]{"\
                + str(self.seq_info_dict["right_junction_len"])\
                + "})([atcgATCG]{"\
                + str(self.seq_info_dict["right_out_indicator_len"])\
                + "})"
            perfect_ki_regex = re.compile(perfect_ki_regex_str)
            identifying_ki_factors = perfect_ki_regex.search(self.seq_info_dict["expected_ki_amplicon_seq"])
            if identifying_ki_factors is None:
                raise ReferenceRegexException('expected_ki_amplicon_seq does not exist in expected perfect knock-in sequence :%s' % ' '.join([perfect_ki_regex_str, self.seq_info_dict["expected_ki_amplicon_seq"]]))
            
            # OUT-LEFT INDICATOR
            self.left_out_indicator_seq = identifying_ki_factors.group(1)
            # LEFTSIDE JUNCTION
            self.left_junction_seq = identifying_ki_factors.group(2)
            # RIGHTSIDE JUNCTION
            self.right_junction_seq = identifying_ki_factors.group(3)
            # OUT-RIGHT INDICATOR
            self.right_out_indicator_seq = identifying_ki_factors.group(4)
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            # identify homology arms
            homologybased_ki_regex_str = "("\
                + self.left_out_indicator_seq\
                + ")([atcgATCG]*)"\
                + "([atcgATCG]{"\
                + str(self.seq_info_dict["left_homology_arm_len"])\
                + "})(" + self.seq_info_dict["donor_seq"]\
                + ")([atcgATCG]{"\
                + str(self.seq_info_dict["right_homology_arm_len"])\
                + "})([atcgATCG]*)("\
                + self.right_out_indicator_seq\
                + ")"
            homologybased_ki_regex = re.compile(homologybased_ki_regex_str)
            identifying_homologybased_ki_factors = homologybased_ki_regex.search(self.seq_info_dict["expected_ki_amplicon_seq"])
            if identifying_homologybased_ki_factors is None:
                raise ReferenceRegexException('expected_ki_amplicon_seq does not exist in expected homology-based knock-in sequence :%s' % ' '.join([homologybased_ki_regex_str, self.seq_info_dict["expected_ki_amplicon_seq"]]))
            
            # LEFTSIDE MARGIN
            self.left_margin_seq = identifying_homologybased_ki_factors.group(2)
            # LEFTSIDE HOMOLOGY ARM
            self.left_homologyarm_seq = identifying_homologybased_ki_factors.group(3)
            # RIGHTSIDE HOMOLOGY ARM
            self.right_homologyarm_seq = identifying_homologybased_ki_factors.group(5)
            # RIGHTSIDE MARGIN
            self.right_margin_seq = identifying_homologybased_ki_factors.group(6)

            # identify untreated reference, which is surrounded by indicators
            untreated_regex_str = "("\
                + self.left_out_indicator_seq \
                + ")" + "([atcgATCG]+)" \
                + "(" + self.right_out_indicator_seq \
                + ")"
            untreated_regex = re.compile(untreated_regex_str)
            identifying_untreated_factors = untreated_regex.search(self.seq_info_dict["amplicon_seq"])
            if identifying_untreated_factors is None:
                raise ReferenceRegexException('amplicon_seq does not exist in expected untreated sequence :%s' % ' '.join([untreated_regex_str, self.seq_info_dict["amplicon_seq"]]))

            # UNTREATED SEQUENCE
            self.untreated_seq = identifying_untreated_factors.group(2)

            # identify expected cutout sequence, which is surrounded by indicators
            cutout_regex_str = "("\
                + self.left_out_indicator_seq \
                + "[atcgATCG]+" \
                + self.seq_info_dict["guide_seq"][:17] \
                + ")" \
                + "(" + self.seq_info_dict["guide_seq"][17:] \
                + "[atcgATCG]+" \
                + self.right_out_indicator_seq \
                + ")"
            cutout_regex = re.compile(cutout_regex_str)
            identifying_cutout_factors = cutout_regex.search(self.seq_info_dict["amplicon_seq"])

            if identifying_cutout_factors is None:
                # Try reverse complement version if sequence is not found.
                self.seq_info_dict["guide_seq"] = reverse_complement(self.seq_info_dict["guide_seq"])
                cutout_regex_str = "("\
                    + self.left_out_indicator_seq \
                    + "[atcgATCG]+" \
                    + self.seq_info_dict["guide_seq"][:17] \
                    + ")" \
                    + "(" + self.seq_info_dict["guide_seq"][17:] \
                    + "[atcgATCG]+" \
                    + self.right_out_indicator_seq \
                    + ")"
                cutout_regex = re.compile(cutout_regex_str)
                identifying_cutout_factors = cutout_regex.search(self.seq_info_dict["amplicon_seq"])
                if identifying_cutout_factors is None:
                    raise ReferenceRegexException('amplicon_seq does not exist in expected cutout sequence :%s' % ' '.join([cutout_regex_str, self.seq_info_dict["amplicon_seq"]]))
            
            # CUTOUT SEQUENCE
            self.left_cutout_seq = identifying_cutout_factors.group(1)
            self.right_cutout_seq = identifying_cutout_factors.group(2)
            # MAXIMUM DELETION CUTOUT SEQUENCE
            self.left_cutout_maxdel_seq = self.left_cutout_seq[ : len(self.left_cutout_seq) - self.oneside_deletion_length]
            self.right_cutout_maxdel_seq = self.right_cutout_seq[self.oneside_deletion_length : ]


            # identify inside indicators, which are only in knock-in alleles
            ki_insert_regex_str = "(^[atcgATCG]{" \
                + str(self.seq_info_dict["left_in_indicator_len"]) \
                + "})" + "([atcgATCG]+)" \
                + "([atcgATCG]{" \
                + str(self.seq_info_dict["right_in_indicator_len"]) \
                + "}$)"
            ki_insert_regex = re.compile(ki_insert_regex_str)
            identifying_ki_insert_factors = ki_insert_regex.search(self.seq_info_dict["donor_seq"])
            if identifying_ki_insert_factors is None:
                raise ReferenceRegexException('donor_seq does not exist in expected donor (insert) sequence :%s' % ' '.join([ki_insert_regex_str, self.seq_info_dict["donor_seq"]]))

            # IN-LEFT INDICATOR
            self.left_in_indicator_seq = identifying_ki_insert_factors.group(1)
            # IN-RIGHT INDICATOR
            self.right_in_indicator_seq = identifying_ki_insert_factors.group(3)
        except ReferenceRegexException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Reference error, please check your input.\n\nERROR: %s' % e)
            sys.exit(1)

    def _make_references(self):
        # PERFECT KNOCK-IN REFERENCE SEQUENCE
        self.indicator_untreated_ref_seq = (self.left_out_indicator_seq \
            + self.untreated_seq \
            + self.right_out_indicator_seq)

        # PERFECT KNOCK-IN REFERENCE SEQUENCE
        self.knockin_ref_seq = (self.left_out_indicator_seq \
            + self.left_junction_seq \
            + self.seq_info_dict["donor_seq"] \
            + self.right_junction_seq \
            + self.right_out_indicator_seq)

        # The reverse complement sequence of specific sequences
        self.rev_left_out_indicator_seq = reverse_complement(self.left_out_indicator_seq)
        self.rev_right_out_indicator_seq = reverse_complement(self.right_out_indicator_seq)
        self.rev_left_junction_seq = reverse_complement(self.left_junction_seq)
        self.rev_right_junction_seq = reverse_complement(self.right_junction_seq)
        self.rev_untreated_seq = reverse_complement(self.untreated_seq)
        self.rev_left_homologyarm_seq  = reverse_complement(self.left_homologyarm_seq)
        self.rev_right_homologyarm_seq = reverse_complement(self.right_homologyarm_seq)
        self.rev_left_margin_seq = reverse_complement(self.left_margin_seq)
        self.rev_right_margin_seq = reverse_complement(self.right_margin_seq)
        self.rev_left_in_indicator_seq = reverse_complement(self.left_in_indicator_seq)
        self.rev_right_in_indicator_seq = reverse_complement(self.right_in_indicator_seq)

        # Note : A system limitation
        # You might think 20 bp deletion is too long.
        # I take outer sequence flanking homology-arm (e.g. the rest sequence of sgRNA targeting site) into consideration.
        # NHEJ-assisted deletion insert
        left_nhej_deletion_size = len(self.left_homologyarm_seq)
        right_nhej_deletion_size = len(self.right_homologyarm_seq)
        self.left_nhej_del_homologyarm_seq  = self.left_homologyarm_seq[left_nhej_deletion_size:]
        self.right_nhej_del_homologyarm_seq = self.right_homologyarm_seq[:len(self.right_homologyarm_seq) - right_nhej_deletion_size]
        self.opposite_left_nhej_del_homologyarm_seq  = self.rev_right_homologyarm_seq[right_nhej_deletion_size:]
        self.opposite_right_nhej_del_homologyarm_seq = self.rev_left_homologyarm_seq[:len(self.rev_left_homologyarm_seq) - left_nhej_deletion_size]

        # LEFTSIDE PERFECT KNOCK-IN REFERENCE SEQUENCE
        self.left_knockin_ref_seq = self.left_out_indicator_seq \
            + self.left_junction_seq \
            + self.left_in_indicator_seq
        # RIGHTSIDE PERFECT KNOCK-IN REFERENCE SEQUENCE
        self.right_knockin_ref_seq = self.right_in_indicator_seq \
            + self.right_junction_seq \
            + self.right_out_indicator_seq
        # OPPOSITE LEFTSIDE PERFECT KNOCK-IN REFERENCE SEQUENCE
        self.opposite_left_knockin_ref_seq = self.left_out_indicator_seq \
            + self.left_margin_seq \
            + self.rev_right_homologyarm_seq \
            + self.rev_right_in_indicator_seq
        # OPPOSITE RIGHTSIDE PERFECT KNOCK-IN REFERENCE SEQUENCE
        self.opposite_right_knockin_ref_seq = self.rev_left_in_indicator_seq \
            + self.rev_left_homologyarm_seq \
            + self.right_margin_seq \
            + self.right_out_indicator_seq
        
        # LEFTSIDE NHEJ-ASSISTED KNOCK-IN REFERENCE SEQUENCE
        self.left_nhej_knockin_ref_seq = self.left_out_indicator_seq \
            + self.left_margin_seq \
            + self.left_homologyarm_seq \
            + self.left_homologyarm_seq \
            + self.left_in_indicator_seq
        # RIGHTSIDE NHEJ-ASSISTED KNOCK-IN REFERENCE SEQUENCE
        self.right_nhej_knockin_ref_seq = self.right_in_indicator_seq \
            + self.right_homologyarm_seq \
            + self.right_homologyarm_seq \
            + self.right_margin_seq \
            + self.right_out_indicator_seq
        # OPPOSITE LEFTSIDE NHEJ-ASSISTED KNOCK-IN REFERENCE SEQUENCE
        self.opposite_left_nhej_knockin_ref_seq = self.left_out_indicator_seq \
            + self.left_margin_seq \
            + self.rev_right_homologyarm_seq \
            + self.rev_right_homologyarm_seq \
            + self.rev_right_in_indicator_seq
        # OPPOSITE RIGHTSIDE NHEJ-ASSISTED KNOCK-IN REFERENCE SEQUENCE
        self.opposite_right_nhej_knockin_ref_seq = self.rev_left_in_indicator_seq \
            + self.rev_left_homologyarm_seq \
            + self.rev_left_homologyarm_seq \
            + self.right_margin_seq \
            + self.right_out_indicator_seq
        
        # LEFTSIDE NHEJ-ASSISTED DELETION KNOCK-IN REFERENCE SEQUENCE
        self.left_nhej_del_knockin_ref_seq = self.left_out_indicator_seq \
            + self.left_margin_seq \
            + self.left_nhej_del_homologyarm_seq \
            + self.left_homologyarm_seq \
            + self.left_in_indicator_seq
        # RIGHTSIDE NHEJ-ASSISTED DELETION KNOCK-IN REFERENCE SEQUENCE
        self.right_nhej_del_knockin_ref_seq = self.right_in_indicator_seq \
            + self.right_homologyarm_seq \
            + self.right_nhej_del_homologyarm_seq \
            + self.right_margin_seq \
            + self.right_out_indicator_seq
        # OPPOSITE LEFTSIDE NHEJ-ASSISTED DELETION KNOCK-IN REFERENCE SEQUENCE
        self.opposite_left_nhej_del_knockin_ref_seq = self.left_out_indicator_seq \
            + self.left_margin_seq \
            + self.opposite_left_nhej_del_homologyarm_seq \
            + self.rev_right_homologyarm_seq \
            + self.rev_right_in_indicator_seq
        # OPPOSITE RIGHTSIDE NHEJ-ASSISTED DELETION KNOCK-IN REFERENCE SEQUENCE
        self.opposite_right_nhej_del_knockin_ref_seq = self.rev_left_in_indicator_seq \
            + self.rev_left_homologyarm_seq \
            + self.opposite_right_nhej_del_homologyarm_seq \
            + self.right_margin_seq \
            + self.right_out_indicator_seq
	
    def _make_regex(self):
        # regex of out-out sequence
        self.targetseq_regex = re.compile(self.left_out_indicator_seq \
            + "[ATCG]*" \
            + self.right_out_indicator_seq)
        self.rc_targetseq_regex = re.compile(self.rev_right_out_indicator_seq \
            + "[ATCG]*" \
            + self.rev_left_out_indicator_seq)
        self.indicator_untreated_regex = re.compile(self.left_out_indicator_seq \
            + self.untreated_seq \
            + self.right_out_indicator_seq)
        self.rc_indicator_untreated_regex = re.compile(self.rev_right_out_indicator_seq \
            + self.rev_untreated_seq \
            + self.rev_left_out_indicator_seq)
        self.all_seq_regex = re.compile("[ATCG]+")

        # regex of out-in sequence
        self.insertseq_regex = re.compile(self.left_out_indicator_seq \
            + "([ATCG]*)" + self.right_out_indicator_seq)

        self.left_knockin_regex = re.compile(self.left_out_indicator_seq \
            + "([ATCG]*)" + self.left_in_indicator_seq)# greedy matching
        self.opp_left_knockin_regex = re.compile(reverse_complement(self.left_in_indicator_seq) \
            + "([ATCG]*)" + reverse_complement(self.left_out_indicator_seq))# greedy matching

        self.left_rc_knockin_regex = re.compile(self.left_out_indicator_seq \
            + "([ATCG]*)" + reverse_complement(self.right_in_indicator_seq))# greedy matching
        self.opp_left_rc_knockin_regex = re.compile(self.right_in_indicator_seq \
            + "([ATCG]*)" + reverse_complement(self.left_out_indicator_seq))# greedy matching


        self.right_knockin_regex = re.compile(self.right_in_indicator_seq \
            + "([ATCG]*?)" + self.right_out_indicator_seq)# non-greedy matching
        self.opp_right_knockin_regex = re.compile(reverse_complement(self.right_out_indicator_seq) \
            + "([ATCG]*?)" + reverse_complement(self.right_in_indicator_seq))# non-greedy matching

        self.right_rc_knockin_regex = re.compile(reverse_complement(self.left_in_indicator_seq) \
            + "([ATCG]*?)" + self.right_out_indicator_seq)# non-greedy matching
        self.opp_right_rc_knockin_regex = re.compile(reverse_complement(self.right_out_indicator_seq) \
            + "([ATCG]*?)" + self.left_in_indicator_seq)# non-greedy matching
    
     #--- helper methods ------------------------------------------------------------------------------------#

    def _register_count(self, key, pddf):
        self.cnt_dict[key] = len(pddf)
        self.read_cnt_dict[key] = pddf["#Reads"].sum()
    
    def _show_count(self, key):
         info('%s(%s) / %s(%s) allele types (reads) have "%s" flag.' % \
                (str(self.cnt_dict[key]), str(self.read_cnt_dict[key]), str(self.cnt_dict["total"]), str(self.read_cnt_dict["total"]), str(key)))


    #--- Bit-by-Bit Alignment ------------------------------------------------------------------------------------#

    def _align_bit_by_bit(self, readseq, refseq, left_align_threshold, right_align_threshold):
        """
        Find a point matching the reference sequence

        The case of out-out sequence
        --------------------------//#########&&&&&&Insert&&&############//-------------------------- sequence
        --------------------------//###############Insert###############//-------------------------- reference sequence
        <-----------left_match_len----------><-------------><---------right_match_len---------------> range
                            insert_start_point^               ^insert_end_point						  points
        
        The case of out-in sequence
        -------------------------------------&&&&&&&&&&&&&&&//###############Insert################# sequence
        ----------------------------------------------------//###############Insert################# reference sequence
        <-----------left_match_len----------><-------------><---------right_match_len---------------> range
                            insert_start_point^               ^insert_end_point						  points

        '&' indicating indel or substitution
        '//' indicating boundary between genome and Knock-In insert

        Args:
            readseq : String represents aligned sequence.
            refseq : String represents reference sequence.

        Returns:
            left_match_len : The length matching left-side(5') reference sequence
            right_match_len	: The length matching right-side(3') reference sequence
            insert_start_point : The relative positon starting insert in left-side (start of left indicator : 0)
            insert_end_point : The relative positon starting insert in right-side (start of left indicator : 0)

        Raises:
            UnexpectedException: An error occurred processing under impossible condition.
        """
        # Margin bp for judgement
        # 1 <= MARGIN
        MARGIN = 5

        readlen = len(readseq)
        reflen = len(refseq)
        minlen = min(readlen, reflen)
        maxlen = max(readlen, reflen)

        left_match_len=0
        right_match_len=0
        insert_start_point=0 # reference-based position
        insert_end_point=reflen # reference-based position
        left_align_score=1.0
        right_align_score=1.0
        is_perfect = False # flag indicating perfect sequence under the setting

        #Search for the maximum point that matches refseq in 5 '-> 3' direction
        for bothind in range(1000):
                
            if bothind >= minlen:
                insert_start_point = reflen # force it to be maxlen due to (if  insert_start_point > insert_end_point:)
                break

            if readseq[:bothind] in refseq[:bothind]:
                left_match_len = len(refseq[:bothind])
                insert_start_point = bothind
            else:
                left_align_score = pairwise2.align.globalms(readseq[:bothind], refseq[:bothind], 2, -1, -10, -.5)[0][2] \
                / float(len(refseq[:bothind]) * 2)
                # I wish to reduce times using pairwise2.align.globalms().

                if left_align_score < left_align_threshold:
                    break
                else:
                    left_match_len = len(refseq[:bothind])
                    insert_start_point = bothind

        #Search for the maximum point that matches refseq in 3 '-> 5' direction
        for bothind in range(1000):
            readind = readlen - bothind
            refind = reflen - bothind

            if readind < 1 or refind < 1:
                break

            if readseq[readind:] in refseq[refind:]:
                right_match_len = len(refseq[refind:])
                if readlen <= reflen:
                    insert_end_point = refind
                if readlen > reflen:
                    insert_end_point = readind
            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                right_align_score = pairwise2.align.globalms(readseq[readind:], refseq[refind:], 2, -1, -10, -.5)[0][2] \
                    / float(len(refseq[refind:]) * 2)
                # I wish to reduce times using pairwise2.align.globalms().
                if right_align_score < right_align_threshold:
                    break
                else:
                    right_match_len = len(refseq[refind:])
                    if readlen <= reflen:
                        insert_end_point = refind
                    if readlen > reflen:
                        insert_end_point = readind
            
            #If the matching position exceeds the maximum point of the 5 '-> 3' direction search, it ends on the spot
            if  insert_start_point > (insert_end_point + 1) + MARGIN:
                right_match_len = len(refseq[refind + 1:])
                insert_end_point = insert_start_point
                is_perfect = True
                break
            

        #if not is_perfect:
        #	import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        return left_match_len,\
                right_match_len,\
                insert_start_point,\
                insert_end_point,\
                left_align_score,\
                right_align_score,\
                is_perfect

    #--- Generator for Optimization ------------------------------------------------------------------------------#

    def generator_first_optimization(self, left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len):
        """perform loop for first optimization. It is useful for reducing optimization time."""
        try:
            for i in range(10000): # This is considered as infinite loop.
                
                self.seq_info_dict["left_out_indicator_len"] = left_out_indicator_len
                self.seq_info_dict["right_out_indicator_len"] = right_out_indicator_len
                self.seq_info_dict["left_junction_len"] = left_junction_len
                self.seq_info_dict["right_junction_len"] = right_junction_len

                # reset parameter
                self._define_sequence_indentifers()
                self._make_references()
                self._make_regex()

                self.out_pddf = self.temp_pddf.copy()
                
                # OUT-OUT Indicator
                self.out_pddf = self.out_pddf.join(self.out_pddf.Aligned_Sequence.apply(self._has_out_indicators)).copy()
                
                self.cnt_dict["total"] = len(self.out_pddf)
                self.cnt_dict["has_outout_indicators"] = len(self.out_pddf.query('has_outout_indicators == True and has_outout_rc_indicators == False'))
                self.cnt_dict["has_outout_rc_indicators"] = len(self.out_pddf.query('has_outout_indicators == False and has_outout_rc_indicators == True'))
                self.cnt_dict["has_outout_reversible_indicators"] = len(self.out_pddf.query('has_outout_indicators == True and has_outout_rc_indicators == True'))
                self.cnt_dict["has_outout_no_indicators"] = len(self.out_pddf.query('has_outout_indicators == False and has_outout_rc_indicators == False'))
                self.cnt_dict["is_nonspecific_outout_indicators"] = len(self.out_pddf.query('is_nonspecific_outout_indicators == True'))
                self.cnt_dict["has_outout_any_indicators"] = self.cnt_dict["has_outout_indicators"] \
                    + self.cnt_dict["has_outout_rc_indicators"] \
                    + self.cnt_dict["has_outout_reversible_indicators"]

                yield

        except Exception as e:
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def generator_second_optimization(self, outleft_align_threshold, outright_align_threshold):
        """perform loop for second optimization. It is useful for reducing optimization time."""
        try:
            for i in range(10000): # This is considered as infinite loop.

                self.seq_info_dict["outleft_align_threshold"] = outleft_align_threshold
                self.seq_info_dict["outright_align_threshold"] = outright_align_threshold

                # reset parameter
                self._define_sequence_indentifers()
                self._make_references()
                self._make_regex()

                self.out_pddf = self.temp_pddf.copy()

                # Detecting Untreated
                self.out_pddf = self.out_pddf.join(self.out_pddf.apply(self._is_untreated_seq, axis = 1)).copy()

                self.cnt_dict["is_untreated"] = len(self.out_pddf.query('is_untreated == True'))
                self.cnt_dict["is_treated"] = self.cnt_dict["has_outout_any_indicators"] - self.cnt_dict["is_untreated"]

                # OUT-OUT BbB Alignment
                self.out_pddf = self.out_pddf.join(self.out_pddf.apply(self._scan_outout_junction, axis = 1)).copy()
                
                self.cnt_dict["is_perfect_outout_ki"] = len(self.out_pddf.query('is_perfect_outout_ki == True'))

                yield

        except Exception as e:
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def generator_third_optimization(self, left_in_indicator_len, right_in_indicator_len):
        """perform loop for third optimization. It is useful for reducing optimization time."""
        try:
            for i in range(10000): # This is considered as infinite loop.

                self.seq_info_dict["left_in_indicator_len"] = left_in_indicator_len
                self.seq_info_dict["right_in_indicator_len"] = right_in_indicator_len

                # reset parameter
                self._define_sequence_indentifers()
                self._make_references()
                self._make_regex()

                self.out_pddf = self.temp_pddf.copy()

                # IN-OUT Indicator
                self.out_pddf = self.out_pddf.join(self.out_pddf.apply(self._has_rightinout_indicators, axis = 1)).copy()
                self.out_pddf = self.out_pddf.join(self.out_pddf.apply(self._has_leftinout_indicators, axis = 1)).copy()

                self.cnt_dict["has_leftinout_indicators"] = len(self.out_pddf.query('has_leftinout_indicators == True'))
                self.cnt_dict["has_multi_leftinout_indicators"] = len(self.out_pddf.query('has_multi_leftinout_indicators == True'))
                self.cnt_dict["has_rightinout_indicators"] = len(self.out_pddf.query('has_rightinout_indicators == True'))
                self.cnt_dict["has_multi_rightinout_indicators"] = len(self.out_pddf.query('has_multi_rightinout_indicators == True'))

                yield

        except Exception as e:
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def generator_forth_optimization(self, inleft_align_threshold, inright_align_threshold):
        """perform loop for forth optimization. It is useful for reducing optimization time."""
        try:
            for i in range(10000): # This is considered as infinite loop.
                
                self.seq_info_dict["inleft_align_threshold"] = inleft_align_threshold
                self.seq_info_dict["inright_align_threshold"] = inright_align_threshold

                # reset parameter
                self._define_sequence_indentifers()
                self._make_references()
                self._make_regex()

                self.out_pddf = self.temp_pddf.copy()

                # IN-OUT Threshold
                self.out_pddf = self.out_pddf.join(self.out_pddf.apply(self._scan_leftinout_junction, axis = 1)).copy()
                self.out_pddf = self.out_pddf.join(self.out_pddf.apply(self._scan_rightinout_junction, axis = 1)).copy()

                self.cnt_dict["is_nearly_leftinout_ki"] = len(self.out_pddf.query('is_nearly_leftinout_ki == True'))
                self.cnt_dict["is_nearly_rightinout_ki"] = len(self.out_pddf.query('is_nearly_rightinout_ki == True'))
                
                yield

        except Exception as e:
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    #--- Whole Classification ------------------------------------------------------------------------------------#

    def _has_out_indicators(self, seq):

        match_target_seq = findall_both(self.targetseq_regex, self.rc_targetseq_regex, seq)
        match_rc_target_seq = findall_both(self.rc_targetseq_regex, self.targetseq_regex, seq)

        # debug
        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        if len(match_target_seq) > 0:# <Dose the seqeuence have both indicators ?>
                
            if len(match_rc_target_seq) > 0:# <Dose the seqeuence have both indicators in opposite direction ?>

                if self.seq_info_dict["optimization_mode"] is None:
                    warn('The indicators are unspecific. Try it with longer length of indicators')
                    should_stop()
                else:
                    # If above command run in the optimization mode, the optimization will be disruputed.
                    pass
                
                return(pd.Series({'has_outout_indicators' : True, 'has_outout_rc_indicators' : True, 'is_nonspecific_outout_indicators' : True, 'OutOut_Aligned_Sequence' : match_target_seq[0]}))

            if len(match_target_seq) > 1:# <Dose the seqeuence have both indicators in multiple sites ?>
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                return(pd.Series({'has_outout_indicators' : True, 'has_outout_rc_indicators' : False, 'is_nonspecific_outout_indicators' : True, 'OutOut_Aligned_Sequence' : match_target_seq[0]}))

            return(pd.Series({'has_outout_indicators' : True, 'has_outout_rc_indicators' : False, 'is_nonspecific_outout_indicators' : False, 'OutOut_Aligned_Sequence' : match_target_seq[0]}))

        else:

            if len(match_rc_target_seq) > 0:# <Dose the seqeuence have both indicators in opposite direction ?>

                if len(match_rc_target_seq) > 1:# <Dose the seqeuence have both indicators in multiple sites ?>
                    return(pd.Series({'has_outout_indicators' : False, 'has_outout_rc_indicators' : True, 'is_nonspecific_outout_indicators' : True, 'OutOut_Aligned_Sequence' : match_rc_target_seq[0]}))


                return(pd.Series({'has_outout_indicators' : False, 'has_outout_rc_indicators' : True, 'is_nonspecific_outout_indicators' : False, 'OutOut_Aligned_Sequence' : reverse_complement(match_rc_target_seq[0])}))

            else:
                    
                return(pd.Series({'has_outout_indicators' : False, 'has_outout_rc_indicators' : False, 'is_nonspecific_outout_indicators' : False, 'OutOut_Aligned_Sequence' : None}))

    def _is_untreated_seq(self, pddf_row):
        """ask whether Aligned_Sequence is untreated sequence"""
        try:
            
            seq = pddf_row["OutOut_Aligned_Sequence"]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            if seq is None:

                return pd.Series({'is_untreated' : None})

            else:

                # Length is equal ?
                if(len(seq) != len(self.indicator_untreated_ref_seq)):
                    return pd.Series({'is_untreated' : False})

                # Find points matching prefect knock-in reference sequence
                temp1,\
                temp2,\
                temp3,\
                temp4,\
                temp5,\
                temp6,\
                is_untreated\
                = self._align_bit_by_bit(seq, self.indicator_untreated_ref_seq,\
                    self.seq_info_dict["outleft_align_threshold"],\
                    self.seq_info_dict["outright_align_threshold"])
                return pd.Series({'is_untreated' : is_untreated})

        except UnexpectedException as e:
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)
	
    def _scan_outout_junction(self, pddf_row):
        """ Find points matching prefect knock-in reference sequence """
        try:
            seq = pddf_row["OutOut_Aligned_Sequence"]
            is_untreated = pddf_row["is_untreated"]

            if ((seq is None) or (is_untreated)):

                # Even if you change the following values to None, they will be forced to be NaN.
                # I hate such an implicit type conversion. I will put 0 in it, instead.
                return pd.Series({'outout_left_match_len' : 0,\
                        'outout_right_match_len' : 0,\
                        'outout_insert_start_point' : 0,\
                        'outout_insert_end_point': 0,\
                        'outout_outleft_align_score': 0,\
                        'outout_outright_align_score': 0,\
                        'is_len_perfect_outout_ki' : False,\
                        'is_perfect_outout_ki': False})
            else:
                
                # Find points matching prefect knock-in reference sequence
                outout_left_match_len,\
                outout_right_match_len,\
                outout_insert_start_point,\
                outout_insert_end_point,\
                left_align_score,\
                right_align_score,\
                is_perfect_outout_ki\
                = self._align_bit_by_bit(seq, self.knockin_ref_seq,\
                    self.seq_info_dict["outleft_align_threshold"],\
                    self.seq_info_dict["outright_align_threshold"])

                is_len_perfect_outout_ki = True
                # Length is equal ?
                if len(seq) != len(self.knockin_ref_seq):
                    is_perfect_outout_ki = False
                    is_len_perfect_outout_ki = False


                return pd.Series({'outout_left_match_len' : outout_left_match_len,\
                        'outout_right_match_len' : outout_right_match_len,\
                        'outout_insert_start_point' : outout_insert_start_point,\
                        'outout_insert_end_point': outout_insert_end_point,\
                        'outout_outleft_align_score': left_align_score,\
                        'outout_outright_align_score': right_align_score,\
                        'is_len_perfect_outout_ki' : is_len_perfect_outout_ki,\
                        'is_perfect_outout_ki' : is_perfect_outout_ki})

        except UnexpectedException as e:
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _has_leftinout_indicators(self, pddf_row):
        
        seq = pddf_row["OutOut_Aligned_Sequence"]

        if seq is None:
                
            return(pd.Series({'has_leftinout_indicators' : False, 'Leftside_InOut_Aligned_Sequence' : None , 'has_multi_leftinout_indicators' : False}))

        match_target_seq = findall_both(self.left_knockin_regex, self.opp_left_knockin_regex, seq)

        if len(match_target_seq) < 1:# <Dose the seqeuence have both indicators ?>
                
            return(pd.Series({'has_leftinout_indicators' : False, 'Leftside_InOut_Aligned_Sequence' : None, 'has_multi_leftinout_indicators' : False}))
        
        elif len(match_target_seq) > 1:
            
            return(pd.Series({'has_leftinout_indicators' : True, 'Leftside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_leftinout_indicators' : True}))

        else:

            return(pd.Series({'has_leftinout_indicators' : True, 'Leftside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_leftinout_indicators' : False}))

    def _has_rc_leftinout_indicators(self, pddf_row):
        
        seq = pddf_row["OutOut_Aligned_Sequence"]

        if seq is None:
                
            return(pd.Series({'has_rc_leftinout_indicators' : False, 'RC_Leftside_InOut_Aligned_Sequence' : None, 'has_multi_rc_leftinout_indicators' : False}))

        match_target_seq = findall_both(self.left_rc_knockin_regex, self.opp_left_rc_knockin_regex, seq)

        if len(match_target_seq) < 1:# <Dose the seqeuence have both indicators ?>
                
            return(pd.Series({'has_rc_leftinout_indicators' : False, 'RC_Leftside_InOut_Aligned_Sequence' : None, 'has_multi_rc_leftinout_indicators' : False}))

        elif len(match_target_seq) > 1:

            return(pd.Series({'has_rc_leftinout_indicators' : True, 'RC_Leftside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_rc_leftinout_indicators' : True}))

        else:

            return(pd.Series({'has_rc_leftinout_indicators' : True, 'RC_Leftside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_rc_leftinout_indicators' : False}))

    def _scan_leftinout_junction(self, pddf_row):
        """ask whether Leftside_InOut_Aligned_Sequence has leftside in and out indicators"""
        try:
            seq = pddf_row["Leftside_InOut_Aligned_Sequence"]
            has_leftinout_indicators = pddf_row["has_leftinout_indicators"]
            is_perfect_outout_ki = pddf_row["is_perfect_outout_ki"]

            # Focusing on treated sequences only
            if not has_leftinout_indicators or is_perfect_outout_ki:

                # Even if you change the following values to None, they will be forced to be NaN.
                # I hate such an implicit type conversion. I will put 0 in it, instead.
                return pd.Series({'leftinout_left_match_len' : 0,\
                        'leftinout_right_match_len' : 0,\
                        'leftinout_insert_start_point' : 0,\
                        'leftinout_insert_end_point': 0,\
                        'leftinout_outleft_align_score': 0,\
                        'leftinout_outright_align_score': 0,\
                        'is_nearly_leftinout_ki': False})

            else:

                # Length is equal ?
                if len(seq) != len(self.left_knockin_ref_seq):
                    return pd.Series({'leftinout_left_match_len' : 0,\
                        'leftinout_right_match_len' : 0,\
                        'leftinout_insert_start_point' : 0,\
                        'leftinout_insert_end_point': 0,\
                        'leftinout_outleft_align_score': 0,\
                        'leftinout_outright_align_score': 0,\
                        'is_nearly_leftinout_ki': False})

                leftinout_left_match_len,\
                leftinout_right_match_len,\
                leftinout_insert_start_point,\
                leftinout_insert_end_point,\
                leftinout_outleft_align_score,\
                leftinout_outright_align_score,\
                is_nearly_leftinout_ki\
                = self._align_bit_by_bit(seq, self.left_knockin_ref_seq,\
                    self.seq_info_dict["outleft_align_threshold"],\
                    self.seq_info_dict["inleft_align_threshold"])

                return pd.Series({'leftinout_left_match_len' : leftinout_left_match_len,\
                        'leftinout_right_match_len' : leftinout_right_match_len,\
                        'leftinout_insert_start_point' : leftinout_insert_start_point,\
                        'leftinout_insert_end_point': leftinout_insert_end_point,\
                        'leftinout_outleft_align_score': leftinout_outleft_align_score,\
                        'leftinout_outright_align_score': leftinout_outright_align_score,\
                        'is_nearly_leftinout_ki' : is_nearly_leftinout_ki})

        except UnexpectedException as e:
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _has_rightinout_indicators(self, pddf_row):

        seq = pddf_row["OutOut_Aligned_Sequence"]

        if seq is None:
                
            return(pd.Series({'has_rightinout_indicators' : False, 'Rightside_InOut_Aligned_Sequence' : None, 'has_multi_rightinout_indicators' : False}))
        
        match_target_seq = findall_both(self.right_knockin_regex, self.opp_right_knockin_regex, seq)

        if len(match_target_seq) < 1:# <Dose the seqeuence have both indicators ?>
                
            return(pd.Series({'has_rightinout_indicators' : False, 'Rightside_InOut_Aligned_Sequence' : None, 'has_multi_rightinout_indicators' : False}))

        elif len(match_target_seq) > 1:
            
            return(pd.Series({'has_rightinout_indicators' : True, 'Rightside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_rightinout_indicators' : True}))

        else:

            return(pd.Series({'has_rightinout_indicators' : True, 'Rightside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_rightinout_indicators' : False}))

    def _has_rc_rightinout_indicators(self, pddf_row):
        
        seq = pddf_row["OutOut_Aligned_Sequence"]

        if seq is None:
                
            return(pd.Series({'has_rc_rightinout_indicators' : False, 'RC_Rightside_InOut_Aligned_Sequence' : None, 'has_multi_rc_rightinout_indicators' : False}))
        
        match_target_seq = findall_both(self.right_rc_knockin_regex, self.opp_right_rc_knockin_regex, seq)

        if len(match_target_seq) < 1:# <Dose the seqeuence have both indicators ?>
                
            return(pd.Series({'has_rc_rightinout_indicators' : False, 'RC_Rightside_InOut_Aligned_Sequence' : None, 'has_multi_rc_rightinout_indicators' : False}))

        elif len(match_target_seq) > 1:

            return(pd.Series({'has_rc_rightinout_indicators' : True, 'RC_Rightside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_rc_rightinout_indicators' : True}))

        else:

            return(pd.Series({'has_rc_rightinout_indicators' : True, 'RC_Rightside_InOut_Aligned_Sequence' : match_target_seq[0], 'has_multi_rc_rightinout_indicators' : False}))

    def _scan_rightinout_junction(self, pddf_row):
        """ask whether Rightside_InOut_Aligned_Sequence has rightside in and out indicators"""
        try:

            seq = pddf_row["Rightside_InOut_Aligned_Sequence"]
            has_rightinout_indicators = pddf_row["has_rightinout_indicators"]
            is_perfect_outout_ki = pddf_row["is_perfect_outout_ki"]

            # Focusing on treated sequences only
            if not has_rightinout_indicators or is_perfect_outout_ki:
                # Even if you change the following values to None, they will be forced to be NaN.
                # I hate such an implicit type conversion. I will put 0 in it, instead.
                return pd.Series({'rightinout_left_match_len' : 0,\
                        'rightinout_right_match_len' : 0,\
                        'rightinout_insert_start_point' : 0,\
                        'rightinout_insert_end_point': 0,\
                        'rightinout_outleft_align_score': 0,\
                        'rightinout_outright_align_score': 0,\
                        'is_nearly_rightinout_ki': False})
            else:

                # Length is equal ?
                if len(seq) != len(self.right_knockin_ref_seq):
                    return pd.Series({'rightinout_left_match_len' : 0,\
                        'rightinout_right_match_len' : 0,\
                        'rightinout_insert_start_point' : 0,\
                        'rightinout_insert_end_point': 0,\
                        'rightinout_outleft_align_score': 0,\
                        'rightinout_outright_align_score': 0,\
                        'is_nearly_rightinout_ki': False})

                rightinout_left_match_len,\
                rightinout_right_match_len,\
                rightinout_insert_start_point,\
                rightinout_insert_end_point,\
                rightinout_outleft_align_score,\
                rightinout_outright_align_score,\
                is_nearly_rightinout_ki\
                = self._align_bit_by_bit(seq, self.right_knockin_ref_seq,\
                    self.seq_info_dict["inright_align_threshold"],\
                    self.seq_info_dict["outright_align_threshold"])

                return pd.Series({'rightinout_left_match_len' : rightinout_left_match_len,\
                        'rightinout_right_match_len' : rightinout_right_match_len,\
                        'rightinout_insert_start_point' : rightinout_insert_start_point,\
                        'rightinout_insert_end_point': rightinout_insert_end_point,\
                        'rightinout_outleft_align_score': rightinout_outleft_align_score,\
                        'rightinout_outright_align_score': rightinout_outright_align_score,\
                        'is_nearly_rightinout_ki' : is_nearly_rightinout_ki})

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _add_whole_classification_label_df(self, pddf_row):
        """add whole_classification_label per row."""

        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        if pddf_row["has_outout_indicators"] or pddf_row["has_outout_rc_indicators"]:
            pass
        else:
            return pd.Series({'whole_classification_label' : "no_outout_indicators"})

        if pddf_row["is_nonspecific_outout_indicators"]:
            return pd.Series({'whole_classification_label' : "multiple_outout_indicators"})

        has_leftinout_indicators = pddf_row["has_leftinout_indicators"]
        has_rightinout_indicators = pddf_row["has_rightinout_indicators"]
        has_rc_leftinout_indicators = pddf_row["has_rc_leftinout_indicators"]
        has_rc_rightinout_indicators =pddf_row["has_rc_rightinout_indicators"]

        if has_leftinout_indicators and has_rightinout_indicators:
            return pd.Series({'whole_classification_label' : "both_inout_indicators"})
        elif has_leftinout_indicators:
            return pd.Series({'whole_classification_label' : "leftonly_inout_indicators"})
        elif has_rightinout_indicators:
            return pd.Series({'whole_classification_label' : "rightonly_inout_indicators"})
        elif has_rc_leftinout_indicators and has_rc_rightinout_indicators:
            return pd.Series({'whole_classification_label' : "both_rc_inout_indicators"})
        elif has_rc_leftinout_indicators:
            return pd.Series({'whole_classification_label' : "leftonly_rc_inout_indicators"})
        elif has_rc_rightinout_indicators:
            return pd.Series({'whole_classification_label' : "rightonly_rc_inout_indicators"})
        else:
            return pd.Series({'whole_classification_label' : "no_inout_indicators"})

    def _scan_whole_seq(self):
        """scan whole sequence to add whole_classification_label"""
        """
        Map of Whole Classification
                
        <> : condition
        () : process
        [] : classification label
        {} : sequence set

        ----------------------------------------------------------------------------------------------- <has_outout_indicators, has_outout_rc_indicators, is_nonspecific_outout_indicators>
        |                  |                     |                             |                       |
        outout_indicators  rc outout_indicators  reversible outout_indicators  [no_outout_indicators]  [multiple_outout_indicators]
        |                  |(converting into rc) |(warning)
        |-----------------------------------------(joining)
        |
        ------------- (asking whether is_untreated is True or False)
        |           |
        untreated   |
        |           |
        --------------- (asking whether is_perfect_outout_ki is True or False)
        |             |
        perfect_ki    |
        |	          |
        ----------------------------------------------------------------------------------------- <has_leftinout_indicators, has_rightinout_indicators, has_rc_leftinout_indicators, has_rc_rightinout_indicators>
        |                       |                            |                                  |
        [both_inout_indicators] [leftonly_inout_indicators] [rightonly_inout_indicators]        |
                                                                                                |
        --------------------------------------------------------------------------------------------
        |                       |                            |                                     |
        [both_rc_inout_indicators] [leftonly_rc_inout_indicators] [rightonly_rc_inout_indicators]  [no_inout_indicators]
        
        Whole_classfication_label
        1.no_outout_indicators
        2.multiple_outout_indicators
        3.both_inout_indicators
        4.leftonly_inout_indicators
        5.rightonly_inout_indicators
        6.both_rc_inout_indicators
        7.leftonly_rc_inout_indicators
        8.rightonly_rc_inout_indicators
        9.no_inout_indicators
        """

        # <How many sequences have out-out indicators ?>

        # Out-Out Indicator

        if self.seq_info_dict["optimization_mode"] in "outoutindicators":
            self.temp_pddf = self.out_pddf.copy()
            return 0
        else:
            tqdm.pandas(desc = "Searching Out-Out Indicators")
            self.out_pddf = self.out_pddf.join(self.out_pddf.Aligned_Sequence.progress_apply(self._has_out_indicators)).copy()

        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        self._register_count("total", self.out_pddf)
        self._register_count("has_outout_indicators", self.out_pddf.query('has_outout_indicators == True and has_outout_rc_indicators == False'))
        self._register_count("has_outout_rc_indicators", self.out_pddf.query('has_outout_indicators == False and has_outout_rc_indicators == True'))
        self._register_count("has_outout_reversible_indicators", self.out_pddf.query('has_outout_indicators == True and has_outout_rc_indicators == True'))
        self._register_count("has_outout_no_indicators", self.out_pddf.query('has_outout_indicators == False and has_outout_rc_indicators == False'))
        self._register_count("is_nonspecific_outout_indicators", self.out_pddf.query('is_nonspecific_outout_indicators == True'))
        self.cnt_dict["has_outout_any_indicators"] = self.cnt_dict["has_outout_indicators"] \
            + self.cnt_dict["has_outout_rc_indicators"] \
            + self.cnt_dict["has_outout_reversible_indicators"]
        self.read_cnt_dict["has_outout_any_indicators"] = self.read_cnt_dict["has_outout_indicators"] \
            + self.read_cnt_dict["has_outout_rc_indicators"] \
            + self.read_cnt_dict["has_outout_reversible_indicators"]
        


        if self.seq_info_dict["optimization_mode"] in "OFF":
            self._show_count("has_outout_indicators")
            self._show_count("has_outout_rc_indicators")
            self._show_count("has_outout_reversible_indicators")
            self._show_count("has_outout_any_indicators")
            self._show_count("has_outout_no_indicators")
        else:
            pass

        # <How many sequences with indicators are untreated sequence ?>

        if self.seq_info_dict["optimization_mode"] in "outout_align_thresholds":
            self.temp_pddf = self.out_pddf.copy()
            return 0
        else:
            tqdm.pandas(desc = "Searching Untreated Sequences")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_untreated_seq, axis = 1)).copy()

        self._register_count("is_untreated", self.out_pddf.query('is_untreated == True'))
        self.cnt_dict["is_treated"] = self.cnt_dict["has_outout_any_indicators"] - self.cnt_dict["is_untreated"]
        self.read_cnt_dict["is_treated"] = self.read_cnt_dict["has_outout_any_indicators"] - self.read_cnt_dict["is_untreated"]

        if self.seq_info_dict["optimization_mode"] in "OFF":
            self._show_count("is_untreated")
            self._show_count("is_treated")
        else:
            pass

        # <What type of sequence are they ?>

        # OUT-OUT BbB Alignment

        tqdm.pandas(desc = "Scanning Whole Junctions")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._scan_outout_junction, axis = 1)).copy()
        
        self._register_count("is_perfect_outout_ki", self.out_pddf.query('is_perfect_outout_ki == True'))

        if self.seq_info_dict["optimization_mode"] in "OFF":
            self._show_count("is_perfect_outout_ki")
        else:
            pass

        # Get In-Out Indicator
        if self.seq_info_dict["optimization_mode"] in "inoutindicators":
            self.temp_pddf = self.out_pddf.copy()
            return 0
        else:
            tqdm.pandas(desc = "Searching Leftside In-Out Indicators")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._has_leftinout_indicators, axis = 1)).copy()
            tqdm.pandas(desc = "Searching Rightside In-Out Indicators")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._has_rightinout_indicators, axis = 1)).copy()
        
        # Get reverse-complement In-Out Indicator
        if self.seq_info_dict["optimization_mode"] in "OFF":
            tqdm.pandas(desc = "Searching Reverse-Complement Leftside In-Out Indicators")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._has_rc_leftinout_indicators, axis = 1)).copy()
            tqdm.pandas(desc = "Searching Reverse-Complement Rightside In-Out Indicators")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._has_rc_rightinout_indicators, axis = 1)).copy()

        self._register_count("has_leftinout_indicators", self.out_pddf.query('has_leftinout_indicators == True'))
        self._register_count("has_multi_leftinout_indicators", self.out_pddf.query('has_multi_leftinout_indicators == True'))
        self._register_count("has_rightinout_indicators", self.out_pddf.query('has_rightinout_indicators == True'))
        self._register_count("has_multi_rightinout_indicators", self.out_pddf.query('has_multi_rightinout_indicators == True'))
        if self.seq_info_dict["optimization_mode"] in "OFF":
            self._show_count("has_leftinout_indicators")
            self._show_count("has_rightinout_indicators")
        else:
            pass

        # Right & Left In-Out BbB Alignment

        if self.seq_info_dict["optimization_mode"] in "inout_align_thresholds":
            self.temp_pddf = self.out_pddf.copy()
            return 0
        else:
            tqdm.pandas(desc = "Scanning Left Junctions")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._scan_leftinout_junction, axis = 1)).copy()
            tqdm.pandas(desc = "Scanning Right Junctions")
            self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._scan_rightinout_junction, axis = 1)).copy()
        
        self._register_count("is_perfect_ki_has_leftinout_indicators", \
            self.out_pddf.query('is_perfect_outout_ki == True and has_leftinout_indicators == True'))
        self._register_count("is_nearly_leftinout_ki", \
            self.out_pddf.query('is_nearly_leftinout_ki == True'))
        self._register_count("is_perfect_ki_has_rightinout_indicators", \
            self.out_pddf.query('is_perfect_outout_ki == True and has_rightinout_indicators == True'))
        self._register_count("is_nearly_rightinout_ki", \
            self.out_pddf.query('is_nearly_rightinout_ki == True'))
        

        # Add whole_classification_label
        tqdm.pandas(desc = "Adding whole_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_whole_classification_label_df, axis = 1)).copy()

        self._register_count("no_outout_indicators", self.out_pddf.query('whole_classification_label == "no_outout_indicators"'))
        self._register_count("multiple_outout_indicators", self.out_pddf.query('whole_classification_label == "multiple_outout_indicators"'))
        self._register_count("both_inout_indicators", self.out_pddf.query('whole_classification_label == "both_inout_indicators"'))
        self._register_count("leftonly_inout_indicators", self.out_pddf.query('whole_classification_label == "leftonly_inout_indicators"'))
        self._register_count("rightonly_inout_indicators", self.out_pddf.query('whole_classification_label == "rightonly_inout_indicators"'))
        self._register_count("both_rc_inout_indicators", self.out_pddf.query('whole_classification_label == "both_rc_inout_indicators"'))
        self._register_count("leftonly_rc_inout_indicators", self.out_pddf.query('whole_classification_label == "leftonly_rc_inout_indicators"'))
        self._register_count("rightonly_rc_inout_indicators", self.out_pddf.query('whole_classification_label == "rightonly_rc_inout_indicators"'))
        self._register_count("no_inout_indicators", self.out_pddf.query('whole_classification_label == "no_inout_indicators"'))

        self._show_count("no_outout_indicators")
        self._show_count("multiple_outout_indicators")
        self._show_count("both_inout_indicators")
        self._show_count("leftonly_inout_indicators")
        self._show_count("rightonly_inout_indicators")
        self._show_count("both_rc_inout_indicators")
        self._show_count("leftonly_rc_inout_indicators")
        self._show_count("rightonly_rc_inout_indicators")
        self._show_count("no_inout_indicators")

    #-------------------------------------------------------------------------------------------------------------#
    #--- Left Classification -------------------------------------------------------------------------------------#

    def _can_scan_leftjunction(self, pddf_row):
        """ask whether Leftside_InOut_Aligned_Sequence can be scanned"""
        try:
            
            #import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            whole_classification_label = pddf_row["whole_classification_label"]

            if is_in(whole_classification_label, ("both_inout_indicators", "leftonly_inout_indicators")):
                return pd.Series({'can_scan_leftjunction' : True})

            else:
                return pd.Series({'can_scan_leftjunction' : False})

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _is_left_nhej_ki(self, pddf_row):
        """ask whether Left InOut_Aligned_Sequence is nhej-assisted knock-in"""
        seq = pddf_row["Leftside_InOut_Aligned_Sequence"]

        if not any((not pddf_row["can_scan_leftjunction"],\
            pddf_row["is_perfect_outout_ki"],\
            pddf_row["is_nearly_leftinout_ki"])\
            ):

            nhej_limit_alignment_score = pairwise2.align.globalms(self.left_nhej_knockin_ref_seq, self.left_nhej_del_knockin_ref_seq, 2, -1, -10, -0.5)[0][2]
            seq_alignment_score = pairwise2.align.globalms(self.left_nhej_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0][2]
            if seq_alignment_score >= nhej_limit_alignment_score:
                return pd.Series({'is_left_nhej_ki' : True})
            else:
                return pd.Series({'is_left_nhej_ki' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.left_nhej_knockin_ref_seq, self.left_out_indicator_seq + self.left_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        else:
            return pd.Series({'is_left_nhej_ki' : False})

    def _is_left_imperfect_ki_complex(self, pddf_row):
        """ask whether Left InOut_Aligned_Sequence is imperfect knock-in complex"""
        seq = pddf_row["Leftside_InOut_Aligned_Sequence"]

        if not any((not pddf_row["can_scan_leftjunction"],\
            pddf_row["is_perfect_outout_ki"],\
            pddf_row["is_nearly_leftinout_ki"],\
            pddf_row["is_left_nhej_ki"])\
            ):

            if len(seq) >= len(self.left_knockin_ref_seq):
                # print pairwise2.align.globalms(self.left_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0]
                return pd.Series({'is_left_imperfect_ki_complex' : True})
            else:
                return pd.Series({'is_left_imperfect_ki_complex' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.left_nhej_knockin_ref_seq, self.left_out_indicator_seq + self.left_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        else:
            return pd.Series({'is_left_imperfect_ki_complex' : False})

    def _add_left_classification_label_df(self, pddf_row):
        """add left_classification_label per row."""

        if not pddf_row["can_scan_leftjunction"]:
            return pd.Series({'left_classification_label' : "no_leftinout_indicators"})
        elif pddf_row["is_perfect_outout_ki"]:
            return pd.Series({'left_classification_label' : "left_perfect_ki"})
        elif pddf_row["is_nearly_leftinout_ki"]:
            return pd.Series({'left_classification_label' : "left_nearly_ki"})
        elif pddf_row["is_left_nhej_ki"]:
            return pd.Series({'left_classification_label' : "left_nhej_ki"})
        elif pddf_row["is_left_imperfect_ki_complex"]:
            return pd.Series({'left_classification_label' : "left_imperfect_ki_complex"})
        else:
            return pd.Series({'left_classification_label' : "other_left_mutations"})

    def _scan_left_seq(self):
        """scan left junction site to add left_classification_label"""
        """
        Map of Left InOut Classification
                
        <> : condition
        () : process
        [] : classification label
        {} : sequence set

        {both_inout_indicators} or {leftonly_inout_indicators}
        ----------------------- <can_scan_leftjunction>
        |                    |
        |                    [no_leftinout_indicators]
        |
        -------------------- <is_perfect_outout_ki>
        |                  |
        [left_perfect_ki]  |
                           |
        -------------------- <is_nearly_leftinout_ki>
        |                  |
        [left_nearly_ki]   |
                           |
        -------------------- <is_left_nhej_ki>
        |                  |
        [left_nhej_ki]     |
                           |
        ------------------------------------ <is_left_imperfect_ki_complex>
        |                                  |
        [left_imperfect_ki_complex]  [other_left_mutations]

        Left_classfication_label
        1.no_leftinout_indicators
        2.left_perfect_ki
        3.left_nearly_ki
        4.left_nhej_ki
        5.left_imperfect_ki_complex
        6.other_left_mutations
        """

        tqdm.pandas(desc = "Searching Left In-Out Indicators")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._can_scan_leftjunction, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Left NHEJ-assisted knock-in")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_left_nhej_ki, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Left imperfect knock-in complex")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_left_imperfect_ki_complex, axis = 1)).copy()

        # Add left_classification_label
        tqdm.pandas(desc = "Adding left_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_left_classification_label_df, axis = 1)).copy()

        self._register_count("no_leftinout_indicators", self.out_pddf.query('left_classification_label == "no_leftinout_indicators"'))
        self._register_count("left_perfect_ki", self.out_pddf.query('left_classification_label == "left_perfect_ki"'))
        self._register_count("left_nearly_ki", self.out_pddf.query('left_classification_label == "left_nearly_ki"'))
        self._register_count("left_nhej_ki", self.out_pddf.query('left_classification_label == "left_nhej_ki"'))
        self._register_count("left_imperfect_ki_complex", self.out_pddf.query('left_classification_label == "left_imperfect_ki_complex"'))
        self._register_count("other_left_mutations", self.out_pddf.query('left_classification_label == "other_left_mutations"'))
        self._show_count("no_leftinout_indicators")
        self._show_count("left_perfect_ki")
        self._show_count("left_nearly_ki")
        self._show_count("left_nhej_ki")
        self._show_count("left_imperfect_ki_complex")
        self._show_count("other_left_mutations")

    #-------------------------------------------------------------------------------------------------------------#
    #--- Right Classification -------------------------------------------------------------------------------------#

    def _can_scan_rightjunction(self, pddf_row):
        """ask whether rightside_InOut_Aligned_Sequence can be scanned"""
        try:
            
            whole_classification_label = pddf_row["whole_classification_label"]

            if is_in(whole_classification_label, ("both_inout_indicators", "rightonly_inout_indicators")):

                return pd.Series({'can_scan_rightjunction' : True})

            else:

                return pd.Series({'can_scan_rightjunction' : False})

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _is_right_nhej_ki(self, pddf_row):
        """ask whether right InOut_Aligned_Sequence is nhej-assisted knock-in"""
        seq = pddf_row["Rightside_InOut_Aligned_Sequence"]

        if not any((not pddf_row["can_scan_rightjunction"],\
            pddf_row["is_perfect_outout_ki"],\
            pddf_row["is_nearly_rightinout_ki"])\
            ):

            nhej_limit_alignment_score = pairwise2.align.globalms(self.right_nhej_knockin_ref_seq, self.right_nhej_del_knockin_ref_seq, 2, -1, -10, -0.5)[0][2]
            seq_alignment_score = pairwise2.align.globalms(self.right_nhej_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0][2]
            if seq_alignment_score >= nhej_limit_alignment_score:

                return pd.Series({'is_right_nhej_ki' : True})

            else:

                return pd.Series({'is_right_nhej_ki' : False})

            # checking score of blank seq
            # print pairwise2.align.globalms(self.right_nhej_knockin_ref_seq, self.right_out_indicator_seq + self.right_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        else:
            return pd.Series({'is_right_nhej_ki' : False})

    def _is_right_imperfect_ki_complex(self, pddf_row):
        """ask whether Right InOut_Aligned_Sequence is imperfect knock-in complex"""
        seq = pddf_row["Rightside_InOut_Aligned_Sequence"]

        if not any((not pddf_row["can_scan_rightjunction"],\
            pddf_row["is_perfect_outout_ki"],\
            pddf_row["is_nearly_rightinout_ki"],\
            pddf_row["is_right_nhej_ki"])\
            ):

            if len(seq) >= len(self.right_knockin_ref_seq):
                # print pairwise2.align.globalms(self.right_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0]
                return pd.Series({'is_right_imperfect_ki_complex' : True})
            else:
                return pd.Series({'is_right_imperfect_ki_complex' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.right_nhej_knockin_ref_seq, self.right_out_indicator_seq + self.right_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        else:
            return pd.Series({'is_right_imperfect_ki_complex' : False})

    def _add_right_classification_label_df(self, pddf_row):
        """add right_classification_label per row."""

        if not pddf_row["can_scan_rightjunction"]:

            return pd.Series({'right_classification_label' : "no_rightinout_indicators"})

        elif pddf_row["is_perfect_outout_ki"]:

            return pd.Series({'right_classification_label' : "right_perfect_ki"})

        elif pddf_row["is_nearly_rightinout_ki"]:

            return pd.Series({'right_classification_label' : "right_nearly_ki"})

        elif pddf_row["is_right_nhej_ki"]:

            return pd.Series({'right_classification_label' : "right_nhej_ki"})

        elif pddf_row["is_right_imperfect_ki_complex"]:

            return pd.Series({'right_classification_label' : "right_imperfect_ki_complex"})

        else:

            return pd.Series({'right_classification_label' : "other_right_mutations"})

    def _scan_right_seq(self):
        """scan right junction site to add right_classification_label"""
        """
        Map of Right InOut Classification
                
        <> : condition
        () : process
        [] : classification label
        {} : sequence set

        {both_inout_indicators} or {rightonly_inout_indicators}
        ----------------------- <can_scan_rightjunction>
        |                    |
        |                    [no_rightinout_indicators]
        |
        -------------------- <is_perfect_outout_ki>
        |                  |
        [right_perfect_ki]  |
                           |
        -------------------- <is_nearly_rightinout_ki>
        |                  |
        [right_nearly_ki]   |
                           |
        -------------------- <is_right_nhej_ki>
        |                  |
        [right_nhej_ki]     |
                           |
        ------------------------------------ <is_right_imperfect_ki_complex>
        |                                  |
        [right_imperfect_ki_complex]  [other_right_mutations]

        Right_classfication_label
        1.no_rightinout_indicators
        2.right_perfect_ki
        3.right_nearly_ki
        4.right_nhej_ki
        5.right_imperfect_ki_complex
        6.other_right_mutations
        """

        tqdm.pandas(desc = "Searching Right In-Out Indicators")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._can_scan_rightjunction, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Right NHEJ-assisted knock-in")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_right_nhej_ki, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Right imperfect knock-in complex")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_right_imperfect_ki_complex, axis = 1)).copy()

        # Add right_classification_label
        tqdm.pandas(desc = "Adding right_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_right_classification_label_df, axis = 1)).copy()

        self._register_count("no_rightinout_indicators", self.out_pddf.query('right_classification_label == "no_rightinout_indicators"'))
        self._register_count("right_perfect_ki", self.out_pddf.query('right_classification_label == "right_perfect_ki"'))
        self._register_count("right_nearly_ki", self.out_pddf.query('right_classification_label == "right_nearly_ki"'))
        self._register_count("right_nhej_ki", self.out_pddf.query('right_classification_label == "right_nhej_ki"'))
        self._register_count("right_imperfect_ki_complex", self.out_pddf.query('right_classification_label == "right_imperfect_ki_complex"'))
        self._register_count("other_right_mutations", self.out_pddf.query('right_classification_label == "other_right_mutations"'))
        
        self._show_count("no_rightinout_indicators")
        self._show_count("right_perfect_ki")
        self._show_count("right_nearly_ki")
        self._show_count("right_nhej_ki")
        self._show_count("right_imperfect_ki_complex")
        self._show_count("other_right_mutations")
    
    #-------------------------------------------------------------------------------------------------------------#
    #--- Reverse-Complement Left Classification -------------------------------------------------------------------------------------#

    def _can_scan_left_rc_junction(self, pddf_row):
        """ask whether RC leftside_InOut_Aligned_Sequence can be scanned using Left Reverse-Complement In-Out indicator"""
        try:
            
            whole_classification_label = pddf_row["whole_classification_label"]

            if is_in(whole_classification_label, ("both_rc_inout_indicators", "leftonly_rc_inout_indicators")):

                return pd.Series({'can_scan_left_rc_junction' : True})

            else:
                return pd.Series({'can_scan_left_rc_junction' : False})

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _is_left_rc_nhej_ki(self, pddf_row):
        """ask whether RC Left InOut_Aligned_Sequence is nhej-assisted knock-in"""
        seq = pddf_row["RC_Leftside_InOut_Aligned_Sequence"]

        if pddf_row["can_scan_left_rc_junction"]:

            nhej_limit_alignment_score = pairwise2.align.globalms(self.opposite_left_nhej_knockin_ref_seq, self.opposite_left_nhej_del_knockin_ref_seq, 2, -1, -10, -0.5)[0][2]
            seq_alignment_score = pairwise2.align.globalms(self.opposite_left_nhej_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0][2]

            if seq_alignment_score >= nhej_limit_alignment_score:

                # print pairwise2.align.globalms(self.opposite_left_nhej_knockin_ref_seq, self.opposite_left_nhej_del_knockin_ref_seq, 2, -1, -10, -0.5)[0]
                return pd.Series({'is_left_rc_nhej_ki' : True})

            else:

                return pd.Series({'is_left_rc_nhej_ki' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.opposite_left_nhej_knockin_ref_seq, self.left_out_indicator_seq + self.left_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        else:

            return pd.Series({'is_left_rc_nhej_ki' : False})

    def _is_left_rc_imperfect_ki_complex(self, pddf_row):
        """ask whether RC Left InOut_Aligned_Sequence is imperfect knock-in complex"""
        seq = pddf_row["RC_Leftside_InOut_Aligned_Sequence"]

        if not any((not pddf_row["can_scan_left_rc_junction"],\
            pddf_row["is_left_rc_nhej_ki"])\
            ):

            if len(seq) >= len(self.opposite_left_knockin_ref_seq):

                # print pairwise2.align.globalms(self.opposite_left_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0]
                return pd.Series({'is_left_rc_imperfect_ki_complex' : True})

            else:

                return pd.Series({'is_left_rc_imperfect_ki_complex' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.opposite_left_nhej_knockin_ref_seq, self.left_out_indicator_seq + self.left_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        else:
            return pd.Series({'is_left_rc_imperfect_ki_complex' : False})

    def _add_left_rc_classification_label_df(self, pddf_row):
        """add left_rc_classification_label per row."""

        if not pddf_row["can_scan_left_rc_junction"]:

            return pd.Series({'left_rc_classification_label' : "no_leftinout_rc_indicators"})

        elif pddf_row["is_left_rc_nhej_ki"]:

            return pd.Series({'left_rc_classification_label' : "left_rc_nhej_ki"})

        elif pddf_row["is_left_rc_imperfect_ki_complex"]:

            return pd.Series({'left_rc_classification_label' : "left_rc_imperfect_ki_complex"})

        else:

            return pd.Series({'left_rc_classification_label' : "other_left_rc_mutations"})


    def _scan_left_rc_seq(self):
        """
        Map of Left Reverse-Complement InOut Classification
                
        <> : condition
        () : process
        [] : classification label
        {} : sequence set

        {both_rc_inout_indicators} or {leftonly_rc_inout_indicators}
        ----------------------- <can_scan_left_rc_junction>
        |                    |
        |                    [no_leftinout_rc_indicators]
        |
        -------------------- <is_left_rc_nhej_ki>
        |                  |
        [left_rc_nhej_ki]  |
                            |
        ------------------------------------ <is_left_rc_imperfect_ki_complex>
        |                                  |
        [left_rc_imperfect_ki_complex]  [other_left_rc_mutations]

        Left_rc_classfication_label
        1.no_leftinout_rc_indicators
        2.left_rc_nhej_ki
        3.left_rc_imperfect_ki_complex
        4.other_left_rc_mutations
        """

        tqdm.pandas(desc = "Searching Left Reverse-Complement In-Out Indicators")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._can_scan_left_rc_junction, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Left Reverse-Complement NHEJ-assisted knock-in")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_left_rc_nhej_ki, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Left Reverse-Complement imperfect knock-in complex")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_left_rc_imperfect_ki_complex, axis = 1)).copy()

        # Add left_rc_classification_label
        tqdm.pandas(desc = "Adding left_rc_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_left_rc_classification_label_df, axis = 1)).copy()

        self._register_count("no_leftinout_rc_indicators", self.out_pddf.query('left_rc_classification_label == "no_leftinout_rc_indicators"'))
        self._register_count("left_rc_nhej_ki", self.out_pddf.query('left_rc_classification_label == "left_rc_nhej_ki"'))
        self._register_count("left_rc_imperfect_ki_complex", self.out_pddf.query('left_rc_classification_label == "left_rc_imperfect_ki_complex"'))
        self._register_count("other_left_rc_mutations", self.out_pddf.query('left_rc_classification_label == "other_left_rc_mutations"'))
        
        self._show_count("no_leftinout_rc_indicators")
        self._show_count("left_rc_nhej_ki")
        self._show_count("left_rc_imperfect_ki_complex")
        self._show_count("other_left_rc_mutations")

    #-------------------------------------------------------------------------------------------------------------#
    #--- Reverse-Complement Right Classification -------------------------------------------------------------------------------------#

    def _can_scan_right_rc_junction(self, pddf_row):
        """ask whether RC rightside_InOut_Aligned_Sequence can be scanned using Right Reverse-Complement In-Out indicator"""
        try:
            
            whole_classification_label = pddf_row["whole_classification_label"]

            if is_in(whole_classification_label, ("both_rc_inout_indicators", "rightonly_rc_inout_indicators")):

                return pd.Series({'can_scan_right_rc_junction' : True})

            else:
                return pd.Series({'can_scan_right_rc_junction' : False})

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _is_right_rc_nhej_ki(self, pddf_row):
        """ask whether RC Right InOut_Aligned_Sequence is nhej-assisted knock-in"""
        seq = pddf_row["RC_Rightside_InOut_Aligned_Sequence"]
        
        if pddf_row["can_scan_right_rc_junction"]:

            nhej_limit_alignment_score = pairwise2.align.globalms(self.opposite_right_nhej_knockin_ref_seq, self.opposite_right_nhej_del_knockin_ref_seq, 2, -1, -10, -0.5)[0][2]
            seq_alignment_score = pairwise2.align.globalms(self.opposite_right_nhej_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0][2]

            if seq_alignment_score >= nhej_limit_alignment_score:

                # print pairwise2.align.globalms(self.opposite_right_nhej_knockin_ref_seq, self.opposite_right_nhej_del_knockin_ref_seq, 2, -1, -10, -0.5)[0]
                return pd.Series({'is_right_rc_nhej_ki' : True})

            else:

                return pd.Series({'is_right_rc_nhej_ki' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.opposite_right_nhej_knockin_ref_seq, self.right_out_indicator_seq + self.right_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        else:

            return pd.Series({'is_right_rc_nhej_ki' : False})

    def _is_right_rc_imperfect_ki_complex(self, pddf_row):
        """ask whether RC Right InOut_Aligned_Sequence is imperfect knock-in complex"""
        seq = pddf_row["RC_Rightside_InOut_Aligned_Sequence"]

        if not any((not pddf_row["can_scan_right_rc_junction"],\
            pddf_row["is_right_rc_nhej_ki"])\
            ):

            if len(seq) >= len(self.opposite_right_knockin_ref_seq):

                # print pairwise2.align.globalms(self.opposite_right_knockin_ref_seq, seq, 2, -1, -10, -0.5)[0]
                return pd.Series({'is_right_rc_imperfect_ki_complex' : True})

            else:

                return pd.Series({'is_right_rc_imperfect_ki_complex' : False})
            # checking score of blank seq
            # print pairwise2.align.globalms(self.opposite_right_nhej_knockin_ref_seq, self.right_out_indicator_seq + self.right_in_indicator_seq, 2, -1, -10, -0.5)[0]
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        else:
            return pd.Series({'is_right_rc_imperfect_ki_complex' : False})

    def _add_right_rc_classification_label_df(self, pddf_row):
        """add right_rc_classification_label per row."""

        if not pddf_row["can_scan_right_rc_junction"]:

            return pd.Series({'right_rc_classification_label' : "no_rightinout_rc_indicators"})

        elif pddf_row["is_right_rc_nhej_ki"]:

            return pd.Series({'right_rc_classification_label' : "right_rc_nhej_ki"})

        elif pddf_row["is_right_rc_imperfect_ki_complex"]:

            return pd.Series({'right_rc_classification_label' : "right_rc_imperfect_ki_complex"})

        else:

            return pd.Series({'right_rc_classification_label' : "other_right_rc_mutations"})


    def _scan_right_rc_seq(self):
        """
        Map of Right Reverse-Complement InOut Classification
                
        <> : condition
        () : process
        [] : classification label
        {} : sequence set

        {both_rc_inout_indicators} or {rightonly_rc_inout_indicators}
        ----------------------- <can_scan_right_rc_junction>
        |                    |
        |                    [no_rightinout_rc_indicators]
        |
        -------------------- <is_right_rc_nhej_ki>
        |                  |
        [right_rc_nhej_ki]  |
                            |
        ------------------------------------ <is_right_rc_imperfect_ki_complex>
        |                                  |
        [right_rc_imperfect_ki_complex]  [other_right_rc_mutations]

        Right_rc_classfication_label
        1.no_rightinout_rc_indicators
        2.right_rc_nhej_ki
        3.right_rc_imperfect_ki_complex
        4.other_right_rc_mutations
        """

        tqdm.pandas(desc = "Searching Right Reverse-Complement In-Out Indicators")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._can_scan_right_rc_junction, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Right Reverse-Complement NHEJ-assisted knock-in")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_right_rc_nhej_ki, axis = 1)).copy()

        tqdm.pandas(desc = "Searching Right Reverse-Complement imperfect knock-in complex")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_right_rc_imperfect_ki_complex, axis = 1)).copy()

        # Add right_rc_classification_label
        tqdm.pandas(desc = "Adding right_rc_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_right_rc_classification_label_df, axis = 1)).copy()

        self._register_count("no_rightinout_rc_indicators", self.out_pddf.query('right_rc_classification_label == "no_rightinout_rc_indicators"'))
        self._register_count("right_rc_nhej_ki", self.out_pddf.query('right_rc_classification_label == "right_rc_nhej_ki"'))
        self._register_count("right_rc_imperfect_ki_complex", self.out_pddf.query('right_rc_classification_label == "right_rc_imperfect_ki_complex"'))
        self._register_count("other_right_rc_mutations", self.out_pddf.query('right_rc_classification_label == "other_right_rc_mutations"'))
        
        self._show_count("no_rightinout_rc_indicators")
        self._show_count("right_rc_nhej_ki")
        self._show_count("right_rc_imperfect_ki_complex")
        self._show_count("other_right_rc_mutations")
    
    #-------------------------------------------------------------------------------------------------------------#
    #--- Mutation Classification -------------------------------------------------------------------------------------#

    def _is_targetsite_edited(self, pddf_row):
        """ask whether OutOut_Aligned_Sequence is imperfect knock-in complex"""
        seq = pddf_row["OutOut_Aligned_Sequence"]
        left_match_len = pddf_row["outout_left_match_len"]
        right_match_len = pddf_row["outout_right_match_len"]

        if not any((\
            (pddf_row["whole_classification_label"] not in "no_inout_indicators"),\
            (pddf_row["whole_classification_label"] in "multiple_outout_indicators"),\
            pddf_row["is_untreated"])\
            ):

            left_cutout_maxdel_len = len(self.left_cutout_maxdel_seq)
            right_cutout_maxdel_len = len(self.right_cutout_maxdel_seq)
            if left_cutout_maxdel_len < left_match_len < left_cutout_maxdel_len + 2 * self.oneside_deletion_length\
                or right_cutout_maxdel_len < right_match_len < right_cutout_maxdel_len + 2 * self.oneside_deletion_length:

                return pd.Series({'is_targetsite_edited' : True})

            else:

                return pd.Series({'is_targetsite_edited' : False})

            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        
        else:

            return pd.Series({'is_targetsite_edited' : False})

    def _is_overall_edited_insertion_or_substitution(self, pddf_row):
        """ask whether OutOut_Aligned_Sequence is imperfect knock-in complex"""
        seq = pddf_row["OutOut_Aligned_Sequence"]

        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        if not any((\
            (pddf_row["whole_classification_label"] not in "no_inout_indicators"),\
            (pddf_row["whole_classification_label"] in "multiple_outout_indicators"),\
            pddf_row["is_untreated"],\
            not pddf_row["is_targetsite_edited"])\
            ):

            if len(seq) == len(self.indicator_untreated_ref_seq):

                return pd.Series({'is_overall_edited_substitution' : True, 'is_overall_edited_insertion' : False, 'is_overall_edited_deletion' : False})

            elif len(seq) > len(self.indicator_untreated_ref_seq):

                return pd.Series({'is_overall_edited_substitution' : False, 'is_overall_edited_insertion' : True, 'is_overall_edited_deletion' : False})

            else:

                return pd.Series({'is_overall_edited_substitution' : False, 'is_overall_edited_insertion' : False, 'is_overall_edited_deletion' : True})

            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        
        else:

            return pd.Series({'is_overall_edited_substitution' : False, 'is_overall_edited_insertion' : False, 'is_overall_edited_deletion' : False})
    
    def _add_mutation_classification_label_df(self, pddf_row):
        """add mutation_classification_label per row."""

        try:

            if is_in(pddf_row["whole_classification_label"], ("no_outout_indicators", "multiple_outout_indicators")):

                return pd.Series({'mutation_classification_label' : "unidentified"})

            if pddf_row["is_untreated"]:
    
                return pd.Series({'mutation_classification_label' : "untreated"})

            if pddf_row["whole_classification_label"] not in "no_inout_indicators":

                return pd.Series({'mutation_classification_label' : "donor_insertion"})

            elif not pddf_row["is_targetsite_edited"]:

                return pd.Series({'mutation_classification_label' : "unexpected_mutations"})

            elif pddf_row["is_overall_edited_substitution"]:
                
                return pd.Series({'mutation_classification_label' : "overall_edited_substitution"})

            elif pddf_row["is_overall_edited_insertion"]:
                
                return pd.Series({'mutation_classification_label' : "overall_edited_insertion"})

            elif pddf_row["is_overall_edited_deletion"]:

                return pd.Series({'mutation_classification_label' : "overall_edited_deletion"})

            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                raise UnexpectedException('Unexpected categorizing occurs.')

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)
            

    def _scan_mutation_seq(self):
        """
        Map of Mutation Classification

        <> : condition
        () : process
        [] : classification label
        {} : sequence set

        ----------------------- <no_outout_indicators, multiple_outout_indicators>
        |                     |                     
        [unidentified]        |            
                              |
        ----------------------- <is_untreated>
        |                     |                     
        [untreated]           |        
                              |
        ---------------------- {no_inout_indicators}
        |                     | 
        |                     [donor_insertion]
        |
        ----------------------- 
        |                     |
        [untreated]           |
                              |
        ----------------------- <is_targetsite_edited, multiple_outout_indicators>
        |                     |
        |                     [unexpected_mutations]
        |                     
        ------------------------------ <is_overall_edited_substitution>
        |                            |
        [overall_edited_substitution]|
                                     |
        ------------------------------ <is_overall_edited_insertion>
        |                            |
        [overall_edited_insertion]   [overall_edited_deletion]

        Mutation_Classification_label
        1.unidentified
        2.donor_insertion
        3.untreated
        4.unexpected_mutations
        5.overall_edited_substitution
        6.overall_edited_insertion
        7.overall_edited_deletion
        """

        tqdm.pandas(desc = "Searching targeted sequence")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_targetsite_edited, axis = 1)).copy()

        tqdm.pandas(desc = "Categorizing targeted mutation sequence")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._is_overall_edited_insertion_or_substitution, axis = 1)).copy()

        # Add right_classification_label
        tqdm.pandas(desc = "Adding mutation_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_mutation_classification_label_df, axis = 1)).copy()

        self._register_count("donor_insertion", self.out_pddf.query('mutation_classification_label == "donor_insertion"'))
        self._register_count("untreated", self.out_pddf.query('mutation_classification_label == "untreated"'))
        self._register_count("unexpected_mutations", self.out_pddf.query('mutation_classification_label == "unexpected_mutations"'))
        self._register_count("overall_edited_substitution", self.out_pddf.query('mutation_classification_label == "overall_edited_substitution"'))
        self._register_count("overall_edited_insertion", self.out_pddf.query('mutation_classification_label == "overall_edited_insertion"'))
        self._register_count("overall_edited_deletion", self.out_pddf.query('mutation_classification_label == "overall_edited_deletion"'))

        self._show_count("donor_insertion")
        self._show_count("untreated")
        self._show_count("unexpected_mutations")
        self._show_count("overall_edited_substitution")
        self._show_count("overall_edited_insertion")
        self._show_count("overall_edited_deletion")

    #-------------------------------------------------------------------------------------------------------------#
    #--- Global Classification -------------------------------------------------------------------------------------#

    def _add_global_left_classification_label_df(self, pddf_row):
        """
        Global Left In-Out Classification
        LEFT_PRECISE_KNOCK_IN {left_perfect_ki, left_nearly_ki}
        LEFT_IMPRECISE_KNOCK_IN {left_nhej_ki, left_rc_nhej_ki}
        LEFT_IMPRECISE_KNOCK_IN_COMPLEX {left_imperfect_ki_complex, left_rc_imperfect_ki_complex}
        LEFT_OTHER_MUATIONS {other_left_mutations, other_left_rc_mutations}
        LEFT_NO_DETECTED {no_leftinout_indicators and no_leftinout_rc_indicators}
        """
        
        try:

            left_classification_label = pddf_row["left_classification_label"]
            left_rc_classification_label = pddf_row["left_rc_classification_label"]

            if is_in(left_classification_label, ("left_perfect_ki", "left_nearly_ki")):

                return pd.Series({'global_left_classification_label' : "LEFT_PRECISE_KNOCK_IN"})

            elif (left_classification_label in "left_nhej_ki") or (left_rc_classification_label in "left_rc_nhej_ki"):

                return pd.Series({'global_left_classification_label' : "LEFT_IMPRECISE_KNOCK_IN"})

            elif (left_classification_label in "left_imperfect_ki_complex") or (left_rc_classification_label in "left_rc_imperfect_ki_complex"):

                return pd.Series({'global_left_classification_label' : "LEFT_IMPRECISE_KNOCK_IN_COMPLEX"})

            elif (left_classification_label in "other_left_mutations") or (left_rc_classification_label in "other_left_rc_mutations"):
                
                return pd.Series({'global_left_classification_label' : "LEFT_OTHER_MUATIONS"})

            elif (left_classification_label in "no_leftinout_indicators") and (left_rc_classification_label in "no_leftinout_rc_indicators"):
                
                return pd.Series({'global_left_classification_label' : "LEFT_NO_DETECTED"})

            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                raise UnexpectedException('Unexpected categorizing occurs.')

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)


    def _add_global_right_classification_label_df(self, pddf_row):
        """
        Global Right In-Out Classification
        RIGHT_PRECISE_KNOCK_IN {right_perfect_ki, right_nearly_ki}
        RIGHT_IMPRECISE_KNOCK_IN {right_nhej_ki, right_rc_nhej_ki}
        RIGHT_IMPRECISE KNOCK_IN_COMPLEX {right_imperfect_ki_complex, right_rc_imperfect_ki_complex}
        RIGHT_OTHER_MUATIONS {other_right_mutations, other_right_rc_mutations}
        RIGHT_NO_DETECTED {no_rightinout_indicators and no_rightinout_rc_indicators}
        """
        
        try:

            right_classification_label = pddf_row["right_classification_label"]
            right_rc_classification_label = pddf_row["right_rc_classification_label"]

            if is_in(right_classification_label, ("right_perfect_ki", "right_nearly_ki")):

                return pd.Series({'global_right_classification_label' : "RIGHT_PRECISE_KNOCK_IN"})

            elif (right_classification_label in "right_nhej_ki") or (right_rc_classification_label in "right_rc_nhej_ki"):

                return pd.Series({'global_right_classification_label' : "RIGHT_IMPRECISE_KNOCK_IN"})

            elif (right_classification_label in "right_imperfect_ki_complex") or (right_rc_classification_label in "right_rc_imperfect_ki_complex"):

                return pd.Series({'global_right_classification_label' : "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX"})

            elif (right_classification_label in "other_right_mutations") or (right_rc_classification_label in "other_right_rc_mutations"):
                
                return pd.Series({'global_right_classification_label' : "RIGHT_OTHER_MUATIONS"})

            elif (right_classification_label in "no_rightinout_indicators") and (right_rc_classification_label in "no_rightinout_rc_indicators"):
                
                return pd.Series({'global_right_classification_label' : "RIGHT_NO_DETECTED"})

            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                raise UnexpectedException('Unexpected categorizing occurs.')

        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)

    def _add_global_outout_classification_label_df(self, pddf_row):
        """
        Global Out-Out Classification
        OTHERS {unexpected_mutations, no_outout_indicators, multiple_outout_indicators}
        UNTREATED {untreated}
        PRECISE_KNOCK_IN {LEFT_PRECISE_KNOCK_IN and RIGHT_PRECISE_KNOCK_IN}
        IMPRECISE_KNOCK_IN {^(LEFT_PRECISE_KNOCK_IN and RIGHT_PRECISE_KNOCK_IN) and ^(LEFT_NO_DETECTED and RIGHT_NO_DETECTED)}
        INSERT_EDITING {overall_edited_insertion}
        DELETION_EDITING {overall_edited_deletion}
        SUBSTITUTION_EDITING {overall_edited_substitution}
        """

        try:

            global_left_classification_label = pddf_row["global_left_classification_label"]
            global_right_classification_label = pddf_row["global_right_classification_label"]
            mutation_classification_label = pddf_row["mutation_classification_label"]
            whole_classification_label = pddf_row["whole_classification_label"]

            is_both_precise = (global_left_classification_label in "LEFT_PRECISE_KNOCK_IN") and (global_right_classification_label in "RIGHT_PRECISE_KNOCK_IN")
            is_both_no_detected = (global_left_classification_label in "LEFT_NO_DETECTED") and (global_right_classification_label in "RIGHT_NO_DETECTED")

            # The Ordering of condition is important for classification because each condition is not independent.
            if is_in(mutation_classification_label, ("unexpected_mutations", "unidentified")):
                
                return pd.Series({'global_outout_classification_label' : "OTHERS"})

            elif mutation_classification_label in "untreated":
                
                return pd.Series({'global_outout_classification_label' : "UNTREATED"})

            elif is_both_precise:

                return pd.Series({'global_outout_classification_label' : "PRECISE_KNOCK_IN"})

            elif (not is_both_precise) and (not is_both_no_detected):

                return pd.Series({'global_outout_classification_label' : "IMPRECISE_KNOCK_IN"})

            elif mutation_classification_label in "overall_edited_insertion":

                return pd.Series({'global_outout_classification_label' : "INSERT_EDITING"})

            elif mutation_classification_label in "overall_edited_deletion":

                return pd.Series({'global_outout_classification_label' : "DELETION_EDITING"})

            elif mutation_classification_label in "overall_edited_substitution":
                
                return pd.Series({'global_outout_classification_label' : "SUBSTITUTION_EDITING"})

            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                raise ReferenceRegexException('Unexpected categorizing occurs.')
        
        except UnexpectedException as e:
            #traceback.print_exc(file = sys.stdout)
            error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
            sys.exit(1)


    def _add_global_classification_labels(self):
        tqdm.pandas(desc = "Adding global_left_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_global_left_classification_label_df, axis = 1)).copy()

        self._register_count("LEFT_PRECISE_KNOCK_IN", self.out_pddf.query('global_left_classification_label == "LEFT_PRECISE_KNOCK_IN"'))
        self._register_count("LEFT_IMPRECISE_KNOCK_IN", self.out_pddf.query('global_left_classification_label == "LEFT_IMPRECISE_KNOCK_IN"'))
        self._register_count("LEFT_IMPRECISE KNOCK_IN_COMPLEX", self.out_pddf.query('global_left_classification_label == "LEFT_IMPRECISE KNOCK_IN_COMPLEX"'))
        self._register_count("LEFT_OTHER_MUATIONS", self.out_pddf.query('global_left_classification_label == "LEFT_OTHER_MUATIONS"'))
        self._register_count("LEFT_NO_DETECTED", self.out_pddf.query('global_left_classification_label == "LEFT_NO_DETECTED"'))

        self._show_count("LEFT_PRECISE_KNOCK_IN")
        self._show_count("LEFT_IMPRECISE_KNOCK_IN")
        self._show_count("LEFT_IMPRECISE KNOCK_IN_COMPLEX")
        self._show_count("LEFT_OTHER_MUATIONS")
        self._show_count("LEFT_NO_DETECTED")

        tqdm.pandas(desc = "Adding global_right_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_global_right_classification_label_df, axis = 1)).copy()

        self._register_count("RIGHT_PRECISE_KNOCK_IN", self.out_pddf.query('global_right_classification_label == "RIGHT_PRECISE_KNOCK_IN"'))
        self._register_count("RIGHT_IMPRECISE_KNOCK_IN", self.out_pddf.query('global_right_classification_label == "RIGHT_IMPRECISE_KNOCK_IN"'))
        self._register_count("RIGHT_IMPRECISE KNOCK_IN_COMPLEX", self.out_pddf.query('global_right_classification_label == "RIGHT_IMPRECISE KNOCK_IN_COMPLEX"'))
        self._register_count("RIGHT_OTHER_MUATIONS", self.out_pddf.query('global_right_classification_label == "RIGHT_OTHER_MUATIONS"'))
        self._register_count("RIGHT_NO_DETECTED", self.out_pddf.query('global_right_classification_label == "RIGHT_NO_DETECTED"'))

        self._show_count("RIGHT_PRECISE_KNOCK_IN")
        self._show_count("RIGHT_IMPRECISE_KNOCK_IN")
        self._show_count("RIGHT_IMPRECISE KNOCK_IN_COMPLEX")
        self._show_count("RIGHT_OTHER_MUATIONS")
        self._show_count("RIGHT_NO_DETECTED")

        tqdm.pandas(desc = "Adding global_outout_classification_label")
        self.out_pddf = self.out_pddf.join(self.out_pddf.progress_apply(self._add_global_outout_classification_label_df, axis = 1)).copy()

        self._register_count("PRECISE_KNOCK_IN", self.out_pddf.query('global_outout_classification_label == "PRECISE_KNOCK_IN"'))
        self._register_count("IMPRECISE_KNOCK_IN", self.out_pddf.query('global_outout_classification_label == "IMPRECISE_KNOCK_IN"'))
        self._register_count("INSERT_EDITING", self.out_pddf.query('global_outout_classification_label == "INSERT_EDITING"'))
        self._register_count("DELETION_EDITING", self.out_pddf.query('global_outout_classification_label == "DELETION_EDITING"'))
        self._register_count("SUBSTITUTION_EDITING", self.out_pddf.query('global_outout_classification_label == "SUBSTITUTION_EDITING"'))
        self._register_count("UNTREATED", self.out_pddf.query('global_outout_classification_label == "UNTREATED"'))
        self._register_count("OTHERS", self.out_pddf.query('global_outout_classification_label == "OTHERS"'))

        self._show_count("PRECISE_KNOCK_IN")
        self._show_count("IMPRECISE_KNOCK_IN")
        self._show_count("INSERT_EDITING")
        self._show_count("DELETION_EDITING")
        self._show_count("SUBSTITUTION_EDITING")
        self._show_count("UNTREATED")
        self._show_count("OTHERS")
