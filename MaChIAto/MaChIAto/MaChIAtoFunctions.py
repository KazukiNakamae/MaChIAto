# -*- coding: utf8 -*-
"""
MaChIAto - Kazuki Nakamae 2019
Software pipeline for the analysis of outcomes of MMEJ-assisted gene knock-in using CRISPR-Cas9 from deep sequencing data
https://github.com/Kazuki-Nakamae/MaChIAto
"""

### FUNCTIONS ############################

def convert_crispresso2row_to_crispresso(crispresso2_row):

    err_msg = 'This input file is not CRISPResso2 format. Check input and setting.'

    if crispresso2_row['Reference_Name'] == "Reference":
        
        if crispresso2_row['Read_Status'] == "UNMODIFIED":

            is_nhej, is_unmodified, is_hdr = False, True, False
        elif crispresso2_row['Read_Status'] == "MODIFIED":

            is_nhej, is_unmodified, is_hdr = True, False, False
        else:
            error(err_msg)
            sys.exit(1)

    elif crispresso2_row['Reference_Name'] == "HDR":

        if crispresso2_row['Read_Status'] == "UNMODIFIED":

            is_nhej, is_unmodified, is_hdr = False, False, True
        elif crispresso2_row['Read_Status'] == "MODIFIED":

            is_nhej, is_unmodified, is_hdr = False, False, False
        else:
            error(err_msg)
            sys.exit(1)
    elif crispresso2_row['Reference_Name'] == "AMBIGUOUS_Reference":

        is_nhej, is_unmodified, is_hdr = True, False, False
    else:
        error(err_msg)
        sys.exit(1)
    
    return crispresso2_row['Aligned_Sequence'],\
        crispresso2_row['Reference_Sequence'],\
        is_nhej,\
        is_unmodified,\
        is_hdr,\
        crispresso2_row['n_deleted'],\
        crispresso2_row['n_inserted'],\
        crispresso2_row['n_mutated'],\
        crispresso2_row['#Reads'],\
        crispresso2_row['%Reads']

# NOTE : dependent on "sys"
def should_stop():
    while True:
        choice = raw_input("Do you wish to stop this program ? (Press y|Y for Yes, any other key for No)").lower()
        if choice in ['y', 'Y']:
            sys.exit(0)
        else:
            break

def reverse_complement(sequence):
    """
    @function   reverse_complement();
    make reverse complement of a DNA sequence
    @param  {string} sequence :  DNA sequence
    @return {string} reverse complement of a DNA sequence
    """
    nt_complement = dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
    return "".join([nt_complement[c] for c in sequence.upper()[-1::-1]])

def is_in(str, tuple_str):
    f = lambda x:x in str
    if any(map(f, tuple_str)):
        return True
    else:
        return False

def nonzero_float_division(a, b):
    if b == 0:
        return 0.
    else:
        return float(a) / float(b)

def findall_both(regex, rc_regex, seq):
    match_target_seq = regex.findall(seq, overlapped=True)
    opp_match_target_seq = rc_regex.findall(reverse_complement(seq), overlapped=True)
    match_target_seq.extend([reverse_complement(opp_seq) for opp_seq in opp_match_target_seq])
    match_target_seq_set = set(match_target_seq)
    match_target_seq = list(match_target_seq_set)
    # remove null string
    match_target_seq = [x for x in match_target_seq if x]
    return(match_target_seq)

def calc_MaChIAto_rates_from_CRISPResso(read_cnt_dict_sum):

    return {
        "Left accuracy" : nonzero_float_division(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"]),
        "Left MMEJ-tendency" : nonzero_float_division(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"]),
        "Left NHEJ-tendency" : nonzero_float_division(read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"]),
        "Left unpredictability" : nonzero_float_division(read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"]
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"]),
        "Left unpredictable knock-in insert" : nonzero_float_division(read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"]),
        "Left unpredictable knock-in deletion" : nonzero_float_division(read_cnt_dict_sum["LEFT_OTHER_MUATIONS"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"]),
        "Left efficiency" : nonzero_float_division(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"] \
            + read_cnt_dict_sum["LEFT_NO_DETECTED"]),
        "Left occupancy" : nonzero_float_division(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"] \
            + read_cnt_dict_sum["LEFT_NO_DETECTED"]),
        "Left MMEJ/NHEJ joinability" : nonzero_float_division(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["LEFT_OTHER_MUATIONS"] \
            + read_cnt_dict_sum["LEFT_NO_DETECTED"]),
        "Right accuracy" : nonzero_float_division(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"]),
        "Right MMEJ-tendency" : nonzero_float_division(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"]),
        "Right NHEJ-tendency" : nonzero_float_division(read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"]),
        "Right unpredictability" : nonzero_float_division(read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"] \
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"]),
        "Right unpredictable knock-in insert" : nonzero_float_division(read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"]
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"]),
        "Right unpredictable knock-in deletion" : nonzero_float_division(read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"]
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"]),
        "Right efficiency" : nonzero_float_division(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"] \
            + read_cnt_dict_sum["RIGHT_NO_DETECTED"]),
        "Right occupancy" : nonzero_float_division(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"] \
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"] \
            + read_cnt_dict_sum["RIGHT_NO_DETECTED"]),
        "Right MMEJ/NHEJ joinability" : nonzero_float_division(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"] \
            + read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"] \
            + read_cnt_dict_sum["RIGHT_NO_DETECTED"]),
        "Precise Knock-in Efficiency" : nonzero_float_division(read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"] \
            + read_cnt_dict_sum["UNTREATED"] \
            + read_cnt_dict_sum["OTHERS"]),
        "Knock-in Efficiency" : nonzero_float_division(read_cnt_dict_sum["PRECISE_KNOCK_IN"]
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"] \
            + read_cnt_dict_sum["UNTREATED"] \
            + read_cnt_dict_sum["OTHERS"]),
        "Indels Editing Efficiency" : nonzero_float_division(read_cnt_dict_sum["PRECISE_KNOCK_IN"]
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"] \
            + read_cnt_dict_sum["UNTREATED"] \
            + read_cnt_dict_sum["OTHERS"]),
        "Editing Efficiency" : nonzero_float_division(read_cnt_dict_sum["PRECISE_KNOCK_IN"]
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"] \
            + read_cnt_dict_sum["UNTREATED"] \
            + read_cnt_dict_sum["OTHERS"]),
        "Accuracy" : nonzero_float_division(read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"]),
        "Occupancy" : nonzero_float_division(read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"]),
        "Conservative" : nonzero_float_division(read_cnt_dict_sum["UNTREATED"] \
            , read_cnt_dict_sum["PRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["IMPRECISE_KNOCK_IN"] \
            + read_cnt_dict_sum["INSERT_EDITING"] \
            + read_cnt_dict_sum["DELETION_EDITING"] \
            + read_cnt_dict_sum["SUBSTITUTION_EDITING"] \
            + read_cnt_dict_sum["UNTREATED"] \
            + read_cnt_dict_sum["OTHERS"]),
        "Precise Knock-in zygomorphy" : -int(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"]) \
            + int(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"]),
        "NHEJ-assisted Knock-in zygomorphy" : -int(read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"]) \
            + int(read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"]),
        "Integration zygomorphy" : -int(read_cnt_dict_sum["LEFT_PRECISE_KNOCK_IN"]) \
            - int(read_cnt_dict_sum["LEFT_IMPRECISE_KNOCK_IN"]) \
            - int(read_cnt_dict_sum["LEFT_IMPRECISE KNOCK_IN_COMPLEX"]) \
            - int(read_cnt_dict_sum["LEFT_OTHER_MUATIONS"]) \
            + int(read_cnt_dict_sum["RIGHT_PRECISE_KNOCK_IN"]) \
            + int(read_cnt_dict_sum["RIGHT_IMPRECISE_KNOCK_IN"]) \
            + int(read_cnt_dict_sum["RIGHT_IMPRECISE KNOCK_IN_COMPLEX"]) \
            + int(read_cnt_dict_sum["RIGHT_OTHER_MUATIONS"]),
        "Left reversibility" : -int(read_cnt_dict_sum["left_rc_nhej_ki"]) \
            - int(read_cnt_dict_sum["left_rc_imperfect_ki_complex"]) \
            - int(read_cnt_dict_sum["other_left_rc_mutations"]) \
            + int(read_cnt_dict_sum["left_nhej_ki"]) \
            + int(read_cnt_dict_sum["left_imperfect_ki_complex"]) \
            + int(read_cnt_dict_sum["other_left_mutations"]),
        "Right reversibility" : -int(read_cnt_dict_sum["right_rc_nhej_ki"]) \
            - int(read_cnt_dict_sum["right_rc_imperfect_ki_complex"]) \
            - int(read_cnt_dict_sum["other_right_rc_mutations"]) \
            + int(read_cnt_dict_sum["right_nhej_ki"]) \
            + int(read_cnt_dict_sum["right_imperfect_ki_complex"]) \
            + int(read_cnt_dict_sum["other_right_mutations"])}
"""
def show_count(MaChIAtoClassifier, label):
    logging.info('%s(%s) / %s(%s) allele types (reads) have "%s" flag.' % \
        (str(MaChIAtoClassifier.cnt_dict[key]), \
            str(MaChIAtoClassifier.read_cnt_dict[key]), \
            str(MaChIAtoClassifier.cnt_dict["total"]), \
            str(MaChIAtoClassifier.read_cnt_dict["total"]), str(key)))
"""

