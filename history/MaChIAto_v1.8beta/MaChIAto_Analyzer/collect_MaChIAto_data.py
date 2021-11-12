#!/usr/local/bin/python3.7
# -*- coding: utf-8 -*-

# @TODO: 全く同一の標的とサンプルラベルをもつものは平均を算出して、一つの列にまとめる

__author__ = "Kazuki Nakamae <kazkinakamae@gmail.com>"
__version__ = "3.00"
__date__ = "4 March 2019"

import argparse
import sys
from os import listdir
from os.path import join
from os.path import abspath
import subprocess
import regex as re
from tqdm import tqdm
import pandas as pd

def _put_indicator(label, untreated_label, knock_out_label_list, knock_in_label_list):
    """Convert label into sample labels

    Convert label into sample labels according to the commandline argument.

    Args:
        label: converted label
        untreated_label: untreated label of commandline argument
        knock_out_label_list: knock_out label list of commandline argument
        knock_in_label_list: knock_in label list of commandline argument

    Returns:
        string object of a matched label that {untreated_label, knock_out_label_list, knock_in_label_list} include.
        If the matched label is not found, the function returns "other" label.

    """
    if label in untreated_label:
        return "untreated"
    elif any(label in indicator for indicator in knock_out_label_list):
        for index, indicator in enumerate(knock_out_label_list):
            if label in indicator:
                return "knock-out_" + str(index + 1)
    elif any(label in indicator for indicator in knock_in_label_list):
        for index, indicator in enumerate(knock_in_label_list):
            if label in indicator:
                return "knock-in_" + str(index + 1)
    else:
        return "other"

###EXCEPTIONS############################

class LabelingException(Exception):
    pass

class FormatException(Exception):
    pass

class NoLabelingException(Exception):
    pass

#########################################

def main():
    try:
        parser = argparse.ArgumentParser(description='collect_MaChIAto_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # arguments
        parser.add_argument('-i','--indir', type=str,  help='input directory which includes the results of MaChIAto Classifier', required=True)
        parser.add_argument('-o','--outdir', type=str,  help='output directory', required=True)
        parser.add_argument('-ul','--untreated_label', type=str,  help='untreated label (This must be one.)', default="machiato_dummy_sample", required=False)
        parser.add_argument('-ol','--knock_out_label', type=str,  help='negative control label', default=[], required=False, nargs='+')
        parser.add_argument('-il','--knock_in_label', type=str,  help='knock-in_label', default=[], required=False, nargs='+')
        parser.add_argument('-sc','--scaffold_seq', type=str,  help='scaffold sequence of sgRNA', default='gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc', required=False)
        parser.add_argument('-t','--target_type', type=str,  help='target sequence type <lmh | rmh | bmh |elmh | ermh | ebmh| protospacer>', default='bmh', required=False)
        parser.add_argument('--ignore_list', type=str,  help='The list of ignore target set contains target names which are not desired to analyze for some reasons. The data (e.g. DBF4B-A, DBF4B-B, DBF4B-C, DBF4B-D) including target name (e.g. DBF4B) shown in the list is skipped through the process of MaChIAto Analyzer. The format should be comma-separated like “TargetA, TargetB, …” “example_data” directory has “ignore_list.csv” as example.', default='', required=False)
        # parser.add_argument('--fake_untreated', help='When you use the flag, An untreated sample is automatically generated. The sample has zero value in each item.', action='store_true', required=False)
        args = parser.parse_args()

        # set input/output directory path
        OUTPUT_DIRECTORY = abspath(args.outdir)
        INPUT_DIRECTORY = abspath(args.indir)
        indir_list = [filename for filename in listdir(INPUT_DIRECTORY) if not filename.startswith('.')]
        if args.ignore_list != "":
            ignore_pddf = pd.read_csv(args.ignore_list)
        else:
            ignore_pddf = pd.Series([])
        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

        rate_temp_pddf = pd.DataFrame()
        consistency_temp_pddf = pd.DataFrame()

        # Insert fake untreated sample

        # untreated_label = args.untreated_label
        # if (args.fake_untreated) or (untreated_label is None):
        #     if "dummy_untreated" == args.untreated_label:
        #         untreated_label = "_dummy_untreated"
        #     else:
        #         untreated_label = "dummy_untreated"

        # Check labeling pattrn

        # import pdb; pdb.set_trace()
        has_untreated_label = (type(args.untreated_label) is str)
        has_knock_out_label = (len(args.knock_out_label) != 0)
        is_multi_knock_out_label = (len(args.knock_out_label) > 1)
        has_knock_in_label = (len(args.knock_in_label) != 0)
        is_multi_knock_in_label = (len(args.knock_in_label) > 1)

        # There are available labeling pattrns.
        # 1) Double knock-in analysis : -ul A -ol B -il C D
        # 2) Single knock-in analysis : -ul A -ol B -il C
        # 3) Simple knock-in analysis : -il C
        # 4) Double knock-out analysis : -ul A -ol B C
        # 5) Single knock-out analysis : -ul A -ol B

        if (not has_untreated_label) and (not has_knock_out_label) and (not has_knock_in_label):# X
            print('Label dose not exist.')
            raise LabelingException('Labeling Error')
        elif (has_knock_out_label) and (has_knock_in_label):# 1, 2) Double knock-in analysis / Single knock-in analysis
            if is_multi_knock_in_label:# 1) Double knock-in analysis
                print("ANALYSIS TYPE: Double knock-in analysis")
                mode_num = 1
            else:# 2) Single knock-in analysis
                print("ANALYSIS TYPE: Single knock-in analysis")
                mode_num = 2
        elif (not has_knock_out_label) and (has_knock_in_label):# 3) Simple knock-in analysis
            print("ANALYSIS TYPE: Simple knock-in analysis")
            mode_num = 3
        elif (has_knock_out_label) and (not has_knock_in_label):# 4, 5) Double knock-out analysis / Single knock-out analysis
            if is_multi_knock_out_label:# 4) Double knock-out analysis
                print("ANALYSIS TYPE: Double knock-out analysis")
                mode_num = 4
            else:# 5) Single knock-out analysis
                print("ANALYSIS TYPE: Single knock-out analysis")
                mode_num = 5
        else: # X
            print('Available label pattrn is not found.')
            raise LabelingException('Labeling Error')
        
        label_type_list = pd.DataFrame({"label" : pd.Series([args.untreated_label]), "sample.type": pd.Series(["untreated"])})
        if mode_num != 3:
            for index, indicator in enumerate(args.knock_out_label):
                temp_label_type_list = pd.DataFrame({"label" : pd.Series([indicator]), "sample.type": pd.Series(["knock-out_" + str(index + 1)])})
                label_type_list = pd.concat([label_type_list, temp_label_type_list])
        if mode_num < 4:
            for index, indicator in enumerate(args.knock_in_label):
                temp_label_type_list = pd.DataFrame({"label" : pd.Series([indicator]), "sample.type": pd.Series(["knock-in_" + str(index + 1)])})
                label_type_list = pd.concat([label_type_list, temp_label_type_list])
        print(label_type_list)
        
        # Get output of MaChIAto Classifier
        i = 1
        for dirname in tqdm(indir_list):
            dirname_words = dirname.split('_')
            if len(dirname_words) == 8:
                MaChIAto, fr, CRISPResso, at, date1, on, date2, samplename = dirname.split('_')
            elif len(dirname_words) == 7:
                MaChIAto, fr, CRISPResso, at, date, on, samplename = dirname.split('_')
            else:
                raise FormatException("Unexpected number of words:" + str(len(dirname_words)))
            
            # split label with sample name
            if len(samplename.split('-')) == 2:
                gene, label = samplename.split('-')
            elif len(samplename.split('-')) == 1:
                if (mode_num != 3) and (mode_num != 5):
                    raise NoLabelingException(samplename)
                else:
                    gene = samplename
                    # If there is no label in samplename, register label according to commandline argument.
                    if mode_num == 5: # knock-out
                        label = args.knock_out_label[1]
                    elif mode_num == 3: # knock-in
                        label = args.knock_in_label[1]
                    else:
                        raise NoLabelingException(samplename)
            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                raise NoLabelingException(samplename)
            
            # exclude gene that is listed in ignore_list
            # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
            if gene in ignore_pddf:
                print('Ignored ' + gene)
                continue
            
            # collect data
            rate_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'ALL_read_rates.csv'), sep=',')
            consistency_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'CRISPResso_MaChIAto_consistency_persentage.csv'), sep=',')
            read_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'ALL_read_countdata.csv'), sep=',')
            seq_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'MaChIAto_specific_sequences.csv'), sep=',')
            param_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'MaChIAto_optimized_param.csv'), sep=',')

            # search forty base surrounding sgRNA
            fourtybase_regex_str = "[atcgATCG]{10}"\
                + param_temp_pddf['guide_seq'].values[0]\
                + "[atcgATCG]{10}"
            fourtybase_regex = re.compile(fourtybase_regex_str)
            fourtybase_factors = fourtybase_regex.search(param_temp_pddf['amplicon_seq'].values[0])
            if fourtybase_factors is None:
                raise FormatException('40bp sgRNA surround sequence is not found.')
            fourtybase = fourtybase_factors.group(0)

            # add sequence info to rate table
            # If there is untreated/knock-out sample, 
            rate_temp_pddf['lmh'] = pd.Series(seq_temp_pddf['Leftside homology arm'], index=rate_temp_pddf.index)
            rate_temp_pddf['rmh'] = pd.Series(seq_temp_pddf['Rightside homology arm'], index=rate_temp_pddf.index)
            rate_temp_pddf['bmh'] = rate_temp_pddf['lmh'] + rate_temp_pddf['rmh']
            rate_temp_pddf['elmh'] = pd.Series(seq_temp_pddf['Leftside edited homology arm'], index=rate_temp_pddf.index)
            rate_temp_pddf['ermh'] = pd.Series(seq_temp_pddf['Rightside edited homology arm'], index=rate_temp_pddf.index)
            rate_temp_pddf['ebmh'] = rate_temp_pddf['elmh'] + rate_temp_pddf['ermh']
            rate_temp_pddf['protospacer'] = pd.Series(param_temp_pddf['guide_seq'], index=rate_temp_pddf.index)
            rate_temp_pddf['scaffold'] = pd.Series([args.scaffold_seq]).repeat(len(param_temp_pddf['guide_seq']))
            rate_temp_pddf['fourtybase.region.surrounding.targetseq'] = pd.Series(fourtybase, index=rate_temp_pddf.index)
            rate_temp_pddf['sample.name'] = pd.Series(samplename, index=rate_temp_pddf.index)
            rate_temp_pddf['group.name'] = pd.Series(gene, index=rate_temp_pddf.index)
            rate_temp_pddf['sample.type'] = pd.Series(_put_indicator(label, args.untreated_label, args.knock_out_label, args.knock_in_label), index=rate_temp_pddf.index)
            rate_temp_pddf['name'] = pd.Series([samplename], index=rate_temp_pddf.index)
            consistency_temp_pddf['name'] = pd.Series([samplename], index=consistency_temp_pddf.index)
            read_temp_pddf['name'] = pd.Series([samplename], index=read_temp_pddf.index)
            seq_temp_pddf['name'] = pd.Series([samplename], index=seq_temp_pddf.index)
            if i == 1:
                rate_sum_pddf = rate_temp_pddf.copy()
                consistency_sum_pddf = consistency_temp_pddf.copy()
                read_sum_pddf = read_temp_pddf.copy()
                seq_sum_pddf = seq_temp_pddf.copy()
            else:
                # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
                rate_sum_pddf = pd.concat([rate_sum_pddf, rate_temp_pddf])
                consistency_sum_pddf = pd.concat([consistency_sum_pddf, consistency_temp_pddf])
                read_sum_pddf = pd.concat([read_sum_pddf, read_temp_pddf])
                seq_sum_pddf = pd.concat([seq_sum_pddf, seq_temp_pddf])
            i = i + 1
        
        # sort
        rate_sum_pddf = rate_sum_pddf.sort_values(by='name') 
        consistency_sum_pddf = consistency_sum_pddf.sort_values(by='name')
        read_sum_pddf = read_sum_pddf.sort_values(by='name')
        seq_sum_pddf = seq_sum_pddf.sort_values(by='name')

        # adjust col/row name
        rate_sum_cols = rate_sum_pddf.columns.tolist()
        consistency_sum_cols = consistency_sum_pddf.columns.tolist()
        read_sum_pddf_cols = read_sum_pddf.columns.tolist()
        seq_sum_pddf_cols = seq_sum_pddf.columns.tolist()

        rate_sum_cols = rate_sum_cols[-1:] + rate_sum_cols[:-1]
        consistency_sum_cols = consistency_sum_cols[-1:] + consistency_sum_cols[:-1]
        read_sum_pddf_cols = read_sum_pddf_cols[-1:] + read_sum_pddf_cols[:-1]
        seq_sum_pddf_cols = seq_sum_pddf_cols[-1:] + seq_sum_pddf_cols[:-1]

        rate_sum_pddf = rate_sum_pddf[rate_sum_cols]
        consistency_sum_pddf = consistency_sum_pddf[consistency_sum_cols]
        read_sum_pddf = read_sum_pddf[read_sum_pddf_cols]
        seq_sum_pddf = seq_sum_pddf[seq_sum_pddf_cols]

        ### aggregate by target type # transfered into R script
        # string_cname = ["name", "lmh", "rmh", "bmh", "elmh", "ermh", "ebmh", "protospacer", "scaffold", "fourtybase.region.surrounding.targetseq", "sample.name", "group.name", "sample.type"]
        # average_rate_sum_pddf = rate_sum_pddf.groupby(by=[args.target_type, "sample.type"], axis=0).mean()
        # for row_ind, target_multi_index in enumerate(average_rate_sum_pddf.index):
        #     temp_index_list = (rate_sum_pddf[args.target_type] == target_multi_index[0]) & (rate_sum_pddf["sample.type"] == target_multi_index[1])
        #     temp_pddf = pd.concat([ average_rate_sum_pddf.loc[target_multi_index, ], rate_sum_pddf[temp_index_list].iloc[0, ][string_cname] ])
        #     if row_ind == 0:
        #         average_rate_sum_seq_pddf = temp_pddf
        #     else:
        #         average_rate_sum_seq_pddf = pd.concat([average_rate_sum_seq_pddf, temp_pddf], axis=1)
        # average_rate_sum_seq_pddf = average_rate_sum_seq_pddf.transpose()
        # import pdb; pdb.set_trace()
        # sort
        # average_rate_sum_seq_pddf = average_rate_sum_seq_pddf.sort_values(by='name') 
        ###

        # save
        rate_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_ALL_read_rates.csv'), index=False)
        # average_rate_sum_seq_pddf.to_csv(join(OUTPUT_DIRECTORY, 'average_read_rates.csv'), index=False) # transfered into R script
        consistency_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_CRISPResso_MaChIAto_consistency_persentage.csv'), index=False)
        read_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_ALL_read_countdata.csv'), index=False)
        seq_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_MaChIAto_specific_sequences.csv'), index=False)
        label_type_list.to_csv(join(OUTPUT_DIRECTORY, 'label_sample_type.csv'), index=False)
        pd.Series([args.target_type]).to_csv(join(OUTPUT_DIRECTORY, 'target_type.csv'), header=False, index=False)

        sys.exit(0)
    except LabelingException as e:
        print("""
    There are available labeling pattrns.
    1) Double knock-in analysis : -ul A -ol B -il C D
    2) Single knock-in analysis : -ul A -ol B -il C
    3) Simple knock-in analysis : -il C
    4) Double knock-out analysis : -ul A -ol B C
    5) Single knock-out analysis : -ul A -ol B
    You can choose a pattrn from these pattrns.
    If you do not enter untreated sample (-ul), the tool automatically enter "machiato_dummy_sample".
            """)
        print('%s' % e)
        sys.exit(1)
    except FormatException as e:
        print('%s' % e)
        sys.exit(1)
    except NoLabelingException as e:
        print("When you analize non-simple knock-out/knock-in analysis, you must add label to sample name (e.g. [sample name]-[label]).")
        sys.exit(1)

if __name__ == '__main__':
    main()
quit()
