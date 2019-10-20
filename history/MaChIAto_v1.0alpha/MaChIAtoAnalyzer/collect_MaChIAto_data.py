#!/usr/local/bin/python2
# -*- coding: utf-8 -*-

__version__ = "1.00"
__date__ = "20 May 2018"

import argparse
import sys
from os import listdir
from os.path import join
from os.path import abspath
import pandas as pd
import subprocess
from tqdm import tqdm
import regex as re

def _put_indicator(label, ctrl_indicator, cas9_indicator, knockin_indicator):
    if label in ctrl_indicator:
        return "ctrl"
    elif any(label in indicator for indicator in cas9_indicator):
        for index, indicator in enumerate(cas9_indicator):
            if label in indicator:
                return "Cas9_" + str(index + 1)
    elif any(label in indicator for indicator in knockin_indicator):
        for index, indicator in enumerate(knockin_indicator):
            if label in indicator:
                return "KI_" + str(index + 1)
    else:
        raise("Unexpected label.")


def main():
    parser = argparse.ArgumentParser(description='collect_MaChIAto_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # arguments
    parser.add_argument('-i','--indir', type=str,  help='input directory', required=True)
    parser.add_argument('-o','--outdir', type=str,  help='output directory', required=True)
    parser.add_argument('-co','--ctrl_indicator', type=str,  help='control indicator (This must be one.)', required=True)
    parser.add_argument('-ca','--cas9_indicator', type=str,  help='Cas9 indicator', required=True, nargs='+')
    parser.add_argument('-ki','--knockin_indicator', type=str,  help='Knock-in indicator', required=True, nargs='+')
    parser.add_argument('-sc','--scaffold_seq', type=str,  help='scaffold sequence', required=True, nargs='+')
    args = parser.parse_args()

    OUTPUT_DIRECTORY = abspath(args.outdir)
    INPUT_DIRECTORY = abspath(args.indir)
    indir_list = [filename for filename in listdir(INPUT_DIRECTORY) if not filename.startswith('.')]

    rate_temp_pddf = pd.DataFrame()
    consistency_temp_pddf = pd.DataFrame()

    label_type_list = pd.DataFrame({"label" : pd.Series([args.ctrl_indicator]), "sample.type": pd.Series(["ctrl"])})
    for index, indicator in enumerate(args.cas9_indicator):
        temp_label_type_list = pd.DataFrame({"label" : pd.Series([indicator]), "sample.type": pd.Series(["Cas9_" + str(index + 1)])})
        label_type_list = pd.concat([label_type_list, temp_label_type_list])
    for index, indicator in enumerate(args.knockin_indicator):
        temp_label_type_list = pd.DataFrame({"label" : pd.Series([indicator]), "sample.type": pd.Series(["KI_" + str(index + 1)])})
        label_type_list = pd.concat([label_type_list, temp_label_type_list])
    print(label_type_list)
    
    i = 1
    for dirname in tqdm(indir_list):
        dirname_words = dirname.split('_')
        if len(dirname_words) == 8:
            MaChIAto, fr, CRISPResso, at, date1, on, date2, samplename = dirname.split('_')
        elif len(dirname_words) == 7:
            MaChIAto, fr, CRISPResso, at, date, on, samplename = dirname.split('_')
        else:
            raise Exception('Unexpected number')
        
        # collect data
        rate_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'ALL_read_rates.csv'), sep=',')
        consistency_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'CRISPesso_MaChiIAto_consistency_persentage.csv'), sep=',')
        read_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'ALL_read_countdata.csv'), sep=',')
        seq_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'MaChiIAto_specific_sequences.csv'), sep=',')
        param_temp_pddf = pd.read_csv(join(INPUT_DIRECTORY, dirname, 'MaChiIAto_optimized_param.csv'), sep=',')

        # search forty base surrounding sgRNA
        fourtybase_regex_str = "[atcgATCG]{10}"\
            + param_temp_pddf['guide_seq'].values[0]\
            + "[atcgATCG]{10}"
        fourtybase_regex = re.compile(fourtybase_regex_str)
        fourtybase_factors = fourtybase_regex.search(param_temp_pddf['amplicon_seq'].values[0])
        if fourtybase_factors is None:
            raise Exception('40bp sgRNA surround sequence is not found.')
        fourtybase = fourtybase_factors.group(0)

        gene, label = samplename.split('-')

        # add sequence info to rate table
        rate_temp_pddf['lmh'] = pd.Series(seq_temp_pddf['Leftside homology arm'], index=rate_temp_pddf.index)
        rate_temp_pddf['rmh'] = pd.Series(seq_temp_pddf['Rightside homology arm'], index=rate_temp_pddf.index)
        rate_temp_pddf['bmh'] = rate_temp_pddf['lmh'] + rate_temp_pddf['rmh']
        rate_temp_pddf['protospacer'] = pd.Series(param_temp_pddf['guide_seq'], index=rate_temp_pddf.index)
        rate_temp_pddf['fourtybase.region.surrounding.targetseq'] = pd.Series(fourtybase, index=rate_temp_pddf.index)
        rate_temp_pddf['sample.name'] = pd.Series(samplename, index=rate_temp_pddf.index)
        rate_temp_pddf['group.name'] = pd.Series(gene, index=rate_temp_pddf.index)
        rate_temp_pddf['sample.type'] = pd.Series(_put_indicator(label, args.ctrl_indicator, args.cas9_indicator, args.knockin_indicator), index=rate_temp_pddf.index)
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
    
    # order
    rate_sum_pddf = rate_sum_pddf.sort_values(by='name') 
    consistency_sum_pddf = consistency_sum_pddf.sort_values(by='name')
    read_sum_pddf = read_sum_pddf.sort_values(by='name')
    seq_sum_pddf = seq_sum_pddf.sort_values(by='name')

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

    # save
    rate_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_ALL_read_rates.csv'), index=False)
    consistency_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_CRISPesso_MaChiIAto_consistency_persentage.csv'), index=False)
    read_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_ALL_read_countdata.csv'), index=False)
    seq_sum_pddf.to_csv(join(OUTPUT_DIRECTORY, 'sum_MaChiIAto_specific_sequences.csv'), index=False)
    label_type_list.to_csv(join(OUTPUT_DIRECTORY, 'label_sample.type_table.csv'), index=False)

    return 0

if __name__ == '__main__':
    main()
quit()
