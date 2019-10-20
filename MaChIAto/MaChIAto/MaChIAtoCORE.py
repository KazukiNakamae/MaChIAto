# -*- coding: utf8 -*-
"""
MaChIAto - Kazuki Nakamae 2019
Software pipeline for the analysis of outcomes of MMEJ-assisted gene knock-in using CRISPR-Cas9 from deep sequencing data
https://github.com/Kazuki-Nakamae/MaChIAto
"""

__version__ = "1.5.0"

from MaChIAtoShared import *
from MaChIAtoClassifier import MaChIAtoClassifier
from MaChIAtoFunctions import convert_crispresso2row_to_crispresso, is_in, calc_MaChIAto_rates_from_CRISPResso

def main():
	try:
		print '  \n===\MaChIAto/==='
		print """
-Software pipeline for the analysis of outcomes of MMEJ-assisted gene knock-in using CRISPR-Cas9 from deep sequencing data-
		"""
		print """
			__  __        ____ _     ___    _   _        
			|  \/  | __ _ / ___| |__ |_ _|  / \ | |_ ___  
			| |\/| |/ _` | |   | '_ \ | |  / _ \| __/ _ \ 
			| |  | | (_| | |___| | | || | / ___ \ || (_) |
			|_|  |_|\__,_|\____|_| |_|___/_/   \_\__\___/ 
──────────────────────────────────────────────────────────────────────────────────────────────────────
		"""
        
		parser = argparse.ArgumentParser(description='MaChIAto Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
		# Common arguments
		parser.add_argument('-m','--mode', type=str,  help='This parameter allows for the specification of type of analysis: “CRISPResso” and “CRISPResso2” is allowed in the latest version', required=True)
		parser.add_argument('-a','--amplicon_seq', type=str,  help='This parameter allows the user to enter the amplicon sequence used for the CRISPResso. The length should be >105bp due to setting for the other parameter.', required=True)
		parser.add_argument('-g','--guide_seq', type=str,  help="This parameter allows for the specification of the sgRNA sequence used for the CRISPResso. The length of sequence should be 20nt without PAM. The MaChIAto convention is to depict the expected cleavage position using the value of the parameter 3 nt 3' from the end of the guide.", required=True)
		parser.add_argument('-d','--donor_seq', type=str,  help='This parameter allows for the specification of the expected HDR amplicon used for the CRISPResso. The length of sequence should be >12bp. In knock-out analysis, this parameter should not be entered, and then the value will be given fake parameter (“TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT”) and some internal settings is changed for knock-out analysis. However, fake parameter will be poly-C|G|A if amplicon sequence contains poly-T sequence.', default='', required=False)
		parser.add_argument('-e','--expected_ki_amplicon_seq', type=str,  help='This parameter allows for the specification of the expected knock-in amplicon sequence which used for the CRISPResso after HDR. The length of sequence should be >12bp. In knock-out analysis, this parameter should not be entered, and then the value will be given fake parameter including fake donor sequence and some internal settings is changed for knock-out analysis.', default='', required=False)
		parser.add_argument('-o','--output_folder', type=str,  help='This parameter allows for the specification of the output directory to use for the analysis (default: current directory)', default='./', required=False)
		parser.add_argument('-lh','--length_left_homologyarm', type=int,  help='This parameter allows for the specification of the length of 5’ homology arm (default: 20). The length of sequence should be >17bp, and the flanking sequence with the homology arm needs length of >24bp in the expected amplicon sequence.', default='20', required=False)
		parser.add_argument('-rh','--length_right_homologyarm', type=int,  help='This parameter allows for the specification of the length of 3’ homology arm (default: 20). The length of sequence should be >17bp, and the flanking sequence with the homology arm needs length of >24bp in the expected amplicon sequence.', default='20', required=False)
		# For CRISPResso
		parser.add_argument('-cf','--crispreeso_file', type=str, help='This parameter allows for the specification of the “Alleles_frequency_table.txt” from CRISPResso. When this parameter is used, “CRISPResso” should be entered as -m parameter.', default='./Alleles_frequency_table.txt', required=False)
		# For CRISPResso2
		parser.add_argument('-ccf','--crispreeso2_file', type=str, help='This parameter allows for the specification of the “Alleles_frequency_table.zip” from CRISPResso2. When this parameter is used, “CRISPResso2” should be entered as -m parameter.', default='./Alleles_frequency_table.zip', required=False)
		# Optional
		parser.add_argument('-n','--name',  help='This parameter allows for the specification of the name which will be included output directory (default: “untitled”). If MaChIAto Analyzer and MaChIAto Reviewer will be used in the following analysis, the value should be “<target_name>-<sample label >” (e.g. DBF4B-C)', default='untitled', required=False)
		parser.add_argument('--force_knockout_analysis', help='Usually, MaChIAto re-classify the data as knock-in analysis. This option forces the setting to change for knock-out analysis. Under this mode, the length of indicator on knock-in donor is maximum value, and threshold value for alignment of knock-in sequence is 1.0. If this parameter is not entered, MaChIAto can automatically set up this mode by finding some characteristics of knock-out sample. For example, MaChIAto checks that there is no donor sequence or expected knock-in sequence as input, and there is less 3 kinds HDR variants among input data.', action='store_true', required=False)
		parser.add_argument('--skip_optimization', help='Usually, MaChIAto run Bayesian optimization for finding the optimized setting. This option allows MaChIAto to skip the process of optimization. The option is made for debugging. So, the option should not be used in usual analysis. However, if the optimization process disturbs an accurate analysis, this option might be useful.', action='store_true', required=False)

		args = parser.parse_args()

		### LAMDA FUNCTIONS ############################
		_time = lambda : datetime.datetime.now().strftime("%Y%m%d%H%M%S")

		database_id=_time()

		is_knockout_analysis = False

		if (args.expected_ki_amplicon_seq is "") or (args.donor_seq is ""):
    		# set knock-out analysis flag
			info("MaChIAto runs knock-out analysis.")
			is_knockout_analysis = True
			info("Make fake parameter for knock-out analysis.")
			# cut amplicon seq at sgRNA site to make fake expected_ki_amplicon_seq
			guide_searched = re.search('(^[ATCGatcg]{17})([ATCGatcg]{3}$)', args.guide_seq)
			amplicon_searched = re.search('([ATCGatcg]*' + guide_searched.group(1) + ')(' + guide_searched.group(2) + '[ATCGatcg]*)', args.amplicon_seq)
			# check poly-N sequence in amplicon_seq
			poly_t_matched = re.search('[tT]{5}', args.amplicon_seq)
			poly_c_matched = re.search('[cC]{5}', args.amplicon_seq)
			poly_g_matched = re.search('[gG]{5}', args.amplicon_seq)
			poly_a_matched = re.search('[aA]{5}', args.amplicon_seq)
			if not poly_t_matched:
				info("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT is used as fake donor sequence.")
				fake_donor = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
			elif not poly_c_matched:
				info("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG is used as fake donor sequence.")
				fake_donor = 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
			elif not poly_g_matched:
				info("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG is used as fake donor sequence.")
				fake_donor = 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
			elif not poly_a_matched:
				info("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA is used as fake donor sequence.")
				fake_donor = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
			else:
				info('The amplicon_seq has poly-A|T|C|G. However, AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA is used as fake donor sequence.')
				fake_donor = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
			args.donor_seq = fake_donor
			args.expected_ki_amplicon_seq = amplicon_searched.group(1) + fake_donor + amplicon_searched.group(2)
			info(args.expected_ki_amplicon_seq + 'is used as fake knock-in amplicon sequence.')
		if not int(args.length_left_homologyarm) > 17:
			raise InputException('length_left_homologyarm should be >17bp')
		if not int(args.length_right_homologyarm) > 17:
			raise InputException('length_right_homologyarm should be >17bp')
		if not len(args.donor_seq) > 12:
			raise InputException('length_right_homologyarm should be >12bp')
		if args.force_knockout_analysis:
			if not is_knockout_analysis:
				info("MaChIAto runs knock-out analysis.")
			is_knockout_analysis = True

		if args.mode in "CRISPResso":
			# last version : 20190311142608
			# INPUT_DIRECTORY = "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190218_development_of_MaChAto/output/CRISPResso_on_atp5b-KI-C-genome_S16_L001_R1_001_atp5b-KI-C-genome_S16_L001_R2_001"
			# last version : 20190311193245
			# INPUT_DIRECTORY = "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190218_development_of_MaChAto/output/CRISPResso_on_atp5b-KI-A-genome_S14_L001_R1_001_atp5b-KI-A-genome_S14_L001_R2_001"
			# last version : 20190311210415
			# INPUT_DIRECTORY = "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190218_development_of_MaChAto/output/CRISPResso_on_atp5b-control-A-genome_S13_L001_R1_001_atp5b-control-A-genome_S13_L001_R2_001"
			# OPTIMIZATION_LOG_DIRECTORY = "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190218_development_of_MaChAto/test_log/"
			# filename = "Alleles_frequency_table.txt"

			info("Loading Alleles_frequency_table")
			pddf_alleles_freq = pd.read_csv(args.crispreeso_file, sep='\t')

			"""
			pddf_alleles_freq.columns ==>
			Index([u'Aligned_Sequence', u'Reference_Sequence', u'NHEJ', u'UNMODIFIED',
			u'HDR', u'n_deleted', u'n_inserted', u'n_mutated', u'#Reads',
			u'%Reads'],
			dtype='object')
			pddf_alleles_freq.dtypes ==>
			Aligned_Sequence       object
			Reference_Sequence     object
			NHEJ                     bool
			UNMODIFIED               bool
			HDR                      bool
			n_deleted             float64
			n_inserted            float64
			n_mutated               int64
			#Reads                  int64
			%Reads                float64
			dtype: object
			"""

		elif args.mode in "CRISPResso2":

			if not os.path.exists(args.crispreeso2_file):
				raise InputException(args.crispreeso2_file + ' dosen\'t exist.')

			info("Loading Alleles_frequency_table.zip")
			pddf_alleles_freq = pd.read_csv(args.crispreeso2_file, compression='zip', sep='\t')
			"""
			pddf_alleles_freq.columns ==>
			Index([u'Aligned_Sequence', u'Reference_Sequence', u'Reference_Name',
			u'Read_Status', u'n_deleted', u'n_inserted', u'n_mutated', u'#Reads',
			u'%Reads'],
			dtype='object')
			pddf_alleles_freq.dtypes ==>
			Aligned_Sequence       object # this is same as CRISPResso1
			Reference_Sequence     object # this is not used.
			Reference_Name         object # this is unique in CRISPResso2.
			Read_Status            object # this is unique in CRISPResso2.
			n_deleted               int64 this is same as CRISPResso1
			n_inserted              int64 this is same as CRISPResso1
			n_mutated               int64 this is same as CRISPResso1
			#Reads                  int64 this is same as CRISPResso1
			%Reads                float64 this is same as CRISPResso1
			dtype: object
			"""

			def convert_crispresso2_to_crispresso(row):
    			
				aligned_sequence,\
				reference_sequence,\
				is_nhej,\
				is_unmodified,\
				is_hdr,\
				n_deleted,\
				n_inserted,\
				n_mutated,\
				n_Reads,\
				r_Reads\
				= convert_crispresso2row_to_crispresso(row)

				return(pd.Series({u'Aligned_Sequence' : aligned_sequence
				, u'Reference_Sequence' : reference_sequence
				, u'NHEJ' : is_nhej
				, u'UNMODIFIED' : is_unmodified
				, u'HDR' : is_hdr
				, u'n_deleted' : n_deleted
				, u'n_inserted' : n_inserted
				, u'n_mutated' : n_mutated
				, u'#Reads' : n_Reads
				, u'%Reads' : r_Reads
				}))

			# modify table as table of CRISPResso1.
			tqdm.pandas(desc = "Converting CRISPResso2 data into CRISPResso1 format")
			pddf_alleles_freq = pd.concat([pddf_alleles_freq.progress_apply(convert_crispresso2_to_crispresso, axis=1).copy()])
			pddf_alleles_freq = pddf_alleles_freq[[u'Aligned_Sequence', u'Reference_Sequence', u'NHEJ', u'UNMODIFIED', u'HDR', u'n_deleted', u'n_inserted', u'n_mutated', u'#Reads', u'%Reads']] # ordering
		else:
			pass

		# It is general process below.
		# The dummy table is used if table is empty.
		null_pddf = pd.DataFrame({'Aligned_Sequence' : ["A"]
			, 'Reference_Sequence' : ["A"]
			, 'NHEJ' : [False]
			, 'UNMODIFIED' : [False]
			, 'HDR' : [False]
			, 'n_deleted' : [0]
			, 'n_inserted' : [0]
			, 'n_mutated' : [0]
			, '#Reads' : [0]
			, '%Reads' : [0]
		})

		# import pdb; pdb.set_trace()

		### Make output directory ##################################################################################
		if args.name:
			database_id = database_id + "_on_" + args.name
		OUTPUT_DIRECTORY='MaChIAto_from_' + args.mode + '_at_%s' % database_id
		
		if args.output_folder:
			OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
		else:
			raise InputException("output folder or directory is required.")
		
		_outputjp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename)
		log_filename=_outputjp('MaChIAto_from_' + args.mode + '_RUNNING_LOG.txt')

		try:
			os.makedirs(OUTPUT_DIRECTORY)
			info('Creating Folder %s' % OUTPUT_DIRECTORY)
			info('Done.')

		except:
			warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

		finally:
			logging.getLogger().addHandler(logging.FileHandler(log_filename)) # make log file
			with open(log_filename,'w+') as outfile:
				outfile.write('[Command used]:\nMaChIAto %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
		###########################################################################################################

		pddf_NHEJ_alleles_freq = pddf_alleles_freq.query('NHEJ == True')
		if pddf_NHEJ_alleles_freq.empty:
			pddf_NHEJ_alleles_freq = null_pddf.copy()
		pddf_UNMODIFIED_alleles_freq = pddf_alleles_freq.query('UNMODIFIED == True')
		if pddf_UNMODIFIED_alleles_freq.empty:
			pddf_UNMODIFIED_alleles_freq = null_pddf.copy()
		pddf_HDR_alleles_freq = pddf_alleles_freq.query('HDR == True')
		if pddf_HDR_alleles_freq.empty:
			pddf_HDR_alleles_freq = null_pddf.copy()
		pddf_MIXED_alleles_freq = pddf_alleles_freq.query('NHEJ == False and UNMODIFIED == False and HDR == False')
		if pddf_MIXED_alleles_freq.empty:
			pddf_MIXED_alleles_freq = null_pddf.copy()

		# defalt parameter (it is used for debug. you can ignore that.)
		amplicon_seq = args.amplicon_seq
		guide_seq = args.guide_seq
		donor_seq = args.donor_seq
		expected_ki_amplicon_seq = args.expected_ki_amplicon_seq
		left_homology_arm_len = args.length_left_homologyarm
		right_homology_arm_len = args.length_right_homologyarm
		init_left_junction_len = args.length_left_homologyarm
		init_right_junction_len = args.length_right_homologyarm
		init_left_out_indicator_len = 5
		init_right_out_indicator_len = 5
		init_left_in_indicator_len = 6
		init_right_in_indicator_len = 6
		init_outleft_align_threshold = 0.97772495
		init_outright_align_threshold = 0.97390871
		init_inleft_align_threshold = 0.5488135
		init_inright_align_threshold = 0.97861834
		
		info("Running Optimizations")

		test_dict = {"amplicon_seq" : amplicon_seq.strip().upper(),
			"guide_seq" : guide_seq.strip().upper(),
			"expected_ki_amplicon_seq" : expected_ki_amplicon_seq.strip().upper(),
			"donor_seq" : donor_seq.strip().upper(),
			"left_homology_arm_len" : left_homology_arm_len,
			"right_homology_arm_len" : right_homology_arm_len,
			"left_junction_len" : init_left_junction_len,
			"right_junction_len" : init_right_junction_len,
			"left_out_indicator_len" : init_left_out_indicator_len,
			"right_out_indicator_len" : init_right_out_indicator_len,
			"left_in_indicator_len" : init_left_in_indicator_len,
			"right_in_indicator_len" : init_right_in_indicator_len,
			"outleft_align_threshold" : init_outleft_align_threshold,
			"outright_align_threshold" : init_outright_align_threshold,
			"inleft_align_threshold" : init_inleft_align_threshold,
			"inright_align_threshold" : init_inright_align_threshold,
			"optimization_mode" : "OFF"}

		#########test################
		# mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq, test_dict.copy())
		# mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq.head(100), test_dict.copy())
		# nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq.head(9), test_dict.copy())
		# nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq, test_dict.copy())
		# hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq.head(10), test_dict.copy())
		# unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq, test_dict.copy())
		# unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq.head(11), test_dict.copy())
		# print unmodified_info.cnt_dict["has_outout_indicators"] + unmodified_info.cnt_dict["has_outout_rc_indicators"] - 1000 * unmodified_info.cnt_dict["has_outout_reversible_indicators"] - unmodified_info.cnt_dict["has_outout_no_indicators"]
		# rate_dict = calc_MaChIAto_rates_from_CRISPResso(mixed_info.read_cnt_dict)
		# import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
		
		#########################

		if not args.skip_optimization:
			### BayesianOptimization 1 ##################################################################################
			# optimization of both out-side indicators

			info("Preparing for first optimization ------------------------------------------------------------->")

			test_dict["optimization_mode"] = "outoutindicators"

			info("Making test data using NHEJ data")
			nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq, test_dict.copy())
			info("Making test data using UNMODIFIED data")
			unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq, test_dict.copy())
			info("Making test data using HDR data")
			hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq, test_dict.copy())
			info("Making test data using Mixed HDR-NHEJ data")
			#mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq.head(min(len(pddf_MIXED_alleles_freq), 1000)), test_dict.copy())
			mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq, test_dict.copy())

			info("Running first optimization : Optimizing Out-Out indicators------------------------------------>")

			# np.random.seed(0)

			def optimized_func_outoutindicators(x):
				left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len\
				= str(int(x[:,0])), str(int(x[:,1])), str(int(x[:,2])), str(int(x[:,3]))
				info("<First optimization>\n[parameters] left_out_indicator_len : " + left_out_indicator_len\
					+ " right_out_indicator_len : " + right_out_indicator_len\
					+ " left_junction_len : " + left_junction_len\
					+ " right_junction_len : " + right_junction_len)
				# debug
				# left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len = (10, 5, 40, 45)
				# run generator method
				nhej_info.generator_first_optimization(left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len).next()
				unmodified_info.generator_first_optimization(left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len).next()
				hdr_info.generator_first_optimization(left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len).next()
				mixed_info.generator_first_optimization(left_out_indicator_len, right_out_indicator_len, left_junction_len, right_junction_len).next()
				# debug
				# import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

				if nhej_info.cnt_dict["is_nonspecific_outout_indicators"] > 0 or\
					unmodified_info.cnt_dict["is_nonspecific_outout_indicators"] > 0 or\
					hdr_info.cnt_dict["is_nonspecific_outout_indicators"] > 0:
					# mixed_info.cnt_dict["is_nonspecific_outout_indicators"] > 0:
					# mixed HDR-NHEJ contains chimera-like sequences. so it is ignored in this filter.
					info("Score : 1")
					return 1.
				else:
					
					def calc_indicator_detection_rate(cnt_dict):
						return (float(cnt_dict["has_outout_any_indicators"]) / float(cnt_dict["total"]))

					_fi = lambda x: (x + 1 - 1 / (math.pow(math.e, 4 * math.pi * math.e * x - math.e) + 1)) / 2
					# import pdb; pdb.set_trace()
					if(hdr_info.cnt_dict["total"] > 2 or not is_knockout_analysis):
						score = _fi(1.0 - (calc_indicator_detection_rate(nhej_info.cnt_dict)\
							+ calc_indicator_detection_rate(unmodified_info.cnt_dict)\
							+ calc_indicator_detection_rate(hdr_info.cnt_dict)\
							+ calc_indicator_detection_rate(mixed_info.cnt_dict))/4)
					else:
						score = _fi(1.0 - (calc_indicator_detection_rate(nhej_info.cnt_dict)\
							+ calc_indicator_detection_rate(unmodified_info.cnt_dict))/2)

				info("Score : " + str(score))
				return score
			
			np.random.seed(0)
			bounds_outoutindicators = [{'name': 'left_out_indicator_len', 'type': 'discrete', 'domain': range(5, 25)},\
				{'name': 'right_out_indicator_len', 'type': 'discrete', 'domain': range(5, 25)},\
				{'name': 'left_junction_len', 'type': 'discrete', 'domain': range(test_dict["left_homology_arm_len"], test_dict["left_homology_arm_len"] + 10)},\
				{'name': 'right_junction_len', 'type': 'discrete', 'domain': range(test_dict["right_homology_arm_len"], test_dict["right_homology_arm_len"] + 10)}]
			Bopt_outoutindicators = BayesianOptimization(f=optimized_func_outoutindicators, domain=bounds_outoutindicators, initial_design_numdata=20)
			Bopt_outoutindicators.run_optimization(max_iter = 200,\
			report_file = _outputjp("first_opt_report.txt"),\
			evaluations_file = _outputjp("first_opt_evaluation.txt"))
			info("Minimum Score : " + str(Bopt_outoutindicators.fx_opt))
			test_dict["left_out_indicator_len"] = str(int(Bopt_outoutindicators.x_opt[0]))
			test_dict["right_out_indicator_len"] = str(int(Bopt_outoutindicators.x_opt[1]))
			test_dict["left_junction_len"] = str(int(Bopt_outoutindicators.x_opt[2]))
			test_dict["right_junction_len"] = str(int(Bopt_outoutindicators.x_opt[3]))
			### BayesianOptimization 2 ##################################################################################
			# optimization of out-out threshold

			info("Preparing for second optimization ------------------------------------------------------------>")

			test_dict["optimization_mode"] = "outout_align_thresholds"
			
			info("Making test data using NHEJ data")
			# nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq, test_dict.copy())
			nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq.head(min(len(pddf_NHEJ_alleles_freq), 100)), test_dict.copy())
			info("Making test data using UNMODIFIED data")
			# unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq, test_dict.copy())
			unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq.head(min(len(pddf_UNMODIFIED_alleles_freq), 100)), test_dict.copy())
			info("Making test data using HDR data")
			# hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq, test_dict.copy())
			hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq.head(min(len(pddf_HDR_alleles_freq), 100)), test_dict.copy())
			info("Making test data using Mixed HDR-NHEJ data")
			# mixed_info = MaChIAtoClassifier(pddf_HDR_alleles_freq, test_dict.copy())
			mixed_info = MaChIAtoClassifier(pddf_HDR_alleles_freq.head(min(len(pddf_MIXED_alleles_freq), 100)), test_dict.copy())
			
			'''
			# debug
			nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq.head(1000), test_dict.copy())
			unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq.head(1000), test_dict.copy())
			hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq.head(1000), test_dict.copy())
			mixed_info = MaChIAtoClassifier(pddf_HDR_alleles_freq.head(1000), test_dict.copy())
			#mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq.head(min(len(pddf_MIXED_alleles_freq), 1000)), test_dict.copy())
			'''

			info("Running second optimization : Optimizing Out-Out alignment thresholds ------------------------>")

			def optimized_func_outout_align_thresholds(x):
				outleft_align_threshold, outright_align_threshold = float(x[:,0]), float(x[:,1])
				info("<Second optimization>\n[parameters] outleft_align_threshold : " + str(outleft_align_threshold)\
					+ " outright_align_threshold : " + str(outright_align_threshold))
				
				def calc_untreated_detection_rate(cnt_dict):
					if float(cnt_dict["has_outout_any_indicators"]) == 0.:
						return 0.
					else:
						return (float(cnt_dict["is_untreated"]) / float(cnt_dict["has_outout_any_indicators"]))

				def calc_perfect_outout_ki_detection_rate(cnt_dict):
					if float(cnt_dict["is_treated"]) == 0.:
						return 0.
					else:
						return (float(cnt_dict["is_perfect_outout_ki"]) / float(cnt_dict["is_treated"]))
				
				_fi = lambda x: (x + 1 - 1 / (math.pow(math.e, 4 * math.pi * math.e * x - math.e) + 1)) / 2

				# np.random.seed(0)
				if(hdr_info.cnt_dict["total"] > 2 or not is_knockout_analysis):
					# run generator method
					nhej_info.generator_second_optimization(outleft_align_threshold, outright_align_threshold).next()
					unmodified_info.generator_second_optimization(outleft_align_threshold, outright_align_threshold).next()
					hdr_info.generator_second_optimization(outleft_align_threshold, outright_align_threshold).next()
					mixed_info.generator_second_optimization(outleft_align_threshold, outright_align_threshold).next()
					score = _fi(1.0 \
						- (0.5 \
						+ (calc_untreated_detection_rate(unmodified_info.cnt_dict)
							+ calc_perfect_outout_ki_detection_rate(hdr_info.cnt_dict)) / 4 \
						- (calc_perfect_outout_ki_detection_rate(mixed_info.cnt_dict) \
							+ calc_untreated_detection_rate(nhej_info.cnt_dict)) / 4))
				else:
					# run generator method
					nhej_info.generator_second_optimization(outleft_align_threshold, outright_align_threshold).next()
					unmodified_info.generator_second_optimization(outleft_align_threshold, outright_align_threshold).next()
					score = _fi(1.0 \
						- (0.5 \
						+ (calc_untreated_detection_rate(unmodified_info.cnt_dict)) / 2 \
						- (calc_untreated_detection_rate(nhej_info.cnt_dict)) / 2))

				info("Score : " + str(score))
				return score

			np.random.seed(0)
			bounds_outout_align_thresholds = [{'name': 'outleft_align_threshold', 'type': 'continuous', 'domain': (0.0,1.0)},
				{'name': 'outright_align_threshold', 'type': 'continuous', 'domain': (0.0,1.0)}]
			Bopt_outout_align_thresholds = BayesianOptimization(f=optimized_func_outout_align_thresholds, domain=bounds_outout_align_thresholds, initial_design_numdata=20)
			Bopt_outout_align_thresholds.run_optimization(max_iter = 100,\
				report_file = _outputjp("second_opt_report.txt"),\
				evaluations_file = _outputjp("second_opt_evaluation.txt"))
			info("Minimum Score : " + str(Bopt_outout_align_thresholds.fx_opt))
			test_dict["outleft_align_threshold"] = float(Bopt_outout_align_thresholds.x_opt[0])
			test_dict["outright_align_threshold"] = float(Bopt_outout_align_thresholds.x_opt[1])
			### BayesianOptimization 3 ##################################################################################
			# optimization of in-out indicator

			info("Preparing for third optimization ------------------------------------------------------------->")

			test_dict["optimization_mode"] = "inoutindicators"

			info("Making test data using NHEJ data")
			nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq, test_dict.copy())
			info("Making test data using UNMODIFIED data")
			unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq, test_dict.copy())
			info("Making test data using HDR data")
			hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq, test_dict.copy())
			info("Making test data using Mixed HDR-NHEJ data")
			mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq, test_dict.copy())

			info("Running third optimization : Optimizing In-Out indicators ------------------------------------>")

			# np.random.seed(0)

			def optimized_func_inoutindicators(x):
				left_in_indicator_len, right_in_indicator_len =\
					str(int(x[:,0])), str(int(x[:,1]))
				info("<Third optimization>\n[parameters] left_in_indicator_len : " + str(left_in_indicator_len)\
					+ " right_in_indicator_len : " + str(right_in_indicator_len))
				nhej_info.generator_third_optimization(left_in_indicator_len, right_in_indicator_len).next()
				unmodified_info.generator_third_optimization(left_in_indicator_len, right_in_indicator_len).next()
				hdr_info.generator_third_optimization(left_in_indicator_len, right_in_indicator_len).next()
				mixed_info.generator_third_optimization(left_in_indicator_len, right_in_indicator_len).next()

				_sum_2_cnt_dict = lambda cnt_dict: (cnt_dict["has_leftinout_indicators"] + cnt_dict["has_rightinout_indicators"])
				_sum_multi_cnt_dict = lambda cnt_dict: (cnt_dict["has_multi_leftinout_indicators"] + cnt_dict["has_multi_rightinout_indicators"])

				if (_sum_2_cnt_dict(unmodified_info.cnt_dict) > 0) \
					or (_sum_multi_cnt_dict(hdr_info.cnt_dict) > 0) \
					or (_sum_multi_cnt_dict(mixed_info.cnt_dict) > 0):
					
					info("Score : 1")
					return 1.

				else:

					def calc_sum_both_indicator_detection_rate(cnt_dict):
						if float(cnt_dict["is_treated"]) == 0.:
								return 0.
						else:
							left_detection_rate = float(cnt_dict["has_leftinout_indicators"]) / float(cnt_dict["is_treated"])
							right_detection_rate = float(cnt_dict["has_rightinout_indicators"]) / float(cnt_dict["is_treated"])
							return ((left_detection_rate + right_detection_rate)/ 2)
					
					'''
					score = _fi(1.0 - ((calc_sum_both_indicator_detection_rate(hdr_info.cnt_dict)\
							+ calc_sum_both_indicator_detection_rate(mixed_info.cnt_dict)) / 2))
					'''
					
					_fi = lambda x: (x + 1 - 1 / (math.pow(math.e, 4 * math.pi * math.e * x - math.e) + 1)) / 2
					score = _fi(1.0 \
					- (0.5 \
					+ ((calc_sum_both_indicator_detection_rate(hdr_info.cnt_dict)\
						+ calc_sum_both_indicator_detection_rate(mixed_info.cnt_dict)) / 4) \
					- (calc_sum_both_indicator_detection_rate(nhej_info.cnt_dict)) / 2))

					info("Score : " + str(score))
					return score

			np.random.seed(0)
			if(hdr_info.cnt_dict["total"] > 2 or not is_knockout_analysis):
				bounds_inoutindicators = [{'name': 'left_in_indicator_len', 'type': 'discrete', 'domain': range(5, int(len(donor_seq)/2))},
				{'name': 'right_in_indicator_len', 'type': 'discrete', 'domain': range(5, int(len(donor_seq)/2 - 1))}]
				Bopt_inoutindicators = BayesianOptimization(f=optimized_func_inoutindicators, domain=bounds_inoutindicators, initial_design_numdata=20)
				Bopt_inoutindicators.run_optimization(max_iter = 200,\
					report_file = _outputjp("third_opt_report.txt"),\
					evaluations_file = _outputjp("third_opt_evaluation.txt"))
				info("Minimum Score : " + str(Bopt_inoutindicators.fx_opt))
				test_dict["left_in_indicator_len"] = str(int(Bopt_inoutindicators.x_opt[0]))
				test_dict["right_in_indicator_len"] = str(int(Bopt_inoutindicators.x_opt[1]))
			else:
				info("This is not knock-in sample : left_in_indicator_len " + str(int(len(donor_seq)/2 - 1)) + " right_in_indicator_len : " + str(int(len(donor_seq)/2 - 1)))
				test_dict["left_in_indicator_len"] = str(int(len(donor_seq)/2 - 1))
				test_dict["right_in_indicator_len"] = str(int(len(donor_seq)/2 - 1))
			

			### BayesianOptimization 4 ##################################################################################
			# optimization of in-out threshold

			info("Preparing for last optimization -------------------------------------------------------------->")

			test_dict["optimization_mode"] = "inout_align_thresholds"

			info("Making test data using HDR data")
			hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq, test_dict.copy())
			info("Making test data using Mixed HDR-NHEJ data")
			mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq, test_dict.copy())
			
			info("Running last optimization : Optimizing In-Out alignment thresholds --------------------------->")

			# np.random.seed(0)

			def optimized_func_inout_align_thresholds(x):
				inleft_align_threshold, inright_align_threshold =\
					float(x[:,0]), float(x[:,1])
				info("<Last optimization>\n[parameters] inleft_align_threshold : " + str(inleft_align_threshold)\
					+ " inright_align_threshold : " + str(inright_align_threshold))

				hdr_info.generator_forth_optimization(inleft_align_threshold, inright_align_threshold).next()
				mixed_info.generator_forth_optimization(inleft_align_threshold, inright_align_threshold).next()

				def calc_nearly_left_inout_ki_detection_rate(cnt_dict):
					if float(cnt_dict["has_leftinout_indicators"]) == 0.:
						return 0.
					else:
						return (float(cnt_dict["is_nearly_leftinout_ki"]) / float(cnt_dict["has_leftinout_indicators"]))

				def calc_nearly_right_inout_ki_detection_rate(cnt_dict):
					if float(cnt_dict["has_rightinout_indicators"]) == 0.:
						return 0.
					else:
						return (float(cnt_dict["is_nearly_rightinout_ki"]) / float(cnt_dict["has_rightinout_indicators"]))
				
				def calc_sum_inout_nearly_kis_detection_rate(cnt_dict):
					return (calc_nearly_left_inout_ki_detection_rate(cnt_dict) + calc_nearly_right_inout_ki_detection_rate(cnt_dict)) / 2

				_fi = lambda x: (x + 1 - 1 / (math.pow(math.e, 4 * math.pi * math.e * x - math.e) + 1)) / 2
				score = _fi(1.0 \
					- (0.5 \
					+ (calc_sum_inout_nearly_kis_detection_rate(hdr_info.cnt_dict)) / 2 \
					- (calc_sum_inout_nearly_kis_detection_rate(mixed_info.cnt_dict)) / 2 \
					))

				info("Score : " + str(score))
				return score

			np.random.seed(0)
			if(hdr_info.cnt_dict["total"] > 2 or not is_knockout_analysis):
				bounds_inout_align_thresholds = [{'name': 'inleft_align_threshold', 'type': 'continuous', 'domain': (0.0,1.0)},
				{'name': 'inright_align_threshold', 'type': 'continuous', 'domain': (0.0,1.0)}]
				Bopt_inout_align_thresholds = BayesianOptimization(f=optimized_func_inout_align_thresholds, domain=bounds_inout_align_thresholds, initial_design_numdata=20)
				Bopt_inout_align_thresholds.run_optimization(max_iter = 50,\
					report_file = _outputjp("forth_opt_report.txt"),\
					evaluations_file = _outputjp("forth_opt_evaluation.txt"))
				info("Minimum Score : " + str(Bopt_inout_align_thresholds.fx_opt))
				test_dict["inleft_align_threshold"] = float(Bopt_inout_align_thresholds.x_opt[0])
				test_dict["inright_align_threshold"] = float(Bopt_inout_align_thresholds.x_opt[1])
			else:
				info("This is not knock-in sample : inleft_align_threshold " + str(1.0) + " inright_align_threshold : " + str(1.0))
				test_dict["inleft_align_threshold"] = str(1.0)
				test_dict["inright_align_threshold"] = str(1.0)


			#####################################################################################

			info("Compliting optimization")

		info("Starting whole analysis")

		### Whole analysis ##################################################################################
		test_dict["optimization_mode"] = "OFF"

		opt_dict = test_dict.copy()

		info("Analyzing NHEJ data -------------------------------------------------------------------------->")
		opt_nhej_info = MaChIAtoClassifier(pddf_NHEJ_alleles_freq, opt_dict.copy())
		info("Analyzing UNMODIFIED data -------------------------------------------------------------------->")
		opt_unmodified_info = MaChIAtoClassifier(pddf_UNMODIFIED_alleles_freq, opt_dict.copy())
		info("Analyzing HDR data --------------------------------------------------------------------------->")
		opt_hdr_info = MaChIAtoClassifier(pddf_HDR_alleles_freq, opt_dict.copy())
		info("Analyzing Mixed HDR-NHEJ data ---------------------------------------------------------------->")
		opt_mixed_info = MaChIAtoClassifier(pddf_MIXED_alleles_freq, opt_dict.copy())

		info("Analyzing All data --------------------------------------------------------------------------->")

		### Add label ##################################################################################
		def add_CRISPResso_reclassification_labels(pddf_row):
		
			global_outout_classification_label = pddf_row['global_outout_classification_label']
			
			if global_outout_classification_label in "PRECISE_KNOCK_IN":

				return pd.Series({'CRISPResso_reclassification_labels' : "HDR"})

			elif global_outout_classification_label in "IMPRECISE_KNOCK_IN":

				return pd.Series({'CRISPResso_reclassification_labels' : "Mixed HDR-NHEJ"})

			elif is_in(global_outout_classification_label, ("INSERT_EDITING", "DELETION_EDITING", "SUBSTITUTION_EDITING")):

				return pd.Series({'CRISPResso_reclassification_labels' : "NHEJ"})

			elif global_outout_classification_label in "UNTREATED":
				
				return pd.Series({'CRISPResso_reclassification_labels' : "Unmodified"})

			elif global_outout_classification_label in "OTHERS":
				
				return pd.Series({'CRISPResso_reclassification_labels' : "Unclassified"})

			else:
				# import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
				error('Unexpected categorizing occurs.')
				sys.exit(1)

		tqdm.pandas(desc = "Adding CRISPResso_reclassification_labels to NHEJ data")
		opt_nhej_info.out_pddf = opt_nhej_info.out_pddf.join(opt_nhej_info.out_pddf.progress_apply(add_CRISPResso_reclassification_labels, axis = 1)).copy()
		tqdm.pandas(desc = "Adding CRISPResso_reclassification_labels to Unmodified data")
		opt_unmodified_info.out_pddf = opt_unmodified_info.out_pddf.join(opt_unmodified_info.out_pddf.progress_apply(add_CRISPResso_reclassification_labels, axis = 1)).copy()
		tqdm.pandas(desc = "Adding CRISPResso_reclassification_labels to HDR data")
		opt_hdr_info.out_pddf = opt_hdr_info.out_pddf.join(opt_hdr_info.out_pddf.progress_apply(add_CRISPResso_reclassification_labels, axis = 1)).copy()
		tqdm.pandas(desc = "Adding CRISPResso_reclassification_labels to Mixed HDR-NHEJ data")
		opt_mixed_info.out_pddf = opt_mixed_info.out_pddf.join(opt_mixed_info.out_pddf.progress_apply(add_CRISPResso_reclassification_labels, axis = 1)).copy()
		
		### Add count ##################################################################################
		def register_count(MaChIAtoClassifier, label, pddf):
			MaChIAtoClassifier.cnt_dict[label] = len(pddf)
			MaChIAtoClassifier.read_cnt_dict[label] = pddf["#Reads"].sum()

		register_count(opt_nhej_info, "NHEJ", opt_nhej_info.out_pddf.query('CRISPResso_reclassification_labels == "NHEJ"'))
		register_count(opt_unmodified_info, "Unmodified", opt_unmodified_info.out_pddf.query('CRISPResso_reclassification_labels == "Unmodified"'))
		register_count(opt_hdr_info, "HDR", opt_hdr_info.out_pddf.query('CRISPResso_reclassification_labels == "HDR"'))
		register_count(opt_mixed_info, "Mixed HDR-NHEJ", opt_mixed_info.out_pddf.query('CRISPResso_reclassification_labels == "Mixed HDR-NHEJ"'))
		
		### Calculate consistency ##################################################################################
		def show_persentage(MaChIAtoClassifier, label):
			if float(MaChIAtoClassifier.read_cnt_dict["total"]) == 0.:
				return 100.
			else:
				return float(MaChIAtoClassifier.read_cnt_dict[label]) \
					/ float(MaChIAtoClassifier.read_cnt_dict["total"]) * 100.0
			
		
		consistency_dict = {
			"NHEJ" : show_persentage(opt_nhej_info, "NHEJ"),
			"Unmodified" : show_persentage(opt_unmodified_info, "Unmodified"),
			"HDR" : show_persentage(opt_hdr_info, "HDR"),
			"Mixed HDR-NHEJ" : show_persentage(opt_mixed_info, "Mixed HDR-NHEJ")
		}
		info('Degree of consistency in NHEJ is %s percent.' % consistency_dict["NHEJ"])
		info('Degree of consistency in Unmodified is %s percent.' % consistency_dict["Unmodified"])
		info('Degree of consistency in HDR is %s percent.' % consistency_dict["HDR"])
		info('Degree of consistency in Mixed HDR-NHEJ is %s percent.' % consistency_dict["Mixed HDR-NHEJ"])
		# import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

		### Merge Dataframes ##################################################################################
		merge_pddf = pd.merge(opt_nhej_info.out_pddf, opt_unmodified_info.out_pddf, how="outer").copy()
		merge_pddf = pd.merge(merge_pddf, opt_hdr_info.out_pddf, how="outer").copy()
		merge_pddf = pd.merge(merge_pddf, opt_mixed_info.out_pddf, how="outer").copy()

		### Sum countdata ##################################################################################
		sum_four_count_dict = lambda a, b, c, d: dict((k, a[k] + b[k] + c[k] + d[k]) for k in list(set(a.keys()) & set(b.keys()) & set(c.keys()) & set(d.keys())))
		cnt_dict_sum = sum_four_count_dict(opt_unmodified_info.cnt_dict, opt_nhej_info.cnt_dict, opt_hdr_info.cnt_dict, opt_mixed_info.cnt_dict)
		read_cnt_dict_sum = sum_four_count_dict(opt_unmodified_info.read_cnt_dict, opt_nhej_info.read_cnt_dict, opt_hdr_info.read_cnt_dict, opt_mixed_info.read_cnt_dict)
		
		### Calculate rates ##################################################################################

		rate_dict = calc_MaChIAto_rates_from_CRISPResso(read_cnt_dict_sum)

		### Save files ##################################################################################
		info("Save files...")
		opt_nhej_info.out_pddf.to_csv(path_or_buf=_outputjp("CRISPesso_NHEJ_dataframe.csv"), sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_nhej_info.cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_NHEJ_alleletype_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_nhej_info.read_cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_NHEJ_read_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		opt_unmodified_info.out_pddf.to_csv(path_or_buf=_outputjp("CRISPesso_UNMODIFIED_dataframe.csv"), sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_unmodified_info.cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_UNMODIFIED_alleletype_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_unmodified_info.read_cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_UNMODIFIED_read_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		opt_hdr_info.out_pddf.to_csv(path_or_buf=_outputjp("CRISPesso_HDR_dataframe.csv"), sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_hdr_info.cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_HDR_alleletype_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_hdr_info.read_cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_HDR_read_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		opt_mixed_info.out_pddf.to_csv(path_or_buf=_outputjp("CRISPesso_Mixed_HDR_NHEJ_dataframe.csv"), sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_mixed_info.cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_Mixed_HDR_NHEJ_alleletype_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([opt_mixed_info.read_cnt_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_Mixed_HDR_NHEJ_read_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")

		pd.DataFrame.from_dict([consistency_dict]).to_csv(path_or_buf=_outputjp("CRISPesso_MaChiIAto_consistency_persentage.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")

		pd.DataFrame.from_dict([opt_unmodified_info.seq_info_dict]).to_csv(path_or_buf=_outputjp("MaChiIAto_optimized_param.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		special_seq_dict = {"Out-Left indicator sequence" : opt_unmodified_info.left_out_indicator_seq
			,"Leftside junction" : opt_unmodified_info.left_junction_seq
			,"Rightside junction" : opt_unmodified_info.right_junction_seq
			,"Out-Right indicator sequence" : opt_unmodified_info.right_out_indicator_seq
			,"Leftside homology arm" : opt_unmodified_info.left_homologyarm_seq
			,"Rightside homology arm" : opt_unmodified_info.right_homologyarm_seq
			,"Untreated sequence" : opt_unmodified_info.untreated_seq
			,"In-Left indicator sequence" : opt_unmodified_info.left_in_indicator_seq
			,"In-Right indicator sequence" :opt_unmodified_info.right_in_indicator_seq
			,"Left cutout sequence" :opt_unmodified_info.left_cutout_seq
			,"Right cutout sequence" :opt_unmodified_info.right_cutout_seq
			,"Left cutout maximum deletion sequence" :opt_unmodified_info.left_cutout_maxdel_seq
			,"Right cutout maximum deletion sequence" :opt_unmodified_info.right_cutout_maxdel_seq}
		pd.DataFrame.from_dict([special_seq_dict]).to_csv(path_or_buf=_outputjp("MaChiIAto_specific_sequences.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")

		merge_pddf.to_csv(path_or_buf=_outputjp("ALL_dataframe.csv"), sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([cnt_dict_sum]).to_csv(path_or_buf=_outputjp("ALL_alleletype_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")
		pd.DataFrame.from_dict([read_cnt_dict_sum]).to_csv(path_or_buf=_outputjp("ALL_read_countdata.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")

		pd.DataFrame.from_dict([rate_dict]).to_csv(path_or_buf=_outputjp("ALL_read_rates.csv"), header=True, index=False, mode='a', sep=',', encoding="utf-8")



		info("All Done!")
		print """
───────────────────────────────────────────────
	　 　  （
	　　　　`.
	　　　､　｀ヽ. 　｀ヽ
	　　　ﾉ　 　 ,!　　ﾉ
	　  ,r　　　'´　_　　 '　 _,
	　 !　　　　 r'´　　 r‐'
	　  　, 　- ―― -　、
	　　 !ゝ;;;;;;;;;;;ノiー.、
	　　 | 　 ￣￣￣ 　.|'⌒i}
	　 _.! 　 　 　 　 　 |/／
	　,´.ゝ ＿___＿_ .ノ'´｀
	　ゝ　＿＿￣￣￣＿＿　ノ
	　　｀ ￣￣￣￣￣￣￣ ´　　　
───────────────────────────────────────────────
		"""
		sys.exit(0)
	except Exception as e:
		# import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
		error('Unexpected error, please check your input.\n\nERROR: %s' % e)
		sys.exit(-1)
	except InputException as e:
		#traceback.print_exc(file = sys.stdout)
		error('Input error.\n\nERROR: %s\n\nPlease check MaChIAto manual:XXX.' % e)
		sys.exit(1)
	except UnexpectedException as e:
		#traceback.print_exc(file = sys.stdout)
		error('Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s' % e)
		sys.exit(1)
	