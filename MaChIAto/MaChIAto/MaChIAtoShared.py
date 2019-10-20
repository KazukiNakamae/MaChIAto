# -*- coding: utf8 -*-
"""
MaChIAto - Kazuki Nakamae 2019
Software pipeline for the analysis of outcomes of MMEJ-assisted gene knock-in using CRISPR-Cas9 from deep sequencing data
https://github.com/Kazuki-Nakamae/MaChIAto
"""

### MODULES ############################
import argparse
from Bio import pairwise2
import csv
import datetime
from GPyOpt.methods import BayesianOptimization
import logging
import pandas as pd
import numpy as np
import math
import os
# import re
import regex as re
import sys
from tqdm import tqdm

### MODULE SETTING ############################

# logging
logging.basicConfig(level=logging.INFO,
	format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
	datefmt='%a, %d %b %Y %H:%M:%S',
	stream=sys.stderr,
	filemode="w"
	)
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info

### EXCEPTIONS ############################
class ReferenceRegexException(Exception):
    pass
class UnexpectedException(Exception):
  	pass
class InputException(Exception):
  	pass
