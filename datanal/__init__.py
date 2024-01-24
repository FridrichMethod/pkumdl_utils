"""A package for data analysis

datanal.utils: A module for general data analysis

datanal.enzyme: A module for enzyme data analysis

datanal.spr: A module for SPR data analysis
"""

__version__ = "1.0.1"
__author__ = "Zhaoyang Li"
__email__ = "zhaoyangli@stanford.edu"

import os

MODULE_PATH = os.path.abspath(__file__)
MODULE_DIR = os.path.dirname(MODULE_PATH)
CONFIG_PATH = os.path.join(MODULE_DIR, 'config/config.json')
