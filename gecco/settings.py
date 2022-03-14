"""
Paths and files which are accessed repeatedly throughout the project
"""

# Paths
import os
SRC_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SRC_DIR, os.pardir))
DEFAULT_PARAMETERS = os.path.join(SRC_DIR,'default_param.csv')
DEFAULT_V_PARAM_DIR = os.path.join(SRC_DIR, 'visualization_param')
