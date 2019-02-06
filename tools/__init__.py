#!/usr/bin/env python3

import sys

# Check the Python version. Be sure there is at least 3.4 available!
# Ideally, do use the latest version of Intel Pythin 3 distribution.
    
version = sys.version_info
    
# We require at least Python 3.4 for the successful execution of the code:

ver_req_maj = 3 # Major Python version number -> 3
ver_req_min = 4 # Minor Python version number -> 4

if version.major < ver_req_maj or (version.major == ver_req_maj and version.minor < ver_req_min):

     print('\nRunning this script requires Python version 3.4 or higher.' +\
           'Ideally, use the latest Intel Python 3 distribution.\n')

     sys.exit(1)


