#!/usr/bin/env python
"""
    Command line app for mtbdiff
    Created July 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import os,sys,subprocess
import numpy as np
import pandas as pd
from . import utils, analysis

def run_tests():
    """test run"""

    path = 'test'
    names = analysis.run_genomes(path)
    struct, snp = utils.get_nucdiff_results('results', names)
    struct.to_csv('results/ref_struct.csv')
    print()
    print ('structural differences')
    print (struct.groupby('species').agg({'start':np.size}))
    rds = utils.find_regions(struct)
    analysis.run_RD_checker(rds)
    return

def main():
    "Run the application"

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Do tests")
    opts, remainder = parser.parse_args()
    if opts.test == True:
        run_tests()

if __name__ == '__main__':
    main()