#!/usr/bin/env python
"""
    Analysis routines for mtbdiff
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

import os, sys, io, random, subprocess, glob
import string
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from . import utils

module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
mtbref = os.path.join(datadir, 'MTB-H37Rv.fna')

def run_genomes(path):
    """Run multiple genome files in path"""

    filenames = glob.glob(path+'/*.f*a')
    for f in filenames:
        n = os.path.basename(f)
        print (f, n)
        utils.run_nucdiff(mtbref,f)

    return