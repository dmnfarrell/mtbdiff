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
import seaborn as sns
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from . import utils

module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
mtb_ref = os.path.join(datadir, 'MTB-H37Rv.fna')
mtb_gff = os.path.join(datadir, 'MTB-H37Rv.gff')

def run_genomes(path, outpath='results', ref=None):
    """Run multiple genome files in path.

    Args:
        path: folder with input genome files
        outpath: output folder
        ref: reference genome fasta file, default is MTB-H37Rv
    """

    filenames = glob.glob(path+'/*.f*a')
    if ref is None:
        ref = mtb_ref
    if not os.path.exists(ref):
        print ('no such file %s' %ref)
        return
    names = []
    for f in filenames:
        n = os.path.splitext(os.path.basename(f))[0]
        names.append(n)
        print (f, n)
        utils.run_nucdiff(ref, f, outpath)
    return names

def run_RD_checker(rds):
    """Check for presence of RD"""

    X = pd.pivot_table(rds,index='RD_name',columns=['species'],values='Start')
    X[X.notnull()] = 1
    X = X.fillna(0)
    return X

def plot_RD(df, width=12, row_colors=None):

    h=len(df)/8+5
    g=sns.clustermap(df,figsize=(width,h),lw=.2,linecolor='gray',cmap='gray_r',
                      yticklabels=True, row_colors=row_colors)# cbar=False)
    return g

def sites_matrix(struct, index=['start','end'], freq=0):
    """Pivot by start site"""

    X = pd.pivot_table(struct,index=index,columns=['label'],values='Name',aggfunc='first')
    X[X.notnull()] = 1
    X = X.fillna(0)
    #remove unique?
    X = X[X.sum(1)>freq]
    return X
