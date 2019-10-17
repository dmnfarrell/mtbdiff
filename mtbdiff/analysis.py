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
    Returns:
        list of labels of the genomes run
    """

    filenames = glob.glob(path+'/*.f*a')
    if len(filenames) == 0:
        print ('no fasta files found in %s' %path)
        return
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

def get_summary(struct):
    """Get summary of common variants for all samples"""
    
    S = struct.groupby(['start','end','gene','RD','region_type','Name'],as_index=False).agg({'ID':np.size,'length':np.mean})
    S = S.rename(columns={'ID':'freq'})
    S = S.sort_values('freq',ascending=False)
    return S

def run_RD_checker(rds):
    """Check for presence of RD"""

    X = pd.pivot_table(rds,index='RD_name',columns=['species'],values='Start')
    X[X.notnull()] = 1
    X = X.fillna(0)
    return X

def plot_RD(df, width=14, row_colors=None, **kwargs):
    """Plot sites matrix as clustermap"""

    h=len(df)/8+6
    g=sns.clustermap(df,figsize=(width,h),linecolor='gray',cmap='gray_r',xticklabels=True,
                      yticklabels=True, row_colors=row_colors, **kwargs)
    return g
