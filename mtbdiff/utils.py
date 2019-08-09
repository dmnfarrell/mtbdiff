#!/usr/bin/env python
"""
    Implements utilities for mtbdiff
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

import os, sys, io, random, subprocess
import string
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
RD = pd.read_csv('RD.csv',comment='#')

def features_to_dataframe(features, cds=False):
    """Get features from a biopython seq record object into a dataframe
    Args:
        features: bio seqfeatures
       returns: a dataframe with a row for each cds/entry.
      """

    allfeat = []
    for (item, f) in enumerate(features):
        x = f.__dict__
        qual = f.qualifiers
        #print(qual)
        x.update(qual)
        d = {}
        d['start'] = f.location.start
        d['end'] = f.location.end
        d['strand'] = f.location.strand
        for i in qual:
            if i in x:
                if type(x[i]) is list:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)
    cols = list(qual.keys())+['start','end']
    df = pd.DataFrame(allfeat,columns=cols)
    return df

def gff_to_dataframe(filename):
    feats = utils.gff_to_features(filename)
    return features_to_dataframe(feats)

def run_nucdiff(ref, query, outpath='results'):
    """Run nucfdiff"""

    r = os.path.splitext(os.path.basename(ref))[0]
    q = os.path.splitext(os.path.basename(query))[0]
    out = outpath + f'/{r}_{q}'
    if not os.path.exists(out):
        cmd = f'nucdiff {ref} {query} {out} query'
        print (cmd)
        subprocess.check_output(cmd,shell=True)
    else:
        print ('folder already present')
    return

def find_regions(result):
    """Find known regions overlap in results from nucdiff"""

    #result = result[result.Name.isin(['insertion','deletion'])]

    regions = ['RD'+str(i) for i in range(1,15)]
    RD = RD[RD.RD_name.isin(regions)]
    found=[]
    for i,r in list(result.iterrows()):
        df = RD[ (abs(RD.Start-r.start)<800) | (abs(RD.Stop-r.end)<800)].copy()
        #df = RD[(RD.Start>=r.start) & (RD.Stop<=r.end) ]
        if len(df)>0:
            idcol = r.keys()[0]
            df['name'] = r[idcol]
            df['species'] = r.species
            df['start'] = r.start
            found.append(df)
    found = pd.concat(found)
    return found
