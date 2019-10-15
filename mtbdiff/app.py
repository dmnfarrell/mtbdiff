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

def run(path, outpath, ref=None):
    """Run workflow for mtbdiff"""
    
    print ('running genomes in %s' %path)
    names = analysis.run_genomes(path, outpath=outpath, ref=ref)
    print ('reading results..')
    struct, snp = utils.get_nucdiff_results(outpath, names)    
    print('-------------------')
    print ('found results for %s genomes' %len(names))
    print ('finding overlapping RD regions..')
    struct['RD'] = struct.apply(utils.get_region,1)
    mtb_feat = utils.get_mtb_features()
    print ('finding overlapping annotations..')
    struct['gene'] = struct.apply(lambda x: utils.get_overlapping_annotations(x, mtb_feat), 1)
    struct['region_type'] = struct.apply(utils.get_region_type,1)
    struct.to_csv(os.path.join(outpath,'ref_struct.csv'))
    print('-------------------')
    #print (struct.region_type.value_counts())
        
    S = analysis.get_summary(struct)
    print('most frequent variants:')
    print (S[:20])

    X = pd.pivot_table(S,index=['Name'],values=['freq'],columns=['region_type'],aggfunc='count')
    print('-------------------')
    print('summary:')
    print (X[:10])
    
    S.to_csv(os.path.join(outpath,'summary.csv'))
    rdmat = utils.RD_matrix(struct)
    smat = utils.sites_matrix(struct, freq=2)
        
    #plot
    fig = analysis.plot_RD(rdmat.T)
    fig.savefig(os.path.join(outpath,'RD_matrix.png'),dpi=150)    
    print ('results saved to %s' %outpath)
    return

def fetch_test_data():
    """Download test genome assemblies"""
    
    ids=['GCA_003431725','GCA_003431765','GCA_003431735','GCA_003431775']
    utils.get_assemblies(ids=ids, path='test_genomes')
    return

def run_test():
    """Run test folder"""
    run('test_genomes', outpath='test_results')
    return

def main():
    "Run the application"

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", 
                        default=None, help="Input folder")
    parser.add_option("-o", "--output", dest="outpath", 
                        default='mtbdiff_results', help="Output folder")  
    parser.add_option("-r", "--ref", dest="ref", 
                        default=None, help="Reference fasta file")    
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Do tests")
    
    opts, remainder = parser.parse_args()
    if opts.input != None:
        run(opts.input, outpath=opts.outpath, ref=opts.ref)
    elif opts.test == True:
        run_test()

if __name__ == '__main__':
    main()