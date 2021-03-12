#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from label_classes import SVRecord
import re
import csv
import argparse
import pysam
import os
import pandas as pd
import statistics as st
import os 
from time import time
import logging
import matplotlib.pyplot as plt
from upsetplot import generate_counts, plot

def read_vcf(invcf):
    
    '''
    Function written by L.santuari
    from sv-channels 
    '''
    
    
    sv_list = []
    with pysam.VariantFile(invcf, 'r') as vcf:
        for rec in vcf.fetch():
            var = SVRecord(rec, svcaller = None)
            chrom1 = var.chrom
            pos1_start = var.start + var.cipos[0]
            pos1_end = var.start + var.cipos[1] + 1
            chrom2 = var.chrom2
            pos2_start = var.end + var.ciend[0]
            pos2_end = var.end + var.ciend[1] + 1
            svtype = var.svtype
            sv_list.append((chrom1, pos1_start, pos1_end,
                            chrom2, pos2_start, pos2_end, svtype))
        return sv_list
def get_bplist(svlist):
    
    bplist = []
    
    for SV in svlist:
        chr1, pos1_start, pos1_end, chr2, pos2_start, pos2_end, sv = SV
        if sv in svs:
            bplist.append([pos1_start,pos2_start])
        
    return bplist

def get_bplists(vcf_lists):
    bp_lists = []
    for vcf in vcf_list:
        svlist = read_vcf(svlist)
        bplist = get_bplist(svlist)
        bp_lists.append(bplist)
    
    return bp_lists
def make_UpSetPlot(GT ,bp_lists, vcf_names):
    bools = []
    for bplist in bp_lists:
        boollist = []
        for bp in GT:
            boollist.append(bpset in bp_lists)
            bools.append(boollist)
    
    dic = {}
    for nr, vcf in enumerate(vcf_names):
        dic[vcf] = bools[nr]
    dic['breakpoints'] = GT
    
    df = pd.DataFrame(dic)
    cols = df.columns.difference(['breakpoints']).tolist()
    s = df.groupby(cols).size()
    
    plot(s, show_counts='%d', sort_by="cardinality")
    plt.title('Intersection breakpoints SV callers')
    
    plt.show()

    # save plot
    currentfig = plt.gcf()
    currentfig.savefig('UpSetPlot %s' % ' '.join(vcf_names))

def get_intersection_SVs(GTvcf, vcflist, vcfnames):
    
    bplists = get_bplists(vcf_lists)
    GT = get_data(GTvcf)
    
    # make pot
    make_UpSetPlot(GT, bplists, vcfnames)
    
    # get list with overlapping SVs
    if len(vcfnames) == 3:
        for nr, bp in enumerate(bplists):
            
            dif1vs2 = [[item, vcfnames] for item in bp if item not in bplist[nr-1] and not in bplist[nr-2]]
            dif2vs1 = [[item,'%s_%s' % (vcfnames[nr], vcfnames[nr-1])] for item in bp if item in bplist[nr-1] and not in bplist[nr-2]]
            difference = dif1vs2 + dif2vs1
            fname = 'Uniq_svcalls_%s' % ('_'.join(vcfnames))
            fout = os.path.join(outdir,fname)
            with open('your_file.txt', 'w') as f:
                for item in difference:
                    f.write("%s\n" % item)
    
    # see difference sv-channels compared to gridss and manta
    if dif_svgridman:
        svchan, grid, man = bplist
        
        difference = [item for item in bplst[] if item]
        
        
        

def main():
    parser = argparse.ArgumentParser(description='Extract read positions (discordant/clipped/split)')
    parser.add_argument('-lvcf',
                        '--listvcf',
                        nargs ="+",
                        help="List with vcf files")

    parser.add_argument('-nvcf',
                        '--namesvcf',
                        type=str,
                        help = "Names of vcf files")
    parser.add_argument('-sv',
                        '--svtype',
                        nargs ="+",
                        help = "Specify list of SVs")
 
    parser.add_argument('-p',
                        '--outpath',
                        type=str,
                        default='.',
                        help = 'Specify output path')
    parser.add_argument('-GT',
                        '--GTvcf'
                        type=str,
                        help = 'file with ground truth SV breakpoint')
    parser.add_argument('-DN',
                        '--DataName',
                        type=str,
                        help = "Specify name data/sample")
    
    args = parser.parse_args()  

    cmd_name = 'Intersection_BP'
    outdir = os.path.join(args.outpath, args.DataName, cmd_name)
    os.makedirs(outdir, exist_ok=True)
    
    
    FORMAT = '%(asctime)s %(message)s'
    logfilename = os.path.join(outdir, '_'.join(['Logfile',str(time())]))
    
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)    
    logging.info('Succesfully created logfile')
    
    t0 = time()
    
    get_intersection(GTvcf =  = args.GTvcf, vcflist = args.lvcf, vcfnames = args.fread)
    
    logging.info('Time: Get intersection between sv callers %s: %f' % (args.bam, (time() - t0)))
    
    
if __name__ == '__main__':
       main()    
                    

