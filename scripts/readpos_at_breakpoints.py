#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##### import packages
from label_classes import SVRecord
import pysam
import os
import pandas as pd
import statistics as st
import os 
from time import time
import logging

def av_distance_reads(bam):
    
    distances = []
    for read in bam.fetch():
        
        if read.is_proper_pair and (not read.is_unmapped) \
        and read.is_reverse != read.mate_is_reverse \
        and read.reference_name == read.next_reference_name:
            
            dis = abs(read.next_reference_start - read.reference_start)
            if dis < 10**3:
                
                distances.append(dis)
            
            if len(distances) > 2 * 10 ** 6:
                break
    mean = st.mean(distances)
    stdev = st.stdev(distances)
    return mean, stdev

def is_clipped(read):

    if not read.is_unmapped:
        if left_clipped(read) or right_clipped(read):
            if not read.has_tag('SA'):
                return True
            else:
                return False
def left_clipped(read):
    '''
    Function to determine is read is left clipped
    by L.santuari (sv-channels)
    '''
    
    if read.cigartuples is not None:
            if read.cigartuples[0][0] in [4,5]:
                return True
            else:
                return False
    else:
        return False
    
    
def right_clipped(read):
    '''
    Function to determine is read is left clipped
    by L.santuari (sv-channels)
    '''
    
    if read.cigartuples is not None:
            if read.cigartuples[-1][0] in [4,5]:
                return True
            else:
                return False
    else:
        return
def get_read_info(read):
    
    r_start = read.reference_start
    r_end = read.reference_end
    r_chr = read.reference_name
    
    m_start = read.next_reference_start
    m_chr = read.next_reference_name
    
    
    if left_clipped(read) and not right_clipped(read):
        side = 'LEFT'
    elif not left_clipped(read) and right_clipped(read):
        side = 'RIGHT'
    elif left_clipped(read) and right_clipped(read):
        side = 'BOTH'
    else:
        side = 'NaN'
    
    
    

    if not read.is_reverse and read.mate_is_reverse:
        orien = 'FR'
    elif not read.is_reverse and not read.mate_is_reverse:
        orien = 'FF'

    elif read.is_reverse and not read.mate_is_reverse:
        orien = 'RF'
    else:
        orien = 'RR'

    
    return r_start, r_chr, r_end, m_start, m_chr, side, orien
    
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
def get_sa_info(read, end = True):
    
    def get_sa_end(sa_cigar, sa_start):
        
        for nr, i in enumerate(sa_cigar):
            if i == 'M':
                m = nr
                print(m)
            if i in ['S','H']:
                s = nr
        if m > s:
            bp = sa_cigar[s+1:m]
            sa_side = 'LEFT'
        else:
            bp = sa_cigar[:int(m)]
            sa_side = 'RIGHT'
        sa_end =int(bp) + sa_start
        
        return sa_end, sa_side
    
    satag = read.get_tag('SA')
    if len(satag) > 0:
        sa = satag.split(",")
        sa_chr, sa_start, sa_orien, cigar = sa[0], sa[1], sa[2], sa[3]
        
        # determine orientation
        
        orien = {'+':True,'-':False}
        if not read.is_reverse and not orien[sa_orien]:
            orien = 'FR'
        elif read.is_reverse and orien[sa_orien]:
            orien = 'RF'
        elif read.is_reverse and not orien[sa_orien]:
            orien = 'RR'
        else:
            orien = 'FF'
        sa_end = get_sa_end(sa_start,cigar)
       
    return sa_chr, sa_start, sa_end, orien, sa_side, cigar
    
def get_reads_SVpos(bamfile, winsize, svlist, get_SA = True):
     
    cols = ['breakpoint_pos','read.qname','read_start','read_end','read_chrom', 'orientation', 'side','type']
    cols_SA = ['read.qname','SA_start','SA_end','SA_chrom', 'orientation', 'side','type']

    dfs = [pd.DataFrame(), pd.DataFrame()]
    dfSA = [pd.DataFrame(), pd.DataFrame()]
    for SV in SV_list:
        chr1, pos1_start, pos1_end, chr2, pos2_start, pos2_end, sv = SV
        
        
  
        for nr, breakpoint in enumerate([[pos1_start, pos1_end],[pos2_start, pos2_end]]):
        
            logging.info('start processing reads at %d' % int(breakpoint[0]))
            
            fetchstart = int(breakpoint[0]) - winsize
            fetchend = int(breakpoint[1]) + winsize 

            with pysam.AlignmentFile(bamfile,'r') as bam:
                mean,stdev = av_distance_reads(bam)
                
                n_r = 10 ** 6
                last_t = time()
                
                for i,read in enumerate(bam.fetch(chr1,fetchstart,fetchend)):
                    
                    if not i % n_r:
                        logging.info("%d clipped reads processed (%f alignments / s)" % \
                                     (i, n_r / (time() - last_t)))
                        last_t = time()
                    
                    brkpnt = breakpoint[0]
                    if is_clipped(read):
                        logging.info('start process clipped read %s' % read.qname)
                        type = 'CLIPPED'
                        ID = read.qname
                        r_start, r_chr, r_end, m_start, m_chr, side, orien = get_read_info(read)
                        
                        df = pd.DataFrame([[brkpnt, ID,r_start, r_chr, r_end, side, orien, type]], columns = cols)
                        #print(df)
                        dfs[nr] = dfs[nr].append(df)
                    elif not read.is_unmapped and read.has_tag('SA'):
                        logging.info('start process split read %s' % read.qname)
                        type = 'SPLIT'
                        if len(read.get_tag('SA').split(",")) == 6:
                            ID = read.qname


                            
                            r_start, r_chr, r_end, m_start, m_chr, side, orien = get_read_info(read)
                            
                            df = pd.DataFrame([[brkpnt, ID, r_start, r_chr, r_end, side, orien, type]], columns = cols)
                            dfs[nr] = dfs[nr].append(df)
                            
                            if get_SA:
                                logging.info('start process Supplementary Alignment of read %s' % read.qname)
                                sa_chr, sa_start, sa_end, sa_orien, sa_cigar = get_sa_info(read)
                                SA = pd.DataFrame([[ID, sa_chr, sa_start, sa_end, sa_orien, 'SA']]columns = cols_SA)
                                dfsa = pd.concat([df,SA], axis = 1)
                
                                dfSA[nr] = dfSA[nr].append(dfsa)
                                
                                
                    if not read.is_unmapped and not read.mate_is_unmapped and not right_clipped(read) and not left_clipped(read):
                        dis = abs(read.next_reference_start - read.reference_start)
                        if dis > 10:
                            if dis >= (mean+(stdev)) or dis <= (mean-(stdev)):
                                ID = read.qname
                                r_start, r_chr, r_end, m_start, m_chr, side, orien = get_read_info(read)
                                type = 'DISCOR_1X'

                                if dis >= (mean+(2*stdev)) or dis <= (mean-(2*stdev)):
                                    type = 'DISCOR_2X'

                                if dis >= (mean+(3*stdev)) or dis <= (mean-(3*stdev)):
                                    type = 'DISCOR_3X'

                                lst = [brkpnt, ID, r_start, r_chr, r_end, side, orien,type]
                                df = pd.DataFrame([lst], columns = cols)
                                dfs[nr] = dfs[nr].append(df)
    dfbp1, dfbp2 = dfs
    
    if get_SA: 
        
        dfSA1, dfSA2 = dfSA
        return dfbp1, dfbp2, dfSA1, dfSA2
    
    else:
        return dfbp1, dfbp2

def get_readpair_sv_overlap(df_bp1, df_bp2):
    
    overlap = list(set(list(df_bp1['read.qname'])) & set(list(df_bp2['read.qname'])))
    cols = ['bp_pos_r1','ID_r1','start_r1','end_r1','chrom_r1', 'side_r1', 'orien_r1','type_r1', \
            'bp_pos_r2','ID_r2','start_r2','end_r2','chrom_r2', 'side_r2', 'orien_r2','type_r2']
    
    overlapdf = pd.DataFrame()
    
    for i in overlap:
        df1 = df_bp1[df_bp1['read.qname'] == i]
        df2 = df_bp2[df_bp2['read.qname'] == i]
        df3 = pd.concat([df1, df2], axis=1)

        overlapdf = overlapdf.append(df3)
        #print(df1)
        #print(df2)

    overlapdf.columns = cols 
    
    return overlapdf
def get_df_split_SA(df_SA1,df_SA2, df_bp1,df_bp2):
    dfs = [pd.DataFrame,pd.DataFrame]
    SAs = [df_SA1, df_SA2]
    for nr, breakpoint in enumerate([df_bp1,df_bp2]):
        IDs = list(breakpoint['read.qname'])
        
        for ID in IDs:
            readinfo = breakpoint[breakpoint['read.qname'] == ID]
            SAinfo = SAs[nr][SAs[nr['read.qname'] == ID]
            read_SA = pd.concat([readinfo,SAinfo])
            dfs[nr] = dfs[nr].append(read_SA)
    
    read_SAbp1, read_SAbp2 = dfs
    
    return read_SAbp1, read_SAbp2
        
def get_reads_at_breakpoints(bamfile, vcf, winsize = 100, DataName, returnCSV, get_SA = True):
    # make logfile
    
    FORMAT = '%(asctime)s %(message)s'
    logfilename = os.path.join(outdir, '_'.join(['Logfile',str(time())]))   
    logging.basicConfig(format=FORMAT,
                            filename=logfilename,
                            filemode='w',
                            level=logging.INFO)
    t0 = time()
    
    # make outdir
    parent_dir  = os.getcwd()
    outdir = os.path.join(parent_dir,DataName)
    os.makedirs(outdir, exist_ok=True)
    
    # generate SV list

    assert os.path.isfile(vcf)
    SV_list = read_vcf(vcf)
    
    # get reads overlapping the breakpoints
    df_bp1, df_bp2, dfSA1, dfSA2 = get_reads_SVpos(bamfile, winsize = 100, svlist = SV_list, get_SA = True)
    
    # get readpairs overlapping SV
    overlapdf = get_readpair_sv_overlap(df_bp1, df_bp2)
    
    # record in csv format
    if returnCSV:
        
        df_bp1_fname = 'reads_bp1_{}_overlap.csv'.format(DataName)
        df_bp2_fname = 'reads_bp2_{}_overlap.csv'.format(DataName)
        overlap_fname = 'readpairs_{}_overlap.csv'.format(DataName)
        
        filenames = [df_bp1_fname,df_bp2_fname, overlap_fname]
        for nr,df in enumerate([df_bp1,df_bp2, overlapdf]):
            outfile = os.path.join(outdir,filenames[nr])
            df.to.csv(outfile,index=False)
            logging.info('Processed  %' % df)
        
        if get_SA:
            
            df_bp1_SA, df_bp2_SA = get_df_split_SA(dfSAa, dfSA2, df_bp1,df_bp2)
            
            df_bp1_SA_fname = 'reads_bp1_{}_overlap.csv'.format(DataName)
            df_bp2_SA_fname = 'reads_bp2_{}_overlap.csv'.format(DataName)
            
            for nr,df in enumerate([df_bp1,df_bp2, overlapdf]):
            
            outfile = os.path.join(outdir,filenames[nr])
            df.to.csv(outfile,index=False)
            
            logging.info('Processed  %' % df)
    
def main():
    parser = argparse.ArgumentParser(description='Extract read positions at SV breakpoints (discordant/clipped/split)')
    parser.add_argument('-b',
                        '--bamfile',
                        type=str,
                        help="Specify input file (BAM)")

    parser.add_argument('-DN',
                        '--DataName',
                        type=str,
                        help = "Specify name data/sample")
 
    parser.add_argument('-rcsv',
                        '--returnCSV',
                        type=bool,
                        default=True,
                        help = 'Specify whether you want the output to be in .csv format')
    parser.add_argument('-vcf',
                        '--SV_vcf',
                        type=bool,
                        default=False,
                        help = 'Specify wheter you want .json formatted output')
    parser.add_argument('-ws',
                        '--winsize',
                        type=int,
                        default=100,
                        help = '')
    parser.add_argument('-gsa',
                        '--get_sa',
                        type=bool,
                        default=True,
                        help = '')                      
    args = parser.parse_args()  
    
        
    get_reads_at_breakpoints(bamfile = args.bamfile, vcf = args.SV_vcf, winsize = args.winsize, DataName = args.DataName, rcsv = args.returnCSV, get_SA = args.get_sa)
    
if __name__ == '__main__':
       main()
    

