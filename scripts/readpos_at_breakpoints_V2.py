#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from label_classes import SVRecord
import argparse
import pysam
import os
import pandas as pd
import statistics as st
import os 
from time import time
import logging

def av_distance_reads(bamfile):
    
    distances = []
    with pysam.AlignmentFile(bamfile,'r') as bam:
        for read in bam.fetch():

            if read.is_proper_pair and (not read.is_unmapped)             and read.is_reverse != read.mate_is_reverse             and read.reference_name == read.next_reference_name:

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
        m = []
        s = []

        for nr, i in enumerate(sa_cigar):
        
            if i == 'M':
                m.append(nr)
            if i in ['S','H']:
                s.append(nr)
        if len(s) == 1:
            if m[0] > s[0]:
                bp = sa_cigar[s[0]+1:m[0]]
                sa_side = 'LEFT'
            
            else:
                bp = sa_cigar[:int(m[0])]
                sa_side = 'RIGHT'
            
        elif len(s) > 1:
            
            if m[0] > (s[0] and s[1]):
                bp = cigar[s[1]+1:m[0]]
                sa_side = 'LEFT'
                
            elif m[0] > s[0] and m[0] < s[1]:
                bp = cigar[s[0]+1:m[0]]
                sa_side = 'BOTH'
            
            else:
                bp = cigar[:m[0]]
                sa_side = 'RIGHT'
        
        sa_end =int(bp) + int(sa_start)
        
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
        sa_end, sa_side = get_sa_end(cigar, sa_start)
       
    return sa_chr, sa_start, sa_end, orien, sa_side

def get_read_info(read):
    
    r_start = read.reference_start
    r_end = read.reference_end
    r_chr = read.reference_name
    
    
    
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

    
    return r_start, r_chr, r_end, side, orien

def get_mate(bamfile, read):
    with pysam.AlignmentFile(bamfile,'r') as bam:
        for i in bam.fetch('12',1000,1010):
            if i == read:

                try: 
                    mate = bam.mate(read) 
                    return mate

                except ValueError: 

                    # Invalid mate (usually post-filtered) 

                    return False
            
def get_mate_type(read,mate, mean, stdev):
    
    if mate.has_tag('SA'):
        matetype = 'SPLIT'
    elif is_clipped(mate):
        matetype = 'CLIPPED'
    elif not mate.is_unmapped:
        dis = abs(read.reference_start - mate.reference_start)
        
        if dis >= (mean+(3*stdev)) or dis <= (mean-(3*stdev)):
            matetype = 'DISCOR_3X'
            
        elif dis >= (mean+(stdev)) or dis <= (mean-(stdev)):
            matetype = 'DISCOR_1X'
            
        elif dis >= (mean+(2*stdev)) or dis <= (mean-(2*stdev)):
            matetype = 'DISCOR_2X'
        else:
            matetype = 'NORMAL'
            
    return matetype


def get_reads_SVpos(bamfile, svlist, mean, stdev, winsize, get_SA = True):
     
    cols = ['breakpoint_pos','read_chrom','read_start','read_end','read_side','read_type','mate_chrom','mate_start', 'mate_end', 'mate side','mate_type', 'orientation']
    cols_SA = ['read.qname','SA_start','SA_end','SA_chrom', 'orientation', 'side','type']

    dfs = [pd.DataFrame(), pd.DataFrame()]
    dfSA = [pd.DataFrame(), pd.DataFrame()]
    for SV in svlist:
        
        chr1, pos1_start, pos1_end, chr2, pos2_start, pos2_end, sv = SV
        
        chroms = [chr1, chr2]
        
        for nr, breakpoint in enumerate([[pos1_start, pos1_end],[pos2_start, pos2_end]]):
        
            logging.info('start processing reads at %d' % int(breakpoint[0]))
            
            fetchstart = int(breakpoint[0]) - winsize
            fetchend = int(breakpoint[1]) + winsize 

            with pysam.AlignmentFile(bamfile,'r') as bam:
                
                n_r = 10 ** 6
                last_t = time()
                
                for i,read in enumerate(bam.fetch(chroms[nr],fetchstart,fetchend)):
                    
                    if not i % n_r:
                        logging.info("%d reads processed (%f alignments / s)" %                                      (i, n_r / (time() - last_t)))
                        last_t = time()
                    
                    brkpnt = breakpoint[0]
                    if is_clipped(read):
                        logging.info('start process clipped read %s' % read.qname)
                        r_type = 'CLIPPED'
                        r_start, r_chr, r_end, r_side, orien = get_read_info(read)
                        
                        
                        mate = get_mate(read=read, bamfile= bamfile)
                        if mate:
                            m_type = get_mate_type(read, mate, mean, stdev)
                                
                            m_start, m_chr, m_end, m_side, m_orien = get_read_info(mate)
                        else:
                            m_type = 'UNMAPPED'
                            m_start, m_chr, m_end, m_side, m_orien = ['NaN','NaN','NaN','NaN','NaN']
                        
                        df = pd.DataFrame([[brkpnt,r_chr ,r_start, r_end, r_side, r_type, m_chr, m_start, m_end, m_side, m_type, orien]], columns = cols)
                        #print(df)
                        dfs[nr] = dfs[nr].append(df)
                    elif not read.is_unmapped and read.has_tag('SA'):
                        logging.info('start process split read %s' % read.qname)
                        r_type = 'SPLIT'
                        if len(read.get_tag('SA').split(",")) == 6:
                            ID = read.qname


                            
                            r_start, r_chr, r_end,r_side, orien = get_read_info(read)
                            
                            
                            mate = get_mate(read=read, bamfile= bamfile)
                            if mate:
                                m_type = get_mate_type(read, mate, mean, stdev)

                                m_start, m_chr, m_end, m_side, m_orien = get_read_info(mate)
                            else:
                                m_type = 'UNMAPPED'
                                m_start, m_chr, m_end, m_side, m_orien = ['NaN','NaN','NaN','NaN','NaN']
                            
                            
                            df = pd.DataFrame([[brkpnt,r_chr ,r_start, r_end, r_side, r_type, m_chr ,m_start, m_end, m_side, m_type, orien]], columns = cols)
                            dfs[nr] = dfs[nr].append(df)
                            
                            if get_SA:
                                logging.info('start process Supplementary Alignment of read %s' % read.qname)
                                sa_chr, sa_start, sa_end, sa_orien, sa_side = get_sa_info(read)
                                SA = pd.DataFrame([[ID, sa_chr, sa_start, sa_end, sa_orien, sa_side, 'SA']], columns = cols_SA)
                                dfsa = pd.concat([df,SA], axis = 1)
                
                                dfSA[nr] = dfSA[nr].append(dfsa)
                                
                                
                    elif not read.is_unmapped and not read.mate_is_unmapped and not right_clipped(read) and not left_clipped(read):
                        dis = abs(read.next_reference_start - read.reference_start)
                        if dis > 10:
                            if dis >= (mean+(stdev)) or dis <= (mean-(stdev)):
                                ID = read.qname
                                r_start, r_chr, r_end, r_side, orien = get_read_info(read)
                                r_type = 'DISCOR_1X'
                                
                                mate = get_mate(read=read, bamfile= bamfile)
                                if mate:
                                    m_type = get_mate_type(read, mate, mean, stdev)

                                    m_start, m_chr, m_end, m_side, m_orien = get_read_info(mate)
                                else:
                                    m_type = 'UNMAPPED'
                                    m_start, m_chr, m_end, m_side, m_orien = ['NaN','NaN','NaN','NaN','NaN']
                                

                                if dis >= (mean+(2*stdev)) or dis <= (mean-(2*stdev)):
                                    r_type = 'DISCOR_2X'

                                if dis >= (mean+(3*stdev)) or dis <= (mean-(3*stdev)):
                                    r_type = 'DISCOR_3X'

        
                                df = pd.DataFrame([[brkpnt,r_chr ,r_start, r_end, r_side, r_type, m_chr, m_start, m_end, m_side, m_type, orien]], columns = cols)
                                dfs[nr] = dfs[nr].append(df)
    dfbp1, dfbp2 = dfs
    
    if get_SA: 
        
        dfSA1, dfSA2 = dfSA
        return dfbp1, dfbp2, dfSA1, dfSA2
    
    else:
        return dfbp1, dfbp2

    
def get_reads_at_breakpoints(bamfile, vcf, DataName, returnCSV, get_SA = True):
    
    # make outdir
    
    parent_dir = os.getcwd()
    outdir = os.path.join(parent_dir,DataName)
    os.makedirs(outdir, exist_ok=True)
    
    
    # make logfile
    
    FORMAT = '%(asctime)s %(message)s'
    logfilename = os.path.join(outdir, '_'.join(['Logfile',str(time())]))   
    logging.basicConfig(format=FORMAT,
                            filename=logfilename,
                            filemode='w',
                            level=logging.INFO)
    t0 = time()
    
    
    # generate SV list and start, end list

    SV_list = read_vcf(vcf)
    
    # determine mean and stdev, define winsize (average length insert)
    mean,stdev = av_distance_reads(bamfile)
    winsize = mean
                             
    # get reads overlapping the breakpoints
    df_bp1, df_bp2, dfSA1, dfSA2 = get_reads_SVpos(bamfile, svlist = SV_list, mean= mean, stdev = stdev, winsize = winsize, get_SA = True)
    
    if returnCSV:
        
        df_bp1_fname = 'reads_bp1_{}_overlap.csv'.format(DataName)
        df_bp2_fname = 'reads_bp2_{}_overlap.csv'.format(DataName)
        
        
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
    parser.add_argument('-gsa',
                        '--get_sa',
                        type=bool,
                        default=True,
                        help = '')
 

    args = parser.parse_args()  
    
        
    get_reads_at_breakpoints(bamfile = args.bamfile, vcf = args.SV_vcf, DataName = args.DataName, returnCSV = args.returnCSV, get_SA = args.get_sa)
    
if __name__ == '__main__':
       main()   

