#!/usr/bin/env python
# coding: utf-8

# In[1]:

from label_classes import SVRecord
import argparse
import pysam
import os
import pandas as pd
import statistics as st
import os 
from time import time
import logging
import sys

def av_distance_reads(bamfile):
    
    distances = []
    with pysam.AlignmentFile(bamfile,'r') as bam:
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
def get_svbp_list(svlist, svtype):
    
    svbp_pos = []
    
    for SV in svlist:
        chr1, pos1_start, pos1_end, chr2, pos2_start, pos2_end, sv = SV
        
        if sv in svtype:

            if svtype in ['DEL','INS']:
                if abs(int(pos1_start) - int(pos2_start)) > 49:
                    svbp_pos.append([chr1,pos1_start,pos1_end])
                    svbp_pos.append([chr2, pos2_start, pos2_end])
            else:
                svbp_pos.append([chr1,pos1_start,pos1_end])
                svbp_pos.append([chr2, pos2_start, pos2_end])

    
    return svbp_pos

def get_reads_SVpos(bamfile, bp_pos, mean, stdev, winsize):
    chromdic = {}
    svpos = []
    cols = ['breakpoint_pos','read_chrom','read_start','read_end','read_side','read_type','mate_chrom','mate_start', 'mate_end', 'mate_side','mate_type', 'orientation']
    
    lsread = []
            
    with pysam.AlignmentFile(bamfile,'r') as bam:
        
        # make dictionary for chromosomes with their length
        for i in bam.header['SQ']:
            chromdic[i['SN']] = i['LN']

        for nr, breakpoint in enumerate(bp_pos):
            #print(breakpoint)

            logging.info('start processing reads at %d' % int(breakpoint[1]))
            logging.info('%d of the %d, %d breakpoints left' % (nr, len(bp_pos), (len(bp_pos) - nr)))
            print('%d of the %d, %d breakpoints left' % (nr, len(bp_pos), (len(bp_pos) - nr)))
            
            fetchstart = int(breakpoint[1]) - winsize
            fetchend = int(breakpoint[2]) + winsize 
            chrom = breakpoint[0]
            
            chromlen = chromdic[chrom]
            
            # make sure window falls within range
            if fetchstart < 1:
                logging.info('windowstart falls out of range: %d, adding to start' % fetchstart)
                while fetchstart <1:
                    fetchstart += 1
            if fetchend > chromlen:
                logging.info('windowend falls out of range: %d, substracting from end' % fetchend)
                while fetchend > chromlen:
                    fetchend -= 1
            
            n_r = 10 ** 6
            last_t = time()

            for i,read in enumerate(bam.fetch(chrom, fetchstart, fetchend)):

                if not i % n_r:
                    logging.info("%d reads processed (%f alignments / s)" % \
                                 (i, n_r / (time() - last_t)))
                    last_t = time()

                brkpnt = breakpoint[1]
                if is_clipped(read):
                    logging.info('start process clipped read %s' % read.qname)
                    r_type = 'CLIPPED'
                    r_start, r_chr, r_end, r_side, orien = get_read_info(read)


                    
                    if not read.mate_is_unmapped and read.is_paired:
                        mate = bam.mate(read)
                        m_type = get_mate_type(read, mate, mean, stdev)

                        m_start, m_chr, m_end, m_side, m_orien = get_read_info(mate)
                    else:
                        m_type = 'UNMAPPED'
                        m_start, m_chr, m_end, m_side, m_orien = [0, 0, 0,'NaN','NaN']

        
                    #print(df)
                    lsread.append([brkpnt,r_chr ,r_start, r_end, r_side, r_type, m_chr, m_start, m_end, m_side, m_type, orien])

                elif not read.is_unmapped and read.is_paired and read.has_tag('SA'):
                    logging.info('start process split read %s' % read.qname)
                    r_type = 'SPLIT'
                    if len(read.get_tag('SA').split(",")) == 6:
                        ID = read.qname

                        r_start, r_chr, r_end,r_side, orien = get_read_info(read)


                        if not read.mate_is_unmapped:
                            mate = bam.mate(read)
                            m_type = get_mate_type(read, mate, mean, stdev)

                            m_start, m_chr, m_end, m_side, m_orien = get_read_info(mate)
                        else:
                            m_type = 'UNMAPPED'
                            m_start, m_chr, m_end, m_side, m_orien = [0, 0, 0,'NaN','NaN']



                        lsread.append([brkpnt,r_chr ,r_start, r_end, r_side, r_type, m_chr ,m_start, m_end, m_side, m_type, orien])

                        # get SA alignment
                        logging.info('start process Supplementary Alignment of read %s' % read.qname)
                        sa_chr, sa_start, sa_end, sa_orien, sa_side = get_sa_info(read)
                    

                        lsread.append([brkpnt,r_chr ,r_start, r_end, r_side, r_type, sa_chr, sa_start, sa_end, sa_side,'SA', sa_orien])


                elif not read.is_unmapped and not read.mate_is_unmapped and read.is_paired and not right_clipped(read) and not left_clipped(read):
                    dis = abs(read.next_reference_start - read.reference_start)
                    if dis > 10:
                        if dis >= (mean+(stdev)) or dis <= (mean-(stdev)):
                            ID = read.qname
                            r_start, r_chr, r_end, r_side, orien = get_read_info(read)
                            r_type = 'DISCOR_1X'

                            
                            if not read.mate_is_unmapped:
                                mate = bam.mate(read)
                                m_type = get_mate_type(read, mate, mean, stdev)

                                m_start, m_chr, m_end, m_side, m_orien = get_read_info(mate)
                            else:
                                m_type = 'UNMAPPED'
                                m_start, m_chr, m_end, m_side, m_orien = [0, 0, 0,'NaN','NaN']


                            if dis >= (mean+(2*stdev)) or dis <= (mean-(2*stdev)):
                                r_type = 'DISCOR_2X'

                            if dis >= (mean+(3*stdev)) or dis <= (mean-(3*stdev)):
                                r_type = 'DISCOR_3X'
                            
                            lsread.append([brkpnt,r_chr ,r_start, r_end, r_side, r_type, m_chr, m_start, m_end, m_side, m_type, orien])
                            
    dfread = pd.DataFrame(lsread, columns = cols)
    
    return dfread

def get_reads_at_svbp(bamfile, vcf, DataName, sv, outdir):
    
    # determine mean and stdev
    mean, stdev = av_distance_reads(bamfile)
    winsize = mean
    logging.info('Calculated mean and stdev from file')
    
    #make svlist
    
    svlist = read_vcf(vcf)
    if len(svlist) == 0:
        logging.info('Could not extract sv positions based on given vcf file: %s' % vcf)
        sys.exit('Could not generate SV list ')
    else:
        logging.info('Generated list of SVs')
    
    #make breakpoint pos list
    if type(sv) != list:
        sv = sv.split()
    bp_pos = get_svbp_list(svlist, sv)
    if len(bp_pos) == 0:
        logging.info('Could not generate breakpoint list based on given sv list')
        sys.exit('Could not generate breakpoint list')
    else:
        logging.info('Generated lits of breakpoints of SV: %s' % sv)
    
    # Extract reads at breakpoint positions

    dfread = get_reads_SVpos(bamfile, bp_pos, mean, stdev, winsize)
    logging.info('Extracted all reads')
    
    # to csv
    fname = 'read_at_SVbp_{}.csv'.format(DataName)
   
    outfile = os.path.join(outdir,fname)
    dfread.to_csv(outfile,index=False)
    logging.info('Created csv file for %s' % fname)


def main():
    parser = argparse.ArgumentParser(description='Extract read positions (discordant/clipped/split)')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        help="Specify input file (BAM)")

    parser.add_argument('-DN',
                        '--DataName',
                        type=str,
                        help = "Specify name data/sample")
    parser.add_argument('-sv',
                        '--svtype',
                        nargs ="+",
                        help = "Specify list of SVs")
 
    parser.add_argument('-p',
                        '--outpath',
                        type=str,
                        default='.',
                        help = 'Specify output path')
    parser.add_argument('-vcf',
                        '--svvcf',
                        type=str,
                        help = 'Specify vcf with SVs')
    args = parser.parse_args()  

    cmd_name = 'readpos_SVbp'
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
    
    get_reads_at_svbp(bamfile = args.bam, vcf = args.svvcf, DataName = args.DataName, sv = args.svtype, outdir = outdir)
    
    logging.info('Time: read positions at SV breakpoints on BAM %s: %f' % (args.bam, (time() - t0)))
    
    
if __name__ == '__main__':
       main()

