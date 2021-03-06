#!/usr/bin/env python
# coding: utf-8

# In[17]:


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



def av_distance_reads(bamfile):
    
    distances = []
    with pysam.AlignmentFile(bamfile,'r') as bam:
        for read in bam.fetch():

            if read.is_proper_pair and (not read.is_unmapped) and read.is_reverse != read.mate_is_reverse \
            and read.reference_name == read.next_reference_name:

                dis = abs(read.next_reference_start - read.reference_start)
                if dis < 10**3:

                    distances.append(dis)

                if len(distances) > 2 * 10 ** 6:
                    break
        mean = st.mean(distances)
        stdev = st.stdev(distances)
        return mean, stdev




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
    
def get_coor(file):
    coor_list = []
    with open(file,'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'breakpoint_pos':
                continue
            coor_list.append(row)
    return coor_list

def get_svbp_list(svlist, svtype):
    
    svbp_pos = []
    
    for SV in svlist:
        chr1, pos1_start, pos1_end, chr2, pos2_start, pos2_end, sv = SV
        
        if sv in svtype:

            if svtype in ['DEL','INS']:
                if abs(int(pos1_start) - int(pos2_start)) > 49:
                    svbp_pos.append([pos1_start,pos2_start])
            else:
                svbp_pos.append([pos1_start,pos2_start])
            

    
    return svbp_pos

def subset_twobp_support(df, bplist, winsize):
    
    lsreads = []
    twosup = []
    cols = ['breakpoint1','breakpoint2','read_chrom','read_start','read_end','read_side','read_type',\
            'mate_chrom','mate_start', 'mate_end', 'mate_side','mate_type', 'orientation']

  
    m_chr = list(df['mate_chrom'])
    m_start = list(df['mate_start'])
    m_end = list(df['mate_end'])
    m_type = list(df['mate_type'])
    
    nr = 0
    nrlen = len(m_start)
    for readnr, (mchr, mstart, mend, mtype) in enumerate(zip(m_chr, m_start, m_end, m_type)):
        if nr == 1000:
            logging.info('%d to go' % nrlen)
            print('%d to go' % nrlen)
            nr = 0
        nrlen -= 1
        
        
        if m_type != 'SA':
            r2= df[(df['read_chrom'] == mchr) & (df['read_start'] == mstart) & (df['read_end'] == mend)]

            if (len(r2)) > 0:
                r1= df[(df['mate_chrom'] == mchr) & (df['mate_start'] == mstart) & (df['mate_end'] == mend)]
                bpr1, bpr2 = int(list(r1['breakpoint_pos'])[0]), int(list(r2['breakpoint_pos'])[0])
                
                if bpr1 == bpr2:
                    continue
                r1inf = list(r1[list(r1.columns)[1:]].values[0])
                
                if [bpr1, bpr2] in bplist:
                    r1inf.insert(0,bpr1)
                    r1inf.insert(1, bpr2)
                
                    lsreads.append(r1inf) 
                    twosup.append(readnr)
                    nr +=1

        else:
            split = df[(df['mate_type'] == 'SA'),(df['mate_start'] == m_start) &(df['mate_end'] == m_end)]
            if len(split) != 0:
                for bp in bplist:
                    if int(split['breakpoint_pos'].values[0]) == bp[0]:
                        otherbp = bp[1]
                        bpnr = 1
                    elif int(split['breakpoint_pos']) == bp[1]:
                        otherbp = bp[0]
                        bpnr = 2
                    else:
                        continue
                
                    splitlist = list(split.values[0])
                
                    if int(split['mate_start']) > (otherbp - winsize) and int(split['mate_start']) < (otherbp + winsize):
                        print('Found match')
                        if bpnr == 1:
                            splitlist.insert(0,otherbp)
                        else:
                            splitlist.insert(1,otherbp)
                        twosupport.append(splitlist)
                        print('SA concatenated')
                        nr +=1 
                    else:
                        twosup.append(readnr) 
                            
    dfreads = pd.DataFrame(lsreads, columns = cols)
    return dfreads, twosup

def one_support(lsreads, index_twosupport):
    
    onesup = []
    cols = ['breakpoint','read_chrom','read_start','read_end','read_side', 'read_type',\
            'mate_chrom','mate_start', 'mate_end', 'mate_side','mate_type','orientation']
    
    for nr, read in enumerate(lsreads):
        
        if nr not in index_twosupport:
            onesup.append(read)
            logging.info('append info one_support readpair')
    dfonesup = pd.DataFrame(onesup, columns = cols)
    logging.info('generated DF one breakpoint supporting readpair')
    return dfonesup

def get_specific_readtype_sv_overlap(dfreads, bplist, split = False, dis = True, clip = True, gts = False):
    
    twobp = pd.DataFrame()
    if gts:
        bpref = ['breakpoint1', 'breakpoint2']
    else:
        bpref = ['breakpoint_pos', 'breakpoint_pos']
        onebp = pd.DataFrame()
    
    filters = []
    filterlist = ['SPLIT','CLIP','DISCOR']
    
    for nr, f in enumerate([split,clip,dis]):
        
        if not f:
            
            filters.append(filterlist[nr])
        
    print(filters)

    for bpset in bplist:
        
        data1 = dfreads[dfreads[bpref[0]] == bpset[0]]
        data2 = dfreads[dfreads[bpref[1]] == bpset[1]]
        
        test = []
        
        if data1.empty or data2.empty:
            
            continue
        
        for data in [data1, data2]:
            
            rtype = '_'.join(list(data.read_type.unique()))
            #print(rtype)
            match = 0
            
            for f in filters:
                
                if re.search(f,rtype):
                    match +=1
                
                if match == 0:
                    test.append('YES')
                
                else:
                    test.append('NO')
        
        if not 'NO' in test:
            #print('found match')
            
            twobp.append(data1)
            twobp.append(data2)
            
        elif ('NO' and 'YES') in test:
            #print('found match')
            if not gts:
                if test[0] == 'YES':
                    onebp.append(data1)
                else:
                    onebp.append(data2)
                
    if gts:
        return twobp
    else:
        return onebp, twobp

def subset_reads(bamfile, vcf, fread, sv, outdir, DataName, split, clip, dis, gts = True):
    
    logging.info('Generating list of read pos file')
    lread = get_coor(fread)
    logging.info('Generating dataframe of read file')
    dfread = pd.read_csv(fread)
    
    svlist = read_vcf(vcf)
    mean, stdev = av_distance_reads(bamfile)
    winsize = mean
    
    # make breakpoint list
    
    bplist = get_svbp_list(svlist, sv)
    
    # start subsetting reads
    
    ## subset reads: do readpairs bridge both Breakpoints of SV?
    logging.info('start subsetting reads, (readpairs bridge SV)')
    
    twosup, index = subset_twobp_support(dfread, bplist, winsize)
    if len(twosup) == 0:
        logging.info('Could not identify readpairs supporting both breakpoints of SV')
    
    
    # subset reads: Get dataframe reads only supporting one BP of SV
    logging.info('start subsetting reads (readpairs supporting only one breakpoint of SV)')
    onesup = one_support(lread, index)
    if len(onesup) == 0:
        logging.info('Could not identify readpairs supporting only one breakpoint of SV')
        
    # Filter reads based on type
    n = ['SPLIT','CLIP','DIS']
    rtype = []
    for nr, f in enumerate([split,clip,dis]):
        if not f:
            rtype.append(n[nr])


    logging.info('filtering BP based on readtypes, implemented filters %s' % ' '.join(rtype))
    onebp, twobp = get_specific_readtype_sv_overlap(dfread, bplist, split = split, dis = dis, clip = clip)
    

     # Filter reads that cover both bp of SV on type
        
    if gts:
        logging.info('filtering BP supported by readpairs based on readtypes, implemented filters %s' % ' '.join(rtype))
        twobp_twosup = get_specific_readtype_sv_overlap(twosup, bplist, split = split, dis = dis, clip = clip, gts = True)
        
    # generate fname
    
    filters = ''.join(rtype)
    fnames = [string.format(DataName) for string in ['readpair_1bp_{}.csv','readpair_2bp_{}.csv']] + \
    [string.format(filters,DataName) for string in ['1bp_{}_filtered_{}.csv','2bp_{}_filtered_{}.csv']] 
    
    files = [onesup,twosup,onebp,twobp]
    
    if gts:
        fnames.append('2bp_twosup_{}_filtered_{}.csv'.format(filters,DataName))
        files.append(twobp_twosup)
    # to file

    for nr, f in enumerate(files):
        if not f.empty:
            fout = os.path.join(outdir,fnames[nr])
            f.to_csv(fout,index=False)
            logging.info('Generated file %s' % fnames[nr])
        

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
    parser.add_argument('-fr',
                        '--fread',
                        type=str,
                        help = 'file with read positions at SV breakpoint')
    parser.add_argument('-sp',
                        '--split',
                        type=bool,
                        default=False,
                        help = 'filter split reads, False == breakpoints must not be supported by splitreads')
    parser.add_argument('-cl',
                        '--clip',
                        type=bool,
                        default=True,
                        help = 'filter clipped reads,False == breakpoint must not be supported by clipped reads')
    parser.add_argument('-ds',
                        '--dis',
                        type=bool,
                        default=True,
                        help = 'filter discordant reads,False == breakpoint must not be supported by discordant reads')
    parser.add_argument('-gts',
                        '--gtwosup',
                        type=bool,
                        default=True,
                        help = ' get readpairs supporting two BPs, Filter data of readpairs that bridge SV breakpoints')
    
    args = parser.parse_args()  

    cmd_name = 'readpos_SVbp/subset'
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
    
    subset_reads(bamfile = args.bam, vcf = args.svvcf, fread = args.fread,sv = args.svtype, outdir = outdir, DataName = args.DataName, split = args.split, clip = args.clip, dis = args.dis, gts = args.gtwosup)
    
    logging.info('Time: subset reads at SV breakpoint positions %s: %f' % (args.bam, (time() - t0)))
    
    
if __name__ == '__main__':
       main()

    


# In[ ]:




