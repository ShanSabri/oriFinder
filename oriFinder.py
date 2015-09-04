#!/usr/bin/env python

'''
oriFinder.py

Author: Shan Sabri
Date Created: 09/04/2015
Date Modified: 09/04/2015

Implementation of an algorithm to identify
the origin of replication within a DNA sequence
(working with e.coli at the moment)
'''

from __future__ import print_function
from collections import defaultdict
from itertools import combinations, product, izip
from peaks import detect_peaks
from colorama import Fore

import matplotlib.pyplot as plt, numpy as np, time, textwrap, re


def read_FASTA(f):
    """Read in a file in FASTA format"""
    print ("File: ", f)
    print ("Reading in FASTA...", end=' ')

    in_file = open(f,'r')
    seqDict = {}
    name = None
    for line in in_file:
        line = line.rstrip()
        if line[0]=='>':
            name = line[1:]
            seqDict[name] = ''
        else:
            seqDict[name] = seqDict[name] + line

    print ("DONE!")
    return seqDict


def calc_GC_skew(seq):
    """Calculate the GC-skew for a given sequence"""
    print ("Sliding window.....", end=' ')

    values = []
    window = int(0.000005 * len(seq))
    # window = 25
    for i in xrange(0, len(seq), window):
        s = seq[i: i + window]
        g = s.count('G') + s.count('g')
        c = s.count('C') + s.count('c')
        skew = (g-c)#/float(g+c)
        values.append(skew)

    print ("DONE!")
    return values, window


def smoothListGaussian(myarray, degree=5000):
    """Smooth out the cumulative skew values for re-plotting"""
    myarray  = np.pad(myarray, (degree-1, degree-1), mode='edge')
    window   = degree*2-1
    weight   = np.arange(-degree+1, degree)/window
    weight   = np.exp(-(16*weight**2))
    weight  /= sum(weight)
    smoothed = np.convolve(myarray, weight, mode='valid')
    return smoothed


def plot_skew(skewArray):
    """Plots a skew diagram and returns the respective cumulative values"""
    print ("Plotting...........", end=' ')

    skewCumValues = np.cumsum(skewArray)
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # ax1.plot(range(len(skewArray)), skewArray, 'r', color='blue')  # noise values
    ax1.plot(range(len(skewArray)), skewCumValues, 'r', color='red') # cumulative values

    smoothedY = smoothListGaussian(skewCumValues)
    ax2.plot(range(len(smoothedY)), smoothedY, 'r', color='black')   # smooth values

    fig.savefig('SkewPlot.png')

    print ("DONE!")
    return(skewCumValues, smoothedY)


def identifyPeaks(data, stepsize):
    """Identify peaks/valleys (global min/max) from the smoothed data"""
    ind = detect_peaks(data, mpd=stepsize ,valley=True , show=True)
    return (ind)


def estimate_loc(skewCumValues, window, sequence):
    """Returns the ballpark sequence of the oriC given the skewplot data, specifically the global minimum"""
    print ("Estimate the oriC..", end=' ')

    windowNum = skewCumValues.argmin()
    estLoc = windowNum * window
    print ("\n\tMin point : ", estLoc)
    print ("\n\tidx from  : ", estLoc)
    print ("\tidx to    : ", estLoc+700)
    estSeq = sequence[estLoc-250:estLoc+700]

    print ("DONE!")
    return (estSeq)


def rev_comp(seq):
    """Returns the reverse complement of a given sequence"""
    function = {"A":"T", "T":"A", "C":"G", "G":"C"}
    output = [function[x] for x in seq[::-1]]
    return "".join(output)


def kmer_mismatches(kmer, d):
    """Returns all k-mers that are within d mismatches of the given k-mer"""
    mismatches = [kmer]
    alt_bases = {'A':'CGT', 'C':'AGT', 'G':'ACT', 'T':'ACG'}
    for dist in xrange(1, d+1):
        for change_indices in combinations(xrange(len(kmer)), dist):
            for substitutions in product(*[alt_bases[kmer[i]] for i in change_indices]):
                new_mistmatch = list(kmer)
                for idx, sub in izip(change_indices, substitutions):
                    new_mistmatch[idx] = sub
                mismatches.append(''.join(new_mistmatch))
    return mismatches


def identify_DnaA_boxes(seq, k, d):
    """Identifies DnaA boxes by solving the most frequent works problem with mismatches"""
    print ("\nIdentifying Unknown DnaA boxes...", end=' ')

    kmer_freq = defaultdict(int)
    for i in xrange(len(seq)-k+1):
        kmer_freq[seq[i:i+k]] += 1
        kmer_freq[rev_comp(seq[i:i+k])] += 1

    mismatch_count = defaultdict(int)
    for kmer, freq in kmer_freq.iteritems():
        for mismatch in kmer_mismatches(kmer, d):
            mismatch_count[mismatch] += freq

    max_count = max(mismatch_count.values())

    print ("DONE!")
    return ([kmer for kmer, count in mismatch_count.iteritems() if count == max_count])

def identify_known_DnaA_boxes(seq, pattern, d):
    """Identifies DnaA boxes by solving the most frequent works problem with mismatches given a specific pattern"""
    print ("\nIdentifying Known DnaA boxes...", end=' ')

    revPattern = rev_comp(pattern)

    approx_match = []
    for i in xrange(len(seq)-len(pattern)+1):
        mismatch_count = 0
        for j in xrange(len(pattern)):
            if seq[i:i+len(pattern)][j] != pattern[j]:
                mismatch_count += 1

        if mismatch_count <= d:
            approx_match.append(str(i))

    for i in xrange(len(seq)-len(revPattern)+1):
        mismatch_count = 0
        for j in xrange(len(revPattern)):
            if seq[i:i+len(revPattern)][j] != revPattern[j]:
                mismatch_count += 1

        if mismatch_count <= d:
            approx_match.append(str(i))

    print ("DONE!")
    return approx_match



def main():
    start = time.clock()

    fastaFile  = 'ecoli_genome.fasta'
    # fastaFile = 'vibrio_cholerae_chr1.fasta'
    # fastaFile = 'vibrio_cholerae_chr2.fasta'
    # fastaFile = 'Thermus_thermophilus_HB27.fasta'
    # fastaFile = 'thermotoga.fasta'

    stepsize = 10000 # 10000
    kmerLen  = 9
    mismatch = 1

    seqDict           = read_FASTA(fastaFile)
    sequence          = ''.join(seqDict.values())
    skewArray, window = calc_GC_skew(sequence)
    skewCumValues, smoothValues = plot_skew(skewArray)
    identifiedPeaks = identifyPeaks(smoothValues, stepsize)
    estSeq          = estimate_loc(skewCumValues, window, sequence)

    if fastaFile == 'ecoli_genome.fasta':           # http://tubic.tju.edu.cn/Ori-Finder/out/152186321123.html
        oriLoc = sequence[3925634:3926011]
            DnaAPattern = 'TTATCCACA'
    if fastaFile == 'vibrio_cholerae_chr1.fasta':   # No results -- Nothing significant
        oriLoc = sequence[2555404:2555904]
        DnaAPattern = 'ATGATCAAG'
    if fastaFile == 'vibrio_cholerae_chr2.fasta':   # http://tubic.tju.edu.cn/Ori-Finder/out/152288434166.html
        oriLoc = sequence[773980:774440]
        DnaAPattern = 'ATGATCAAG'
    if fastaFile == 'thermotoga.fasta':             # http://tubic.tju.edu.cn/Ori-Finder/out/15228852792.html
        oriLoc = sequence[788317:788661]
        DnaAPattern = 'CCTACCACC'
    # else:
        # print("There is no known data for this file!")


    print("\tBallpark location of oriC:  ( length: ", len(estSeq), ")\n",textwrap.fill(estSeq, initial_indent='\t', subsequent_indent='\t'))

    try: oriLoc
    except NameError:
      print ("\n\tThere is no experimental data for this sequence! Cannot give an actual location of the oriC!\n")
    else:
        print("\n\tActual location by Ori-Finder:\n", textwrap.fill(oriLoc, initial_indent='\t', subsequent_indent='\t'))


    ########################################
    #### DnaA Box has NOT been verified ####
    ########################################
    DnaA_boxes = identify_DnaA_boxes(estSeq, kmerLen, mismatch)
    highlight = "|".join(DnaA_boxes)
    highlightSeq = re.sub('(' + highlight + ')' , Fore.RED + r'\1' + Fore.RESET, estSeq)
    print ("\t# of DNA boxes : ", len(DnaA_boxes))
    print ("\tDnaA boxes     : ",' '.join(DnaA_boxes))
    print ("\tResults        : ", highlightSeq)


    ########################################
    ###### DnaA Box has been verified ######
    ########################################
    # Has the oriC previously been experimentally identified? If so, search for that specific pattern.
    try: DnaAPattern
    except NameError:
      print ("\nIdentifying Known DnaA boxes...\n\tThere is no experimental data for this sequence!\n")
    else:
        known_boxes = identify_known_DnaA_boxes(estSeq, DnaAPattern, mismatch)
        boxes = []
        [boxes.append(estSeq[int(box):int(box)+kmerLen]) for box in known_boxes]
        highlight = "|".join(boxes)
        highlightSeq = re.sub('(' + highlight + ')' , Fore.RED + r'\1' + Fore.RESET, estSeq)
        print ("\t# of DNA boxes : ", len(known_boxes))
        print ("\tDnaA boxes     : ",' '.join(boxes))
        print ("\tResults        : ", highlightSeq)


    print ("\nFinish:",time.clock() - start,"seconds")

if __name__ == '__main__':
    main()