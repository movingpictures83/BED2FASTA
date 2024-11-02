#!/usr/bin/env python2
"""
This script is used to extract sequences from bed file.

USAGE 
  getseqfrombed.py {OPTIONS} <.bed file|-> <reference .fasta file> 

OPTIONS

  -b/--seqerror [error file]\tSpecify the positional error profile to be used. The file should include at least 100 lines, each containing a positive number. The number at line x is the weight that an error is occured at x% position of the read. If no positional error file specified, uniform weight is assumed.

  -r/--errorrate [error rate]\tSpecify the overall error rate, a real positive number. The number of errors of each read will follow a Poisson distribution with its mean value specified by --errorrate.  Default 0 (no errors).

  -l/--readlen [read length]\tSpecify the read length. Default 75. 

  -f/--fill [seq]\tFill at the end of each read by the sequence seq, if the read is shorter than the read length. Default A (to simulate poly-A tails in RNA-Seq reads).

NOTE

  1. The input .bed file is best to sort according to chromosome names. Use - to input from STDIN.
  2. Biopython and numpy package are required.
  
  3. When applying models, we assume that all sequences are in the same length. The length information is given by the -l parameter. If the sequence length is greater than read length, nucleotides outside the read length will not be simulated for error.

HISTORY

	02/01/2013:
	  Fix a bug with no read errors generated.
	  Fix a bug with error profiles in the minus strand.
	08/25/2011:
    	  Rename makebedseq.py to getseqfrombed.py.
    	  Print results to stdout.
"""

import sys;
import pydoc;
import os;
import random;
import bisect;
import math;
import numpy;
from Bio import SeqIO;
from Bio.SeqRecord import SeqRecord;

# import argparse;
# parser=argparse.ArgumentParser('Extract sequences from bed file');
# parser.add_argument('-b','--seqerror',help='Specify the positional error profile to be used. The file should include at least 100 lines, each containing a positive number. The number at line x is the weight that an error is occured at x% position of the read. If no positional error file specified, uniform weight is assumed.');
# parser.add_argument('-r','--errorrate',type=float,default=0.0,help='Specify the overall error rate, a number between 0 and 1. Default 0 (no errors).');
# parser.add_argument('-l','--readlen',type=int,default=75,help='Specify the read length. Default 75.');
# parser.add_argument('-f','--fill',default='A',help='Fill at the end of each read by the sequence seq, if the read is shorter than the read length. Default A (to simulate poly-A tails in RNA-Seq reads).');





# analyzing parameters
posweight=[];
errrate=0.00;
readlength=75;
forcelength=False;
filledseq='A';




#python2 getseqfrombed.py -b $READERR -f A -r 0.01 -l $READLEN -  $REFERENCE

import PyPluMA
import PyIO
class BED2FASTAPlugin:
 def input(self, inputfile):
  self.parameters = PyIO.readParameters(inputfile)
 def run(self):
     pass
 def output(self, outputfile):
  bline=0;
  tbweight=0;
  for lines in open(PyPluMA.prefix()+"/"+self.parameters["errorfile"]):
        bline=bline+1;
        if bline>100:
          break;
        tbweight=float(lines.strip());
        posweight.append(tbweight);
  if len(posweight)!=100:
        sys.exit();
  forcelength=True
  fillseq = self.parameters["fill"]
  errorrate = float(self.parameters["errorrate"])
  readlength = int(self.parameters["readlen"])
  reffile = PyPluMA.prefix()+"/"+self.parameters["reffile"]
  bedfile = PyPluMA.prefix()+"/"+self.parameters["bedfile"]
  # construct weight probability for read length, if possible
  rlenweight=[];
  if len(posweight)!=0:
    kweight=0;
    for i in range(readlength):
      nfrac=i*100.0/readlength;
      lower=int(math.floor(nfrac));
      higher=int(math.ceil(nfrac));
      if higher==lower: higher=lower+1;
      #print('higher:'+str(higher)+',lower:'+str(lower));
      if higher<100:
        val=posweight[lower]*(nfrac-lower)+posweight[higher]*(higher-nfrac);
      else:
        val=posweight[99];
      kweight+=val;
      rlenweight.append(kweight);

  #bedfile=sys.argv[-2];
  #reffile=sys.argv[-1];
  #ofastafile=sys.argv[-1];

  # build reference
  seqref=SeqIO.index(reffile,'fasta');
  refkeys = list(seqref.keys())

  # read bed file, and ready for writing
  if bedfile!="-":
    fid=open(bedfile);
  else:
    fid=sys.stdin;
  ofid=open(outputfile,'w');
  #ofid=sys.stdout

  nlines=0;

  prevchr='';
  previndex='';

  for lines in fid:
    # update line counter
    nlines=nlines+1;
    # parse lines
    bedfield=lines.strip().split('\t');
    # clustering
    fieldrange=[int(bedfield[1]),int(bedfield[2])];
    # parse all exons
    exonlen=[int(x) for x in bedfield[10][:-1].split(',')];
    exonstart=[int(x)+fieldrange[0] for x in bedfield[11][:-1].split(',')];
    if bedfield[0]!=prevchr:
      prevchr=bedfield[0];
      previndex=seqref[bedfield[0]];
    # extract sequences
    thisseq=SeqRecord('');
    for i in range(len(exonlen)):
      thisseq+=previndex[exonstart[i]:(exonstart[i]+exonlen[i])];
    if forcelength:
      if sum(exonlen)<readlength:
        thisseq+=filledseq*(readlength-sum(exonlen));
    thisseq.id=bedfield[3];
    thisseq.description='';
    # mutation
    nmut=numpy.random.poisson(errrate);
    if nmut>0:
      newseq=thisseq.seq;
      for n in range(nmut):
        if len(posweight)==0:
          # uniform distrib
          modifyposition=random.choice(range(len(newseq)));
        else:
          rchosen=random.random()*kweight;
          modifyposition=bisect.bisect_right(posweight,rchosen);
        # mutate the position
        if len(newseq)>modifyposition:
          topos=random.choice('ATGC');
          while topos==newseq[modifyposition]:
            topos=random.choice('ATGC');
          print('MUTATION at position '+str(modifyposition)+','+newseq[modifyposition]+'->'+topos);
          newseq=newseq[:modifyposition]+topos+newseq[(modifyposition+1):];
      thisseq.seq=newseq;
    # reverse-complement the sequence if it is on the negative strand
    if bedfield[5]=='-':
      thisseq.seq=thisseq.seq.reverse_complement();
    # write to record
    try:
       SeqIO.write(thisseq,ofid,'fasta');
    except ValueError:
       print("Skip")
  
        

  # ofid.close();
  if bedfile!="-":
    fid.close();
    



