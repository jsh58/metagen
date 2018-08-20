#!/usr/bin/python

# JMG 6/2018

# Simulate PE reads from nt sequences for a given taxon.

import sys
import gzip
import random

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def revComp(seq):
  rc = ''
  comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  for nuc in seq[::-1]:
    if nuc.upper() not in comp:
      return None
    rc += comp[nuc.upper()]
  return rc

def printOutput(f1, f2, d, fragLen, readLen, number):
  '''
  Randomly produce reads; print results.
  '''
  acc = d.keys()
  for i in range(number):
    while True:

      # choose an nt sequence and position
      j = random.randint(0, len(acc)-1)
      pos = random.randint(0, len(d[acc[j]])-fragLen)

      # generate reads
      seqR1 = d[acc[j]][pos:pos+readLen].upper()
      seqR2 = revComp(d[acc[j]][pos+fragLen-readLen:pos+fragLen])
      if not seqR2 or not revComp(seqR1):
        # if non-N ambiguous nucleotides, try again
        continue
      if seqR1.count('N') > 2 or seqR2.count('N') > 2:
        # if more than 2 Ns, try again
        continue

      # print output
      f1.write('@read%d %s_%d-%d R1\n' % (i, acc[j], pos, pos+fragLen) \
        + seqR1 + '\n+\n' + 'I' * readLen + '\n')
      f2.write('@read%d %s_%d-%d R2\n' % (i, acc[j], pos, pos+fragLen) \
        + seqR2 + '\n+\n' + 'I' * readLen + '\n')
      break


def parseNT(f, acc, minLen):
  '''
  Parse nt: save seqs of given accessions.
  '''
  d = {}
  head = seq = ''
  save = False
  for line in f:
    if line[0] == '>':
      if save and len(seq) >= minLen:
        d[head] = seq
      head = line.rstrip().split(' ')[0][1:]
      seq = ''
      save = False
      if head in acc:
        save = True
    elif save:
      seq += line.rstrip()
  if save and len(seq) >= minLen:
    d[head] = seq
  return d

def loadAcc(f, taxon):
  '''
  Load accessions for given taxon.
  '''
  acc = {}
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 2:
      sys.stderr.write('Error! Improperly formatted acc2taxid file\n')
      sys.exit(-1)
    if spl[1] == taxon:
      acc[spl[0]] = 1
  return acc

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 8:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '''<acc2taxid>  <taxon>  <fasta>
    <len1>  <len2>  <num>  <outR1>  <outR2>
  <acc2taxid>  File listing accessions and taxonomic IDs
  <taxon>      Taxonomic ID of interest
  <fasta>      Reference fasta file (e.g. nt.fa)
  <len1>       Length of simulated DNA fragments
  <len2>       Length of simulated reads
  <num>        Number of simulated read pairs
  <outR1>      Output file for R1 reads
  <outR2>      Output file for R2 reads
''')
    sys.exit(-1)

  # load accessions for given taxon
  fAcc = openRead(args[0])
  acc = loadAcc(fAcc, args[1])
  if fAcc != sys.stdin:
    fAcc.close()
  sys.stderr.write('Accessions for taxon %s in %s: %d\n' \
    % (args[1], args[0], len(acc)))

  # save nt sequences
  fIn = openRead(args[2])
  d = parseNT(fIn, acc, int(args[3]))
  if fIn != sys.stdin:
    fIn.close()
  sys.stderr.write('Sequences loaded from %s: %d\n' % (args[2], len(d)))

  # print output
  fOut1 = openWrite(args[6])
  fOut2 = openWrite(args[7])
  printOutput(fOut1, fOut2, d, int(args[3]), \
    int(args[4]), int(args[5]))
  if fOut1 != sys.stdout:
    fOut1.close()
  if fOut2 != sys.stdout:
    fOut2.close()

if __name__ == '__main__':
  main()
