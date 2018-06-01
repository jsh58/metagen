#!/usr/bin/python

# JMG 6/2018

# Produce a summary of sequences in nt.

import sys
import gzip

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

class Node:
  '''
  Node: contains parent node,
    count (number of sequences in nt), and
    length (length of sequences in nt).
  '''
  def __init__(self, parent, count = 0, length = 0):
    self.parent = parent
    self.count = count
    self.length = length

def printOutput(f, d):
  '''
  Print results.
  '''
  for n in sorted(d, key=int):
    f.write('\t|\t'.join(map(str, [n, d[n].parent,
      d[n].count, d[n].length])) + '\n')

def addCount(length, node, d):
  '''
  Add counts to node and parents (recursively).
  '''
  node.count += 1
  node.length += length
  if node.parent and node.parent in d:
    addCount(length, d[node.parent], d)

def saveInfo(seq, length, acc2tax, d):
  '''
  Save info about seq.
  '''
  if seq in acc2tax and acc2tax[seq] in d:
    addCount(length, d[acc2tax[seq]], d)
  else:
    d['0'].count += 1
    d['0'].length += length

def parseNT(f, acc2tax, d):
  '''
  Parse nt: save seq info to each taxon.
  '''
  total = totalLen = 0
  seq = ''
  length = 0
  for line in f:
    if line[0] == '>':
      if seq:
        saveInfo(seq, length, acc2tax, d)
        total += 1
        totalLen += length
        length = 0
      seq = line.rstrip().split(' ')[0][1:]
    else:
      length += len(line) - 1
  if seq:
    saveInfo(seq, length, acc2tax, d)
    total += 1
    totalLen += length
  return total, totalLen

def loadTax(f):
  '''
  Load parents of each taxon from tree.
  '''
  d = {}
  for line in f:
    spl = line.split('|')
    if len(spl) < 3:
      sys.stderr.write('Error! Improperly formatted tree file\n')
      sys.exit(-1)
    base = spl[0].strip()
    if base == '1':
      d[base] = Node(None)
    else:
      d[base] = Node(spl[1].strip())
  d['0'] = Node(None)  # node for seqs with unassigned taxonomy
  return d

def loadAcc(f):
  '''
  Load accession -> taxID info.
  '''
  d = {}
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 2:
      sys.stderr.write('Error! Improperly formatted acc2taxid file\n')
      sys.exit(-1)
    d[spl[0]] = spl[1]
  return d

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 4:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '<acc2taxid>  <taxTree>  \ \n' \
      + '  <fasta>  <out>\n')
    sys.exit(-1)

  # load acc2taxid
  fAcc = openRead(args[0])
  acc2tax = loadAcc(fAcc)
  if fAcc != sys.stdin:
    fAcc.close()

  # load tax tree
  fTax = openRead(args[1])
  d = loadTax(fTax)
  if fTax != sys.stdin:
    fTax.close()

  # parse nt.fa
  fIn = openRead(args[2])
  total, totalLen = parseNT(fIn, acc2tax, d)
  if fIn != sys.stdin:
    fIn.close()
  sys.stderr.write('Total seqs in %s: %d\n' % (args[2], total) \
    + '  Total length (bp): %d\n' % totalLen)

  # print output
  fOut = openWrite(args[3])
  printOutput(fOut, d)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
