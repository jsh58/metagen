#!/usr/bin/python

# JMG 12/2017

# Filter a fasta file (e.g. the NCBI nt database):
#   - remove sequences of pure Ns
#   - remove sequences shorter than a specified
#       length (def. 25bp)
#   - mask sequences in a given BED file
#   - remove sequences whose headers are in a
#       given list

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

def loadBed(filename):
  '''Load BED regions to dict.'''
  d = {}
  f = openRead(filename)
  for line in f:
    if line[0] == '#': continue
    spl = line.rstrip().split('\t')
    if len(spl) < 3:
      sys.stderr.write('Error! Improperly formatted BED record:\n' + line)
      sys.exit(-1)
    if spl[0] not in d:
      d[spl[0]] = [(int(spl[1]), int(spl[2]))]
    else:
      d[spl[0]].append((int(spl[1]), int(spl[2])))
  if f != sys.stdin:
    f.close()
  return d

def parseFasta(fIn, fOut, minLen, mask, headers):
  '''
  Parse fasta file, write output on the fly.
  '''
  count = short = pureNs = xReads = total = 0
  head = ''    # header (1st space-delim token)
  read = ''    # full read (header + sequence)
  nseq = True  # sequence is pure Ns
  length = 0   # length of sequence
  inter = []   # intervals to mask

  # analyze fasta reads
  for line in fIn:
    if line[0] == '>':

      # process previous read
      if read:
        count += 1
        if length < minLen:
          short += 1
        elif nseq:
          pureNs += 1
        elif head in headers:
          xReads += 1
        else:
          fOut.write(read)
          total += 1

      # start new read
      read = line
      head = line.rstrip().split(' ')[0][1:]
      inter = []
      if head in mask:
        # load masked sites
        for m in mask[head]:
          inter += range(m[0], m[1])
        inter.sort()
      length = 0
      nseq = True

    elif read:
      # mask sequence
      while inter and length <= inter[0] < length + len(line) - 1:
        line = line[:inter[0]-length] + 'N' + line[inter[0]-length+1:]
        del[inter[0]]

      # save sequence
      read += line
      length += len(line) - 1
      if nseq and line.rstrip() != 'N' * (len(line) - 1):
        nseq = False

  if fIn != sys.stdin:
    fIn.close()

  # process last read
  if read:
    count += 1
    if length < minLen:
      short += 1
    elif nseq:
      pureNs += 1
    elif head in headers:
      xReads += 1
    else:
      fOut.write(read)
      total += 1

  if fOut != sys.stdout:
    fOut.close()

  return count, short, pureNs, xReads, total

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python filterNT2.py  <input>  <output> \ \n' \
      + '    [<minLen>]  [<BED>]  [<headers]\n')
    sys.stderr.write('  <minLen>    Minimum sequence length (def. 25bp)\n')
    sys.stderr.write('  <BED>       BED file of regions to mask\n')
    sys.stderr.write('  <headers>   File listing headers of sequences to exclude\n')
    sys.exit(-1)

  # get CL args
  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  minLen = 25
  if len(args) > 2:
    minLen = int(args[2])

  # load BED regions
  if len(args) > 3:
    mask = loadBed(args[3])

  # load headers of seqs to exclude
  headers = {}
  if len(args) > 4:
    fRead = openRead(args[4])
    for line in fRead:
      headers[line.rstrip()] = 1

  # parse fasta
  count, short, pureNs, xReads, total \
    = parseFasta(fIn, fOut, minLen, mask, headers)

  sys.stderr.write('Total fasta sequences in %s: %d\n' % (args[0], count))
  sys.stderr.write('  Shorter than %dbp: %d\n' % (minLen, short))
  sys.stderr.write('  Pure Ns: %d\n' % pureNs)
  if len(args) > 4:
    sys.stderr.write('  Excluded: %d\n' % xReads)
  sys.stderr.write('  Written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()
