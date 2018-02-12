# JMG 2/2018

# making 1000 copies of fasta sequences,
#   80 chars per line (to match nt)

import sys

if len(sys.argv) < 3:
  sys.stderr.write('Usage: python copy1000.py  <in>  <out>\n')
  sys.exit(-1)

fIn = open(sys.argv[1], 'rU')
fOut = open(sys.argv[2], 'w')

head = ''
seq = ''
for line in fIn:
  if line[0] == '>':
    if seq:
      fOut.write(head)
      seq *= 1000
      for i in range(0, len(seq), 80):
        fOut.write(seq[i:i+80] + '\n')
      head = ''
      seq = ''
    head = line
  else:
    seq += line.rstrip()
if seq:
  fOut.write(head)
  seq *= 1000
  for i in range(0, len(seq), 80):
    fOut.write(seq[i:i+80] + '\n')

fIn.close()
fOut.close()
