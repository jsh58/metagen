#!/usr/bin/python

# JMG 1/2018

# Producing an html summary of the top 20 taxa from
#   centrifuge's kraken-style report.
# Version 3: include nt summary stats

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
  Node: contains taxon name, taxon number,
    list of child nodes, parent node,
    score (percent of the sample),
    count (number of reads),
    ntSeqs (number of seqs in nt with read assignments),
    ntTotal (total number of seqs in nt for this taxon).
  '''
  def __init__(self, parent, name, taxon, score, count,
      ntSeqs, ntTotal):
    self.child = []
    self.parent = parent
    self.name = name
    self.taxon = taxon
    self.score = float(score)
    self.count = int(count)
    self.ntSeqs = int(ntSeqs)
    self.ntTotal = int(ntTotal)

def printFooter(f, num, version, date, count, length):
  '''
  Print centrifuge/nt info, etc. for the output.
  '''
  f.write('<p>Top %d taxa, identified by ' % num \
    + '<a href="http://www.ccb.jhu.edu/software/centrifuge/manual.shtml">' \
    + 'Centrifuge</a>')
  if version:
    f.write(' (version %s)' % version)
  f.write(', querying the NCBI ' \
    + '<a href="https://www.ncbi.nlm.nih.gov/nucleotide">nt</a> database ' \
    + 'of %d sequences spanning %.1fGbp' % (count, length/1.0e9))
  if date:
    f.write(' (downloaded %s)' % date)
  f.write('.</p>\n')
  f.write('''<h4>Caveats:</h4>
<ul>
  <li>The "unclassified" category includes both reads that did not match
    anything in the nt database, plus those that matched a sequence with
    an unspecified or unknown taxonomy.</li>
  <p>
  <li>The value for a given taxon is the percent of all the sequence reads
    assigned to that taxon <strong>or</strong> any lower node in its tree.</li>
  <p>
  <li>Reads derived from one organism may align equally well to other
    related organisms, especially those well-represented in the nt database.
    For example, reads from a sequencing run of a <i>Homo sapiens</i> sample
    may align to <i>Pan</i> or even <i>Mus</i>. Reads that map to multiple
    taxa are counted as a fractional portion to each taxon, rather than to
    the lowest common ancestor (LCA).</li>
  <p>
  <li>Some sequences in the nt database are mislabeled, and some are
    contaminated with miscellaneous DNA (e.g. vectors). Reads may be erroneously
    assigned to sundry taxa (e.g. <strong><i>Cyprinus carpio</i></strong> [in
    class Actinopteri] and <strong><i>Ralstonia solanacearum</i></strong> [in
    class Betaproteobacteria] are frequently observed) due to matching
    contaminated reference sequences.</li>
  <p>
  <li>The depiction of the results above is based on the major levels
    (DKPCOFGS) in the NCBI
    <a href="https://www.ncbi.nlm.nih.gov/taxonomy">taxonomy</a>
    tree. Some branches in that tree skip a major level; hence, the
    columns in the above table should not be interpreted as
    corresponding to a specific taxonomic level.</li>
  <p>
  <li>The nt database was edited prior to querying:
    <ol>
      <li>short sequences (<30bp) were removed</li>
      <li>subsequences matching common Illumina adapters (Nextera, TruSeq)
        were masked</li>
      <li>Illumina adapter sequences (Nextera, TruSeq) were appended</li>
    </ol>
  </li>
</ul>
<p>Questions/concerns/comments/suggestions?
<a href="mailto:jgaspar@fas.harvard.edu">Please let us know.</a>
</p>
''')

def printLevel(f, n, level, cutoff):
  '''
  Print results for a node (if its count meets cutoff).
    Continue printing for children nodes (recursively).
  '''
  if n.count >= cutoff:  # or level == 0: # to include all children of root
    f.write('  <tr>\n' \
      + '    <td align="right">%.2f&emsp;</td>\n' % n.score \
      + '    <td>%s%s</td>\n' % (level * 2 * '&emsp;', n.name) \
      + '    <td align="right">%d</td>\n' % n.ntSeqs \
      + '    <td align="right">%d</td>\n' % n.ntTotal \
      + '  </tr>\n')
  for m in n.child:
    printLevel(f, m, level + 1, cutoff)

def printOutput(f, unclass, root, num, cutoff, version, date,
    count, length):
  '''
  Begin printing results (header and unclassified).
    Start recursive tree printing.
  '''
  f.write('''<h2>Taxonomy Analysis</h2>
<strong><font color="red" size="4">Warning:</font></strong>
<font size="4"> experimental software; not suitable for publication</font>
<p>
<table style="width:100%;border:1px solid;">
  <tr>
    <th align="right" width=10%>Percent&emsp;</th>
    <th align="left">Taxon</th>
    <th align="right">nt sequences with read assignments</th>
    <th align="right">Total nt sequences</th>
  </tr>
''')
  f.write('  <tr>\n' \
    + '    <td align="right">%.2f&emsp;</td>\n' % unclass[0] \
    + '    <td>unclassified</td>\n' \
    + '    <td align="right">-</td>\n' \
    + '    <td align="right">%d</td>\n' % unclass[1] \
    + '  </tr>\n')

  # rearrange children of root
  node = []
  for i in range(len(root.child))[::-1]:
    if root.child[i].taxon in ['12908', '28384']:
      node.append(root.child.pop(i))
  for n in node[::-1]:
    root.child.append(n)

  # print tree
  for n in root.child:
    printLevel(f, n, 0, cutoff)
  f.write('</table>\n')

  printFooter(f, num, version, date, count, length)

def findCutoff(score, x):
  '''
  Determine threshold for top x scores.
  '''
  if x >= len(score):
    return 0
  res = sorted(score, reverse=True)
  return res[x-1]

def checkTree(node, taxon):
  '''
  Search node and its children for given taxon.
  '''
  if node.taxon == taxon:
    return node
  for n in node.child:
    parent = checkTree(n, taxon)
    if parent:
      return parent
  return None

def loadScores(f, d):
  '''
  Create taxonomic tree (including scores) from a
    Centrifuge report.
  '''
  rank = 'DKPCOFGS'
  unclass = 0.0  # 'unclassified' score
  root = Node(None, 'root', '1', -1, -1, -1, -1) # root of tree
  temp = root    # pointer to previous node
  score = []     # list of scores (read counts)

  for line in f:
    spl = line.split('\t')
    if len(spl) < 7:
      sys.stderr.write('Error! Improperly formatted ' \
        + 'centrifuge-kreport file\n')
      sys.exit(-1)

    # skip non-canonical levels
    if spl[3] == '-':
      if spl[4] not in ['12908', '28384']:
        # make exception for these '-' taxa
        #   (see loadTax(), below)
        continue
    elif spl[3] == 'U':
      # save unclassified value automatically
      unclass = (float(spl[0]), d['0'][1], d['0'][2])
      continue
    elif spl[3] not in rank:
      sys.stderr.write('Warning! Unknown taxonomic rank:' \
        + ' %s\n' % spl[3] + '  ' + line)
      continue

    # find parent node in hierarchy
    if spl[4] not in d:
      sys.stderr.write('Warning! Unknown taxon: %s\n' % spl[4])
      continue
    parent = None
    # check current branch first
    while temp != None:
      if temp.taxon == d[spl[4]][0]:
        parent = temp
        break
      temp = temp.parent
    # if not found, check whole tree from root
    if parent == None:
      parent = checkTree(root, d[spl[4]][0])
    if parent == None:
      sys.stderr.write('Warning! Cannot find parent for ' \
        + 'taxon %s\n' % spl[4])
      continue

    # save to tree
    name = spl[6].strip()
    if spl[3] in 'GS':
      name = '<i>' + name + '</i>'  # italicize genus/species
    n = Node(parent, name, spl[4], spl[0], spl[1], spl[5],
      d[spl[4]][1])
    parent.child.append(n)
    temp = n

    # save score
    score.append(int(spl[1]))

  return unclass, root, score

def findParent(d, taxon):
  '''
  Find parent taxon (recursively).
  '''
  if taxon not in d or d[taxon][0] not in d:
    return None
  if d[ d[taxon][0] ][1]:
    return d[taxon][0]
  return findParent(d, d[taxon][0])

def loadTax(f):
  '''
  Load parents of each taxon, keeping only
    canonical (DKPCOFGS) ones.
    (include 0 -> unclassified
             1 -> 'root'
         12908 -> 'unclassified sequences'
         28384 -> 'other sequences')
  '''
  # load immediate parent taxa and counts to dict
  count = length = 0
  temp = {}
  for line in f:
    spl = line.split('|')
    if len(spl) < 3:
      sys.stderr.write('Error! Improperly formatted tree file\n')
      sys.exit(-1)

    level = False
    if spl[2].strip() in ['superkingdom', 'kingdom', 'phylum', \
        'class', 'order', 'family', 'genus', 'species'] \
        or spl[0].strip() in ['0', '1', '12908', '28384']:
      level = True
    temp[spl[0].strip()] = (spl[1].strip(), level, int(spl[3]),
      int(spl[4]))

  # save canonical parent taxa (and stats) to dict
  d = {}
  for taxon in temp:
    # for root nodes, save counts only (parents are 'None')
    if taxon in ['0', '1']:
      d[taxon] = (None, temp[taxon][2], temp[taxon][3])
      count += temp[taxon][2]
      length += temp[taxon][3]
      continue
    # find canonical parent
    if temp[taxon][1]:
      parent = findParent(temp, taxon)
      if parent:
        d[taxon] = (parent, temp[taxon][2], temp[taxon][3])
      else:
        sys.stderr.write('Warning! Cannot find parent ' \
          + 'of taxon %s\n' % taxon)

  return d, count, length

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '<kreport>  <taxTree>  <out> \ \n' \
      + '    [<num>]  [<version>  <date>]\n' \
      + '  <num>      Number of taxa to print (def. 20)\n' \
      + '  <version>  Version of centrifuge\n' \
      + '  <date>     Date of nt download\n')
    sys.exit(-1)

  # load tax tree
  fTax = openRead(args[1])
  d, count, length = loadTax(fTax)
  if fTax != sys.stdin:
    fTax.close()

  # load scores and create taxonomic tree
  fIn = openRead(args[0])
  unclass, root, score = loadScores(fIn, d)
  if fIn != sys.stdin:
    fIn.close()

  # find cutoff score for top n taxa
  num = 20
  if len(args) > 3:
    num = int(args[3])
  cutoff = findCutoff(score, num)

  # load centrifuge version, date of nt download
  version = date = ''
  if len(args) > 5:
    version = args[4]
    date = args[5]

  # print output
  fOut = openWrite(args[2])
  printOutput(fOut, unclass, root, num, cutoff, version, date,
    count, length)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
