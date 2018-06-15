#!/usr/bin/python

# JMG 1/2018

# Producing an html summary of the top N taxa (def. 20)
#   from centrifuge's kraken-style report.

import sys
import gzip
import math
import scipy.stats

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
    nt90 (number of seqs in nt accounting for 90% of read assignments),
    ntTotal (total number of seqs in nt for this taxon).
  '''
  def __init__(self, parent, name, taxon, score, count,
      nt90, ntTotal):
    self.child = []
    self.parent = parent
    self.name = name
    self.taxon = taxon
    self.score = float(score)
    self.count = int(count)
    self.nt90 = int(nt90)
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
  f.write(',\n  querying the NCBI ' \
    + '<a href="https://www.ncbi.nlm.nih.gov/nucleotide">nt</a>\n  database' \
    + ' of %d sequences spanning %.1fGbp' % (count, length/1.0e9))
  if date:
    f.write(' (downloaded %s)' % date)
  f.write('.</p>\n')
  f.write('''<h4>Notes:</h4>
<ul>
  <li>The <strong>Percent</strong> value for a given taxon is the percent of
    all the sequence reads assigned to that taxon <strong>or</strong> any lower
    node in its tree.</li>
  <p>
  <li>The <strong>Enriched</strong> column provides information about the
    relative abundance of each taxon:
    <ul>
      <li><font style="color:green">&#10003;</font> = enriched (either the
        measured abundance is significantly more than expected, given the nt
        representations of the taxon and its parent [<i>p</i>-value &le; 0.05],
        <strong>or</strong> the taxon accounts for &ge; 90% of its parent's read
        assignments)</li>
      <li><font style="color:red">&#10008;</font> = not enriched (neither of the
        above criteria is met)</li>
      <li><font>&#9472;</font> = not measured (for taxa with an unenriched
        ancestor node, as well as unclassified sequences, "other sequences"
        [e.g. adapters], and top level taxa [e.g. Eukaryota])</li>
    </ul>
  </li>
  <p>
  <li>The <strong>nt90</strong> column indicates if the read assignments to the
    given taxon are due to a limited number of nt sequences:
    <ul>
      <li><font>&#10071;</font> = <strong>warning!</strong> nt90 &lt; 5
        (the nt90 value is the number of nt sequences that collectively account
        for 90% of the read assignments to the taxon; for non-leaf nodes whose
        direct read assignments do not reach this level, the given nt90 value is
        the sum of the nt90 values for the subset of nodes [children and itself]
        that do account for 90%; a generic assignment given by centrifuge
        [e.g. "species"] is counted as 5 nt sequences; a low nt90 value may
        indicate a mislabeled nt sequence or a dearth of nt sequences for this
        taxon)</li>
      <li><font>&#9472;</font> = not measured (for unclassified sequences and
        "other sequences" [e.g. adapters])</li>
    </ul>
  </li>
  <p>
  <li>The <strong>Total nt sequences</strong> value is the number of
    sequences in the nt database that are labeled with the given taxon
    <strong>or</strong> any lower node in its tree.</li>
  <p>
  <li>The "unclassified" category includes both reads that did not match
    anything in the nt database, plus those that matched a sequence with
    an unspecified or unknown taxonomy.</li>
</ul>
<h4>Caveats:</h4>
<ul>
  <li>Reads derived from one organism may align equally well to other
    related organisms, especially those well-represented in the nt database.
    For example, reads from a sequencing run of a <i>Homo sapiens</i> sample
    may align to <i>Pan</i> or even <i>Mus</i>. Reads that map to multiple
    taxa are counted as a fractional portion to each taxon, rather than to
    the lowest common ancestor (LCA).</li>
  <p>
  <li>Some sequences in the nt database are mislabeled, and some are
    contaminated with miscellaneous DNA (e.g. vectors). Reads may be erroneously
    assigned to sundry taxa (e.g. <strong><i>Cyprinus carpio</i></strong>
    [in class Actinopteri] is frequently observed) due to matching
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

def calcPval(count, n, p1, p2):
  '''
  Calculate 1-sample proportion p-value
    (one-sided, testing p_hat > p0).
  '''
  p0 = p1 / float(p2)  # null p-value (based on nt counts)
  p_hat = count / float(n)
  div = math.sqrt( p0 * (1-p0) / n )
  if div:
    z = (p_hat - p0) / div
  if not div or z > 100:
    z = 100
  return scipy.stats.norm.sf(z), p_hat

def printLevel(f, n, level, cutoff, signif, nt90Bool):
  '''
  Print results for a node (if its count meets cutoff).
    Continue printing for children nodes (recursively).
  '''
  if n.count >= cutoff:  # or level == 0: # to include all children of root

    # do not include significance/nt90 results for 'other' taxa
    if n.taxon in ['12908', '28384']:
      signif = False
      nt90Bool = False

    # determine significance (p-value <= 0.05 or prop >= 0.9)
    if level > 0 and signif:
      pval, p_hat = calcPval(n.count, n.parent.count, \
          n.ntTotal, n.parent.ntTotal)
      if pval > 0.05 and p_hat < 0.9:
        signif = False
        sigRes = '    <td align="center" style="color:red">&#10008;</td>\n'
      else:
        sigRes = '    <td align="center" style="color:green">&#10003;</td>\n'
    else:
      sigRes = '    <td align="center">&#9472;</td>\n'

    # add nt90 value
    if not nt90Bool:
      nt90 = '    <td align="center">&#9472;</td>\n'
    elif n.nt90 < 5:
      nt90 = '    <td align="center">&#10071;</td>\n'
    else:
      nt90 = '    <td></td>\n'

    # write results
    f.write('  <tr>\n' \
      + '    <td align="right">%.2f&emsp;</td>\n' % n.score \
      + '    <td>%s%s</td>\n' % (level * 2 * '&emsp;', n.name) \
      + sigRes \
      + nt90 \
      + '    <td align="right">%d</td>\n' % n.ntTotal \
      + '  </tr>\n')

  for m in n.child:
    printLevel(f, m, level + 1, cutoff, signif, nt90Bool)

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
    <th align="left" width=50%>Taxon</th>
    <th align="center">Enriched</th>
    <th align="center">nt90</th>
    <th align="right">Total nt sequences</th>
  </tr>
''')
  f.write('  <tr>\n' \
    + '    <td align="right">%.2f&emsp;</td>\n' % unclass[0] \
    + '    <td>unclassified</td>\n' \
    + '    <td align="center">&#9472;</td>\n' \
    + '    <td align="center">&#9472;</td>\n' \
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
    printLevel(f, n, 0, cutoff, True, True)
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
    name = spl[7].strip()
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

  # find cutoff score for top N taxa
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
