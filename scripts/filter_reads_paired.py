import pandas as pd
import pysam
from collections import Counter
import numpy as np


def createTag(d):
    return ''.join([''.join(key) + str(d[key]) + ';' for key in d.keys()])[:-1]
def parseSCTag(read):
    tag = read.get_tag('SC')
    splittag = tag.split(';')
    specific_conversions = {}
    for c in splittag:
        specific_conversions[(c[0], c[1])] = np.int_(c[2:])
    return specific_conversions
def parseTCTag(read):
    tag = read.get_tag('TC')
    splittag = tag.split(';')
    total_content = {}
    for c in splittag:
        total_content[c[0]] = np.int_(c[1:])
    return total_content
def modifyLocationTag(read, snp_pos):
    if read.get_tag('ST') == '+':
        convs = set(read.get_tag('TL'))
    else:
        convs = set(read.get_tag('AL'))
    new_convs = np.array(list(convs - snp_pos))
    return new_convs
def modifySCTag(read, new_convs):
    SC_dict = parseSCTag(read)
    if read.get_tag('ST') == '+':
        SC_dict[('t', 'C')] = len(set(new_convs) - set([0]))
    else:
        SC_dict[('a', 'G')] = len(set(new_convs) - set([0]))
    return SC_dict
def concatConvLocs(new_convs_read, new_convs_mate):
    conv_locs_concat = list(set(new_convs_read).union(set(new_convs_mate)) - set([0]))
    if len(conv_locs_concat) == 0:
        conv_locs_concat.append(0)
    conv_locs_concat = [int(loc) for loc in conv_locs_concat]
    return conv_locs_concat
def modifyTagsPaired(read, mate,snp_pos ):
    TC_dict_read = parseTCTag(read)
    TC_dict_mate = parseTCTag(mate)
    TC_dict_both = dict(Counter(TC_dict_read) + Counter(TC_dict_mate))
    shared_locs = set(read.get_reference_positions()).intersection(set(mate.get_reference_positions()))
    read_refpos = pd.Series(read.get_reference_positions())
    idx = read_refpos[[pos in shared_locs for pos in read_refpos]].index.values
    overlap = ''.join([read.get_reference_sequence()[i] for i in idx]).lower()
    total_content = {'a' : 0, 'c' : 0, 'g' : 0, 't' : 0}
    for base in total_content.keys():
        total_content[base] += overlap.count(base)
    TC_dict_both = Counter(TC_dict_both)
    TC_dict_both.subtract(Counter(total_content))
    TC_dict_both = dict(TC_dict_both)
    SC_dict_read = convInReadPaired(read, snp_list=snp_pos)
    SC_dict_mate = convInReadPaired(mate, shared_locs = shared_locs, snp_list = snp_pos)
    SC_dict_read = Counter(SC_dict_read)
    SC_dict_read.update(Counter(SC_dict_mate))
    SC_dict_both = dict(SC_dict_read)
    return TC_dict_both, SC_dict_both
def convInReadPaired(read, qual = 20, shared_locs = [], snp_list = []):
    specific_conversions = {}
    specific_conversions[('c', 'A')] = 0
    specific_conversions[('g', 'A')] = 0
    specific_conversions[('t', 'A')] = 0
    specific_conversions[('a', 'C')] = 0
    specific_conversions[('g', 'C')] = 0
    specific_conversions[('t', 'C')] = 0
    specific_conversions[('a', 'G')] = 0
    specific_conversions[('c', 'G')] = 0
    specific_conversions[('t', 'G')] = 0
    specific_conversions[('a', 'T')] = 0
    specific_conversions[('c', 'T')] = 0
    specific_conversions[('g', 'T')] = 0
    specific_conversions[('a', 'N')] = 0
    specific_conversions[('c', 'N')] = 0
    specific_conversions[('g', 'N')] = 0
    specific_conversions[('t', 'N')] = 0
    avoid_list = list(set(shared_locs).union(snp_list))
    for pair in read.get_aligned_pairs(with_seq=True):
        try:
            if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                if str(pair[2]).islower() and not read.query_qualities[pair[0]] < qual and pair[1] not in avoid_list:
                    specific_conversions[(pair[2],read.seq[pair[0]])] += 1
        except (UnicodeDecodeError, KeyError):
            continue
    return specific_conversions

if __name__ == "__main__":
    import sys
    bamfilename = sys.argv[1]
    output = sys.argv[2]
    posfile = sys.argv[3]
    read_dict = {}
    conv_dict = {}
    bamfile = pysam.AlignmentFile(bamfilename, 'rb')
    mod_bamfile = pysam.AlignmentFile(output, mode='wb', template=bamfile)
    snp_positions = pd.read_csv(posfile, index_col=0, header=0)
    for read in bamfile.fetch():
      read_snp_positions = set(snp_positions[snp_positions['chrom'] == read.reference_name]['pos2'][snp_positions[snp_positions['chrom'] == read.reference_name]['pos2'].between(read.reference_start,read.reference_end)])
      if read.is_paired:
              if read.query_name in read_dict:
                  conv_locs = concatConvLocs(modifyLocationTag(read, read_snp_positions), conv_dict[read.query_name])
                  TC_tag, SC_tag = modifyTagsPaired(read, read_dict[read.query_name], read_snp_positions)
                  read.set_tag('SC', createTag(SC_tag))
                  read.set_tag('TC', createTag(TC_tag))
                  if read.get_tag('ST') == '+':
                      read.set_tag('TL', conv_locs)
                  else:
                      read.set_tag('AL', conv_locs)
                  mod_bamfile.write(read)
                  del read_dict[read.query_name]
                  del conv_dict[read.query_name]
              else:
                  read_dict[read.query_name] = read
                  conv_dict[read.query_name] = modifyLocationTag(read, read_snp_positions)
