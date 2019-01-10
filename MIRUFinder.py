import os
import sys
import argparse
import pandas as pd
import statistics
import subprocess
from statistics import mode
from collections import Counter


def chooseMode(name, table, mode1, mode2):
    x = 0
    for k,v in table.items():
        if name in k:
            x += 1
    mode1_mm = 0
    mode2_mm = 0
    for i in range(x):
        string = name + '_' + str(i+1)
        if table[string][1] == mode1:
            mode1_mm += table[string][0]
        elif table[string][1] == mode2:
            mode2_mm += table[string][0]
    finalMode = 0
    if mode1_mm > mode2_mm:
        finalMode = mode2
    elif mode1_mm < mode2_mm:
        finalMode = mode2
    else:
        finalMode = str(mode1) + '/' + str(mode2)
    return finalMode

'''
Main function
'''

parser = argparse.ArgumentParser()
parser.add_argument('reads', help='input reads file in fasta format')
optional_group = parser.add_argument_group('Optional argument')
optional_group.add_argument('--amplicons', help='provide output from primersearch and summarize MIRU profile directly', action='store_true')
args = parser.parse_args()

if not os.path.exists(args.reads):
    sys.exit('Error: ' + args.reads + ' is not found!')

sample_prefix = os.path.splitext(os.path.basename(args.reads))[0]
sample_dir = os.path.dirname(os.path.abspath(args.reads))
psearchOut = sample_dir + '/' + sample_prefix + '.18.primersearch.out'
script_dir = os.path.dirname(sys.argv[0])
MIRU_table = script_dir + "/MIRU_table"
MIRU_primers = script_dir + "/MIRU_primers"

df = pd.read_table(MIRU_table)
miru = ['0154','0424','0577','0580','0802','0960','1644','1955','2059','2163b','2165','2347','2401','2461','2531','2687','2996','3007','3171','3192','3690','4052','4156','4348']

if not args.amplicons:
    subprocess_args = ['primersearch', '-seqall', args.reads, '-infile', MIRU_primers, '-mismatchpercent', '18', '-outfile', psearchOut]
    subprocess.call(subprocess_args)

if not os.path.exists(psearchOut):
    sys.exit('Error: ' + psearchOut + ' is not found!')

lookup = {}
repeats = {}
with open(psearchOut, 'r') as infile:
    for line in infile.read().splitlines():
        if line.startswith('Primer'):
            col = line.split(' ')
            loci = str(col[2])
            repeats.setdefault(loci, [])
        elif line.startswith('Amplimer'):
            col = line.split(' ')
            primerID = loci + '_' + str(col[1])
            lookup.setdefault(primerID, [])
            mm = 0
        elif 'mismatches' in line:
            mm += int(line.partition('with ')[2].rstrip(' mismatches'))
        elif 'Amplimer length' in line:
            lookup.setdefault(primerID).append(mm)
            field = line.split(':')
            amplicon = int(field[1].strip(' ').rstrip(' bp'))
            for i in range(16):
                if amplicon < df[loci][i]:
                    if i != 0:
                        first = df[loci][i-1]
                        second = df[loci][i]
                        if abs(amplicon - first) > abs(amplicon - second):
                            repeats.setdefault(loci).append(i)
                            lookup.setdefault(primerID).append(i)
                            break
                        else:
                            repeats.setdefault(loci).append(i-1)
                            lookup.setdefault(primerID).append(i-1)
                            break
                    else:
                        repeats.setdefault(loci).append(0)
                        lookup.setdefault(primerID).append(0)

#example: lookup = {'0154_1':[2,4]} total no. of mismatches, repeat number
miru_repeats = pd.DataFrame(columns = ['sample_prefix'] + miru, index = range(1))
miru_repeats['sample_prefix'] = sample_prefix
for item in miru:
    if repeats[item] != []:
        try:
            repeat = mode(repeats[item])
            miru_repeats[item][0] = repeat
        except statistics.StatisticsError:
            common_mode = Counter(repeats[item]).most_common(2)
            mode1 = common_mode[0][0]
            mode2 = common_mode[1][0]
            repeat = chooseMode(item, lookup, mode1, mode2)
            miru_repeats[item][0] = repeat
    else:
        miru_repeats[item][0] = "nohit"

print(miru_repeats.to_csv(sep='\t', index=False, header=True))
