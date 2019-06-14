#Copyright 2019 NUS pathogen genomics
#Written by Cheng Yee Tang (chengyee.tang@nus.edu.sg)

import os
import sys
import gzip
import argparse
import pandas as pd
import statistics
import subprocess
from statistics import mode
from collections import Counter


#function to determine repeat number based on total number of mismatches in primer sequence
def chooseMode(name, table, CounterList):
    maxcount = max(CounterList.values())
    repeatToCheck = []
    for k, v in CounterList.items():
        if v == maxcount:
            repeatToCheck.append(k)
    x = 0
    for i, j in table.items():
        if name in i:
            x += 1
    mismatchDict = {}
    for rp in repeatToCheck:
        mismatchDict[rp] = 0
    for i in range(x):
        string = name + '_' + str(i+1)
        if table[string][1] in repeatToCheck:
            mismatchDict[table[string][1]] += table[string][0]
    checklist2 = []
    for m, n in mismatchDict.items():
        checklist2.append(n)
    duplicates = ''
    for item in checklist2:
        if checklist2.count(item) > 1:
            duplicates = 'yes'
    finalMode = ''
    if duplicates == 'yes':
        finalMode = '/'.join(str(r) for min_value in (min(mismatchDict.values()),) for r in mismatchDict if mismatchDict[r]==min_value)
    else:
        finalMode = min(mismatchDict.keys(), key=(lambda k: mismatchDict[k]))
    return finalMode    

'''
Main function
'''

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reads', required=True, help='input reads file in fastq/fasta format (required)')
parser.add_argument('-p', '--prefix', required=True, help='sample ID (required)')
optional_group = parser.add_argument_group('Optional argument')
optional_group.add_argument('--amplicons', help='provide output from primersearch and summarize MIRU profile directly', action='store_true')
optional_group.add_argument('--details', help='for inspection', action='store_true')
optional_group.add_argument('--nofasta', help='delete the fasta reads file generated if your reads are in fastq format', action='store_true')
args = parser.parse_args()

if not os.path.exists(args.reads):
    sys.exit('Error: ' + args.reads + ' is not found!')

sample_prefix = args.prefix
sample_dir = os.path.dirname(os.path.abspath(args.reads))
mismatch_allowed = 18
psearchOut = sample_dir + '/' + sample_prefix + '.' + str(mismatch_allowed) + '.primersearch.out'
script_dir = os.path.dirname(sys.argv[0])
MIRU_table = script_dir + "/MIRU_table"
MIRU_table_0580 = script_dir + "/MIRU_table_0580"
MIRU_primers = script_dir + "/MIRU_primers"

df = pd.read_csv(MIRU_table, sep='\t')
df_0580 = pd.read_csv(MIRU_table_0580, sep='\t')
miru = ['0154','0424','0577','0580','0802','0960','1644','1955','2059','2163b','2165','2347','2401','2461','2531','2687','2996','3007','3171','3192','3690','4052','4156','4348']

#auto detect .fastq, .fastq.gz, .fasta, .fasta.gz
#convert fastq to fasta

fastaReads = sample_dir + '/' + sample_prefix + '.fasta'
if '.fastq' in args.reads:
    if '.gz' in args.reads:
        tmpH = open(fastaReads, 'w')
        p1 = subprocess.Popen(['zcat', args.reads], stdout=subprocess.PIPE)
        subprocess_args1 = ['sed', '-n', '1~4s/^@/>/p;2~4p']
        subprocess.call(subprocess_args1, stdin=p1.stdout, stdout=tmpH)
        tmpH.close()
    else:
        tmpH = open(fastaReads, 'w')
        subprocess_args1 = ['sed', '-n', '1~4s/^@/>/p;2~4p', args.reads]
        subprocess.call(subprocess_args1, stdout=tmpH)
        tmpH.close()
elif '.fasta' in args.reads:
    if '.gz' in args.reads:
        with open(fastaReads, 'w') as f:
            for line in gzip.open(args.reads, 'rb').readlines():
                f.write(line)
    else:
        fastaReads = args.reads

if not args.amplicons:
    subprocess_args = ['primersearch', '-seqall', fastaReads, '-infile', MIRU_primers, '-mismatchpercent', str(mismatch_allowed), '-outfile', psearchOut]
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
        elif (line.startswith('Amplimer') and len(line) < 12):
            col = line.split(' ')
            primerID = loci + '_' + str(col[1])
            lookup.setdefault(primerID, [])
            mm = 0
        elif 'mismatches' in line:
            mm += int(line.partition('with ')[2].rstrip(' mismatches'))
        elif 'Amplimer length' in line:
            field = line.split(':')
            amplicon = int(field[1].strip(' ').rstrip(' bp'))
            lookup.setdefault(primerID).append(mm)
            if amplicon > 1828:
                lookup.setdefault(primerID).append('NA')
            elif loci == '0580':
                if amplicon > df_0580[loci][25]:
                    lookup.setdefault(primerID).append('NA')
                else:
                    for i in range(26):
                        if amplicon < df_0580[loci][i]:
                            if i != 0:
                                first = df_0580[loci][i-1]
                                second = df_0580[loci][i]
                                if abs(amplicon - first) > abs(amplicon - second):
                                    repeats.setdefault(loci).append(df_0580['No.'][i])
                                    lookup.setdefault(primerID).append(df_0580['No.'][i])
                                    break
                                else:
                                    repeats.setdefault(loci).append(df_0580['No.'][i-1])
                                    lookup.setdefault(primerID).append(df_0580['No.'][i-1])
                                    break
                            else:
                                repeats.setdefault(loci).append(0)
                                lookup.setdefault(primerID).append(0)
                        
            else:
                if amplicon > df[loci][15]:
                    lookup.setdefault(primerID).append('NA')
                else:
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

if args.details:
    for key, value in lookup.items():
        #example: lookup = {'0154_1':[2,4]} total no. of mismatches, repeat number
        print(key, value)
    for item in miru:
        #array that used to determine repeat number
        print(Counter(repeats[item]))

miru_repeats = pd.DataFrame(columns = ['sample_prefix'] + miru, index = range(1))
miru_repeats['sample_prefix'] = sample_prefix
for item in miru:
    if repeats[item] != []:
        try:
            repeat = mode(repeats[item])
            miru_repeats[item][0] = repeat
        except statistics.StatisticsError:
            repeat = chooseMode(item, lookup, Counter(repeats[item]))
            miru_repeats[item][0] = repeat
    else:
        miru_repeats[item][0] = "ND"

if args.nofasta:
    if ('.fastq' in args.reads) or ('.fasta.gz' in args.reads):
        os.remove(fastaReads)

print(miru_repeats.to_csv(sep='\t', index=False, header=True))
