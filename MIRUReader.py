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


# Function that corrects the mode() function where it does not always return statistical error
def custom_mode(List):
    counts = Counter(List)
    max_count = max(counts.values())

    modes = [key for key, value in counts.items() if value == max_count]

    if len(modes) == 1:
        return modes[0]
    else:
        raise statistics.StatisticsError

# Function to extract multiples modes
def modes(List):
    counts = Counter(List)
    max_count = max(counts.values())

    modes = [key for key, value in counts.items() if value == max_count]
    return modes


# Function to determine repeat number based on total number of mismatches in primer sequence
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

script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
MIRU_table = script_dir + "/MIRU_table"
MIRU_table_0580 = script_dir + "/MIRU_table_0580"
MIRU_primers = script_dir + "/MIRU_primers"

parser = argparse.ArgumentParser()

main_group = parser.add_argument_group('Main options')

main_group.add_argument('-r', '--reads', required=True, help='input reads file in fastq/fasta format (required)')
main_group.add_argument('-p', '--prefix', required=True, help='sample ID (required)')
main_group.add_argument('--table', type=str, default=MIRU_table, help='allele calling table')
main_group.add_argument('--primers', type=str, default=MIRU_primers, help='primers sequences')

optional_group = parser.add_argument_group('Optional options')

optional_group.add_argument('--amplicons', help='provide output from primersearch and summarize MIRU profile directly', action='store_true')
optional_group.add_argument('--mismatch', type=int, dest="mismatch", required=False, default=18, help="Allowed percent mismatch. Default: 18")
optional_group.add_argument('--nofasta', help='delete the fasta reads file generated if your reads are in fastq format', action='store_true')
optional_group.add_argument('--min_amplicons', type=int, dest='min_amplicons', required=False, default=3, help='Minimum number of amplicons required for a reliable result. Below this threshold, the program returns "Warning 1" for low coverage. Default: 3')
optional_group.add_argument('--freq', type=float, dest='freq', required=False, default=0.6, help='Minimum frequency required for a reliable result. Below this threshold, the program returns "Warning 2" for an unfixed allele. Default: 0.6')
# optional_group.add_argument('--amplicon_freq', type=int, dest='amplicon_freq', required=False, default=20, help='Number of amplicons required for reliable results with mixed alleles. Below this threshold, the program returns "Warning 2" for an unfixed allele. Default: 20')
optional_group.add_argument('--amplicon_mode', type=int, dest='amplicon_mode', required=False, default=10, help='Number of amplicons required for reliable results. Below this threshold, the program returns "Warning 3" for a possible polyclonal allele with low coverage. Default: 10')

args = parser.parse_args()

if not os.path.exists(args.reads):
    sys.exit('Error: ' + args.reads + ' is not found!')

sample_prefix = args.prefix
sample_dir = os.path.dirname(os.path.abspath(args.reads))
mismatch_allowed = args.mismatch
psearchOut = sample_dir + '/' + sample_prefix + '.' + str(args.mismatch) + '.primersearch.out'

df = pd.read_csv(MIRU_table, sep='\t')
df_0580 = pd.read_csv(MIRU_table_0580, sep='\t')
miru = []
with open(args.primers) as primerFile:
    for line in primerFile:
        miru.append(line.split()[0])

#auto detect .fastq, .fastq.gz, .fasta, .fasta.gz
#convert fastq to fasta

fastaReads = sample_dir + '/' + sample_prefix + '.fasta'
if not args.amplicons:
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
    try:
        subprocess_args = ['primersearch', '-seqall', fastaReads, '-infile', args.primers, '-mismatchpercent', str(args.mismatch), '-outfile', psearchOut]
        subprocess.call(subprocess_args)
    except OSError:
        print('OSError: primersearch command is not found.')
        sys.exit()
        
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

miru_repeats = pd.DataFrame(columns = ['sample_prefix'] + miru, index = range(1))
miru_repeats['sample_prefix'] = sample_prefix
for item in miru:
    if repeats[item] != []:
        try:
            if len(repeats[item]) < args.min_amplicons:
                repeat = f"{custom_mode(repeats[item])} (Warning 1: Low Coverage)"
            # elif repeats[item].count(mode(repeats[item])) / len(repeats[item]) <= args.freq and len(repeats[item]) <= args.amplicon_freq: ## If you need to put some minimum number of amplicon for those unfixed alleles, uncomment this line and the flag corresponded.
            elif repeats[item].count(mode(repeats[item])) / len(repeats[item]) <= args.freq:
                repeat = f"{custom_mode(repeats[item])} (Warning 2: Unfixed allele)"
            else:
                repeat = custom_mode(repeats[item])
        except statistics.StatisticsError:
            if len(repeats[item]) < args.amplicon_mode:
                repeat = f"{chooseMode(item, lookup, Counter(repeats[item]))} (Warning 3: Possible polyclonal {modes(repeats[item])}, Low Coverage)"
            else:
                repeat = f"{chooseMode(item, lookup, Counter(repeats[item]))} (Warning 4: Possible polyclonal {modes(repeats[item])})"
    else:
        repeat = "ND"

    miru_repeats[item][0] = repeat

if args.nofasta:
    if ('.fastq' in args.reads) or ('.gz' in args.reads):
        os.remove(fastaReads)

print(miru_repeats.to_csv(sep='\t', index=False, header=True))
