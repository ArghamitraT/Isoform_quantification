import sys

file = sys.argv[1]

reads = []

with open(file, 'r') as handle:
    for line in handle:
        reads.append(line.split('\t')[0])

with open('reads.txt', 'w') as handle:
    for r in reads:
        handle.write(r)
        handle.write('\n')