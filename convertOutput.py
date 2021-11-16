# this script converts our output format into AZiele's

cR = open('output_files/crossedResult.csv', 'r')
taxonomy = open('input_files/taxonomy.csv', 'r')
converted = open('output_files/convCR.csv', 'w')


for line in cR:
    line = line.split()
    # lvl percentage org1 org2 gene1 gene2 genus1 genus2
    converted.write(f'{line[-1].split(sep=":")[1]}\t{line[2]}\t{line[0].split(sep="|")[0]}\t{line[1].split(sep="|")[0]}\t{line[0].split(sep="|")[1]}\t{line[1].split(sep="|")[1]}\n')
