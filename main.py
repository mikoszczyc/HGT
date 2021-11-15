# HORIZONTAL GENE TRANSFER

# import argparse
from operator import itemgetter
from pathlib import Path
# import generate_small_alignment as generator
import create_objects

if __name__ == '__main__':

    # identity_file = id_file + '.identity'
    taxonomy_file = 'input_files/taxonomy.csv'
    output_dir = 'output_files'
    fh = open('input_files/alignment.identity')
    taxonomy = create_objects.createTaxObj(taxonomy_file)
    hitPercentage = 0
    levels = ['family', 'order', 'cl', 'phylum', 'sk']
    Path(output_dir).mkdir(parents=True, exist_ok=True)  # creates directory if it didn't exist before

    # number of all sequences
    seq_number = fh.readline()

    # get all ids in one array
    ids = []
    for line in fh:
        ids.append(line.split()[0])
    fh.seek(0)  # go back to beginning of the file
    fh.readline()  # read first line (number of sequences)
    hit = open(f'{output_dir}/hit.txt', 'w')

    for line in fh:
        temp = line.split()
        species = temp[0].split('|')[0]
        protein = temp[0].split('|')[1]
        # array containing percentage identity to other sequences
        percentage_identity = temp[1:]
        identity = []
        for i, percentage in enumerate(percentage_identity):
            pair = (ids[i], float(percentage))
            identity.append(pair)  # append tuple to list

        identity.sort(key=itemgetter(1), reverse=True)  # sort by % DESC
        found = False
        originOrganism = taxonomy[species]

        for organism in identity:
            if found:
                if organism[1] == hitPercentage:
                    tmpOrganism = taxonomy[organism[0].split('|')[0]]
                    for i, level in enumerate(levels):
                        if level == 'family':
                            continue
                        if getattr(tmpOrganism, level) == getattr(originOrganism, level):
                            hitPercentage = organism[1]
                            hit.write(f'{species}|{protein}\t{organism[0]}\t{organism[1]}\tl:{levels[i-1]}\n')
                            break
                    continue
                else:
                    break
            elif not found:
                hitPercentage = 0
                if organism[1] != 100.0:
                    tmpOrganism = taxonomy[organism[0].split('|')[0]]
                    if (tmpOrganism.genus != originOrganism.genus) and (tmpOrganism.family != originOrganism.family):
                        found = True
                        for i, level in enumerate(levels):
                            if level == 'family':
                                continue
                            if getattr(tmpOrganism, level) == getattr(originOrganism, level):
                                hitPercentage = organism[1]
                                hit.write(f'{species}|{protein}\t{organism[0]}\t{organism[1]}\tl:{levels[i-1]}\n')
                                break
    hit.close()

    hit = open(f'{output_dir}/hit.txt', 'r')
    result = open(f'{output_dir}/crossedResult.csv', 'w')
    hitRead = hit.readlines()
    foundItems = []

    for i, line in enumerate(hitRead):
        for line2 in hitRead[(i+1):]:
            temp = sorted(line.split()[:2])
            pair = temp[0] + temp[1]
            temp2 = sorted(line2.split()[:2])
            pair2 = temp2[0] + temp2[1]
            if (pair == pair2) and (pair not in foundItems):
                foundItems.append(pair)
                result.write(f'{line}')
