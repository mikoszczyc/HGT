# HORIZONTAL GENE TRANSFER

# import argparse
from operator import itemgetter
from pathlib import Path
# import generate_small_alignment as generator
import create_objects

if __name__ == '__main__':

    # num_seq = input("Enter number of sequences you want to keep:")
    # # save_f = input("Enter name of a file to save alignment to:")
    # save_f = 'test'
    # id_file = generator.generate_alignment(save_f, num_seq)

    # identity_file = id_file + '.identity'
    taxonomy_file = 'input_files/taxonomy.csv'
    output_dir = 'output_files'
    fh = open('input_files/alignment.identity')
    taxonomy = create_objects.createTaxObj(taxonomy_file)

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

        # if el ∉ genus ∧ el ∈ family
        for organism in identity:
            if organism[1] != 100.0:
                tmpOrganism = taxonomy[organism[0].split('|')[0]]
                if (tmpOrganism.family == originOrganism.family) and (tmpOrganism.genus != originOrganism.genus):
                    found = True
                    hit.write(f'{species}|{protein}\t{organism[0]}\t{organism[1]}\n')  # HIT!
                    break

        # if el ∉ family ∧ el ∈ order
        if not found:
            for organism in identity:
                tmpOrganism = taxonomy[organism[0].split('|')[0]]
                if (tmpOrganism.order == originOrganism.order) and (tmpOrganism.family != originOrganism.family):
                    found = True
                    hit.write(f'{species}|{protein}\t{organism[0]}\t{organism[1]}\n')  # HIT!
                    break

        # if el ∉ order ∧ el ∈ class
        if not found:
            for organism in identity:
                tmpOrganism = taxonomy[organism[0].split('|')[0]]
                if (tmpOrganism.cl == originOrganism.cl) and (tmpOrganism.order != originOrganism.order):
                    found = True
                    hit.write(f'{species}|{protein}\t{organism[0]}\t{organism[1]}\n')  # HIT!
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
