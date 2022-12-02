#!/usr/bin/env python3
# HORIZONTAL GENE TRANSFER

import argparse
import os
from operator import itemgetter
from pathlib import Path
# from ete3 import NCBITaxa
# import generate_small_alignment as generator
import create_objects

parser = argparse.ArgumentParser(description="Find horizontal gene transfer in given sequences.")
parser.add_argument('-i', '--input', dest='input_file', required=True,
                    help='Input file (json format)')
parser.add_argument('-o', '--output', dest='output_file', default='resultHGT.csv',
                    help='Output file consisting results of HTG analysis.')
args = parser.parse_args()


if __name__ == '__main__':

    # probably not gonna use NCBI
    # ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()     # downloading and parsing latest database from NCBI

    # taxonomy_file = 'input_files/taxonomy.csv'
    taxonomy_file = 'taxonomy_from_id.csv'
    output_dir = 'output_files'
    # alignment_file = open('input_files/alignment.identity')
    taxonomy = create_objects.createTaxObj(taxonomy_file)
    hitPercentage = 0
    taxonomy_levels = ['genus', 'family', 'order', 'cl', 'phylum', 'sk']
    Path(output_dir).mkdir(parents=True, exist_ok=True)  # creates directory if it didn't exist before

    fasta_file = 'input_files/proteins.fa'
    clustalo = 'clustalo -i ' + fasta_file + ' -o output_files/alignment.fasta --full --distmat-out ' \
                                             'output_files/identity.txt --force --percent-id'
    # os.system(clustalo)
    alignment_file = open('output_files/identity.txt')

    # number of all sequences
    seq_number = alignment_file.readline()

    # get all ids in one array
    ids = []
    for line in alignment_file:
        ids.append(line.split()[0])
    alignment_file.seek(0)  # go back to beginning of the file
    alignment_file.readline()  # read first line (number of sequences)
    hit = open(f'{output_dir}/hit.txt', 'w')

    for line in alignment_file:
        temp = line.split()
        analysed_species = temp[0].split('|')[0]
        analysed_protein = temp[0].split('|')[1]
        # array containing percentage identity to other sequences
        percentage_identity = temp[1:]
        identity = []
        for i, percentage in enumerate(percentage_identity):
            pair = (ids[i], float(percentage))
            identity.append(pair)  # append tuple to list

        identity.sort(key=itemgetter(1), reverse=True)  # sort by % DESC
        found = False
        originOrganism = taxonomy[analysed_species]
        k = 1
        for i, level in enumerate(taxonomy_levels):
            for j in range(k, len(identity)):
                organism = identity[j]
                other_species = organism[0].split(sep='|')[0]
                other_protein = organism[0].split(sep='|')[1]
                similarityPercentage = organism[1]
                tmpOrganism = taxonomy[other_species]
                if found is True:
                    ori_org_level = getattr(originOrganism, level)
                    tmp_org_level = getattr(tmpOrganism, level)
                    if ori_org_level != tmp_org_level:
                        if level == 'sk':
                            break
                        else:
                            nextLevel = taxonomy_levels[i+1]
                            ori_org_next = getattr(originOrganism, nextLevel)
                            tmp_org_next = getattr(tmpOrganism, nextLevel)
                            if (similarityPercentage == hitPercentage) and ori_org_next != tmp_org_next:
                                # HIT
                                hit.write(f'{hitPercentage}\t{analysed_species}\t{other_species}\t{analysed_protein}\t{other_protein}\t{level}:{getattr(originOrganism,level)}\n')
                            else:
                                k = j + 1
                                break  # excluding next lvl
                elif not found:
                    ori_org_level = getattr(originOrganism, level)
                    tmp_org_level = getattr(tmpOrganism, level)
                    if ori_org_level != tmp_org_level:
                        if level == 'sk':
                            found = True
                        else:
                            nextLevel = taxonomy_levels[i+1]
                            ori_org_next = getattr(originOrganism, nextLevel)
                            tmp_org_next = getattr(tmpOrganism, nextLevel)
                            if (similarityPercentage == hitPercentage) and ori_org_next != tmp_org_next:
                                # HIT
                                hitPercentage = similarityPercentage
                                hit.write(f'{hitPercentage}\t{analysed_species}\t{other_species}\t{analysed_protein}\t{other_protein}\t{level}:{getattr(originOrganism,level)}\n')
                                found = True
                            else:
                                k = j + 1
                                break  # excluding next lvl
    hit.close()

    # save pairs to crossedResult.csv
    inputFile = open(f'{output_dir}/hit.txt', 'r')
    outputFile = open(f'{output_dir}/crossedResult.csv', 'w')
    inputContent = inputFile.readlines()
    foundItems = []

    for i, line in enumerate(inputContent):
        line_organism1 = line.split()[1]
        line_organism1_protein = line.split()[3]

        line_organism2 = line.split()[2]
        line_organism2_protein = line.split()[4]
        for line2 in inputContent[(i+1):]:
            line2_organism1 = line2.split()[1]
            line2_organism1_protein = line2.split()[3]

            line2_organism2 = line2.split()[2]
            line2_organism2_protein = line2.split()[4]

            if line_organism1 == line2_organism2 and line_organism2 == line2_organism1:
                if line_organism1_protein == line2_organism2_protein and line_organism2_protein == line2_organism1_protein:
                    lvl2 = line2.split()[5]
                    lvl_diff = 0
                    lvl_name = ''
                    for idx, level in enumerate(taxonomy_levels):
                        org1 = taxonomy[line_organism1]
                        org2 = taxonomy[line_organism2]
                        if getattr(org1, level) == getattr(org2, level):
                            lvl_diff = str(idx+1)
                            lvl_name = taxonomy_levels[idx-1]
                            break
                    line.strip()

                    line = lvl_diff + '\t' + lvl_name + '\t' + line + '\t' + lvl2 + '\n'

                    foundItems.append(line.split())

    sorted(foundItems, key=itemgetter(1))
    for line in foundItems:
        outputFile.write('\t'.join(line) + '\n')

    inputFile.close()
    outputFile.close()
