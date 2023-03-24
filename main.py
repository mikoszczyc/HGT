#!/usr/bin/env python3
# HORIZONTAL GENE TRANSFER

import argparse
import os
from operator import itemgetter
from pathlib import Path
# from ete3 import NCBITaxa
# import generate_small_alignment as generator
from create_objects import createTaxObj

# parser = argparse.ArgumentParser(description="Find horizontal gene transfer in given sequences.")
# parser.add_argument('-i', '--input', dest='input_file', required=True,
#                     help='Input file (json format)')
# parser.add_argument('-o', '--output', dest='output_file', default='resultHGT.csv',
#                     help='Output file consisting results of HTG analysis.')
# args = parser.parse_args()

OUT_DIR = 'output_files'

def create_alignment(filename, out_alignment=f'{OUT_DIR}/alignment.fasta', out_identity=f'{OUT_DIR}/identity.txt'):
    """
    Performing alignment for the given fasta sequences using clustalo. Two files are created: one with alignment,
    the other with identity matrix calculated after alignment.

    :param filename: file containing fasta sequences
    :param out_alignment: output file containing alignment
    :param out_identity: output file containing identity matrix of the sequences after alignment
    :return: name of the file with identity matrix (out_identity)
    """
    clustalo = 'clustalo -i ' + filename + ' -o ' + out_alignment + ' --full --distmat-out ' + \
               out_identity + ' --force --percent-id'
    # os.system(clustalo)

    return out_identity


def get_ids(aln_identity_filename):
    """
    IDs of the sequences are retrieved from the identity matrix.

    :param aln_identity_filename: file containing identity matrix of the sequences after alignment
    :return: list of IDs of the sequences
    """
    file = open(aln_identity_filename)
    file.readline()  # read first line (number of sequences)

    seq_ids = []
    for line in file:
        seq_ids.append(line.split()[0])
    file.close()

    return seq_ids


def find_possible_HGT(aln_identity_filename, out_hit_filename='hit.txt'):
    """

    :param aln_identity_filename:
    :param out_hit_filename:
    :return:
    """
    # open file with identity percentage after alignment
    alignment_file = open(aln_identity_filename)

    # number of all sequences
    seq_number = alignment_file.readline()

    # open file for saving all found hits
    hit = open(f'{OUT_DIR}/{out_hit_filename}', 'w')

    hitPercentage = 0
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
                                hit.write(f'{hitPercentage}\t{analysed_species}\t{other_species}\t{analysed_protein}\t'
                                          f'{other_protein}\t{level}:{getattr(originOrganism,level)}\n')
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
                            if ori_org_next != tmp_org_next:
                                # HIT
                                hitPercentage = similarityPercentage
                                hit.write(f'{hitPercentage}\t{analysed_species}\t{other_species}\t{analysed_protein}\t'
                                          f'{other_protein}\t{level}:{getattr(originOrganism,level)}\n')
                                found = True
                            else:
                                k = j + 1
                                break  # excluding next lvl
    alignment_file.close()
    hit.close()

    # return name of the file with found hits
    return out_hit_filename


def crossHits(hits_filename, out_results='crossedResult.csv'):
    """

    :param hits_filename:
    :param out_results:
    :return:
    """
    # save pairs to crossedResult.csv
    inputFile = open(f'{OUT_DIR}/{hits_filename}', 'r')
    outputFile = open(f'{OUT_DIR}/{out_results}', 'w')
    inputContent = inputFile.readlines()
    inputFile.close()
    foundItems = []

    for i, line in enumerate(inputContent):
        line_organism1 = line.split()[1]
        line_organism1_protein = line.split()[3]

        line_organism2 = line.split()[2]
        line_organism2_protein = line.split()[4]
        for line2 in inputContent[(i + 1):]:
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
                            lvl_diff = str(idx + 1)
                            lvl_name = taxonomy_levels[idx - 1]
                            break
                    line.strip()

                    line = lvl_diff + '\t' + lvl_name + '\t' + line + '\t' + lvl2 + '\n'

                    foundItems.append(line.split())

    sorted(foundItems, key=itemgetter(1))
    # write results to the file
    for line in foundItems:
        outputFile.write('\t'.join(line) + '\n')

    outputFile.close()

    # return name of the file with results
    return out_results


if __name__ == '__main__':

    # ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()     # downloading and parsing latest database from NCBI

    taxonomy_file = 'input_files/taxonomy.csv'
    # taxonomy_file = 'taxonomy_from_id.csv'
    output_dir = OUT_DIR
    # alignment_file = open('input_files/alignment.identity')
    taxonomy = createTaxObj(taxonomy_file)
    taxonomy_levels = ['genus', 'family', 'order', 'cl', 'phylum', 'sk']
    Path(output_dir).mkdir(parents=True, exist_ok=True)  # creates directory if it didn't exist before

    fasta_file = 'input_files/proteins.fa'

    print("Performing alignment of sequences and calculating identity matrix...")
    align_identity_file = create_alignment(fasta_file)
    print("Alignment is done.")

    # get all ids in one array
    ids = get_ids(align_identity_file)

    print("Looking for horizontal gene transfer...")
    possibleHits_file = find_possible_HGT(align_identity_file)
    results_file = crossHits(possibleHits_file)
    print(f'Results can be found in {results_file} file in directory {output_dir}.')

