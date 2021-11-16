# HORIZONTAL GENE TRANSFER

# import argparse
from operator import itemgetter
from pathlib import Path
# import generate_small_alignment as generator
import create_objects

if __name__ == '__main__':

    taxonomy_file = 'input_files/taxonomy.csv'
    output_dir = 'output_files'
    alignment_file = open('input_files/alignment.identity')
    taxonomy = create_objects.createTaxObj(taxonomy_file)
    hitPercentage = 0
    taxonomy_levels = ['family', 'order', 'cl', 'phylum', 'sk']
    Path(output_dir).mkdir(parents=True, exist_ok=True)  # creates directory if it didn't exist before

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
        analysed_genus = taxonomy[analysed_species].genus
        # array containing percentage identity to other sequences
        percentage_identity = temp[1:]
        identity = []
        for i, percentage in enumerate(percentage_identity):
            pair = (ids[i], float(percentage))
            identity.append(pair)  # append tuple to list

        identity.sort(key=itemgetter(1), reverse=True)  # sort by % DESC
        found = False
        originOrganism = taxonomy[analysed_species]

        for organism in identity:
            other_species = organism[0].split(sep='|')[0]
            other_protein = organism[0].split(sep='|')[1]
            other_genus = taxonomy[other_species].genus
            similarityPercentage = organism[1]
            if found:
                if organism[1] == hitPercentage:
                    tmpOrganism = taxonomy[other_species]
                    for i, level in enumerate(taxonomy_levels):
                        if level == 'family':
                            continue
                        if getattr(tmpOrganism, level) == getattr(originOrganism, level):
                            hitPercentage = similarityPercentage
                            hit.write(f'{i+2}\t{taxonomy_levels[i-1]}\t{hitPercentage}\t{analysed_species}\t{other_species}\t{analysed_protein}\t{other_protein}\n')
                            break
                    continue
                else:
                    break
            elif not found:
                hitPercentage = 0
                if organism[1] != 100.0:
                    tmpOrganism = taxonomy[other_species]
                    if (tmpOrganism.genus != originOrganism.genus) and (tmpOrganism.family != originOrganism.family):
                        found = True
                        for i, level in enumerate(taxonomy_levels):
                            if level == 'family':
                                continue
                            if getattr(tmpOrganism, level) == getattr(originOrganism, level):
                                hitPercentage = similarityPercentage
                                hit.write(f'{i+2}\t{taxonomy_levels[i-1]}\t{hitPercentage}\t{analysed_species}\t{other_species}\t{analysed_protein}\t{other_protein}\n')
                                break
    hit.close()

# zapis do par do pliku (crossed)
    inputFile = open(f'{output_dir}/hit.txt', 'r')
    outputFile = open(f'{output_dir}/crossedResult.csv', 'w')
    inputContent = inputFile.readlines()
    foundItems = []

    for i, line in enumerate(inputContent):
        line_organism1 = line.split(sep='\t')[3]
        line_organism1_protein = line.split(sep='\t')[5]

        line_organism2 = line.split(sep='\t')[4]
        line_organism2_protein = line.split(sep='\t')[6]
        for line2 in inputContent[(i+1):]:
            line2_organism1 = line2.split(sep='\t')[3]
            line2_organism1_protein = line2.split(sep='\t')[5]

            line2_organism2 = line2.split(sep='\t')[4]
            line2_organism2_protein = line2.split(sep='\t')[6]

            if line_organism1 == line2_organism2 and line_organism2 == line2_organism1:
                if line_organism1_protein == line2_organism2_protein and line_organism2_protein == line2_organism1_protein:
                    foundItems.append(line.split(sep='\t'))

    sorted(foundItems, key=itemgetter(3))
    for line in foundItems:
        outputFile.write('\t'.join(line))
