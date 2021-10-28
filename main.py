# HORIZONTAL GENE TRANSFER

# import argparse
from operator import itemgetter
from pathlib import Path
import generate_small_alignment as generator
import create_objects

if __name__ == '__main__':

    num_seq = input("Enter number of sequences you want to keep:")
    # save_f = input("Enter name of a file to save alignment to:")
    save_f = 'test'
    id_file = generator.generate_alignment(save_f, num_seq)

    identity_file = id_file + '.identity'
    taxonomy_file = 'input_files/taxonomy.csv'
    output_dir = 'output_files'
    fh = open(identity_file)

    Path(output_dir).mkdir(parents=True, exist_ok=True)  # creates directory if it didn't exist before

    sf = open(f'{output_dir}/out.txt', 'w')

    # number of all sequences
    seq_number = fh.readline()

    # get all ids in one array
    ids = []
    for line in fh:
        ids.append(line.split()[0])
    fh.seek(0)  # go back to beginning of the file
    fh.readline()  # read first line (number of sequences)

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

        identity.sort(key=itemgetter(1), reverse=True)  # sortowanie wg % w krotce (DESC)
        sf.write(f'{species}|{protein}\n')
        sf.write('\t'.join('{} - {}'.format(x[0], x[1]) for x in identity))
        sf.write('\n\n')

    # ----------------------------------------------------------------
    # Load taxonomy data:
    taxonomy = create_objects.createTaxObj(taxonomy_file)

    print(taxonomy)
