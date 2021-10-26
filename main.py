# HORIZONTAL GENE TRANSFER

# import argparse
import generate_small_alignment as generator

if __name__ == '__main__':

    num_seq = input("Enter number of sequences you want to keep:")
    save_f = input("Enter name of a file to save alignment to:")
    id_file = generator.generate_alignment(save_f, num_seq)

    identity_file = id_file + '.identity'
    taxonomy_file = 'input_files/taxonomy.csv'
    fh = open(identity_file)
    # number of all sequences
    seq_number = fh.readline()

    # get all ids in one array
    ids = []
    for line in fh:
        ids = [line.split()[0]]

    for line in fh:
        temp = line.split()
        species = temp[0].split('|')[0]
        protein = temp[0].split('|')[1]
        # array containing percentage identity to other sequences
        percentage_identity = temp[1:]
        identity = []
        for i, percentage in enumerate(percentage_identity):
            pair = (ids[i], percentage)
            identity = identity.append(pair)
        print(identity)
