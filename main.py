# HORIZONTAL GENE TRANSFER

# import argparse

if __name__ == '__main__':

    identity_file = 'input_files/alignment.identity'
    taxonomy_file = 'input_files/taxonomy.csv'
    fh = open(identity_file)
    # number of all sequences
    seq_number = fh.readline()

    # get all ids in one array
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
            pair = [ids[i], percentage]
