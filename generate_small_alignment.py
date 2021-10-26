# script generates smaller alignment.identity files for testing

def generate_alignment(file, size):
    fh = open('input_files/alignment.identity')
    fh.readline()
    file = 'input_files/' + file
    f = open(file + '.identity', 'w')
    f.write(f'{size}\n')
    for i in range(0, int(size)):
        line = fh.readline()
        line = line.split()
        temp = " ".join(line[:(int(size)+1)])
        f.write(f'{temp}\n')
    f.close()
    return file
