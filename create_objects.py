class Tax:
    # genome_id = ''
    def __init__(self, genome_id, sk, phylum, cl, order, family, genus, species):
        self.genome_id = genome_id
        self.sk = sk
        self.phylum = phylum
        self.cl = cl
        self.order = order
        self.family = family
        self.genus = genus
        self.species = species

    def __str__(self):
        return str(self.__dict__)


def createTaxObj(file): # Creates objects from taxonomy.csv
    # Create obj list
    a = []

    # Open input file
    fh = open(file)
    fh.readline()
    for i, line in enumerate(fh):
        tmp = line.split(sep=',')
        tmp[-1] = tmp[-1][:-1] # removing new line symbol

        a.append(Tax(*tmp))

    return a


if __name__ == '__main__':
    createTaxObj('input_files/taxonomy.csv')
