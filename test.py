from ete3 import NCBITaxa
import os
ncbi = NCBITaxa()
# ncbi.update_taxonomy_database() # wystarczy raz pobrać - polecenie aktualizuje lokalną bazę danych co zajmuje trochę czasu
file_name = 'input_files/taxIDs.txt'
file = open(file_name)

proteins = open('input_files/proteins.fa')
genome_ids = []
for line in proteins:
    if line[0] != '>':
        continue
    else:
        genome_ids.append(line.strip('>').strip().split('|')[0])

# --------------------------------------- TEST 1 -------------------------------------
hier = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]  # hierarchia jaką chcemy otrzymać w pliku - zmiana tej listy nie wystarczy, wymagana jest zmiana kodu poniżej
# file = ['2162078'] # na potrzeby testowania dla 1 ID


d = {}  # słownik zawierajacy taksonomię dla każdego ID
# zapis do pliku
f = open('taxonomy_from_id.csv', 'w')
f.write(f'genome_id\tncbi_id\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n')
# id == NCBI_ID NIE GENOME_ID
for i, id in enumerate(file):
    id = id.strip() # usuwanie znaków białych z id

    lineage = ncbi.get_lineage(id)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)
    #print(f'{ranks.keys()}\n{names.keys()}')
    #print(f'{ranks.values()}\n{names.values()}')

    d[id] = {}
    for key in ranks:
        if ranks[key] == ('superkingdom'):
            d[id]['superkingdom'] = names[key]
        elif ranks[key] == ('phylum'):
            d[id]['phylum'] = names[key]
        elif ranks[key] == 'class':
            d[id]['class'] = names[key]
        elif ranks[key] == ('order'):
            d[id]['order'] = names[key]
        elif ranks[key] == ('family'):
            d[id]['family'] = names[key]
        elif ranks[key] == ('genus'):
            d[id]['genus'] = names[key]
        elif ranks[key] == ('species'):
            d[id]['species'] = names[key]



    print(id, d[id])
    f.write(f'{genome_ids[i]}\t{id}\t')
    line = ''
    prev_level = 'superkingdom'
    for level in hier:
        if level in d[id]:
            line += d[id][level]+'\t'

        else:
            d[id][level] = d[id][prev_level]
            line += d[id][level]+"\t"

        # line += (d[id][level]+'\t' if level in d[id] else d[id][prev_level]+"\t")

        prev_level = level

    line = line[:-1]
    f.write(line)
    f.write('\n')
f.close()
