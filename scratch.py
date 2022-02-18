from ete3 import NCBITaxa

ncbi = NCBITaxa()
# ncbi.update_taxonomy_database() # wystarczy raz pobrać - polecenie aktualizuje lokalną bazę danych co zajmuje trochę czasu
file_name = 'input_files/taxIDs.txt'
file = open(file_name)
# --------------------------------------- TEST 1 -------------------------------------
hier = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]  # hierarchia jaką chcemy otrzymać w pliku - zmiana tej listy nie wystarczy, wymagana jest zmiana kodu poniżej
# file = ['2162078'] # na potrzeby testowania dla 1 ID


d = {}  # słownik zawierajacy taksonomię dla każdego ID
for id in file:
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
        elif ranks[key] == 'class' or ranks[key] == 'clade':
            d[id]['class'] = names[key]
        elif ranks[key] == ('order'):
            d[id]['order'] = names[key]
        elif ranks[key] == ('family'):
            d[id]['family'] = names[key]
        elif ranks[key] == ('genus'):
            d[id]['genus'] = names[key]
        elif ranks[key] == ('species'):
            d[id]['species'] = names[key]

# zapis do pliku
f = open('taxonomy_from_id.csv', 'w')
f.write(f'id\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n')
for id in d:
    print(id, d[id])
    f.write(f'{id}\t')
    for level in hier:
        f.write(d[id][level]+'\t' if level in d[id] else "\t")
    f.write('\n')
f.close()
