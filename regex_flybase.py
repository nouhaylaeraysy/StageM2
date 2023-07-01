import csv
import os.path
import re

# create a list to store symbols
symbols = set()

# open txt file
with open("FlyBase_Fields_download.txt", "r") as fichier:
    # read lines
    lignes = fichier.readlines()

    # extract column "SYMBOL"
    for ligne in lignes:
        # divide column
        colonnes = ligne.strip().split("\t")

        # if we have more than 4 columns
        if len(colonnes) >= 4:
            # extract column "SYMBOL"
            symbol = colonnes[3]

            # add symbol to the set
            symbols.add(symbol)

search_terms = list(symbols)

# Creating an empty dictionary
best_matches = {}

def search_exact_word(mot, texte):
    pattern = r'\b{}\b'.format(re.escape(mot))  # search for a whole word and ignore partial matches {}
    match = re.search(pattern, texte, re.IGNORECASE)
    if match:
        return texte
    else:
        return 0

with open('List_of_genes_Description.csv', 'r') as rf:
    reader = csv.reader(rf, delimiter=',')
    next(reader)  # Skip the first row (header)
    for row in reader:
        gene_id = row[0]
        role = row[1]
        for term in search_terms:
            result = search_exact_word(term, row[1])
            if result != 0:
                best_matches[gene_id] = (row[1], term)

# Check if the file "matched_withsymbol.csv" already exists
if not os.path.isfile('matched_withsymbol.csv'):
    # Create an output csv file
    with open('matched_withsymbol.csv', 'w', newline='') as wf:
        writer = csv.writer(wf)

        # header of csv file
        writer.writerow(['gene_id', 'role', 'term'])

        # write data of best_matches to csv file
        for gene_id, (role, term) in best_matches.items():
            writer.writerow([gene_id, role, term])

    print("The file matched_withsymbol.csv is created.")
else:
    print("The file matched_withsymbol.csv already exists.")
