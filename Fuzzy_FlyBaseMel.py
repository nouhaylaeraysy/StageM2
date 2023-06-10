import csv
from fuzzywuzzy import process
import os.path

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
best_matches = {}

with open('List_of_genes_Description.csv', 'r') as rf:
    reader = csv.reader(rf, delimiter=',')
    next(reader)  # Skip the first row (header)
    for row in reader:
        gene_id = row[0]  # extract gene_id
        sol = process.extractOne(row[1], search_terms)  # return tuple of the best correspondence and the similarity score
        if sol[1] >= 90:
            if gene_id not in best_matches or sol[1] > best_matches[gene_id][1]:
                best_matches[gene_id] = (row[1], sol[1])

output_file = "matched_withsymbolDmel.csv"

# Write best matches to a new CSV file
write_header = not os.path.exists(output_file) or os.path.getsize(output_file) == 0

with open(output_file, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile)
    if write_header:
        writer.writerow(['Gene ID', 'Matched Value', 'Score'])  # Write header only if the file is empty
    for gene_id, (value, score) in best_matches.items():
        writer.writerow([gene_id, value, score])

"""
sol = process.extractOne(row[1], search_terms)
matched_word = sol[0]  # Mot trouvé
score = sol[1]  # Score de similarité

if len(matched_word) == len(search_terms):
    # La taille du mot trouvé correspond à la taille de search_terms
    # Faites quelque chose avec le mot trouvé
    print("Mot trouvé avec la même taille :", matched_word)
else:
    # La taille du mot trouvé ne correspond pas à la taille de search_terms
    print("Mot trouvé avec une taille différente :", matched_word)
"""







