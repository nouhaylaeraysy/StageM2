import csv
import re

search_terms = ['zinc finger', 'DMRT', 'homeobox', 'NK','homeodomain', 'HMG', 'Nuclear', 'DNA-binding', 'POLYCOMB', 'trasncription','HLH' ,'Leucin zipper', 'box']
best_matches = {}
no = ['Socs' , 'nuclear pore' , 'ribonucleoprotein']
def search_exact_word(mot, texte):
    pattern = r'\b{}\b'.format(re.escape(mot))  # search for a whole word and ignore partial matches {}
    match = re.search(pattern, texte, re.IGNORECASE)
    if match:
        return texte
    else:
        return 0

with open('List_of_genes_Description.csv', 'r') as rf:
    reader = csv.reader(rf, delimiter=',')
    for row in reader:
        gene_id = row[0]  
        role = row[1]
        for term in search_terms:
            for x in no: 
                if x not in term:  # Exclude terms containing "Sox"
                    result = search_exact_word(term, row[1])
                    if result != 0:
                        best_matches[gene_id] = (row[1], term)

# Write best matches to a new CSV file
with open("matched_genes.csv", 'w', newline='') as wf:
    writer = csv.writer(wf)

    # Header of the CSV file
    writer.writerow(['gene_id', 'role', 'term'])

    # Write the data of best_matches to the CSV file
    for gene_id, (role, term) in best_matches.items():
        writer.writerow([gene_id, role, term])

print("The file matched_genes.csv has been created.")
