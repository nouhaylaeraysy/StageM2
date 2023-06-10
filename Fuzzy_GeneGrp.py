#pip install fuzzywuzzy
#import modules 
import csv
from fuzzywuzzy import process

search_terms = ['zinc finger', 'DMRT', 'homeobox', 'NK','homeodomain', 'HMG', 'Nuclear', 'DNA-binding', 'POLYCOMB', 'trasncription','HLH']
best_matches = {}

with open('List_of_genes_Description.csv', 'r') as rf:
    reader = csv.reader(rf, delimiter=',')
    for row in reader:
        gene_id = row[0]  #extract gene_id 
        sol = process.extractOne(row[1], search_terms)  #return tuple of the best correspondance and the score of similarity 
        #print(sol)  # ('Nuclear', 51)
        if sol[1] >= 85: 
            if gene_id not in best_matches or sol[1] > best_matches[gene_id][1]:
                best_matches[gene_id] = (row[1], sol[1])
        #print(best_matches)  #'LOC100680201': ('nuclear transcription factor Y subunit B-1', 90),

# Write best matches to a new CSV file 
with open("matched_genes.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Gene ID', 'Matched Value', 'Score'])  #heading 
    for gene_id, (value, score) in best_matches.items():
        writer.writerow([gene_id, value, score])
