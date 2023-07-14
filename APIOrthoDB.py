import requests
import json
import time
import csv

# Read gene IDs from the text file
with open("candidat_TF_Data1.txt", "r") as file:
    gene_ids = file.readlines()
    gene_ids = [gene_id.strip() for gene_id in gene_ids]

gene_data = {}
# Loop to execute the search for each gene ID
for gene_id in gene_ids:
    url = f"https://data.orthodb.org/v11/genesearch?query={gene_id}"
    key = gene_id
    # Send GET request and retrieve data
    response = requests.get(url)
    data = response.content

    # Store data in a file
    with open(f"result_{gene_id}.txt", "wb") as file:
        file.write(data)

    print(f"Data for gene ID {gene_id} has been saved in the file 'result_{gene_id}.txt'.")

    # Read the file
    with open(f"result_{gene_id}.txt", "r") as file:
        data = json.load(file)

    # Extract "CG8159" from "CG8159;Dmel\\CG8159"
    if "orthologs_in_model_organisms" in data:
        orthologs = data["orthologs_in_model_organisms"]
        if orthologs and len(orthologs) > 0:
            gene_id = orthologs[0]["genes"][0]["gene_id"]["id"]
            extracted_value = gene_id.split(";")[0].strip()
            value = extracted_value
            print("Extracted:", extracted_value)
    gene_data[key] = value
    print(gene_data)
    # Wait for 5 seconds before the next request
    time.sleep(1)
# Write gene_data dictionary to a text file
# Ouvrir le fichier CSV en mode écriture
with open("candidat_TF_Data1.csv", "w", newline="") as file:
    writer = csv.writer(file, delimiter="\t")

    # Écrire l'en-tête du fichier CSV
    writer.writerow(["gene_id virilis", "gene_id Dmel"])

    # Écrire chaque paire clé-valeur dans une ligne du fichier CSV
    for key, value in gene_data.items():
        writer.writerow([key, value])
