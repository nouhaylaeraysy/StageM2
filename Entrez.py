from Bio import Entrez
import csv 

gene_data = {}
def fetch_gene_names(input_file, output_file):
    Entrez.email = 'ernouha2@gmail.com'  # email address registered in NCBI

    # Read gene IDs from the input file
    gene_ids = []
    with open(input_file, "r") as file:
        for row in file:
            gene_id = row.strip()
            gene_ids.append(gene_id)

    # Fetch the gene information from Batch Entrez
    handle = Entrez.efetch(db='gene', id=gene_ids, rettype='docsum')

    # Read the response as XML
    response = Entrez.read(handle)

    # Extract the gene names from the response
    gene_names = []
    for record in response['DocumentSummarySet']['DocumentSummary']:
        gene_name = record['Description']
        id = record['Name']
        #gene_names.append(gene_name)
        gene_data[id] = gene_name
    #print(gene_data)

    with open("gene_id_name.csv", "w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")

        # Écrire l'en-tête du fichier CSV
        writer.writerow(["gene_id virilis", "gene_ description virilis "])

        # Écrire chaque paire clé-valeur dans une ligne du fichier CSV
        for key, value in gene_data.items():
            writer.writerow([key, value])


# Usage example
input_file = "list_TF_PheatmapKons1.txt"
output_file = "gene_id_names.txt"
fetch_gene_names(input_file, output_file)
