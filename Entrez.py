from Bio import Entrez

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
        gene_names.append(gene_name)

    # Print the gene names
    for gene_name in gene_names:
        print(gene_name)

    # Store the gene names in the output file
    with open(output_file, 'w') as file:
        for gene_name in gene_names:
            file.write(gene_name + '\n')

    print("The gene names have been stored in the file", output_file)

# Usage example
input_file = "list_ofTF_Data2_fetchData.txt"
output_file = "gene_names_Data2_fetch.txt"
fetch_gene_names(input_file, output_file)

