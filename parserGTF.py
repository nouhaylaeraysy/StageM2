import csv

def extract_gene_data_from_gtf(gtf_file):
    gene_data = {}  # Use a dictionary to store unique gene IDs

    # GTF file parser
    with open(gtf_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')

        for line in reader:
            if len(line) >= 9:
                attributes = line[8].split(';')  # Split attributes field by semicolon

                gene_id = None  # Define a null variable
                role = None

                for attr in attributes:
                    attr = attr.strip()  # Remove leading and trailing whitespaces
                    if attr.startswith("gene_id"):
                        gene_id = attr.split('"')[1]  # Extract gene ID from double quotes
                    elif attr.startswith("product"):
                        role = attr.split('"')[1]  # Extract role from double quotes

                if gene_id and role:
                    if gene_id not in gene_data:  # Add gene ID and role to dictionary if it doesn't exist
                        gene_data[gene_id] = role

    return gene_data

# Example usage:
gtf_file = "GCF_000002335.3_Tcas5.2_genomic.gtf"  # Replace with the path to your GTF file

gene_data = extract_gene_data_from_gtf(gtf_file)

# Print the extracted gene data
for gene_id, role in gene_data.items():
    print(f"Gene ID: {gene_id} | Role: {role}")

# Write gene_data to a CSV file
csv_filename = "List_of_genes_Description.csv"
with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Gene ID', 'Role'])  # Write header row
    for gene_id, role in gene_data.items():
        writer.writerow([gene_id, role])
