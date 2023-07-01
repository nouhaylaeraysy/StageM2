import csv

gtf_file = "GCF_009193385.2_Nvit_psr_1.1_genomic.gtf"
gene_data = {}  # Use a dictionary to store unique gene IDs

# Parser of GTF file
with open(gtf_file, 'r') as file:
    reader = csv.reader(file, delimiter='\t')

    for line in reader:
        if len(line) >= 9:
            attributes = line[8].split(';')  # Split attributes field by semicolon

            gene_id = None # define a null variable
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

# Print the extracted gene data
for gene_id, role in gene_data.items():
    print(f"Gene ID: {gene_id} | Role: {role}")

# Write gene_data to a CSV file
with open("List_of_genes_Description.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Gene ID', 'Role'])  # Write header row
    for gene_id, role in gene_data.items():
        writer.writerow([gene_id, role])