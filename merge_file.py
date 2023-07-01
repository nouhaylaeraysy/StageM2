import csv

def merge_csv_file(file1, file2, merge_file):
    # read the first csv file 
    with open(file1, 'r') as f1:
        reader = csv.reader(f1)
        data1 = list(reader)

    # read the second file
    with open(file2, 'r') as f2:
        reader = csv.reader(f2)
        data2 = list(reader)

    # Creating a data structure for the merged file
    fusionne = []

    # Adding entries from the first file
    for row in data1:
        gene_id = row[0]
        matched_value = row[1]
        fusionne.append([gene_id, matched_value, ''])

    # Adding entries from the second file
    for row in data2:
        gene_id = row[0]
        matched_value = row[1]
        # Checking for duplicates
        if gene_id not in [entry[0] for entry in fusionne]:
            fusionne.append([gene_id, matched_value, ''])

    #Writing merged data to a new CSV file
    with open(merge_file, 'w', newline='') as ff:
        writer = csv.writer(ff)
        writer.writerows(fusionne)

    print("merge file is done.")

# Appel de la fonction de fusion
merge_csv_file("matched_genes.csv", "matched_withsymbol.csv", "merge_file.csv")
