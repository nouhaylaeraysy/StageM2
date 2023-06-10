from thefuzz import fuzz
# Check the similarity score
name = "Kurtis Pykes"
full_name = "Kurtis K D Pykes"

print(f"Similarity score: {fuzz.ratio(name, full_name)}")

"""
Similarity Score: 86
"""
# Check the similarity score
name = "Kurtis Pykes"
full_name = "Kurtis K D Pykes"

print(f"Token sort ratio similarity score: {fuzz.token_set_ratio(name, full_name)}")

"""
Token sort ratio similarity score: 100
"""
# extraire du texte d'une collection à l'aide d'une correspondance de chaîne floue.
from thefuzz import process

collection = ["AFC Barcelona", "Barcelona AFC", "barcelona fc", "afc barcalona"]
print(process.extract("barcelona", collection, scorer=fuzz.ratio))

"""
[('barcelona fc', 86), ('AFC Barcelona', 82), ('Barcelona AFC', 82), ('afc barcalona', 73)]
"""

import pandas as pd

# Creating a dataframe
dict_one = {
  "country": ["England", "Scotland", "Wales", "United Kingdom", "Northern Ireland"],
  "population_in_millions": [55.98, 5.45, 3.14, 67.33, 1.89]
}

dict_two = {
  "country": ["Northern Iland", "Wles", "Scotlnd", "Englnd", "United K."],
  "GDP_per_capita": [24900, 23882, 37460, 45101, 46510.28]
}

existing_data = pd.DataFrame(dict_one)
exported_data = pd.DataFrame(dict_two)

print(existing_data, exported_data, sep="\n")

