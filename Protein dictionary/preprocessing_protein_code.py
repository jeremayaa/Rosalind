import pickle

with open('protein_code.txt') as f:
    protein_code = f.read().splitlines()

protein_dict = {}

for line in protein_code:
    pairs = line.split()
    for i in range(0, len(pairs), 2):
        code = pairs[i]
        amino_acid = pairs[i+1]
        protein_dict[code] = amino_acid

with open('protein_dict.pickle', 'wb') as f:
    pickle.dump(protein_dict, f)


# It would be way faster to write that dictionary manually 
