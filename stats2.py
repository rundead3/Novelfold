# Column 1: Class assignment - B for buried or E for Exposed - Threshold: 25% exposure, but not based on RSA
# Column 2: Amino acid
# Column 3: Sequence name
# Column 4: Amino acid number
# Column 5: Relative Surface Accessibility - RSA
# Column 6: Absolute Surface Accessibility
# Column 7: Not used
# Column 8: Probability for Alpha-Helix
# Column 9: Probability for Beta-strand
# Column 10: Probability for Coil

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import os

import pandas as pd
import statsmodels.api as sm


# open csv file
file = open("C:\\Users\\Rundead\\net1.csv", "r")
# save the last column and start at the second line of the csv file removing `\n` at the end of each line
dis = []
phi = []
psi = []
rsa = []
asa = []

for i in file.readlines()[1:]:
    dis.append(float(i.split(',')[-1].replace('\n', '')))
    phi.append(float(i.split(',')[-3].replace('\n', '')))
    psi.append(float(i.split(',')[-2].replace('\n', '')))
    rsa.append(float(i.split(',')[3].replace('\n', '')))
    asa.append(float(i.split(',')[4].replace('\n', '')))





helyx, beta, coil = [], [], []
for line in file.readlines()[1:]:
    elements = line.split(',')
    helyx.append(float(elements[6]))
    beta.append(float(elements[7]))
    coil.append(float(elements[8]))

results = []
for i in range(len(helyx)):
    largest = max(helyx[i], beta[i], coil[i])
    if largest == helyx[i]:
        results.append('H')
    elif largest == beta[i]:
        results.append('E')
    else:
        results.append('-')
# close the file
file.close()


results_string = ''.join(results)


# I divided the results_string into 71 strings of 144 characters each by adding a new line every 144 characters
#print(len(results_string))

#for i in range(0,72):
    #print(results_string[144*i:144*i+144])

#open a new file
file = open("C:\\Users\\Rundead\\Novelfold\dssps\\netsurfdssp.txt", "r")

netsurfresults = file.read()
file.close()

#Compare the two strings, if equal write 1 to a new file, if not write 0

omega_string = open("C:\\Users\\Rundead\\Novelfold\\dssps\\dssps.txt", "r")

omega_string = omega_string.read()
y = []

for i in range(len(omega_string)):
    if omega_string[i] == netsurfresults[i]:
        y.append(1)
    else:
        y.append(0)



# Path to the folder containing the .pdb files
folder_path = 'C:\\pdbs\\'

# List to store the values in the specified columns
column_values = []

# Iterate through the files in the folder in sorted order
for file in sorted(os.listdir(folder_path)):
    # Check if the file is a .pdb file
    if file.endswith('.pdb'):
        # Open the file and read the lines
        with open(os.path.join(folder_path, file), 'r') as f:
            lines = f.readlines()
        # Extract the specified columns from each line and add them to the list
        for line in lines[:-2]:
            value = line[61:67].strip()
            if column_values == [] or value != column_values[-1]:
                column_values.append(value)
# Print the list of values in the specified columns

bi = y[0:10213]
dis = dis[0:10213]
phi = phi[0:10213]
psi = psi[0:10213]
rsa = rsa[0:10213]
asa = asa[0:10213]
pred = column_values[0:10213]

#########################################################################################

df = pd.DataFrame({'outcome': bi, "prediction": pred, "dis": dis, "phi": phi, "psi": psi, "rsa": rsa, "asa": asa})


x = pred,dis,phi,psi,rsa,asa  # predictor variables
y = bi # response variables (binary variables, 0 or 1)

formula = 'outcome ~ prediction+dis+phi+psi+rsa+asa'

df.to_csv('data.csv', index=False)
print("bananana")
