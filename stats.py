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

dis = dis[0:10223]
phi = phi[0:10223]
psi = psi[0:10223]
rsa = rsa[0:10223]
asa = asa[0:10223]


#########################################################################################

df = pd.DataFrame({'outcome': y, 'disorder': dis, 'phi_angle': phi, 'psi_angle': psi, 'relative_solvent_accessibility': rsa, 'accessible_surface_area': asa})

formula = 'outcome ~ disorder + phi_angle + psi_angle + relative_solvent_accessibility + accessible_surface_area'

x = dis,phi,psi,rsa,asa  # predictor variables
y = y # response variables (binary variables, 0 or 1)

# fit the model
model = sm.GLM.from_formula(formula, data=df, family=sm.families.Binomial())

result = model.fit()

print(result.summary())


pca = PCA(n_components=5)
pca.fit(x)



explained_variance = pca.explained_variance_ratio_

plt.bar(range(1, 6), explained_variance)
plt.xlabel('Principal component')
plt.ylabel('Explained variance')
plt.show()

print(explained_variance)

loadings = pca.components_

# Plot the loadings
plt.plot(dis, loadings[0, :], '-o', label='Principal component 1')
plt.xlabel('Original variable')
plt.ylabel('Loading')
plt.legend()
plt.show()
plt.plot(phi, loadings[1, :], '-o', label='Principal component 2')
plt.xlabel('Original variable')
plt.ylabel('Loading')
plt.legend()
plt.show()
plt.plot(psi, loadings[2, :], '-o', label='Principal component 3')
plt.xlabel('Original variable')
plt.ylabel('Loading')
plt.legend()
plt.show()
plt.plot(rsa, loadings[3, :], '-o', label='Principal component 4')
plt.xlabel('Original variable')
plt.ylabel('Loading')
plt.legend()
plt.show()
plt.plot(asa, loadings[4, :], '-o', label='Principal component 5')
plt.xlabel('Original variable')
plt.ylabel('Loading')
plt.legend()
plt.show()

# Plot the loadings
plt.plot(loadings[1, :], loadings[3, :], '-o', label='Principal component 1')
plt.xlabel('Loading2')
plt.ylabel('Loading')
plt.legend()
plt.show()

