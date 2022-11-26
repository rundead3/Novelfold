import pandas as pd
import numpy as np
import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from random import random

import matplotlib as mpl
from matplotlib.lines import Line2D

import os 


## COLOR MAPS
cmap = LinearSegmentedColormap.from_list(
    "mycmap", [(0.0 / 1, "red"), ((0.5) / 1, "whitesmoke"), (1.0, "blue")]
)

vmax=2.5
cmap_enrich = LinearSegmentedColormap.from_list('mycmap', [(0/ vmax, 'red'), (1./vmax, 'whitesmoke'), (vmax / vmax, 'blue')])

cNorm_enrich = matplotlib.colors.Normalize(vmin=0, vmax=2) #re-wrapping normalization
scalarMap_enrich = matplotlib.cm.ScalarMappable(norm=cNorm_enrich, cmap=cmap)

cmap_disorder = plt.get_cmap('PuOr')
cmap_u = plt.get_cmap('PuOr')
#This is when you want to change the scale of colormap
cNorm = matplotlib.colors.Normalize(vmin=-0.3, vmax=0.3) #re-wrapping normalization
scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cmap_u)
cval = scalarMap.to_rgba(0)

blob_length_cutoff_enrichment = pd.read_csv('./Table_S1.csv')
dict_enrich = dict(zip(zip(blob_length_cutoff_enrichment['Hydrophobicity cutoff'], blob_length_cutoff_enrichment['Blob length']), blob_length_cutoff_enrichment['Enrichment ']))


def uversky_color(x):
    ncpr = x[0]
    m_color = scalarMap.to_rgba(ncpr)
    return "rgb" + str(tuple([255 * x for x in m_color[:-1]]))

def NCPR_color(x):
    ncpr = x[0]
    #Setting the limts of the ncpr graph to 0.5 and -0.5, any number higher or lower will be considered max
    if ncpr > 0.5:
        ncpr = 1
    if ncpr < -0.5:
        ncpr = -1
    ncpr_normalized = (ncpr + 1.0) / 2
    m_color = cmap(ncpr_normalized)

    return "rgb" + str(tuple([255 * x for x in m_color[:-1]]))

def disorder_color(x):
    ncpr = x[0]
    m_color = cmap_disorder(ncpr)
    return "rgb" + str(tuple([255 * x for x in m_color[:-1]]))

def enrichment_color(enrich_value):
    m_color = scalarMap_enrich.to_rgba(enrich_value)
    return "rgb" + str(tuple([255 * x for x in m_color[:-1]]))




def writeCMap(func, fname, vmin, vmax, res):
    colorResolution = 2*(10**res)+1  #this defines how many colors will be defined.

    values = [] 
    keys = np.linspace(vmin,vmax,colorResolution).round(res)
    for i in keys:
        values.append(func([i]))
    myDict = dict(zip(keys, values))
    myDF = pd.DataFrame.from_dict(myDict, orient='index')
    #file = open(fname, "wb")
    #pickle.dump(myDict, file)
    #file.close()
    myDF.to_csv(fname)

res = 2 #this determines the resolution in terms of the number of decimal places

writeCMap(uversky_color, "uverskyCMap.csv", -1, 1, res)

writeCMap(NCPR_color, "ncprCMap.csv", -1, 1, res)

writeCMap(disorder_color,"disorderCMap.csv",-1, 1, res)




# Generating the custom color maps for the plots
cmap = LinearSegmentedColormap.from_list(
    "mycmap", [(0.0 / 1, "red"), ((0.5) / 1, "whitesmoke"), (1.0, "blue")]
)

vmax=2.5
cmap_enrich = LinearSegmentedColormap.from_list('mycmap', [(0/ vmax, 'red'), (1./vmax, 'whitesmoke'), (vmax / vmax, 'blue')])

cNorm_enrich = matplotlib.colors.Normalize(vmin=0, vmax=2) #re-wrapping normalization
scalarMap_enrich = matplotlib.cm.ScalarMappable(norm=cNorm_enrich, cmap=cmap)

def enrichment_color(enrich_value):
    m_color = scalarMap_enrich.to_rgba(enrich_value)
    return "rgb" + str(tuple([255 * x for x in m_color[:-1]]))


df = pd.read_csv('./Table_S1.csv')

dfMI = df.set_index(["Hydrophobicity cutoff", "Blob length"])

#Special function for enrichment prediction color which requires the cutoff
dfMI["color"] = ""
errCount = 0
for i in range(df.shape[0]):
    try:
        enrich_value = df.loc[i, "Enrichment "]
        m_color = scalarMap_enrich.to_rgba(enrich_value)
        theColor = "rgb" + str(tuple([255 * x for x in m_color[:-1]]))
        dfMI["color"].iloc[i] = theColor
    except KeyError:
        errCount = errCount+1
print(f'Num errors: {errCount}')

if len(dfMI.index) != len(dfMI.index.unique()):
    print("WARNING: non-unique blob-size/cutoff enrichment predictions!!!")

dfMI.to_csv('./enrichCMap.csv')
#fname = "enrichCMap.pkl"
#file = open(fname, "wb")
#pickle.dump(dfMI, file)
#file.close()