Novelfold
Overview
Novelfold is a Python-based genetic algorithm that utilizes MAP-Elites for simulating the generation, transformation, and evaluation of proteins. The protein sequences (Fastas) serve as genotypes and are transformed into 3D structures (PDBs) using OmegaFold. This transformation, often referred to as the genotype-to-phenotype mapping, is a critical aspect of the simulation.

The proteins go through mutations and crossovers, fitness and survivability of the proteins are calculated based on arbitrary criteria. The MAP-Elites algorithm is employed for population management, ensuring a diverse range of solutions across different dimensions of performance.

Files
The codebase consists of several Python scripts, each with a specific role:

main.py - This is the main program file that runs the simulation. It handles the generation of populations of protein sequences, performs transformations on the sequences, and calculates the fitness and survivability of proteins based on the MAP-Elites algorithm.

chainLogic.py - Defines a class, chainLogic, that represents a protein sequence. It includes attributes and methods for calculating various features of the protein's 3D structure and loading necessary data.

grid.py - Contains a grid class which is a simple grid data structure with functionalities to get, set, and replace values in the grid, as well as to find the neighbours of a cell.

nicheSpace.py - Defines a nicheSpace class representing the performance space (also known as a niche space) that contains the proteins. It includes methods for adding entries to the space, getting random entries, adjusting the range of the space, and writing the archive of the space into fasta format.

triplets.py - Contains functions related to amino acid triplets, such as generating a random triplet or a random protein of triplets, and converting a triplet into a one-letter code (fasta format).

config.py - Contains functions that return different configuration parameters such as paths to directories or files.

How to Run
To run Novelfold, execute the main.py file. Make sure that the required data files and tools are available at the paths specified in config.py.

Dependencies
Novelfold makes use of several libraries and tools including:

NumPy
BioPython
OmegaFold (external tool for genotype-to-phenotype transformation)
DSSP (external tool for secondary structure determination)
FreeSASA (external tool for Solvent Accessible Surface Area (SASA) calculation)
Pythia

