Novelfold
Overview
Novelfold is a Python program used to simulate the generation and transformation of amino acid chains. This simulation allows for the assessment of different features of the chains such as the bounding sphere, root mean square deviation, gyration radius, and start-end distance, among others. The chains are also folded and mutated as part of the simulation.

Files
The codebase consists of several Python scripts, each with a specific role:

main.py - This is the main program file that runs the simulation. It handles the generation of populations of amino acid chains, performs transformations on the chains, and calculates the fitness and survivability of chains.

chainLogic.py - Defines a class, chainLogic, that represents an amino acid chain. It includes attributes and methods for calculating various features of the chain and loading necessary data.

grid.py - Contains a grid class which is a simple grid data structure with functionalities to get, set, and replace values in the grid, as well as to find the neighbours of a cell.

nicheSpace.py - Defines a nicheSpace class representing a niche space that contains the chains. It includes methods for adding entries to the space, getting random entries, adjusting the range of the space, and writing the archive of the space into fasta format.

triplets.py - Contains functions related to amino acid triplets, such as generating a random triplet or a random chain of triplets, and converting a triplet into a one-letter code (fasta format).

config.py - Contains functions that return different configuration parameters such as paths to directories or files.

How to Run
To run Novelfold, execute the main.py file. Make sure that the required data files and tools are available at the paths specified in config.py.

Dependencies
Novelfold makes use of several libraries and tools including:

NumPy
BioPython
OmegaFold (external tool)
DSSP (external tool)
FreeSASA (external tool)
Please ensure these are installed and correctly configured before running the program.

