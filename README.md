# Adsorption_Clus

This repository provides the algorithms used in the paper "_Energy Decomposition to Access the Stability Changes Induced by CO Adsorption on Transition-metal 13-atom Clusters_" by Krys E. A. Batista, Marinalva D. Soares, Marcos G. Quiles, Maurício J.Piotrowski, and Juarez L. F. Da Silva. 

## Code for selecting representative (prototype) molecules:

# Usage:
* Option 1: Only structural information (xyz)
	* $ python script1.py #samples folder_xyz_files
* Option 2: Structural information (xyz) + extra data
	* $ python script1.py #samples folder_xyz_files extra_data.txt
		Obs. when loading extra data, the selection of representatives can be biased with an external variable (e.g. energy) to deliver the results. In this case, the variable must be informed. See Option 3
* Option 3: Structural information (xyz) + extra data + biased selection
	* python script1.py #samples folder_xyz_files extra_data.txt id_variable min_max
		id_variable indicates which column of extra_data is considered
		min_max: 0 - select samples with the lowest value, 1 - the largest




### dasdads

- 1
- 2
- 3
