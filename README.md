# Adsorption_Clus

This repository provides the algorithms used in the paper "_Energy Decomposition to Access the Stability Changes Induced by CO Adsorption on Transition-metal 13-atom Clusters_" by Krys E. A. Batista, Marinalva D. Soares, Marcos G. Quiles, Maurício J.Piotrowski, and Juarez L. F. Da Silva. 

## Code for selecting representative (prototype) molecules:

This algorithm uses the K-means clustering algorithm for selecting representative samples from a dataset of compounds. Both geometric-based features (eigen values of the Coulomb matrix) and external variables, such as energy, are considered as encoding features.

**Usage:**

* Option 1: Only structural information (xyz)
	* $ python representative.py #samples folder_xyz_files
* Option 2: Structural information (xyz) + extra data
	* $ python representative.py #samples folder_xyz_files extra_data.txt
		
		Obs. when loading extra data, the selection of representatives can be biased with an external variable (e.g. energy) to deliver the results. In this case, the variable must be informed. See Option 3

* Option 3: Structural information (xyz) + extra data + biased selection
	* python representative.py #samples folder_xyz_files extra_data.txt id_variable min_max
		
		id_variable indicates which column of extra_data is considered
		
		min_max: 0 - select samples with the lowest value, 1 - the largest

**For more details, see the manuscript and support information**


## Code for adsorption of two-atoms molecules on clusters

The adsorption algorithm consistis of two main steps: 1) distribution of the particles in a spherical surface; 2) deformation of the positions to fit to the cluster surface structure. The first part is responsible for positioning _n_ particles on the surface of a sphere. We start with _n_ particles randomly set over the surface and, by using a force-field-based approach, we optimize the position of all particles to maximize their shortest distances.

**Usage**

* $ python adsorption.py xyz_cluster xyz_CO num_CO min_dist num_samples conection_type heterogeneity show3D
    * xyz_cluster: XYZ file of the cluster
    * xyz_CO: XYZ file of the CO molecule
    * num_CO: number of molecules of CO
    * min_dist: distance from CO to the cluster surface
    * num_samples: number of generated structures
    * heterogeneity: value between 0 and 1. [0 -> homogeneous distribution / 1 - random]
    * show3D: 0 (no) / 1(yes)
    conection_type: select which atom will interact with the cluster (two atoms only, e.g. CO)
    	* 0 - atom closest to the surface (i.e. C)
    	* 1 - second atom (i.e. O)

**For more details, see the manuscript and support information**


