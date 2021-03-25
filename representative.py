#!/usr/bin/python
# Selection of Representative Molecules
# Marcos Quiles (quiles@gmail.com) / Marinalva Soares

"""
INPUT 1 - target number of representatives
INPUT 2 - XYZ files (folder)
INPUT 3 [Optional] - Additional information (Energy_T, Delta_E, etc)
    obs. First column must contain the ids of the files of INPUT 2 (check example)
INPUT 4 [Optional] - Informe if the selection must be biased by any external variable
    obs. By default, samples are selected by their proximity with the cluster centroid
         Range 0..#, # -> number of external variables
INPUT 5 [Optional] - Indicates if the min or max of the external variable should be used
    obs. By default, the code will select the samples with the lowest value in each cluster
"""

import sys
import pandas as pd
from glob import glob
from os import makedirs, getcwd, path
from tools import *
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

def load_args():
    """Checking args."""
    if len(sys.argv) < 3 or len(sys.argv) > 6:
        print("\nSelection of representative Molecules\n")
        print("\tUsage: \n")
        print("\tOption 1: Only structural information (xyz) \n")
        print("\t\t$ python script1.py #samples folder_xyz_files \n\n")
        print("\tOption 2: Structural information (xyz) + extra data\n")
        print("\t\t$ python script1.py #samples folder_xyz_files extra_data.txt \n\n")
        print("\tObs. when loading extra data, the selection of representatives")
        print("\t     can be biased with an external variable (e.g. energy) to ")
        print("\t     deliver the results. In this case, the variable must be")
        print("\t     informed. See Option 3\n")
        print("\tOption 3: Structural information (xyz) + extra data + biased selection\n")
        print("\t\t$ python script1.py #samples folder_xyz_files extra_data.txt id_variable min_max\n")
        print("\t\t\t id_variable indicates which column of extra_data is considered\n")
        print("\t\t\t min_max: 0 - select samples with the lowest value, 1 - the largest\n")
        exit();

    print("Checking parameters.")
    arglen = len(sys.argv)

    nrep = int(sys.argv[1])
    if nrep < 2:
        print('\nERROR: invalid number of representatives. Must be > 1')
        exit()

    dfolder = glob(str(sys.argv[2])+'/*.xyz')
    if len(dfolder) == 0:
        print('\nERROR: parameter failure')
        exit()

    if len(dfolder) < nrep:
        print('\nERROR: number of representatives must be > number of samples')
        exit()

    dfolder.sort()
    if arglen>3: extra_file = sys.argv[3]
    else: extra_file = ""
    if arglen>4: bias_id = int(sys.argv[4])
    else: bias_id = -1 
    if arglen>5: bias_op = int(sys.argv[5])
    else: bias_op = -1

    return nrep, dfolder, extra_file, bias_id, bias_op


def load_data(dfolder, extra):
    """Load, normalize, and return the dataset."""
    print("Loading data.")
    dataxyz = []
    ids = []
    num_files = 1
    dim_extra = 0

    for fin in dfolder:
        natoms, atomtypes, coords = xyzRead(fin)
        mat = eigenCoulomb(fin,natoms)
        dataxyz.append(mat)
        ids.append(num_files)
        num_files+=1

    if extra.strip():
        print("Loading extra data.")
        try:
            datae = pd.read_csv(extra, sep="\t")
            columns = datae.columns.tolist()
            datae.sort_values(columns[0], axis=0, ascending=True, inplace=True)
            datae = datae.iloc[:,1:].values
            dim_extra = datae.shape[1]
            # print("Dimension extra file", dim_extra)
            data = np.concatenate((datae,np.array(dataxyz)), axis=1)
        except IOError:
            print("ERROR: invalid extra data file.")
            exit()
    else:
        data = np.array(dataxyz)

    # print(data.shape)
    data = StandardScaler().fit_transform(data)
    return data, ids, dim_extra

def run_clustering(data, nrep):
    """Clustering the dataset."""
    print("Clustering data.")
    # kmeans = KMeans(init='k-means++', 
    kmeans = KMeans(init='random', 
        n_clusters=nrep, n_init=10, random_state=0)
    kmeans.fit(data)
    return kmeans

def get_representatives(kmeans, data, ids, bias_id, bias_op):
    """Selection of the representative molecules."""
    print("Selecting representative samples.")
    selected = []
    res_clustering = []
    centroids = kmeans.cluster_centers_

    if bias_id == -1: #standard procedure
        # Molecules closest to the cluster centroids are selected as representatives
        for clus in range(len(centroids)):
            idsIn = np.where(kmeans.labels_==clus)
            sel = idsIn[0][0]
            dMin = np.linalg.norm(centroids[clus,:] - data[idsIn[0][0],:])
            for sample in idsIn[0][1:]:
                dist = np.linalg.norm(centroids[clus,:] - data[sample,:])
                if dist < dMin:
                    dMin = dist
                    sel = sample
            selected.append(int(ids[sel]))
            res_clustering.append(idsIn)
    else:
        # A representative molecule is selected from each cluster. 
        # The selection considers the value of the informed external variable
        print("Selecting biased samples ", bias_id, " ",bias_op)
        for clus in range(len(centroids)):
            idsIn = np.where(kmeans.labels_==clus)
            sel = idsIn[0][0]
            val = data[sel,bias_id]
            for sample in idsIn[0][1:]:
                sample_val = data[sample,bias_id]
                if bias_op == 0:
                    if sample_val < val:
                        val = sample_val
                        sel = sample
                else:
                    if sample_val > val:
                        val = sample_val
                        sel = sample
            selected.append(int(ids[sel]))
            res_clustering.append(idsIn)

    return selected, res_clustering


def save_results(selected, clusters, ids, data):
    """Save results into a txt file."""
    print("Saving results [output.txt].")
    with open("output.txt", 'w') as f:
        cc = 0
        for clus in clusters:
            f.write("Cluster  %d - " % cc)
            f.write("[%d] " % selected[cc])
            for item in clus[0]:
                if ids[item] != selected[cc]:
                    f.write("%d " % ids[item])
            f.write("\n")
            cc += 1


def main():
    nrep, dfolder, extra, bias_id, bias_op = load_args()
    data, ids, dim_extra = load_data(dfolder, extra)
    if bias_id >= dim_extra:
        bias_id = -1
        print("\tWARNING: extra file does not contain the choosen variable")
        print("\t         the standard selection will be taken - centroid proximity")
    kmeans = run_clustering(data, nrep)
    selected, res_clustering = get_representatives(kmeans, data, ids, bias_id, bias_op)
    save_results(selected, res_clustering, ids, data)

if __name__ == '__main__':
    main()

