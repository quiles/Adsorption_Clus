#!/usr/bin/python
# Script for adsorption of CO molecules in clusters
# Marcos Quiles (quiles@gmail.com) / Marinalva Soares

# -*- Coding: UTF-8 -*-
#coding: utf-8

from tools import *
import numpy as np
import sys

if len(sys.argv) != 9:
    print("Adsorption of CO molecules on metal clusters");
    print("\nUsage: \n")
    print("\t$ ./adsorption xyz_cluster xyz_CO num_CO min_dist num_samples conection_type heterogeneity show3D\n");
    print("\t\txyz_cluster: XYZ file of the cluster");
    print("\t\txyz_CO: XYZ file of the CO molecule");
    print("\t\tnum_CO: number of molecules of CO");
    print("\t\tmin_dist: distance from CO to the cluster surface");
    print("\t\tnum_samples: number of generated structures");
    print("\t\theterogeneity: value between 0 and 1. [0 -> homogeneous distribution / 1 - random]");
    print("\t\tshow3D: 0 (no) / 1(yes)");
    print("\t\tconection_type: select which atom will interact with the cluster (two atoms only, e.g. CO)");
    print("\t\t\t 0 - atom closest to the surface (i.e. C)");
    print("\t\t\t 1 - second atom (i.e. O) \n\n");

    print("\t i.e. $ python adsorption.py Ru13_001.xyz CO.xyz 10 2 5 1 0.1 1");
    print("\t\truns the script with following inputs")
    print("\t\t\t->Ru13_001.xyz (file containing the cluster configuration)")
    print("\t\t\t->CO.xyz (file containing the ??? molecule) ")
    print("\t\t\t->10 (number of CO molecules) ")
    print("\t\t\t->2 (distance from CO to the Cluster)")
    print("\t\t\t->5 (number of random configurations)")
    print("\t\t\t->1 (informs that atom O will be positioned close to the cluster)")
    print("\t\t\t->0.1 (target heterogeneity)")
    print("\t\t\t->1 (plot the 3D result)\n\n")    
    print("\t\tOutput: set of files: coords_#.txt\n\n\n")

    exit()


inputMol1 = sys.argv[1]
inputMol2 = sys.argv[2]
num_p = int(sys.argv[3])
min_dist = float(sys.argv[4])
num_samples = int(sys.argv[5])
conection_type = int(sys.argv[6])
heterog = float(sys.argv[7])

if (int(sys.argv[8]) == 0):
	showimg = False
else:
	showimg = True

print("Loading data.")
natoms, atomtypes, coords = xyzRead(inputMol1)
natoms_par, atomtypes_par, coords_par = xyzRead(inputMol2)
dist_atom_par = np.linalg.norm(coords_par[0,:]-coords_par[1,:])

print("Evolving particles (optimizing).")
xyz = evolveParticles6(num_p, maxDist(coords)+min_dist, coords, heterog)

print("Generating samples.")

if showimg:
	num_p2 = num_p*2
	import matplotlib
	matplotlib.use("TkAgg")
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	# import matplotlib.tri as mtri

if conection_type==0:
	idxC = 0
	idxO = 1
else:
	idxC = 1
	idxO = 0

for s in range(num_samples):
	if heterog==1:
		xyz = evolveParticles6(num_p, maxDist(coords), coords, heterog)
	else:
		xyz = rotMatrix(xyz)

	# Deformation of the sphere (particles) to approximate the cluster surface
	sample = getSample4(xyz, min_dist, coords) 
	sample, atomtypes_parS = adjustMol(sample,atomtypes_par,conection_type,dist_atom_par)
	fileName = 'coords_'+str(s+1)+'.xyz'
	print("Sphere deformation %d - [%s]" %(s+1, fileName))

	with open(fileName, 'w') as f:
		f.write("  %d X_CG, Y_CG, Z_CG 0.000000 0.000000 0.000000\n\n" %(natoms+(num_p*natoms_par)))
		for coor, ato in zip(coords, atomtypes):
			f.write("%s \t%.18g \t%.18g \t%.18g\n" %(ato, coor[0], coor[1], coor[2]))
		for coor, ato in zip(sample, atomtypes_parS):
			f.write("%s \t%.18g \t%.18g \t%.18g\n" %(ato, coor[0], coor[1], coor[2]))

	if showimg:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(sample[idxC:num_p2:2,0], sample[idxC:num_p2:2,1], sample[idxC:num_p2:2,2], s=100)
		ax.scatter(sample[idxO:num_p2:2,0], sample[idxO:num_p2:2,1], sample[idxO:num_p2:2,2], s=30)
		ax.scatter(coords[:,0], coords[:,1], coords[:,2], s=200)
		ax.set_xlabel('X Label')
		ax.set_ylabel('Y Label')
		ax.set_zlabel('Z Label')
		plt.title('Illustration: #'+str(s+1))

		plt.show()

