# -*- Coding: UTF-8 -*-
#coding: utf-8

import numpy as np
import sys

# __all__ = ["xyzRead", "xyzWrite"]

def xyzRead(fname):
    try: 
        fin = open(fname, "r")
        line1 = fin.readline().split()
        natoms = int(line1[0])
        comments = fin.readline()[:-1]
        coords = np.zeros([natoms, 3], dtype="float64")
        atomtypes = []
        for x in coords:
            line = fin.readline().split()
            atomtypes.append(line[0])
            x[:] = list(map(float, line[1:4]))
    except:
        print("ERROR: failed to read the xyz file.")
        exit()            

    return natoms, atomtypes, coords;


def getCharge(element):
    f = open("mol.txt")
    atomicnum = [line.split()[1] for line in f if line.split()[0] == element]
    f.close()
    return int(atomicnum[0])


def coulombMatrix(fname):
    natoms, atomtypes, coords = xyzRead(fname)
    i=0 ; j=0    
    colM = np.zeros((natoms,natoms))
    chargearray = np.zeros((natoms,1))
    charge = [getCharge(symbol)  for symbol in atomtypes]
    for i in range(0,natoms):
        colM[i,i]=0.5*charge[i]**2.4   # Diagonal term described by Potential energy of isolated atom
        for j in range(i+1,natoms):
            dist= np.linalg.norm(coords[i,:] - coords[j,:])   
            colM[j,i] = charge[i]*charge[j]/dist   #Pair-wise repulsion 
            colM[i,j] = colM[j,i]
    return colM


def eigenCoulomb(fname, num):
    sCoulomb = coulombMatrix(fname)
    eigValues = -np.sort(-np.linalg.eigvals(sCoulomb)).real
    return eigValues[0:num]


def maxDist(coords):
    norm = np.linalg.norm(coords,axis=1)
    return np.amax(norm)


def evolveParticles6(num_p, maxdist, coords, heterog):
    particles = np.random.rand(num_p,3)-0.5
    norm = np.linalg.norm(particles,axis=1)
    norm = norm.reshape((num_p,1))
    norm = np.repeat(norm,3,axis=1)
    norm /= maxdist
    particles = np.divide(particles,norm)
    t = 0
    
    if heterog == 1:
        return particles

    maxstd_max = 0
    while True:
        step = np.zeros((num_p,3))
        mdistv = []
        mindistV = []
        stddistV = []
        stepmean = []
        for i in range(num_p):
            diff = particles[i,:] - particles
            norm = np.linalg.norm(diff,axis=1)/maxdist
            stddistV.append(np.std(norm)) 
            norm[i] = sys.maxsize
            mindistV.append(np.amin(norm))
            sim = np.exp(-norm).reshape((num_p,1))
            step -= sim*diff / norm.reshape((num_p,1))
            stepmean.append(np.linalg.norm(step))

        if num_p!=2:
            varstd = np.std(stddistV)
            if varstd > maxstd_max:
                maxstd_max = varstd
            varstd /= maxstd_max
        else:
            ddiff = 2 - np.max(mindistV)
            varstd = ddiff

        print(str(t)+" - [ "+str(varstd)+" ]")

        if (varstd < heterog or t>1000): 
            print("\n\nThe shortest distances between 1-neighbor particles (before deformation): ")
            print(mindistV)
            break
        t += 1
        particles += step
        for i in range(num_p):
            norm = np.linalg.norm(particles[i,:])
            norm /= maxdist
            particles[i,:] /= norm

    return particles


def getSample4(xyz, maxdist, coords):
    num_p = np.shape(xyz)[0]
    num_a = np.shape(coords)[0]
    
    # expanding the cluster with radius
    norm = np.linalg.norm(coords,axis=1).reshape((num_a,1))
    inc = np.divide(coords,norm.repeat(3,axis=1))
    new_coords = coords# + np.multiply(inc,radius)

    particles = xyz
    step = np.zeros((num_p,3))
    learning = -np.ones(num_p)*0.1
    signal = np.ones(num_p)
    t = 1
    while True:
        diferenca = np.zeros(num_p)
        for i in range(num_p):
            diff = particles[i,:] - new_coords
            norm = np.linalg.norm(diff,axis=1)

            mindist1 = np.amin(norm)
            j1 = np.where(norm == mindist1)
            norm_j1 = np.linalg.norm(new_coords[j1,:])
            norm_i = np.linalg.norm(particles[i,:])
            diff_ij = norm_j1 - norm_i

            diferenca[i] = norm[j1]
            
            if (maxdist - diferenca[i])*signal[i] > 0:
                signal[i] = -signal[i]
                learning[i] *= -0.7

            normalization = norm_i+learning[i]*diferenca[i]
            if normalization < norm_j1:
                normalization = norm_j1
            particles[i,:] /= norm_i
            particles[i,:] *= normalization

        t += 1
        if (t>200): break

    return particles


def distMin(xyz, num_p):
    dist = np.ones((num_p))*10
    for i in range(num_p):
        for j in range(num_p):
            if j==i: continue
            norm = np.linalg.norm(xyz[i,:] - xyz[j,:])
            if norm < dist[i]:
                dist[i] = norm
    return dist

import math
from random import random


def rotMatrix(coords):
    theta = random()*math.pi*2
    co = math.cos(theta)
    si = math.sin(theta)
    rx = np.array([[1,0,0],[0,co,-si],[0, si, co]])
    ry = np.array([[co,0,si],[0,1,0],[-si, 0, co]])
    rz = np.array([[co,-si,0],[si,co,0],[0, 0, 1]])
    R = np.matmul(rz, np.matmul(ry,rx))
    newcoords = np.matmul(coords,R)
    return newcoords


def adjustMol(xyz,atomtypes_par,conection_type,dist_atom_par):
    num_p = np.shape(xyz)[0]
    num_a = np.shape(atomtypes_par)[0]
    xyz_new = []
    atomTypes = []
    ato = conection_type
    for i in range(num_p):
        norm = np.linalg.norm(xyz[i,:])

        # to be extended for molecules with more than 2 atoms
        xyz_new.append(xyz[i,:])
        atomTypes.append(atomtypes_par[ato%2])
        ato += 1

        pxyz = (norm+dist_atom_par)*(xyz[i,:]/norm)
        xyz_new.append(pxyz)
        atomTypes.append(atomtypes_par[ato%2])
        ato += 1
    # print(atomTypes)
    return np.array(xyz_new), atomTypes


