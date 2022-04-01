#! /usr/bin/env python
# coding: utf-8

import pyemma
import numpy as np
import glob
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import pyemma.coordinates as coor
import sys


npz_file = ['trajectory_adapt0.npz','trajectory_adapt1.npz','trajectory_native12.npz']
save_file = ['Allfeatures_adapt0.npz','Allfeatures_adapt1.npz', 'Allfeatures_native12.npz']
for index,npz in enumerate(npz_file):
    traj_npz = np.load(npz)
    trajfile = traj_npz['files'].tolist()
    topfile = str(traj_npz['top'])
    trajfile = trajfile
     
    feat = coor.featurizer(topfile)
    residues=[i for i in range(8)]+[i for i in range(12,20)]+[i for i in range(21,36)]+[i for i in range(41,46)]+[i for i in range(50,55)]
    pairs = []
    for i,r1 in enumerate(residues):
        for r2 in residues[i+1::2]:
            pairs.append([r1,r2])
    pairs = np.array(pairs)
    print(len(pairs))
    feat.add_residue_mindist(pairs, scheme='closest-heavy',periodic=False)
    X = coor.load(trajfile, feat)
    np.savez(save_file[index],X)                                                                                                                         
