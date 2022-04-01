#! /usr/bin/env python
# coding: utf-8

import pyemma
import numpy as np
import glob
import os
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import pyemma.coordinates as coor
import sys
import mdtraj as md
import itertools


npz_file = ['trajectory_adapt0.npz', 'trajectory_adapt1.npz', 'trajectory_beta2143.npz', 'trajectory_fold_unfold.npz', 'trajectory_native12.npz', 'trajectory_native3.npz']
save_file = ['checkshifted12Q_adapt0.npz', 'checkshifted12Q_adapt1.npz', 'checkshifted12Q_beta2143.npz', 'checkshifted12Q_fold_unfold.npz', 'checkshifted12Q_native12.npz', 'checkshifted12Q_native3.npz']                           

traj_npz = np.load(npz_file[0])
trajfile = traj_npz['files'].tolist()
topfile = str(traj_npz['top'])

def best_hummer_q(trajfile, native,selected=False):
    BETA_CONST = 50  # 1/nm    
    LAMBDA_CONST = 1.8         
    NATIVE_CUTOFF = 0.45  # nanometers

    # get the indices of all of the heavy atoms
    trajs=md.load(trajfile[0],top=native)
    feat = coor.featurizer(native)
    native=md.load_pdb(native) 
    if not selected:
        heavy = trajs.top.select_atom_indices('heavy')
    else:
        heavy = trajs.top.select_atom_indices('heavy')[selected]
    # get the pairs of heavy atoms which are farther than 3
    # residues apart           
    heavy_pairs = np.array(    
     [(i,j) for (i,j) in itertools.combinations(heavy, 2)
         if abs(trajs.top.atom(i).residue.index - \
                trajs.top.atom(j).residue.index) > 3])

    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts: ", len(native_contacts))

    # now compute these distances for the whole trajectory
    feat.add_distances(native_contacts,periodic=False)

    r = pyemma.coordinates.load(trajfile, features=feat)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)
    q = [np.mean(1.0 / (1 + np.exp(BETA_CONST * (i- LAMBDA_CONST * r0))), axis=1) for i in r]
    return q                   
for i in range(len(npz_file)):
  traj_npz = np.load(npz_file[i])
  trajfile = traj_npz['files'].tolist()    
  topfile = str(traj_npz['top'])
  topfile = './check_20317.pdb'  
  q=best_hummer_q(trajfile, topfile,selected=list(range(12,64))+list(range(100,152)))
  np.savez(save_file[i],q)
#[0,160];[167,282];[-122:-1]
