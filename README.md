# DeepComplex2
A gradient descent based approach to generate complex structures using inter chain residue-residue contacts

**(1) Download the code (short path is recommended)**

```
git clone git@github.com:jianlin-cheng/DeepComplex2.git

(If fail, try username) git clone https://github.com/jianlin-cheng/DeepComplex2.git

```

**(2) Install PyRosetta4: http://www.pyrosetta.org/dow**

**(3) Run the code for generating protein complexes using gradient descent method**

```
   Usage:
   $ sh scripts/run_dock_gd.sh <target id> <initial pdb file>.pdb <path of restraint file> <output folder> <path of weight file>

   Example:
   $ sh scripts/run_dock_gd.sh 1D3Y  /data/esdft/initial_starts/1D3Y_modified.pdb  /data/esdft/restraints/1D3Y_AB.rr  /data/esdft/output  /data/esdft/weight_files/talaris2013.wts
``` 

**(3) Run the code for generating protein complexes using Markov chain Monte Carlo**

```
   Usage:
   $ sh scripts/run_dock_mc.sh <target id> <initial pdb file>.pdb <path of restraint file> <output folder> <path of weight file>

   Example:
   $ sh scripts/run_dock_mc.sh 1D3Y  /data/esdft/initial_starts/1D3Y_modified.pdb  /data/esdft/restraints/1D3Y_AB.rr  /data/esdft/output  /data/esdft/weight_files/talaris2013.wts
``` 
