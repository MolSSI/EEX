# HOOMD-Blue

The HOOMD-Blue software is capable of performing molecular dynamics and hard particle Monte Carlo simulation and is
available as python modules (installation with conda)

Simulation trajectories and snapshots are saved in gsd ("general simulation data") format, which was created by the 
Glotzer group.

To install HOOMD and GSD:

```
conda config --add channel glotzer
conda install hoomd
conda install gsd
```

The GSD module contains:
* gsd.fl: file layer - "provides low level API for directly accessing gsd files."
* gsd.hoomd: Reference implementation for reading/writing hoomd schema GSD files
* gsd.pygsd: Provides a GSD reader written in pure python

## Examples of use:

### Create empty HOOMD snapshot and populate with data

The following creates a system with 4 particles of type 'A' and 'B' and sets their positions. The configuration is saved
in the file 'test.gsd'
```
import hoomd  
import gsd
import gsd.hoomd

s = gsd.hoomd.Snapshot()  
s.particles.N = 4  
s.particles.types = ['A', 'B']  
s.particles.typeid = [0,0,1,1]  
s.particles.position = [[0,0,0],[1,1,1], [-1,-1,-1], [1,-1,-1]]  
s.configuration.box = [3, 3, 3, 0, 0, 0] 
traj = gsd.hoomd.open(name='test.gsd', mode='wb')
traj.append(s)
```

Alternatively,

```
import hoomd
import numpy

# Initialize snapshot
snap = hoomd.data.make_snapshot(N=4, box=hoomd.data.boxdim(L=10), particle_types=['A', 'B'])

# Assign coordinates
snap.particles.position[0] = [1,2,3]
snap.particles.position[1] = [-1,-2,-3]
snap.particles.position[2] = [3,2,1]
snap.particles.position[3] = [-3,-2,-1]

# Define bonds
snap.bonds.resize(2);
snap.bonds.group[:] = [[0,1], [1, 2]]

# Initialize hoomd blue (can't set forces without doing this)
hoomd.init.read_snapshot(snap);

# Write snapshot
hoomd.dump.gsd("init.gsd", period=None, group=all, dynamic=None, overwrite=True);

# Set LJ parameters
nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=3.0, nlist=nl);
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=2.0);
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0);
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0);
all = hoomd.group.all();

# Define bond types
bond1 = hoomd.md.bond.harmonic(name="polymer")
bond1.bond_coeff.set('polymer', k=330.0, r0=0.84)
```

It seems like the gsd file does not contain force field information and this is set in Python script.
We will need to write a gsd containing positions of particles, bonds, angles, dihedrals, impropers, then write file
Python file with ff information.

### Capture snapshot from simulation
`snapshot = system.take_snapshot(all=True)`

### Load conditions from snapshot
`hoomd.init.read_snapshot(snap)`


## Information about snapshots:  
input: `dir(gsd.hoomd.Snapshot())`  

output: ```['__class__',
 '__delattr__',
 '__dict__',
 '__dir__',
 '__doc__',
 '__eq__',
 '__format__',
 '__ge__',
 '__getattribute__',
 '__gt__',
 '__hash__',
 '__init__',
 '__init_subclass__',
 '__le__',
 '__lt__',
 '__module__',
 '__ne__',
 '__new__',
 '__reduce__',
 '__reduce_ex__',
 '__repr__',
 '__setattr__',
 '__sizeof__',
 '__str__',
 '__subclasshook__',
 '__weakref__',
 'angles',
 'bonds',
 'configuration',
 'constraints',
 'dihedrals',
 'impropers',
 'pairs',
 'particles',
 'validate']```

