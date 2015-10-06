# Lsqrsolve

it solves *very large* (seismological) sparse system of linear equations, using lsqr (from C. C. Paige and M. A. Saunders), allowing positive *damping*.

The core algorithms (lsqr.c, lsqr.h) are provided by SOL (System Optimisation Laboratory), from Stanford University. The source code from C. C. Paige and M. A. Saunders is avalaible here :
	http://www.stanford.edu/group/SOL/software.html

I have added and wrapped around a *[sparse](https://github.com/marcopovitch/sparse)* matrix implementation using linked list (row and column) of non nul elements, to speed up things. Some time ago, I have developed and used it for global sesimic travel time tomography.

lsqrsolve code is a very specific tool (may works with *ray2mesh*) but 
I provide also some sample codes (rse, rse2, rse3) which implement generic
lsqr solver Ax=B (A could be sparse or not) in command line.

Use also *[mesh](https://github.com/marcopovitch/mesh)* and *[sparse](https://github.com/marcopovitch/sparse)* libraries.

# Beware
 
Still usable but written 10 years ago ! Use it only if you know what you are doing ... no support will be provided !

So far, no code parallelization (should be easy).


# Compilation

## prerequisite

Install the *[mesh](https://github.com/marcopovitch/mesh)* (trying to remove that dependency ) and *[sparse](https://github.com/marcopovitch/sparse)* libraries following the steps decribed in each relating README.md files. 

## lsqrsolve compilation  

To make all the needed *autotools* files, go into the `lsqrsolve` directory run:

`glibtoolize -i`

To configure and install the lsqrsolve files in your installation directory $INSTALL_PATH, run:

`./autogen.sh --with-mesh-prefix=$MESH_INSTALL_PATH --prefix=$INSTALL_PATH`

`make install`
 

# Links

* [http://web.stanford.edu/group/SOL/software/lsqr/](http://web.stanford.edu/group/SOL/software/lsqr/)
* [ http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c](http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c)

