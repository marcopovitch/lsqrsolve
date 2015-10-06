# Lsqrsolve

Solve *very large* (seismological) sparse system of linear equations, using lsqr from C. C. Paige and M. A. Saunders.

Use also *[mesh](https://github.com/marcopovitch/mesh)* and *[sparse](https://github.com/marcopovitch/sparse)* libraries.

# Beware
 
Still usable but written 10 years ago ! Use it only if you know what you are doing ... no support will be provided !


# Compilation

## prerequisite

Install the *[mesh](https://github.com/marcopovitch/mesh)* and *[sparse](https://github.com/marcopovitch/sparse)* libraries following the steps decribed in each relating README.md files. 

## lsqrsolve compilation  

To make all the needed *autotools* files, go into the `lsqrsolve` directory run:

`glibtoolize -i`

To configure and install the lsqrsolve files in your installation directory $INSTALL_PATH, run:

`./autogen.sh --with-mesh-prefix=$MESH_INSTALL_PATH --prefix=$INSTALL_PATH`

`make install`
 

# Links

* [http://web.stanford.edu/group/SOL/software/lsqr/](http://web.stanford.edu/group/SOL/software/lsqr/)
* [ http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c](http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c)

