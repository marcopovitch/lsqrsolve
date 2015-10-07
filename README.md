# Lsqrsolve

it solves *very large* (seismological) sparse system of linear equations, using lsqr (from C. C. Paige and M. A. Saunders), allowing positive *damping*.

The core algorithms (lsqr.c, lsqr.h) are provided by SOL (System Optimisation Laboratory), from Stanford University. The source code from C. C. Paige and M. A. Saunders is avalaible here :
	http://www.stanford.edu/group/SOL/software.html

I have added and wrapped around a *[sparse](https://github.com/marcopovitch/sparse)* matrix implementation using linked list (row and column) of non nul elements, to speed up things. Some time ago, I have developed and used it for global sesimic travel time tomography.

lsqrsolve code is a very specific tool (may works with *ray2mesh*) but 
I provide also some sample codes (rse, rse2, rse3) which implement generic
lsqr solver :

`A x = B`

(A could be sparse or not) in command line.

Use also *[mesh](https://github.com/marcopovitch/mesh)* and *[sparse](https://github.com/marcopovitch/sparse)* libraries.

# Beware
 
Still usable but written 10 years ago ! Use it only if you know what you are doing ... no support will be provided !

So far, no code parallelization (should be easy).


# Compilation

### prerequisite

Install the *[mesh](https://github.com/marcopovitch/mesh)* (trying to remove that dependency ) and *[sparse](https://github.com/marcopovitch/sparse)* libraries following the steps decribed in each relating README.md files. 

### lsqrsolve compilation  

To make all the needed *autotools* files, go into the `lsqrsolve` directory run:

`glibtoolize -i`

To configure and install the lsqrsolve files in your installation directory $INSTALL_PATH, run:

`./autogen.sh --with-mesh-prefix=$MESH_INSTALL_PATH --prefix=$INSTALL_PATH`

`make install`
 
 
# Files format
<pre>
A x = B
</pre>

All indices *i* and *j* start with *0* and not *1* !

### Matrice A

File for sparse matrix (A) should only contains non null *value* such as:
<pre>
ni nj
i j value
...
</pre>

### Vector B

File format for vector (B) is :
<pre>
nj
value[0]
...
value[j-1]
</pre>

# Example
`rse3 sparseA.txt B.txt x.txt nbiter_max damping [dump_sol_iter]`

where :

- `sparseA.txt`   : sparse matrix (*input*)
- `B.txt`         : standard vector (*input*)
- `x.txt`         : solution (*output*)
- `nbiter_max`   : nb maximun iteration to perform
- `damping`      : damping value to use
- `dump_sol_iter` : dump solution for each `dump_sol_iter` iteration in iter-%d.txt file (optional)

# Tips
When using sparse matrix please sort your data first. It speeds up things (a lot) !
You can use the sort unix command line utility :

`sort -n -k1 -k2.1n sparseA.txt > sparseA-sorted.txt`




# Links

* [http://web.stanford.edu/group/SOL/software/lsqr/](http://web.stanford.edu/group/SOL/software/lsqr/)
* [ http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c](http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c)

