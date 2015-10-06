## Process this file with automake to produce Makefile.in
 
#SUBDIRS = src doc share
SUBDIRS = src
 
lsqrsolvedocdir = ${prefix}/doc/lsqrsolve
lsqrsolvedoc_DATA = \
	           README\
	           COPYING\
	           AUTHORS\
	           ChangeLog\
	           INSTALL\
	           NEWS\
	           TODO
 
EXTRA_DIST = lsqrsolve.spec.in lsqrsolve.spec
 
rpms : distcheck
	        rpm -ta @PACKAGE@-@VERSION@.tar.gz

