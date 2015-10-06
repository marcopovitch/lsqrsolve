dnl Process this file with autoconf to produce a configure script.
AC_INIT([lsqrsolve], [0.9])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])

dnl ------------------------------------------------------------------------------------- 
dnl                mesh library and headers information  
dnl ------------------------------------------------------------------------------------- 
dnl ---- allows to specify a different prefix PATH for mesh libray
AC_ARG_WITH(mesh-prefix,
            [  --with-mesh-prefix=PREFIX Where your libmesh was installed (e.g. /usr/local/)])
	if test -n "${with_mesh_prefix}"; then
		extra_ld="${extra_ld} -L${with_mesh_prefix}/lib"
		extra_inc="${extra_inc} -I${with_mesh_prefix}/include"
	fi

dnl ----- now assert if we must check for libs and headers. 
if test "x${with_mesh_prefix}" = "x"; then
dnl ---------------------------------------------------------------------
dnl                 Checks for mesh presence 
dnl ---------------------------------------------------------------------

AC_CHECK_LIB(mesh,mesh_init_from_file,
       echo "+------------------------------------------------------------+";
      echo "| libmesh found. Good. using it.                             |";
       echo "+------------------------------------------------------------+",
       echo "+------------------------------------------------------------+";
       echo "| libmesh is badly or not installed. Will not build !        |";
       echo "+------------------------------------------------------------+",
	 )

dnl ---- check presence of headers 
AC_CHECK_HEADERS(mesh/mesh.h,
       echo "+------------------------------------------------------------+";
	 echo "| headers files for libmesh found.                           |";
       echo "+------------------------------------------------------------+",
       echo "+------------------------------------------------------------+";
	 echo "| headers files for libmesh not found. Wont build.           |";
       echo "+------------------------------------------------------------+")
fi


dnl ------------------------------------------------------------------------------------- 
dnl                popt location  
dnl ------------------------------------------------------------------------------------- 
dnl ---- allows to specify a different prefix PATH for libray
AC_ARG_WITH(popt-prefix,
      [--with-popt-prefix=PREFIX Where your popt library was installed (e.g. /usr/local/)])
      if test -n "${with_popt_prefix}"; then
	              extra_ld="-L${with_popt_prefix}/lib"
	              extra_inc="-I${with_popt_prefix}/include"
      fi


pkgconfig_invoke="pkg-config"
case "${host}" in
   i[[3456]]86-*-linux-gnu*)
          syst="Linux"
          extra_inc="${extra_inc} -Wall"
          pkgconfig_invoke="pkg-config"
          
   ;;
   mip*-*-irix*)
          syst="Irix"
          extra_inc="${extra_inc} -I/usr/freeware/include"
          extra_ld="${extra_ld} -L/usr/freeware/lib32"
          pkgconfig_invoke="pkg-config --define-variable=ROOT=/ --define-variable=ABILIB-lib32=lib32"
   ;;
esac


dnl ---------------- Documentation generation ----------------------------------
DOC_GENERATOR=true ## later write a macro to see if doxygen is there

dnl Let people disable the doxygen doc generation.
AC_ARG_ENABLE(doc, [  --enable-doc  Use doc to build doxygen documentation [default=no]], enable_doc="$enableval", enable_doc=no)

if test x$enable_doc = xyes ; then
   if test x$DOC_GENERATOR = xtrue ; then
      enable_doc=yes
   else
      enable_doc=no
   fi
fi

AM_CONDITIONAL(ENABLE_DOC, test x$enable_doc = xyes)
AC_SUBST(extra_inc)
AC_SUBST(extra_ld)
AC_SUBST(pkgconfig_invoke)

echo "--------------------------Summary-----------------------------"
echo " System                : ${syst}"
echo " C compiler invocation : ${CC}                       "
echo " C Flags               : ${CFLAGS} ${extra_inc}      "
echo " LD flags              : ${LD_ADD} ${extra_ld}       "
echo " defines               : ${extra_defines}                 "
echo " pk-config invocation  : ${pkgconfig_invoke}              "
echo " Documentation build   : ${enable_doc}                    "
echo "--------------------------------------------------------------"

AC_OUTPUT([
Makefile
src/Makefile
lsqrsolve.spec
])

dnl share/Makefile
dnl doc/Makefile
dnl doc/html/Makefile
