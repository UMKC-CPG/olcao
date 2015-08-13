AC_INIT([OLCAO],[0.6],[rulisp@umkc.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])


dnl Link to the appropriate library.  The dtsevx and zherk are just test cases
dnl   to make sure the libraries exist.  If they are found, then they are
dnl   included in the link stage.
dnl AC_CHECK_LIB(lapack, dstevx)
dnl AC_CHECK_LIB([blas], [zherk])
dnl AC_CHECK_LIB([veclib], [zherk])
dnl AC_CHECK_LIB([xt-libsci], [zherk])


dnl ######################################################################
dnl
dnl Require the user to specify the location of the lapack and blas
dnl
dnl ######################################################################



dnl ######################################################################
dnl
dnl File:	hdf5.m4
dnl
dnl Purpose:	Determine the locations of hdf5 includes and libraries
dnl
dnl Version: $Id: hdf5.m4,v 1.26 2003/09/15 20:36:26 cary Exp $
dnl
dnl Tech-X configure system
dnl
dnl Copyright Tech-X Corporation
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Allow the user to specify an overall hdf5 directory.  If specified,
dnl we assume that libs and include will be under this.  If not specified,
dnl then we will look for it in some typical places.  If it cannot be found
dnl then we complain.
dnl
dnl ######################################################################

dnl AC_ARG_WITH(hdf5-dir,[  --with-hdf5-dir=<location of hdf5 installation> ],HDF5_DIR="$withval",HDF5_DIR="")

dnl if test -n "HDF5_DIR"; then
dnl    AC_SUBST([HDF5_DIR])
dnl    HDF5_LIBPATH=${HDF5_DIR}/lib
dnl    HDF5_LDFLAGS='-I${HDF5_DIR}/include -L${HDF5_DIR}/lib -I${HDF5_DIR}/lib'
dnl fi

dnl ######################################################################
dnl
dnl Confirm that we have access to all the necessary files.
dnl
dnl ######################################################################
dnl AC_CHECK_FILE(${HDF5_LIBPATH}/libhdf5_fortran.so,LIBHDF5_A=${HDF5_LIBPATH}/hdf5_fortran.h,LIBHDF5_A="")
dnl if test -n "LIBHDF5_A"; then
dnl    ac_cv_have_hdf5=yes
dnl else
dnl    ac_cv_have_hdf5=no
dnl fi

dnl ######################################################################
dnl
dnl Allow perl modules to be installed as desired.
dnl
dnl ######################################################################

AC_ARG_ENABLE([perllibdir],
   [AS_HELP_STRING([--enable-perllibdir=DIR],
   [perl module install directory ($(bindir)/perl5)])],
   [perllibdir=$enableval],
   [perllibdir='$(bindir)/perl5'])
AC_SUBST([perllibdir])

dnl ######################################################################
dnl
dnl Set up the fortran compile flags.
dnl
dnl ######################################################################

dnl AC_ARG_ENABLE(debug,[  --enable-debug Turns on debugging features],\
dnl       DEBUGFLAGS="-g -C",DEBUGFLAGS="")
dnl FCFLAGS=${FCFLAGS}
dnl AC_SUBST([FCFLAGS],${FCFLAGS})

AC_ARG_VAR(FC,Fortran compiler)
AC_ARG_VAR(FCFLAGS,Fortran compiler flags)

dnl ######################################################################
dnl
dnl Define for whether hdf5 found
dnl
dnl ######################################################################

dnl if test $ac_cv_have_hdf5 = yes; then
dnl   AC_DEFINE(HAVE_HDF5,1,Define to 1 if you have the hdf5 library (-lhdf5).)
  FC=h5fc
  AC_SUBST(FC)
dnl fi


dnl ######################################################################
dnl
dnl Set up correct use of the pre-processor (C or Fortran).
dnl
dnl ######################################################################

dnl if test -z "$CPP"; then
dnl    AC_PATH_PROG(CPP,fpp)
dnl else
   AC_PROG_CPP
dnl fi

dnl if test -z "$CPPFLAGS"; then
dnl   CPPFLAGS="-Wp,-C,-P"
dnl fi

AC_SUBST(CPP)
dnl AC_SUBST(CPPFLAGS)



dnl ######################################################################
dnl
dnl All Done!   Output the results.
dnl
dnl ######################################################################
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
