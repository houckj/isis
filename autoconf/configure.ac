AC_INIT
AC_CONFIG_SRCDIR([src/isis.c])
AC_PREFIX_DEFAULT(/usr/local)
JH_CONFIG_DIR_SANS_AMD
JD_INIT

JH_WITH_SUBDIR

AC_CONFIG_HEADER(\
 src/config.h:src/config.hin \
 test/config.h:test/config.hin \
 modules/xspec/src/config.h:modules/xspec/src/config.hin\
 modules/cfitsio/src/config.h:modules/cfitsio/src/config.hin\
 modules/pgplot/src/config.h:modules/pgplot/src/config.hin\
 modules/maplib/src/config.h:modules/maplib/src/config.hin\
)

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_PROG_FC

JD_ANSI_CC
JD_ELF_COMPILER
JD_IEEE_CFLAGS
ELF_CFLAGS="$ELF_CFLAGS $IEEE_CFLAGS"
CFLAGS="$CFLAGS $IEEE_CFLAGS"

# Extra Headers for xlibs
AC_PATH_XTRA

# Fixup X_LIBS definition produced by AC_PATH_XTRA
case "$host_os" in
  # Fix rpath issue on Macs:
  # The Apple/Xcode gcc understands the -R option but the
  # hpc.sourceforge.net compilers (gcc, g++, gfortran) do not.
  # When using the Apple gcc, AC_PATH_XTRA produces e.g.
  #     X_LIBS=-L/usr/X11/lib -R/usr/X11/lib
  # Later, when linking with the sourceforge gfortran, the
  # '-R' will will not be recognized, causing an error.
  # To avoid the link error, we replace "-R" with "-Wl,-rpath,".
  # The latter option is currently accepted by both Apple
  # (gcc-4.2.1) and hpc.sourceforge.net (gcc-4.6.1) compilers.
  *darwin* )
     X_LIBS=`echo $X_LIBS | sed s/-R/-Wl,-rpath,/`
     AC_SUBST(X_LIBS)
  ;;
esac

AC_FC_WRAPPERS
AC_FC_LIBRARY_LDFLAGS
AC_FC_DUMMY_MAIN(,
  AC_DEFINE(FC_DUMMY_MAIN,MAIN__)
)
SAFE_AC_FC_MANGLE

AC_ARG_WITH(fclibs,
  [  --with-fclibs=libs            Provide extra Fortran libraries, e.g. "-lg2c" ],
  [jh_fclibs_extra=$withval], [jh_fclibs_extra=no])
FC_EXTRA_LIBS=""
if test "x$jh_fclibs_extra" != "xno"
then
   FC_EXTRA_LIBS="$jh_fclibs_extra"
fi
AC_SUBST(FC_EXTRA_LIBS)

JD_LARGE_FILE_SUPPORT

AC_CHECK_TYPES(unsigned long long)
AC_CHECK_SIZEOF(unsigned long long)

AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)

JH_INIT_MODULE_LIST

# initialization
SLANG_LIB=no
CFITSIO_LIB=no
PGPLOT_LIB=no
CCFITS_LIB=no

# reasonable pgplot default
PGPLOT_LIBDIR=/usr/local/pgplot
if test -d $PGPLOT_LIBDIR
then
  PGPLOT_LIB="-L$PGPLOT_LIBDIR"
  PGPLOT_INC="-I$PGPLOT_LIBDIR"
fi

ATOMDB_DIR=""
HEADAS_DIR=""
XSPEC_MODELSONLY=""

dnl This must precede the JD_WITH_LIBRARY() macros
dnl for slang, cfitsio and pgplot
JH_HANDLE_PACKAGE_OPTIONS

if test "x$PGPLOT_LIB" = "xno" \
   || test "${with_pgplot+set}" = set \
   || test "${with_pgplotlib+set}" = set \
   || test "${with_pgplotinc+set}" = set
then
  JD_WITH_LIBRARY(pgplot, cpgplot.h)
  JH_CHECK_PGPLOT_LIBNAME($jd_pgplot_library_dir)
  PGPLOT_LIBDIR="$jh_pgplot_library_dir"
  AC_SUBST(PGPLOT_LIBDIR)
fi

if test "x$CFITSIO_LIB" = "xno" \
   || test "${with_cfitsio+set}" = set \
   || test "${with_cfitsiolib+set}" = set \
   || test "${with_cfitsioinc+set}" = set
then
  JD_WITH_LIBRARY(cfitsio, fitsio.h)
  JH_CHECK_CFITSIO_LIBNAME($jd_cfitsio_library_dir)
fi

if test "x$SLANG_LIB" = "xno"
then
  JD_WITH_LIBRARY(slang, slang.h)
  SLANG_LIBDIR="$jd_slang_library_dir"
  AC_SUBST(SLANG_LIBDIR)
fi

AC_ARG_WITH(local-setup,
  [  --with-local-setup=FILE      Load site-local customizations from FILE ],
  [ jh_local_setup=$withval ], [jh_local_setup=""])
SITE_LOCAL_SETUP=$jh_local_setup
AC_SUBST(SITE_LOCAL_SETUP)

AC_ARG_WITH(atomdb,
  [  --with-atomdb[=DIR]            location of spectroscopy database (ATOMDB)],
  [jd_use_atomdb=$withval], [jd_use_atomdb=no])
if test "x$jd_use_atomdb" != "xno"
then
   ATOMDB_DIR=$jd_use_atomdb
fi
AC_SUBST(ATOMDB_DIR)

AC_ARG_WITH(djbfft,
[  --with-djbfft=DIR            Use djbfft for pileup FFTs
                                (from DIR/lib and DIR/include)],
  [ jh_use_djbfft=$withval ], [jh_use_djbfft=no])
if test "x$jh_use_djbfft" != "xno"
then
   EXTRA_LIB="$EXTRA_LIB $jh_use_djbfft/lib/djbfft.a"
   EXTRA_INC="$EXTRA_INC -I$jh_use_djbfft/include"
   AC_SUBST(EXTRA_LIB)
   AC_SUBST(EXTRA_INC)
   AC_DEFINE(HAVE_DJBFFT)
fi

AC_ARG_ENABLE(extras,
  [  --disable-extras        Don't install extras ],
  [ jh_extras_enable=$enableval ], [jh_extras_enable=yes])
BUILD_ISIS_EXTRAS=no
if test X"$jh_extras_enable" = X"yes" ; then
   if test -d extras ; then
     BUILD_ISIS_EXTRAS=yes
     AC_DEFINE(HAVE_ISIS_EXTRAS)
     ac_configure_args="$ac_configure_args --with-slanginc=$jd_slang_include_dir --with-slanglib=$jd_slang_library_dir"
     AC_CONFIG_COMMANDS_POST([test -d extras && ./config.status --file=extras/Makefile:extras/Makefile.in])
   else
     echo 'No extras directory - ok.'
   fi
fi
AC_SUBST(BUILD_ISIS_EXTRAS)

JH_TERMCAP_LIB
JH_WITH_READLINE

JH_SYS_EXTRA_LIBS

MODULE_INSTALL_DIR='${prefix}/lib/modules'
SL_FILES_INSTALL_DIR='${prefix}/share'
AC_SUBST(MODULE_INSTALL_DIR)
AC_SUBST(SL_FILES_INSTALL_DIR)

dnl Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
sys/stat.h \
sys/wait.h \
sys/types.h \
dlfcn.h \
ieeefp.h \
)

AC_CHECK_FUNCS(\
stat \
sigaction \
sigemptyset \
sigprocmask \
sigaddset \
vsnprintf \
isinf \
isnan \
finite \
)

JD_SET_OBJ_SRC_DIR(src)
JD_GCC_WARNINGS
JH_PURIFY

JD_GET_MODULES(src/modules.lis,src)
JH_GET_FC_MODULES(src/fc_modules.lis,src)
JH_GET_SHAREFILES(share/modules.lis)
JH_GET_ETCFILES(etc/modules.lis)

AC_CONFIG_FILES([\
 Makefile:autoconf/Makefile.in \
 src/Makefile \
 etc/local.sl:etc/local.sl.in \
 test/Makefile \
 test/uncon/Makefile \
 modules/Makefile \
 modules/cfitsio/Makefile:modules/cfitsio/autoconf/Makefile.in \
 modules/cfitsio/src/Makefile:modules/cfitsio/src/Makefile.in \
 modules/pgplot/Makefile:modules/pgplot/autoconf/Makefile.in \
 modules/pgplot/src/Makefile:modules/pgplot/src/Makefile.in \
 modules/xspec/Makefile:modules/xspec/autoconf/Makefile.in \
 modules/xspec/src/Makefile:modules/xspec/src/Makefile.in \
 modules/xspec/src/xspec_libs:modules/xspec/src/xspec_libs.in \
 modules/maplib/Makefile:modules/maplib/autoconf/Makefile.in \
 modules/maplib/src/Makefile:modules/maplib/src/Makefile.in \
])
AC_OUTPUT

dnl JD_CREATE_EXEC_RULE(isis)
JD_CREATE_FLINK_EXEC_RULE(isis)
JD_CREATE_ELFORULE(isis)

JD_CREATE_MODULE_ORULES(src)
JD_APPEND_RULES(src/Makefile)
JD_APPEND_ELFRULES(src/Makefile)

JH_CREATE_FC_MODULE_ORULES(src)
JH_APPEND_FC_RULES(src/Makefile)
JH_APPEND_FC_ELFRULES(src/Makefile)

echo ""
echo "You are compiling with this compiler configuration:"
echo "       CC =" "$CC"
echo "   CFLAGS =" "$CFLAGS"
echo "  LDFLAGS =" "$LDFLAGS" "$DYNAMIC_LINK_FLAGS"
echo "       FC =" "$FC"
echo "  FCFLAGS =" "$FCFLAGS"
echo ""
if test "x$HEADAS_DIR" != "x" ; then
    echo "XSPEC module configuration:"
    echo "    --with-xspec-version = ${with_xspec_version}"
    if test x"${XSPEC_MODELSONLY}" = xyes ; then
       headas_arg="xspec-modelsonly"
    else
       headas_arg="headas"
    fi
    echo "    --with-${headas_arg} = $HEADAS_DIR"
    ./autoconf/compiler_check.sh "$CC" "$FC" "$HEADAS_DIR"
else
    echo "The XSPEC module will not be compiled."
    ./autoconf/compiler_check.sh "$CC" "$FC"
fi
    echo ""
echo "isis will be installed in" "$prefix"
echo ""
echo "If $prefix_input/bin exists,"
echo "a symbolic link to the isis executable will be placed there."
echo ""
if test "x$SITE_LOCAL_SETUP" != "x" ; then
    echo "Site-local setup script"
    echo "    --with-local-setup = $SITE_LOCAL_SETUP"
else
    echo "No site-local setup script was specified."
fi
echo ""
echo "To see the list of configuration options:"
echo "         ./configure --help"
echo ""
