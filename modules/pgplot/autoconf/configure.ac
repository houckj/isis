dnl -*- sh -*-
AC_INIT(src/pgplot-module.c)
AC_PREFIX_DEFAULT(/usr/local)

AC_CONFIG_AUX_DIR(autoconf)
AC_CANONICAL_HOST
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

JD_INIT
JD_ANSI_CC
JD_ELF_COMPILER
JD_IEEE_CFLAGS

dnl JD_FC_COMPILER
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_FC_DUMMY_MAIN(,
  AC_DEFINE(FC_DUMMY_MAIN,MAIN__)
)

AC_PATH_XTRA

JD_WITH_LIBRARY(slang)

dnl# Add libraries here
JD_WITH_LIBRARY(pgplot, cpgplot.h)

dnl# JD_WITH_LIBRARY(png)

dnl# This idea was stolen from John Houck
OTHER_PGPLOT_LIBS=
pgplot_libs_check="passed"
AC_MSG_CHECKING(for addition libraries required by pgplot)
AC_CHECK_LIB(png, cpgvstd, 
    [OTHER_PGPLOT_LIBS="$OTHER_PGPLOT_LIBS -lpng"],
    [pgplot_libs_check="failed"],
    [$PGPLOT_LIB -lcpgplot -lpgplot $OTHER_PGPLOT_LIBS $X_LIBS -lX11 $X_EXTRA_LIBS $FLIBS])
AC_SUBST(OTHER_PGPLOT_LIBS)

dnl# This macro inits the module installation dir
JD_SLANG_MODULE_INSTALL_DIR

dnl Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)

JD_SET_RPATH($libdir)
ELF_CFLAGS="$ELF_CFLAGS $IEEE_CFLAGS"
CFLAGS="$CFLAGS $IEEE_CFLAGS"

AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_OUTPUT(Makefile:autoconf/Makefile.in src/Makefile)

echo ""
echo "You are compiling with the following compiler configuration:"
echo "        CC =" "$CC"
echo " CC_SHARED =" "$CC_SHARED"
echo "    CFLAGS =" "$CFLAGS"
echo "   LDFLAGS =" "$LDFLAGS" "$DYNAMIC_LINK_FLAGS"
echo ""
echo "The modules will be installed in $MODULE_INSTALL_DIR."
echo "Any associated .sl files will be install in $SL_FILES_INSTALL_DIR"
echo ""
echo "If any of these quantities are incorrect, edit src/Makefile accordingly."
echo ""
if test "$pgplot_libs_check" = "failed"
then
   echo "*************************************************************************"
   echo "Additional libraries may be required.  You may need to edit src/Makefile."
   echo "*************************************************************************"
fi
