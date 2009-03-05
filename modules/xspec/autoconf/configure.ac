AC_INIT(src/xspec-module.c)
AC_PREFIX_DEFAULT(..)
JH_CONFIG_DIR_SANS_AMD
JD_INIT

AC_CONFIG_HEADER(src/config.h:src/config.hin)

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

AC_FC_LIBRARY_LDFLAGS

JH_FC_DEFS
JH_FORCE_UNDEFINED_SYMBOLS

JH_SYS_EXTRA_LIBS

# Extra Headers for xlibs
AC_PATH_XTRA

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
sys/stat.h \
sys/types.h \
)

ISIS_SRCDIR=`cd ../..;pwd`
AC_SUBST(ISIS_SRCDIR)

dnl This must precede the JD_WITH_LIBRARY() macros
dnl for slang, cfitsio and pgplot
JH_WITH_CIAO

if test "x$jh_use_ciao" = "xno"
then

dnl Uncomment this to try finding slang automagically
dnl JD_FIND_SLANG
JD_WITH_LIBRARY(slang)

fi

JD_SET_OBJ_SRC_DIR(src)
JD_GCC_WARNINGS
JH_PURIFY

JH_WITH_XSPEC_STATIC

AC_ARG_ENABLE(xspec_tables,
  [  --disable-xspec_tables  compile without XSPEC table models],
  [ ],
  [AC_DEFINE(HAVE_XSPEC_TABLE_MODELS)])

JD_GET_MODULES(src/modules.lis,src)

AC_OUTPUT(Makefile:autoconf/Makefile.in src/Makefile)

dnl JD_CREATE_MODULE_ORULES(src)
dnl JD_APPEND_RULES(src/Makefile)
dnl JD_APPEND_ELFRULES(src/Makefile)
