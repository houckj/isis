dnl# -*- mode: sh; mode: fold -*-

AC_DEFUN(JD_INIT,     dnl#{{{
[
#These variable are initialized by JD init function
CONFIG_DIR=`pwd`
cd $srcdir
if test "`pwd`" != "$CONFIG_DIR"
then 
  AC_MSG_ERROR("This software does not support configuring from another directory.   See the INSTALL file")
fi
dnl# if test "X$PWD" != "X"
dnl# then
dnl#  CONFIG_DIR="$PWD"
dnl# fi
AC_SUBST(CONFIG_DIR)dnl
# Note: these will differ if one is a symbolic link
if test -f /usr/bin/dirname; then
  JD_Above_Dir=`dirname $CONFIG_DIR`
else
# system is a loser
  JD_Above_Dir=`cd ..;pwd`
fi
JD_Above_Dir2=`cd ..;pwd`
])
dnl#}}}

AC_DEFUN(JD_SET_OBJ_SRC_DIR, dnl#{{{
[
#---------------------------------------------------------------------------
# Set the source directory and object directory.   The makefile assumes an
# absolute path name.  This is because src/Makefile cds to OBJDIR and compiles
# the src file which is in SRCDIR
#---------------------------------------------------------------------------
SRCDIR=$CONFIG_DIR
if test "$1" != "."
then
  if test -z "$1"
  then
    SRCDIR=$SRCDIR/src
  else
    SRCDIR=$SRCDIR/$1
  fi
fi

OBJDIR=$SRCDIR/"$ARCH"objs
ELFDIR=$SRCDIR/elf"$ARCH"objs
AC_SUBST(SRCDIR)dnl
AC_SUBST(OBJDIR)dnl
AC_SUBST(ELFDIR)dnl
])
dnl#}}}

RPATH=""
AC_DEFUN(JD_INIT_RPATH, dnl#{{{
[
dnl# determine whether or not -R or -rpath can be used
case "$host_os" in
  *linux*|*solaris* )
    if test "X$GCC" = Xyes
    then
      if test "X$ac_R_nospace" = "Xno"
      then
        RPATH="-Wl,-R,"
      else
        RPATH="-Wl,-R"
      fi
    else
      if test "X$ac_R_nospace" = "Xno"
      then
        RPATH="-R "
      else
	RPATH="-R"
      fi
    fi
  ;;
  *osf*)
    if test "X$GCC" = Xyes
    then
      RPATH="-Wl,-rpath,"
    else
      RPATH="-rpath "
    fi
  ;;
esac
])

dnl#}}}

AC_DEFUN(JD_SET_RPATH, dnl#{{{
[
if test "X$1" != "X"
then
  if test "X$RPATH" = "X"
  then 
    JD_INIT_RPATH
    if test "X$RPATH" != "X"
    then
      RPATH="$RPATH$1"
    fi
  else
    RPATH="$RPATH:$1"
  fi
fi
])
AC_SUBST(RPATH)dnl

dnl#}}}

AC_DEFUN(JD_UPPERCASE, dnl#{{{
[
changequote(<<, >>)dnl
define(<<$2>>, translit($1, [a-z], [A-Z]))dnl
changequote([, ])dnl
])
#}}}

AC_DEFUN(JD_SIMPLE_LIB_DIR, dnl#{{{
[
JD_UPPERCASE($1,JD_UP_NAME)
JD_UP_NAME[]_LIB_DIR=$JD_Above_Dir/$1/libsrc/"$ARCH"objs
JD_UP_NAME[]_INCLUDE=$JD_Above_Dir/$1/libsrc

if test ! -d "[$]JD_UP_NAME[]_INCLUDE"
then
   JD_UP_NAME[]_LIB_DIR=$JD_Above_Dir/$1/src/"$ARCH"objs
   JD_UP_NAME[]_INCLUDE=$JD_Above_Dir/$1/src
   if test ! -d "[$]JD_UP_NAME[]_INCLUDE"
   then
     echo ""
     echo WARNING------Unable to find the JD_UP_NAME directory
     echo You may have to edit $CONFIG_DIR/src/Makefile.
     echo ""
   fi
fi

AC_SUBST(JD_UP_NAME[]_LIB_DIR)dnl
AC_SUBST(JD_UP_NAME[]_INCLUDE)dnl
undefine([JD_UP_NAME])dnl
])


dnl#}}}

AC_DEFUN(JD_FIND_GENERIC, dnl#{{{
[
changequote(<<, >>)dnl
define(<<JD_UP_NAME>>, translit($1, [a-z], [A-Z]))dnl
changequote([, ])dnl
# Look for the JD_UP_NAME package
#JD_UP_NAME[]_INCLUDE=""
#JD_UP_NAME[]_LIB_DIR=""

# This list consists of "include,lib include,lib ..."
JD_Search_Dirs="$JD_Above_Dir2/$1/libsrc,$JD_Above_Dir2/$1/libsrc/"$ARCH"objs \
                $JD_Above_Dir/$1/libsrc,$JD_Above_Dir/$1/libsrc/"$ARCH"objs \
		$JD_Above_Dir2/$1/src,$JD_Above_Dir2/$1/src/"$ARCH"objs \
                $JD_Above_Dir/$1/src,$JD_Above_Dir/$1/src/"$ARCH"objs"

test "x$prefix" = "xNONE" && prefix="$ac_default_prefix"
test "x$exec_prefix" = "xNONE" && exec_prefix="$prefix"
JD_Search_Dirs="$JD_Search_Dirs \
                $includedir,$libdir \
		$prefix/include,$exec_prefix/lib \
		$HOME/include,$HOME/lib"

if test -n "$ARCH"
then
 JD_Search_Dirs="$JD_Search_Dirs $HOME/include,$HOME/$ARCH/lib"
 JD_Search_Dirs="$JD_Search_Dirs $HOME/include,$HOME/sys/$ARCH/lib"
fi

# Now add the standard system includes.  The reason for doing this is that 
# the other directories may have a better chance of containing a more recent
# version.

test "x$exec" = "xNONE" && exec="$ac_default_prefix"
test "x$exec_prefix" = "xNONE" && exec_prefix="$prefix"
JD_Search_Dirs="$JD_Search_Dirs \
                /usr/local/include,/usr/local/lib \
		/usr/include,/usr/lib \
		/usr/include/$1,/usr/lib \
		/usr/include/$1,/usr/lib/$1"

echo looking for the JD_UP_NAME library

for include_and_lib in $JD_Search_Dirs
do
  # Yuk.  Is there a better way to set these variables??
  generic_include=`echo $include_and_lib | tr ',' ' ' | awk '{print [$]1}'`
  generic_lib=`echo $include_and_lib | tr ',' ' ' | awk '{print [$]2}'`
  echo Looking for $1.h in $generic_include
  echo and lib$1.a in $generic_lib
  if test -r $generic_include/$1.h && test -r $generic_lib/lib$1.a
  then
    echo Found it.
    JD_UP_NAME[]_LIB_DIR="$generic_lib"
    JD_UP_NAME[]_INCLUDE="$generic_include"
    break
  else
    if test -r $generic_include/$1.h && test -r $generic_lib/lib$1.so
    then
      echo Found it.
      JD_UP_NAME[]_LIB_DIR="$generic_lib"
      JD_UP_NAME[]_INCLUDE="$generic_include"
      break
    fi
  fi
done

if test -n "[$]JD_UP_NAME[]_LIB_DIR"
then
    jd_have_$1="yes"
else
    echo Unable to find the $JD_UP_NAME library.  
    echo You may have to edit $CONFIG_DIR/src/Makefile.
    JD_UP_NAME[]_INCLUDE=$JD_Above_Dir/$1/src
    JD_UP_NAME[]_LIB_DIR=$JD_Above_Dir/$1/src/"$ARCH"objs
    jd_have_$1="no"
fi

JD_UP_NAME[]_INC="-I[$]JD_UP_NAME[]_INCLUDE"
JD_UP_NAME[]_LIB="-L[$]JD_UP_NAME[]_LIB_DIR"
JD_SET_RPATH([$]JD_UP_NAME[]_LIB_DIR)
dnl if test "X$GCC" = Xyes
dnl then
dnl    RPATH_[]JD_UP_NAME="-Wl,-R[$]JD_UP_NAME[]_LIB_DIR"
dnl else
dnl    RPATH_[]JD_UP_NAME="-R[$]JD_UP_NAME[]_LIB_DIR"
dnl fi

# gcc under solaris is often not installed correctly.  Avoid specifying
# -I/usr/include.
if test "[$]JD_UP_NAME[]_INC" = "-I/usr/include"
then
    JD_UP_NAME[]_INC=""
fi

if test "[$]JD_UP_NAME[]_LIB" = "-L/usr/lib"
then
    JD_UP_NAME[]_LIB=""
    RPATH_[]JD_UP_NAME=""
fi

AC_SUBST(JD_UP_NAME[]_LIB)dnl
AC_SUBST(JD_UP_NAME[]_INC)dnl
AC_SUBST(JD_UP_NAME[]_LIB_DIR)dnl
AC_SUBST(JD_UP_NAME[]_INCLUDE)dnl
dnl AC_SUBST(RPATH_[]JD_UP_NAME)dnl
undefine([JD_UP_NAME])dnl
])


dnl#}}}

AC_DEFUN(JD_FIND_SLANG, dnl#{{{
[
JD_FIND_GENERIC(slang)
])

dnl#}}}

AC_DEFUN(JD_GCC_WARNINGS, dnl#{{{
[
AC_ARG_ENABLE(warnings,
	      [  --enable-warnings       turn on GCC compiler warnings],
	      [gcc_warnings=$enableval])
if test -n "$GCC"
then
  CFLAGS="$CFLAGS -fno-strength-reduce"
  if test -n "$gcc_warnings"
  then
    CFLAGS="$CFLAGS -Wall -W -pedantic -Winline -Wmissing-prototypes \
 -Wnested-externs -Wpointer-arith -Wcast-align -Wshadow -Wstrict-prototypes"
    # Now trim excess whitespace
    CFLAGS=`echo $CFLAGS`
  fi
fi
])


dnl#}}}

IEEE_CFLAGS=""
AC_DEFUN(JD_IEEE_CFLAGS, dnl#{{{
[
case "$host_cpu" in
  *alpha* )
    if test "$GCC" = yes
    then
      IEEE_CFLAGS="-mieee"
    else
      IEEE_CFLAGS="-ieee_with_no_inexact"
    fi
    ;;
  * )
    IEEE_CFLAGS=""
esac
])


dnl#}}}

AC_DEFUN(JD_CREATE_ORULE, dnl#{{{
[
PROGRAM_OBJECT_RULES="$PROGRAM_OBJECT_RULES
\$(OBJDIR)/$1.o : \$(SRCDIR)/$1.c \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(OBJDIR); \$(COMPILE_CMD) \$("$1"_C_FLAGS) \$(SRCDIR)/$1.c
"
])

dnl#}}}

AC_DEFUN(JD_CREATE_ELFORULE, dnl#{{{
[
PROGRAM_ELF_ORULES="$PROGRAM_ELF_ORULES
\$(ELFDIR)/$1.o : \$(SRCDIR)/$1.c \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(ELFDIR); \$(ELFCOMPILE_CMD) \$("$1"_C_FLAGS) \$(SRCDIR)/$1.c
"
])


dnl#}}}

AC_DEFUN(JD_CREATE_EXEC_RULE, dnl#{{{
[  
PROGRAM_OBJECT_RULES="$PROGRAM_OBJECT_RULES
$1 : \$(OBJDIR)/$1
	@echo $1 created in \$(OBJDIR)
\$(OBJDIR)/$1 : \$(OBJDIR)/$1.o \$("$1"_DEPS) \$(EXECDEPS)
	\$(CC) -o \$(OBJDIR)/$1 \$(LDFLAGS) \$(OBJDIR)/$1.o \$("$1"_LIBS) \$(EXECLIBS)
\$(OBJDIR)/$1.o : \$(SRCDIR)/$1.c \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(OBJDIR); \$(COMPILE_CMD) \$("$1"_INC) \$(EXECINC) \$(SRCDIR)/$1.c
"
])


dnl#}}}

AC_DEFUN(JD_CREATE_MODULE_ORULES, dnl#{{{
[
 for program_module in $Program_Modules; do
   JD_CREATE_ORULE($program_module)
   JD_CREATE_ELFORULE($program_module)
 done
])

dnl#}}}

AC_DEFUN(JD_GET_MODULES, dnl#{{{
[
 PROGRAM_HFILES=""
 PROGRAM_OFILES=""
 PROGRAM_CFILES=""
 PROGRAM_OBJECTS=""
 PROGRAM_ELFOBJECTS=""
 PROGRAM_OBJECT_RULES=""
 PROGRAM_ELF_ORULES=""
 if test -z "$1"
 then
   Program_Modules=""
 else
   comment_re="^#"
   Program_Modules=`grep -v '$comment_re' $1 | awk '{print [$]1}'`
   Program_H_Modules=`grep -v '$comment_re' $1 | awk '{print [$]2}'`
   for program_module in $Program_H_Modules; do
     PROGRAM_HFILES="$PROGRAM_HFILES $program_module"
   done
 fi
 for program_module in $Program_Modules; do
   PROGRAM_OFILES="$PROGRAM_OFILES $program_module.o"
   PROGRAM_CFILES="$PROGRAM_CFILES $program_module.c"
   PROGRAM_OBJECTS="$PROGRAM_OBJECTS \$(OBJDIR)/$program_module.o"
   PROGRAM_ELFOBJECTS="$PROGRAM_ELFOBJECTS \$(ELFDIR)/$program_module.o"
 done
dnl echo $PROGRAM_OFILES
dnl echo $PROGRAM_HFILES
AC_SUBST(PROGRAM_OFILES)dnl
AC_SUBST(PROGRAM_CFILES)dnl
AC_SUBST(PROGRAM_HFILES)dnl
AC_SUBST(PROGRAM_OBJECTS)dnl
AC_SUBST(PROGRAM_ELFOBJECTS)dnl
])


dnl#}}}

AC_DEFUN(JD_APPEND_RULES, dnl#{{{
[ 
 echo "$PROGRAM_OBJECT_RULES" >> $1
])


dnl#}}}

AC_DEFUN(JD_APPEND_ELFRULES, dnl#{{{
[ 
 echo "$PROGRAM_ELF_ORULES" >> $1
])

dnl#}}}

AC_DEFUN(JD_CREATE_MODULE_EXEC_RULES, dnl#{{{
[
 for program_module in $Program_Modules; do
   JD_CREATE_EXEC_RULE($program_module)
 done
])

dnl#}}}

AC_DEFUN(JD_TERMCAP, dnl#{{{
[
AC_MSG_CHECKING(for Terminfo)
MISC_TERMINFO_DIRS="$FINKPREFIX/share/terminfo"
if test ! -d $MISC_TERMINFO_DIRS
then
   MISC_TERMINFO_DIRS=""
fi

JD_Terminfo_Dirs="/usr/lib/terminfo \
                 /usr/share/terminfo \
                 /usr/share/lib/terminfo \
		 /usr/local/lib/terminfo \
		 $MISC_TERMINFO_DIRS"

TERMCAP=-ltermcap

for terminfo_dir in $JD_Terminfo_Dirs
do
   if test -d $terminfo_dir 
   then
      AC_MSG_RESULT(yes)
      TERMCAP=""
      break
   fi
done
if test "$TERMCAP"; then
  AC_MSG_RESULT(no)
  AC_DEFINE(USE_TERMCAP)
fi
AC_SUBST(TERMCAP)dnl
AC_SUBST(MISC_TERMINFO_DIRS)dnl
])


dnl#}}}

AC_DEFUN(JD_ANSI_CC, dnl#{{{
[
AC_PROG_CC
AC_PROG_CPP
AC_PROG_GCC_TRADITIONAL
AC_ISC_POSIX
AC_AIX

dnl #This stuff came from Yorick config script
dnl
dnl # HPUX needs special stuff
dnl
AC_EGREP_CPP(yes,
[#ifdef hpux
  yes
#endif
], [
AC_DEFINE(_HPUX_SOURCE)
if test "$CC" = cc; then CC="cc -Ae"; fi
])dnl
dnl
dnl #Be sure we've found compiler that understands prototypes
dnl
AC_MSG_CHECKING(C compiler that understands ANSI prototypes)
AC_TRY_COMPILE([ ],[
 extern int silly (int);], [
 AC_MSG_RESULT($CC looks ok.  Good.)], [
 AC_MSG_RESULT($CC is not a good enough compiler)
 AC_MSG_ERROR(Set env variable CC to your ANSI compiler and rerun configure.)
 ])dnl
])dnl

dnl#}}}

AC_DEFUN(JD_ELF_COMPILER, dnl#{{{
[
dnl #-------------------------------------------------------------------------
dnl # Check for dynamic linker
dnl #-------------------------------------------------------------------------
DYNAMIC_LINK_LIB=""
AC_CHECK_HEADER(dlfcn.h,[
  AC_DEFINE(HAVE_DLFCN_H)
  AC_CHECK_LIB(dl,dlopen,[
    DYNAMIC_LINK_LIB="-ldl"
    AC_DEFINE(HAVE_DLOPEN)
   ],[
    AC_CHECK_FUNC(dlopen,AC_DEFINE(HAVE_DLOPEN))
    if test "$ac_cv_func_dlopen" != yes
    then
      AC_MSG_WARN(cannot perform dynamic linking)
    fi
   ])])
AC_SUBST(DYNAMIC_LINK_LIB)

ELFLIB="lib\$(THIS_LIB).so"
ELFLIB_MAJOR="\$(ELFLIB).\$(ELF_MAJOR_VERSION)"
ELFLIB_MAJOR_MINOR="\$(ELFLIB).\$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION)"

case "$host_os" in
  *linux* )
    DYNAMIC_LINK_FLAGS="-Wl,-export-dynamic"
    ELF_CC="gcc"
    ELF_CFLAGS="-O2 -fno-strength-reduce -fPIC"
    ELF_LINK="gcc -shared -Wl,-soname#"
    ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED="gcc \$(CFLAGS) -shared -fPIC"
    ;;
  *solaris* )
    if test "$GCC" = yes
    then
      DYNAMIC_LINK_FLAGS=""
      ELF_CC="gcc"
      ELF_CFLAGS="-O2 -fno-strength-reduce -fPIC"
      ELF_LINK="gcc -shared -Wl,-ztext -Wl,-h#"
      ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
      ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
      CC_SHARED="gcc \$(CFLAGS) -G -fPIC"
    else
      DYNAMIC_LINK_FLAGS=""
      ELF_CC="cc"
      ELF_CFLAGS="-K pic"
      ELF_LINK="cc -G -h#"
      ELF_LINK_CMD="\$(ELF_LINK)\$(ELFLIB_MAJOR)"
      ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
      CC_SHARED="cc \$(CFLAGS) -G -K pic"
    fi
    ;;
   # osr5 or unixware7 with current or late autoconf
  *sco3.2v5* | *unixware-5* | *sco-sysv5uw7*)
     if test "$GCC" = yes
     then
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="gcc"
       ELF_CFLAGS="-O2 -fno-strength-reduce -fPIC"
       ELF_LINK="gcc -shared -Wl,-h#"
       ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="gcc \$(CFLAGS) -G -fPIC"
     else
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="cc"
       ELF_CFLAGS="-K pic"
       # ELF_LINK="ld -G -z text -h#"
       ELF_LINK="cc -G -z text -h#"
       ELF_LINK_CMD="\$(ELF_LINK)\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="cc \$(CFLAGS) -G -K pic"
     fi
     ;;
  *irix6.5* )
     echo "Note: ELF compiler for host_os=$host_os may not be correct"
     echo "double-check: 'mode_t', 'pid_t' may be wrong!"
     if test "$GCC" = yes
     then
       # not tested
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="gcc"
       ELF_CFLAGS="-O2 -fno-strength-reduce -fPIC"
       ELF_LINK="gcc -shared -Wl,-h#"
       ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="gcc \$(CFLAGS) -shared -fPIC"
     else
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="cc"
       ELF_CFLAGS="-K pic"     # default anyhow
       ELF_LINK="cc -shared -o #"
       ELF_LINK_CMD="\$(ELF_LINK)\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="cc \$(CFLAGS) -shared -K pic"
     fi
     ;;
  *darwin* )
     DYNAMIC_LINK_FLAGS=""
     ELF_CC="cc"
     ELF_CFLAGS="$CFLAGS -O2 -fno-strength-reduce -fno-common"
     ELF_LINK="cc -dynamiclib"
     ELF_LINK_CMD="\$(ELF_LINK) -install_name \$(install_lib_dir)/\$(ELFLIB_MAJOR) -compatibility_version \$(ELF_MAJOR_VERSION) -current_version \$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION)"
     ELF_DEP_LIBS="$LDFLAGS \$(DL_LIB)"
     CC_SHARED="cc -bundle -flat_namespace -undefined suppress \$(CFLAGS) -fno-common"
     ELFLIB="lib\$(THIS_LIB).dylib"
     ELFLIB_MAJOR="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).dylib"
     ELFLIB_MAJOR_MINOR="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION).dylib"
     ;;
  * )
    echo "Note: ELF compiler for host_os=$host_os may be wrong"
    ELF_CC="$CC"
    ELF_CFLAGS="$CFLAGS -fPIC"
    ELF_LINK="$CC -shared"
    ELF_LINK_CMD="\$(ELF_LINK)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED="$CC $CFLAGS -shared -fPIC"
esac

AC_SUBST(ELF_CC)
AC_SUBST(ELF_CFLAGS)
AC_SUBST(ELF_LINK)
AC_SUBST(ELF_LINK_CMD)
AC_SUBST(ELF_DEP_LIBS)
AC_SUBST(DYNAMIC_LINK_FLAGS)
AC_SUBST(CC_SHARED)
AC_SUBST(ELFLIB)
AC_SUBST(ELFLIB_MAJOR)
AC_SUBST(ELFLIB_MAJOR_MINOR)
])


dnl#}}}

AC_DEFUN(JD_F77_COMPILER, dnl#{{{
[
case "$host_os" in
 *linux* )
   F77="g77"
   F77_LIBS="-lg2c"
 ;;
 *solaris*)
   F77=f77
   #F77_LIBS="-lF77 -lM77 -L/opt/SUNWspro/SC4.0/lib -lsunmath"
   F77_LIBS="-lF77 -lM77 -lsunmath"
   ;;
 *)
   echo ""
   echo "WARNING: Assuming f77 as your FORTRAN compiler"
   echo ""
   F77=f77
   F77_LIBS=""
esac
AC_SUBST(F77)
AC_SUBST(F77_LIBS)
])



dnl#}}}

AC_DEFUN(JD_WITH_LIBRARY, dnl#{{{
[
JD_UPPERCASE($1,JD_ARG1)
AC_ARG_WITH($1,
  [  --with-$1=DIR      Use DIR/lib and DIR/include for $1],
  [jd_with_$1_arg=$withval], [jd_with_$1_arg=no])
case "x$jd_with_$1_arg" in
xno)
  ;;
x)
  AC_MSG_ERROR(--with-$1 requres a value)
  ;;
*)
   JD_ARG1[]_INC=-I$jd_with_$1_arg/include
   JD_ARG1[]_LIB=-L$jd_with_$1_arg/lib
   JD_SET_RPATH($jd_with_$1_arg/lib)
  ;;
esac

AC_ARG_WITH($1lib,
  [  --with-$1lib=DIR   $1 library in DIR],
  [jd_with_$1lib_arg=$withval], [jd_with_$1lib_arg=no])
case "x$jd_with_$1lib_arg" in
xno)
  ;;
x)
  AC_MSG_ERROR(--with-$1lib requres a value)
  ;;
*)
   JD_ARG1[]_INC=-I$jd_with_$1lib_arg
   JD_ARG1[]_LIB=-L$jd_with_$1lib_arg
   JD_SET_RPATH($jd_with_$1lib_arg)
  ;;
esac

AC_ARG_WITH($1inc, 
  [  --with-$1inc=DIR   $1 include files in DIR],
  [jd_with_$1inc_arg=$withval], [jd_with_$1inc_arg=no])
case "x$jd_with_$1inc_arg" in
x)
  AC_MSG_ERROR(--with-$1inc requres a value)
  ;;
xno)
  ;;
*)
   JD_ARG1[]_INC=-I$jd_with_$1inc_arg
  ;;
esac
AC_SUBST(JD_ARG1[]_INC)
AC_SUBST(JD_ARG1[]_LIB)
])
dnl#}}}

AC_DEFUN(JD_SLANG_MODULE_INSTALL_DIR, dnl#{{{
[
MODULE_INSTALL_DIR=$libdir/slang/modules
SL_FILES_INSTALL_DIR=$datadir/slsh/local-packages
AC_SUBST(MODULE_INSTALL_DIR)
AC_SUBST(SL_FILES_INSTALL_DIR)
])
#}}
