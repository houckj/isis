dnl# -*- mode: sh; mode: fold -*-
dnl Here are some global variables that need initialized.

AC_DEFUN(JH_CONFIG_DIR_SANS_AMD, dnl#{{{
[
CONFIG_DIR_SANS_AMD=$PWD
AC_SUBST(CONFIG_DIR_SANS_AMD)dnl
])

dnl#}}}

AC_DEFUN(JD_INIT, dnl#{{{
[
#These variable are initialized by JD init function
CONFIG_DIR=`pwd`
cd $srcdir
if test "`pwd`" != "$CONFIG_DIR"
  then
  AC_MSG_ERROR("This software does not support configuring from another directory.   See the INSTALL file")
  fi
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

dnl# This function expand the "prefix variables.  For example, it will expand
dnl# values such as ${exec_prefix}/foo when ${exec_prefix} itself has a
dnl# of ${prefix}.  This function produces the shell variables:
dnl# jd_prefix_libdir, jd_prefix_incdir
AC_DEFUN(JD_EXPAND_PREFIX, dnl#{{{
[
  if test "X$jd_prefix" = "X"
  then
    jd_prefix=$ac_default_prefix
    if test "X$prefix" != "XNONE"
    then
      jd_prefix="$prefix"
    fi
    jd_exec_prefix="$jd_prefix"
    if test "X$exec_prefix" != "XNONE"
    then
      jd_exec_prefix="$exec_prefix"
    fi

    dnl#Unfortunately, exec_prefix may have a value like ${prefix}, etc.
    dnl#Let the shell expand those.  Yuk.
    eval `sh <<EOF
      prefix=$jd_prefix
      exec_prefix=$jd_exec_prefix
      libdir=$libdir
      includedir=$includedir
      echo jd_prefix_libdir="\$libdir" jd_prefix_incdir="\$includedir"
EOF
`
  fi
])
#}}}

AC_DEFUN(JD_SET_OBJ_SRC_DIR, dnl#{{{
[
#---------------------------------------------------------------------------
# Set the source directory and object directory.   The makefile assumes an
# abcolute path name.  This is because src/Makefile cds to OBJDIR and compiles
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

OBJDIR=$SRCDIR/objs
ELFDIR=$SRCDIR/elfobjs
AC_SUBST(SRCDIR)dnl
AC_SUBST(OBJDIR)dnl
AC_SUBST(ELFDIR)dnl
])

dnl#}}}

RPATH=""
AC_SUBST(RPATH)dnl

AC_DEFUN(JD_INIT_RPATH, dnl#{{{
[
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
JD_UP_NAME[]_LIB_DIR=$JD_Above_Dir/$1/src/objs
JD_UP_NAME[]_INCLUDE=$JD_Above_Dir/$1/src

if test -z "[$]JD_UP_NAME[]_INCLUDE"
then
    echo Unable to find the JD_UP_NAME directory
    echo You may have to edit $CONFIG_DIR/src/Makefile.
fi

AC_SUBST(JD_UP_NAME[]_LIB_DIR)dnl
AC_SUBST(JD_UP_NAME[]_INCLUDE)dnl
undefine([JD_UP_NAME])dnl
])

dnl#}}}

AC_DEFUN(JD_FIND_GENERIC, dnl#{{{
[
  AC_REQUIRE([JD_EXPAND_PREFIX])dnl

  changequote(<<, >>)dnl
  define(<<JD_UP_NAME>>, translit($1, [a-z], [A-Z]))dnl
  changequote([, ])dnl
# Look for the JD_UP_NAME package
#JD_UP_NAME[]_INCLUDE=""
#JD_UP_NAME[]_LIB_DIR=""

# This list consists of "include,lib include,lib ..."
JD_Search_Dirs="$JD_Above_Dir2/$1/libsrc,$JD_Above_Dir2/$1/libsrc/objs \
                $JD_Above_Dir/$1/libsrc,$JD_Above_Dir/$1/libsrc/objs \
		$JD_Above_Dir2/$1/src,$JD_Above_Dir2/$1/src/objs \
                $JD_Above_Dir/$1/src,$JD_Above_Dir/$1/src/objs"

JD_Search_Dirs="$JD_Search_Dirs \
                $jd_prefix_incdir,$jd_prefix_libdir \
		$HOME/include,$HOME/lib"

dnl if test -n "$ARCH"
dnl then
dnl  JD_Search_Dirs="$JD_Search_Dirs $HOME/include,$HOME/$ARCH/lib"
dnl  JD_Search_Dirs="$JD_Search_Dirs $HOME/include,$HOME/sys/$ARCH/lib"
dnl fi

# Now add the standard system includes.  The reason for doing this is that
# the other directories may have a better chance of containing a more recent
# version.

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
    JD_UP_NAME[]_LIB_DIR=$JD_Above_Dir/$1/src/objs
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

AC_DEFUN(JH_PURIFY, dnl#{{{
[
AC_ARG_ENABLE(purify,
	      [  --enable-purify         compile with Purify],
	      [use_purify=$enableval])
if test -n "$use_purify"
then
  CC="purify -windows=no -log-file=purify.log $CC"
  FC="purify -windows=no -log-file=purify.log $FC"
fi
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
    CFLAGS="$CFLAGS \
            -ansi -pedantic -D__EXTENSIONS__ -D_XOPEN_SOURCE -D_POSIX_C_SOURCE \
            -W -Wall \
            -Wpointer-arith -Wcast-align -Wcast-qual \
	    -Wstrict-prototypes \
	    -Wshadow \
	    -Waggregate-return \
	    -Wmissing-prototypes \
	    -Wnested-externs"
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

AC_DEFUN(JD_APPEND_RULES, dnl#{{{
[
 echo "$PROGRAM_OBJECT_RULES" >> $1
 PROGRAM_OBJECT_RULES=""
])
dnl#}}}

AC_DEFUN(JH_APPEND_FC_RULES, dnl#{{{
[
 echo "$PROGRAM_FC_OBJECT_RULES" >> $1
 PROGRAM_FC_OBJECT_RULES=""
])
dnl#}}}

AC_DEFUN(JD_APPEND_ELFRULES, dnl#{{{
[
 echo "$PROGRAM_ELF_ORULES" >> $1
 PROGRAM_ELF_ORULES=""
])

dnl#}}}

AC_DEFUN(JH_APPEND_FC_ELFRULES, dnl#{{{
[
 echo "$PROGRAM_FC_ELF_ORULES" >> $1
 PROGRAM_FC_ELF_ORULES=""
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

AC_DEFUN(JH_CREATE_FC_ORULE, dnl#{{{
[
PROGRAM_FC_OBJECT_RULES="$PROGRAM_FC_OBJECT_RULES
\$(OBJDIR)/$1.o : \$(SRCDIR)/$1.f \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(OBJDIR); \$(FC_COMPILE_CMD) \$("$1"_FC_FLAGS) \$(SRCDIR)/$1.f
"
])

dnl#}}}

AC_DEFUN(JD_CREATE_ELFORULE, dnl#{{{
[
PROGRAM_ELF_ORULES="$PROGRAM_ELF_ORULES
\$(ELFDIR)/$1.o : \$(SRCDIR)/$1.c \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(ELFDIR); \$(ELFCOMPILE_CMD) \$("$1"_C_FLAGS) \$("$1"_ELFC_FLAGS) \$(SRCDIR)/$1.c
"
])

dnl#}}}

AC_DEFUN(JH_CREATE_FC_ELFORULE, dnl#{{{
[
PROGRAM_FC_ELF_ORULES="$PROGRAM_FC_ELF_ORULES
\$(ELFDIR)/$1.o : \$(SRCDIR)/$1.f \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(ELFDIR); \$(FC_ELFCOMPILE_CMD) \$("$1"_FC_FLAGS) \$("$1"_ELF_FC_FLAGS) \$(SRCDIR)/$1.f
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
	cd \$(OBJDIR); \$(COMPILE_CMD) \$("$1"_INC) \$("$1"_C_FLAGS) \$(EXECINC) \$(SRCDIR)/$1.c
"
])

dnl#}}}

AC_DEFUN(JD_CREATE_FLINK_EXEC_RULE, dnl#{{{
[
PROGRAM_OBJECT_RULES="$PROGRAM_OBJECT_RULES
$1 : \$(OBJDIR)/$1
	@echo $1 created in \$(OBJDIR)
\$(OBJDIR)/$1 : \$(OBJDIR)/$1.o \$("$1"_DEPS) \$(EXECDEPS)
	\$(FC) \$(FCFLAGS) -o \$(OBJDIR)/$1 \$(LDFLAGS) \$(OBJDIR)/$1.o \$("$1"_LIBS) \$(EXECLIBS)
\$(OBJDIR)/$1.o : \$(SRCDIR)/$1.c \$(DOT_O_DEPS) \$("$1"_O_DEP)
	cd \$(OBJDIR); \$(COMPILE_CMD) \$("$1"_INC) \$("$1"_C_FLAGS) \$(EXECINC) \$(SRCDIR)/$1.c
"
])

dnl#}}}

AC_DEFUN(JD_CREATE_MODULE_ORULES, dnl#{{{
[
 Program_Modules="Program_Modules_$1"
 cmd="eval echo \$$Program_Modules"
 list=`$cmd`
 for program_module in $list; do
   JD_CREATE_ORULE($program_module)
   JD_CREATE_ELFORULE($program_module)
 done
])

dnl#}}}

AC_DEFUN(JH_CREATE_FC_MODULE_ORULES, dnl#{{{
[
 Program_FC_Modules="Program_FC_Modules_$1"
 cmd="eval echo \$$Program_FC_Modules"
 list=`$cmd`
 for program_module in $list; do
   JH_CREATE_FC_ORULE($program_module)
   JH_CREATE_FC_ELFORULE($program_module)
 done
])

dnl#}}}

AC_DEFUN(JD_CREATE_MODULE_EXEC_RULES, dnl#{{{
[
 Program_Modules="Program_Modules_$1"
 cmd="eval echo \$$Program_Modules"
 list=`$cmd`
 for program_module in $list; do
   JD_CREATE_EXEC_RULE($program_module)
 done
])

dnl#}}}

dnl# This macro process the --with-xxx, --with-xxxinc, and --with-xxxlib
dnl# command line arguments and returns the values as shell variables
dnl# jd_xxx_include_dir and jd_xxx_library_dir.  It does not perform any
dnl# substitutions, nor check for the existence of the supplied values.
AC_DEFUN(JD_WITH_LIBRARY_PATHS, dnl#{{{
[
 JD_UPPERCASE($1,JD_ARG1)
 jd_$1_include_dir=""
 jd_$1_library_dir=""
 jd_with_$1_library=""

 AC_ARG_WITH($1,
  [  --with-$1=DIR      Use DIR/lib and DIR/include for $1],
  [jd_with_$1_arg=$withval], [jd_with_$1_arg=unspecified])

 case "x$jd_with_$1_arg" in
   xno)
     jd_with_$1_library="no"
    ;;
   x)
    AC_MSG_ERROR(--with-$1 requires a value-- try yes or no)
    ;;
   xunspecified)
    ;;
   xyes)
    ;;
   *)
    jd_$1_include_dir="$jd_with_$1_arg"/include
    jd_$1_library_dir="$jd_with_$1_arg"/lib
    ;;
 esac

 AC_ARG_WITH($1lib,
  [  --with-$1lib=DIR   $1 library in DIR],
  [jd_with_$1lib_arg=$withval], [jd_with_$1lib_arg=unspecified])
 case "x$jd_with_$1lib_arg" in
   xunspecified)
    ;;
   xno)
    ;;
   x)
    AC_MSG_ERROR(--with-$1lib requres a value)
    ;;
   *)
    jd_$1_library_dir="$jd_with_$1lib_arg"
    ;;
 esac

 AC_ARG_WITH($1inc,
  [  --with-$1inc=DIR   $1 include files in DIR],
  [jd_with_$1inc_arg=$withval], [jd_with_$1inc_arg=unspecified])
 case "x$jd_with_$1inc_arg" in
   x)
     AC_MSG_ERROR(--with-$1inc requres a value)
     ;;
   xunspecified)
     ;;
   xno)
     ;;
   *)
    jd_$1_include_dir="$jd_with_$1inc_arg"
   ;;
 esac
])
dnl#}}}

dnl# This function checks for the existence of the specified library $1 with
dnl# header file $2.  If the library exists, then the shell variables will
dnl# be created:
dnl#  jd_with_$1_library=yes/no,
dnl#  jd_$1_inc_file
dnl#  jd_$1_include_dir
dnl#  jd_$1_library_dir
AC_DEFUN(JD_CHECK_FOR_LIBRARY, dnl#{{{
[
  AC_REQUIRE([JD_EXPAND_PREFIX])dnl
  dnl JD_UPPERCASE($1,JD_ARG1)
  JD_WITH_LIBRARY_PATHS($1)
  AC_MSG_CHECKING(for the $1 library and header files $2)
  if test X"$jd_with_$1_library" != Xno
  then
    jd_$1_inc_file=$2
    dnl# jd_with_$1_library="yes"

    if test "X$jd_$1_inc_file" = "X"
    then
       jd_$1_inc_file=$1.h
    fi
    if test X"$jd_$1_include_dir" = X
    then
       lib_include_dirs="\
            $jd_prefix_incdir \
            $jh_depdir/include \
            $jh_depdir/pgplot \
            /usr/local/$1/include \
            /usr/local/include/$1 \
  	  /usr/local/include \
  	  /usr/include/$1 \
  	  /usr/$1/include \
  	  /usr/include \
  	  /opt/include/$1 \
  	  /opt/$1/include \
  	  /opt/include"

       for X in $lib_include_dirs
       do
          if test -r "$X/$jd_$1_inc_file"
	  then
  	  jd_$1_include_dir="$X"
            break
          fi
       done
       if test X"$jd_$1_include_dir" = X
       then
         jd_with_$1_library="no"
       fi
    fi

    if test X"$jd_$1_library_dir" = X
    then
       lib_library_dirs="\
            $jd_prefix_libdir \
            $jh_depdir/lib \
            $jh_depdir/pgplot \
            /usr/local/lib \
            /usr/local/lib/$1 \
            /usr/local/$1/lib \
  	  /usr/lib \
  	  /usr/lib/$1 \
  	  /usr/$1/lib \
  	  /opt/lib \
  	  /opt/lib/$1 \
  	  /opt/$1/lib"

       case "$host_os" in
         *darwin* )
	   exts="dylib so a"
	   ;;
	 *cygwin* )
	   exts="dll.a so a"
	   ;;
	 * )
	   exts="so a"
       esac

       found=0
       for X in $lib_library_dirs
       do
         for E in $exts
	 do
           if test -r "$X/lib$1.$E"
	   then
  	     jd_$1_library_dir="$X"
	     found=1
	     break
           fi
         done
	 if test $found -eq 1
	 then
	   break
	 fi
       done
       if test X"$jd_$1_library_dir" = X
       then
         jd_with_$1_library="no"
       fi
    fi
  fi

  if test X"$jd_$1_include_dir" != X -a "$jd_$1_library_dir" != X
  then
    AC_MSG_RESULT(yes: $jd_$1_library_dir and $jd_$1_include_dir)
    jd_with_$1_library="yes"
    dnl#  Avoid using /usr/lib and /usr/include because of problems with
    dnl#  gcc on some solaris systems.
    JD_ARG1[]_LIB=-L$jd_$1_library_dir
    if test "X$jd_$1_library_dir" = "X/usr/lib"
    then
      JD_ARG1[]_LIB=""
    else
      JD_SET_RPATH($jd_$1_library_dir)
    fi

    JD_ARG1[]_INC=-I$jd_$1_include_dir
    if test "X$jd_$1_include_dir" = "X/usr/include"
    then
      JD_ARG1[]_INC=""
    fi
  else
    AC_MSG_RESULT(no)
    jd_with_$1_library="no"
    JD_ARG1[]_INC=""
    JD_ARG1[]_LIB=""
  fi
  AC_SUBST(JD_ARG1[]_LIB)
  AC_SUBST(JD_ARG1[]_INC)
])
dnl#}}}

AC_DEFUN(JD_WITH_LIBRARY, dnl#{{{
[
  JD_CHECK_FOR_LIBRARY($1, $2)
  if test "$jd_with_$1_library" = "no"
  then
    AC_MSG_ERROR(unable to find the $1 library and header file $jd_$1_inc_file)
  fi
])
dnl#}}}

AC_DEFUN(JD_SLANG_VERSION, dnl#{{{
[
 slang_h=$jd_slang_include_dir/slang.h
 AC_MSG_CHECKING(SLANG_VERSION in $slang_h)
slang_version=`grep "^#define  *SLANG_VERSION " $slang_h |
               awk '{ print [$]3 }'`
slang_major_version=`echo $slang_version |
 awk '{ print int([$]1/10000) }'`
slang_minor_version=`echo $slang_version $slang_major_version |
 awk '{ print int(([$]1 - [$]2*10000)/100) }'`
slang_patchlevel_version=`echo $slang_version $slang_major_version $slang_minor_version |
 awk '{ print ([$]1 - [$]2*10000 - [$]3*100) }'`

AC_MSG_RESULT($slang_major_version.$slang_minor_version.$slang_patchlevel_version)
AC_SUBST(slang_version)
AC_SUBST(slang_major_version)
AC_SUBST(slang_minor_version)
AC_SUBST(slang_patchlevel_version)
])
#}}}

AC_DEFUN(JD_SLANG_MODULE_INSTALL_DIR, dnl#{{{
[
  AC_REQUIRE([JD_SLANG_VERSION])
  if test "X$slang_major_version" = "X1"
  then
    MODULE_INSTALL_DIR="$libdir/slang/modules"
  else
    MODULE_INSTALL_DIR="$libdir/slang/v$slang_major_version/modules"
  fi
  SL_FILES_INSTALL_DIR=$datadir/slsh/local-packages
  AC_SUBST(MODULE_INSTALL_DIR)
  AC_SUBST(SL_FILES_INSTALL_DIR)
])
#}}}

AC_DEFUN(JD_GET_MODULES, dnl#{{{
[
 if test -z "$1"
 then
   Program_Modules=""
 else
   Program_Modules=`cat $1`
 fi
 PROGRAM_OFILES=""
 PROGRAM_OBJECTS=""
 PROGRAM_ELFOBJECTS=""
 for program_module in $Program_Modules; do
   PROGRAM_OFILES="$PROGRAM_OFILES $program_module.o"
   PROGRAM_OBJECTS="$PROGRAM_OBJECTS \$(OBJDIR)/$program_module.o"
   PROGRAM_ELFOBJECTS="$PROGRAM_ELFOBJECTS \$(ELFDIR)/$program_module.o"
 done
Program_Modules_$2=$Program_Modules
PROGRAM_OFILES_$2=$PROGRAM_OFILES
PROGRAM_OBJECTS_$2=$PROGRAM_OBJECTS
PROGRAM_ELFOBJECTS_$2=$PROGRAM_ELFOBJECTS
AC_SUBST(PROGRAM_OFILES_$2)dnl
AC_SUBST(PROGRAM_OBJECTS_$2)dnl
AC_SUBST(PROGRAM_ELFOBJECTS_$2)dnl
])

dnl#}}}

AC_DEFUN(JH_GET_FC_MODULES, dnl#{{{
[
 if test -z "$1"
 then
   Program_FC_Modules=""
 else
   Program_FC_Modules=`cat $1`
 fi
 PROGRAM_FC_OFILES=""
 PROGRAM_FC_OBJECTS=""
 PROGRAM_FC_ELFOBJECTS=""
 for program_module in $Program_FC_Modules; do
   PROGRAM_FC_OFILES="$PROGRAM_FC_OFILES $program_module.o"
   PROGRAM_FC_OBJECTS="$PROGRAM_FC_OBJECTS \$(OBJDIR)/$program_module.o"
   PROGRAM_FC_ELFOBJECTS="$PROGRAM_FC_ELFOBJECTS \$(ELFDIR)/$program_module.o"
 done
Program_FC_Modules_$2=$Program_FC_Modules
PROGRAM_FC_OFILES_$2=$PROGRAM_FC_OFILES
PROGRAM_FC_OBJECTS_$2=$PROGRAM_FC_OBJECTS
PROGRAM_FC_ELFOBJECTS_$2=$PROGRAM_FC_ELFOBJECTS
AC_SUBST(PROGRAM_FC_OFILES_$2)dnl
AC_SUBST(PROGRAM_FC_OBJECTS_$2)dnl
AC_SUBST(PROGRAM_FC_ELFOBJECTS_$2)dnl
])

dnl#}}}

AC_DEFUN(JH_GET_SHAREFILES, dnl#{{{
[
 if test -z "$1"
 then
   Program_Sharefiles=""
 else
   Program_Sharefiles=`cat $1`
 fi
 PROGRAM_SHAREFILES=""
 for program_share in $Program_Sharefiles; do
   PROGRAM_SHAREFILES="$PROGRAM_SHAREFILES $program_share"
 done
AC_SUBST(PROGRAM_SHAREFILES)dnl
])

dnl#}}}

AC_DEFUN(JH_GET_ETCFILES, dnl#{{{
[
 if test -z "$1"
 then
   Program_Etcfiles=""
 else
   Program_Etcfiles=`cat $1`
 fi
 PROGRAM_ETCFILES=""
 for program_etc in $Program_Etcfiles; do
   PROGRAM_ETCFILES="$PROGRAM_ETCFILES $program_etc"
 done
AC_SUBST(PROGRAM_ETCFILES)dnl
])

dnl#}}}

AC_DEFUN(JH_LIST_ELFOBJECTS, dnl#{{{
[
 Other_Modules=`cat $1/modules.lis`
 other_elfobjects=""
 for other_module in $Other_Modules; do
   other_elfobjects="$other_elfobjects \$(ELFDIR_$2)/$other_module.o"
 done
ELFOBJECTS_$2="$other_elfobjects"
ELFDIR_$2="\$(config_dir)/$1/elfobjs"
AC_SUBST(ELFOBJECTS_$2)dnl
AC_SUBST(ELFDIR_$2)dnl
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
      case "$host_os" in
        *darwin* )
	  dnl # assuming dlcompat is from Fink
	   AC_CHECK_FILE(/sw/lib/libdl.a,[
             DYNAMIC_LINK_LIB="-L/sw/lib -ldl"
	     AC_DEFINE(HAVE_DLOPEN)
	   ],[
             AC_MSG_WARN(cannot perform dynamic linking)
	   ])
 	;;
        * )
          AC_MSG_WARN(cannot perform dynamic linking)
	;;
      esac
    fi
   ])])
AC_SUBST(DYNAMIC_LINK_LIB)

ELFLIB="lib\$(THIS_LIB).so"
ELFLIB_MAJOR="\$(ELFLIB).\$(ELF_MAJOR_VERSION)"
ELFLIB_MAJOR_MINOR="\$(ELFLIB).\$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION)"

SLANG_LIB_FOR_MODULES=#
CFITSIO_LIB_FOR_MODULES=#
ISIS_LIB_FOR_MODULES=#
INSTALL_MODULE="\$(INSTALL_DATA)"
SLANG_DLL_CFLAGS=""

if test "$GCC" = yes
then
  if test X"$CFLAGS" = X
  then
     CFLAGS="-O2"
  fi
fi

case "$host_os" in
  *linux* )
    DYNAMIC_LINK_FLAGS="-Wl,-export-dynamic"
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-soname#"
    ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED="\$(CC) \$(CFLAGS) -shared -fPIC"
    ELF_FC_FCFLAGS="-fPIC"
    ;;
  *solaris* )
    if test "$GCC" = yes
    then
      DYNAMIC_LINK_FLAGS=""
      ELF_CC="\$(CC)"
      ELF_CFLAGS="\$(CFLAGS) -fPIC"
      ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-ztext -Wl,-h#"
      ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
      ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
      CC_SHARED="\$(CC) \$(CFLAGS) -G -fPIC"
      ELF_FC_FCFLAGS="-fPIC"
    else
      DYNAMIC_LINK_FLAGS=""
      ELF_CC="\$(CC)"
      ELF_CFLAGS="\$(CFLAGS) -KPIC"
      ELF_LINK="\$(CC) \$(LDFLAGS) -G -h#"
      ELF_LINK_CMD="\$(ELF_LINK)\$(ELFLIB_MAJOR)"
      ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
      CC_SHARED="\$(CC) \$(CFLAGS) -G -KPIC"
      ELF_FC_FCFLAGS="-KPIC"
    fi
    ;;
   # osr5 or unixware7 with current or late autoconf
  *sco3.2v5* | *unixware-5* | *sco-sysv5uw7*)
     if test "$GCC" = yes
     then
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -fPIC"
       ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-h#"
       ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="\$(CC) \$(CFLAGS) -G -fPIC"
       ELF_FC_FCFLAGS="-fPIC"
     else
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -Kpic"
       # ELF_LINK="ld -G -z text -h#"
       ELF_LINK="\$(CC) \$(LDFLAGS) -G -z text -h#"
       ELF_LINK_CMD="\$(ELF_LINK)\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="\$(CC) \$(CFLAGS) -G -Kpic"
       ELF_FC_FCFLAGS="-Kpic"
     fi
     ;;
  *irix6.5* )
     echo "Note: ELF compiler for host_os=$host_os may not be correct"
     echo "double-check: 'mode_t', 'pid_t' may be wrong!"
     if test "$GCC" = yes
     then
       # not tested
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -fPIC"
       ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-h#"
       ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="\$(CC) \$(CFLAGS) -shared -fPIC"
       ELF_FC_FCFLAGS="-fPIC"
     else
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -Kpic"     # default anyhow
       ELF_LINK="\$(CC) \$(LDFLAGS) -shared -o #"
       ELF_LINK_CMD="\$(ELF_LINK)\$(ELFLIB_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED="\$(CC) \$(CFLAGS) -shared -Kpic"
       ELF_FC_FCFLAGS="-Kpic"
     fi
     ;;
  *darwin* )
     DYNAMIC_LINK_FLAGS=""
     if test "X$DYNAMIC_LINK_LIB" != "X"
     then
         dnl # assume dlcompat from Fink
         dnl CFLAGS="$CFLAGS -Ddlsym=dlsym_prepend_underscore"
	 dnl # hint for pgplot from Fink
	 EXTRA_LIB="$EXTRA_LIB #-framework Foundation -framework AppKit -lpng"
     fi
     ELF_CC="\$(CC)"
     ELF_CFLAGS="\$(CFLAGS) -fno-common"
     ELF_LINK="\$(CC) \$(LDFLAGS) -dynamiclib"
     ELF_LINK_CMD="\$(ELF_LINK) -install_name \$(install_lib_dir)/\$(ELFLIB_MAJOR) -compatibility_version \$(ELF_MAJOR_VERSION) -current_version \$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION)"
     ELF_DEP_LIBS="\$(DL_LIB)"
     CC_SHARED="\$(CC) -bundle -flat_namespace -undefined suppress \$(CFLAGS) -fno-common"
     ELFLIB="lib\$(THIS_LIB).dylib"
     ELFLIB_MAJOR="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).dylib"
     ELFLIB_MAJOR_MINOR="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION).dylib"
     ELF_FC_FCFLAGS=""
     ;;
  *freebsd* )
    ELFLIB_MAJOR_MINOR="\$(ELFLIB).\$(ELF_MAJOR_VERSION)"
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    if test "X$PORTOBJFORMAT" = "Xelf" ; then
      ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-soname,\$(ELFLIB_MAJOR)"
    else
      ELF_LINK="ld -Bshareable -x"
    fi
    ELF_LINK_CMD="\$(ELF_LINK)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm"
    CC_SHARED="\$(CC) \$(CFLAGS) -shared -fPIC"
    ELF_FC_FCFLAGS="-fPIC"
    ;;

  *cygwin* )
    DYNAMIC_LINK_FLAGS=""
    ELF_CC="\$(CC)"
    SLANG_DLL_CFLAGS="-DSLANG_DLL=1"
    ELF_CFLAGS="\$(CFLAGS)"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-soname#"
    ELF_LINK_CMD="\$(ELF_LINK),\$(ELFLIB_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED="\$(CC) \$(CFLAGS) -shared"
    ELF_FC_FCFLAGS=""
    ISIS_LIB_FOR_MODULES="-L\$(config_dir)/src/objs -lisis"
    SLANG_LIB_FOR_MODULES="\$(SLANG_LIB)"
    CFITSIO_LIB_FOR_MODULES="\$(CFITSIO_LIB)"
    INSTALL_MODULE="\$(INSTALL)"
    ;;

  * )
    echo "Note: ELF compiler for host_os=$host_os may be wrong"
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared"
    ELF_LINK_CMD="\$(ELF_LINK)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED="\$(CC) \$(CFLAGS) -shared -fPIC"
    ELF_FC_FCFLAGS="-fPIC"
esac

AC_SUBST(ELF_CC)
AC_SUBST(ELF_CFLAGS)
AC_SUBST(ELF_LINK)
AC_SUBST(ELF_LINK_CMD)
AC_SUBST(ELF_DEP_LIBS)
AC_SUBST(ELF_FC_FCFLAGS)
AC_SUBST(DYNAMIC_LINK_FLAGS)
AC_SUBST(CC_SHARED)
AC_SUBST(ISIS_LIB_FOR_MODULES)
AC_SUBST(SLANG_LIB_FOR_MODULES)
AC_SUBST(CFITSIO_LIB_FOR_MODULES)
AC_SUBST(INSTALL_MODULE)
AC_SUBST(SLANG_DLL_CFLAGS)
])

dnl#}}}

AC_DEFUN(JH_CHECK_CFITSIO_LIBNAME, dnl#{{{
[
dnl This nonsense is necessary because we cannot rely on
dnl the name of the cfitsio library.  Common choices are:
dnl
dnl   libcfitsio.so
dnl   libcxccfitsio.so
dnl   libcfitsio${version_string}.so
dnl
dnl And you can forget about about binary compatibility...
dnl
  jh_cfitsio_library_dir=$1
  case "$host_os" in
    *darwin* )
        ext="dylib"
        ;;
    *cygwin* )
        ext="dll"
        ;;
    * )
        ext="so"
        ;;
  esac
  if test -f "${jh_cfitsio_library_dir}/libcfitsio.$ext" ; then
     cfitsio_libname=cfitsio
  elif test -f "${jh_cfitsio_library_dir}/libcxccfitsio.$ext" ; then
     cfitsio_libname=cxccfitsio
  else
     cfitsio_libfile_exists=0
     cfitsio_libname=`/bin/ls ${jh_cfitsio_library_dir}/libcfitsio*.$ext`
     if test x"$cfitsio_libname" != "x" ; then
        cfitsio_libname=`basename $cfitsio_libname | sed s/lib// | sed s/\.$ext//`
        if test -f "${jh_cfitsio_library_dir}/lib${cfitsio_libname}.$ext" ; then
           cfitsio_libfile_exists=1
        fi
     fi
     if test "$cfitsio_libfile_exists" == 0 ; then
        echo '*** WARNING:  Cannot find cfitsio library in' "$jh_cfitsio_library_dir"
     fi
  fi
  CFITSIO_LIB="-L${jh_cfitsio_library_dir} -l${cfitsio_libname}#"
  AC_SUBST(CFITSIO_LIB)
])
dnl#}}}

AC_DEFUN(JH_CHECK_PGPLOT_PNG, dnl#{{{
dnl
dnl Try to determine whether libpgplot needs libpng
dnl
  jh_pgplot_libfile=$1
  OTHER_PGPLOT_LIBS=""
  if test -f "$jh_pgplot_libfile" ; then
    jh_png_symbols=`nm $jh_pgplot_libfile | grep -i png`
    if test x"$jh_png_symbols" != x ; then
       OTHER_PGPLOT_LIBS="-lpng"
    fi
  fi
  AC_SUBST(OTHER_PGPLOT_LIBS)
)
dnl#}}}

AC_DEFUN(JH_CHECK_PGPLOT_LIBNAME, dnl#{{{
[
dnl This nonsense is necessary because we cannot rely on
dnl the name of the pgplot library.  Common choices are:
dnl
dnl   libpgplot.so
dnl   libpgplot${version_string}.so
dnl
dnl (and the same for libcpgplot)
dnl
  jh_pgplot_library_dir=$1
  case "$host_os" in
    *darwin* )
        ext_list="dylib so a"
        ;;
    *cygwin* )
        ext_list="dll a"
        ;;
    * )
        ext_list="so a"
        ;;
  esac
  pgplot_libnames=no
  for ext in ${ext_list}
  do
    lib_pgplot_file="${jh_pgplot_library_dir}/libpgplot.$ext"
    if test -f "$lib_pgplot_file" ; then
       pgplot_libnames="-lcpgplot -lpgplot"
       break
    else
       lib_pgplot_file=`/bin/ls ${jh_pgplot_library_dir}/libpgplot*.$ext | grep -v perl`
       if test x"$lib_pgplot_file" != "x" ; then
          lib_pgplot_name=`basename $lib_pgplot_file | sed s/lib// | sed s/\.$ext//`
          if test -f "${jh_pgplot_library_dir}/lib${lib_pgplot_name}.$ext" ; then
             pgplot_libnames="-lc${lib_pgplot_name} -l${lib_pgplot_name}"
             break
          fi
       fi
    fi
  done
  JH_CHECK_PGPLOT_PNG($lib_pgplot_file)
  if test "X${pgplot_libnames}" = Xno
  then
     echo '*** WARNING:  Cannot find pgplot libraries in' "$jh_pgplot_library_dir"
  fi
  PGPLOT_LIB="-L${jh_pgplot_library_dir} ${pgplot_libnames}"
  AC_SUBST(PGPLOT_LIB)
])
dnl#}}}

AC_DEFUN(JH_WITH_XSPEC_STATIC, dnl#{{{
[
AC_ARG_WITH(xspec-static,
  [  --with-xspec-static[=DIR]    Statically link XSPEC models from DIR=$HEADAS],
  [jh_use_xspec_static=$withval], [jh_use_xspec_static=no])
if test "x$jh_use_xspec_static" = "xno"
then
   LINK_XSPEC_STATIC="no"
else
   LINK_XSPEC_STATIC="yes"
   AC_DEFINE(WITH_XSPEC_STATIC_LINKED)

   HEADAS_DIR=$jh_use_xspec_static
   AC_SUBST(HEADAS_DIR)
   HEADAS=$jh_use_xspec_static

   HEADAS_LIBDIR="${HEADAS}/lib"
   AC_SUBST(HEADAS_LIBDIR)

   WITH_HEADAS="-DWITH_HEADAS"
   AC_SUBST(WITH_HEADAS)

   XSPEC_MODULE_LIBS="-L\$(config_dir)/modules/xspec/src/objs -lxspec-module \$(XS_LIBS)"
   AC_SUBST(XSPEC_MODULE_LIBS)

   CFITSIO_INC="-I$HEADAS/include"
   AC_SUBST(CFITSIO_INC)
   JH_CHECK_CFITSIO_LIBNAME("$HEADAS/lib")

   PGPLOT_INC="-I$HEADAS/include"
   AC_SUBST(PGPLOT_INC)
   JH_CHECK_PGPLOT_LIBNAME("$HEADAS/lib")
   PGPLOT_LIBDIR="$jh_pgplot_library_dir"
   AC_SUBST(PGPLOT_LIBDIR)

   JD_SET_RPATH($HEADAS_LIBDIR)
fi
AC_SUBST(LINK_XSPEC_STATIC)
])

dnl#}}}

AC_DEFUN(JH_WITH_HEADAS, dnl#{{{
[
AC_ARG_WITH(headas,
  [  --with-headas[=DIR]      HEADAS install directory],
  [jh_use_headas=$withval], [jh_use_headas=no])
if test "x$jh_use_headas" != "xno"
then
   HEADAS_DIR=$jh_use_headas
   HEADAS_LIBDIR="${HEADAS_DIR}/lib"
   WITH_HEADAS="-DWITH_HEADAS"
   MODULE_LIST="$MODULE_LIST xspec"
   AC_SUBST(HEADAS_DIR)
   AC_SUBST(HEADAS_LIBDIR)
   AC_SUBST(WITH_HEADAS)
   AC_SUBST(MODULE_LIST)

   CFITSIO_INC="-I${HEADAS_DIR}/include"
   AC_SUBST(CFITSIO_INC)
   JH_CHECK_CFITSIO_LIBNAME($HEADAS_DIR/lib)

   PGPLOT_INC="-I${HEADAS_DIR}/include"
   AC_SUBST(PGPLOT_INC)
   JH_CHECK_PGPLOT_LIBNAME($HEADAS_DIR/lib)
   PGPLOT_LIBDIR="$jh_pgplot_library_dir"
   AC_SUBST(PGPLOT_LIBDIR)

   JD_SET_RPATH($HEADAS_LIBDIR)
fi
])

dnl#}}}

AC_DEFUN(JH_HANDLE_PACKAGE_OPTIONS, dnl#{{{
[
  case "$host_os" in
    *linux* )
    # For 64-bit Linux, force the address of the mem_ symbol to be in low
    # memory (0x80 chosen so that the symbol is aligned on a 128 byte boundary)
    #
      if test "x$ac_cv_sizeof_long" = "x8"; then
         LD_FLAGS_FOR_XSPEC="-Wl,--defsym=mem_=0x80"
      fi
    ;;
    * )
    ;;
  esac
  AC_SUBST(LD_FLAGS_FOR_XSPEC)

  JH_WITH_HEADAS
  JH_WITH_XSPEC_STATIC

  AC_ARG_ENABLE(xspec_tables,
    [  --disable-xspec_tables  Compile without XSPEC table models],
    [ ],
    [AC_DEFINE(HAVE_XSPEC_TABLE_MODELS)])

  AC_ARG_WITH(xspec-version,
    [  --with-xspec-version[=nn]    XSPEC version 11 or 12?],
    [jh_with_xspec_version=$withval], [jh_with_xspec_version=no])
  if test x"$jh_with_xspec_version" = "x11" || test x"$jh_with_xspec_version" = "x12" ; then
     with_xspec_version=$jh_with_xspec_version
  elif test x"$jh_with_xspec_version" = "xno" ; then
     with_xspec_version="NONE"
  else
     echo "### Unsupported xspec version $jh_with_xspec_version"
     exit 1
  fi
  AC_SUBST(with_xspec_version)
])

dnl#}}}

AC_DEFUN(JH_SYS_EXTRA_LIBS, dnl#{{{
[
SYS_EXTRA_LIBS=""
case "$host_os" in
  *darwin* )
     if test "$GCC" = yes
     then
        jh_gcc_major_version=`$CC -dumpversion | sed -e 's/egcs-//' | cut -d'.' -f1`
        if test $jh_gcc_major_version -lt 4
        then
           SYS_EXTRA_LIBS="-lcc_dynamic"
        fi
     fi
     ;;
  * )
     ;;
esac
AC_SUBST(SYS_EXTRA_LIBS)
])

dnl#}}}

AC_DEFUN(JH_INIT_MODULE_LIST, dnl#{{{
[
MODULE_LIST="cfitsio pgplot maplib"
AC_SUBST(MODULE_LIST)
])

dnl#}}}

AC_DEFUN(JD_LARGE_FILE_SUPPORT, dnl#{{{
[
  AC_SYS_LARGEFILE
  AC_FUNC_FSEEKO
  AC_TYPE_OFF_T
  AC_CHECK_SIZEOF(off_t)
])
#}}}

AC_DEFUN(JH_TERMCAP_LIB, dnl#{{{
[
  termcap_lib=""
  for lib in termcap ncurses curses ; do
      AC_CHECK_LIB(${lib},tgetent,[termcap_lib="-l${lib}"; break])
  done
  AC_SUBST(termcap_lib)
])
#}}}

AC_DEFUN(JH_WITH_READLINE, dnl#{{{
[
dnl # Establish readline support: GNU or S-Lang?
  AC_ARG_WITH(readline,
[  --with-readline[[=arg]]  arg=gnu means use GNU readline
                         arg=DIR means use GNU readline
                             (from DIR/lib and DIR/include)
                         arg=slang means use S-Lang readline only
                         ],
     [ with_readline_arg=$withval ], [with_readline_arg=""])
  AC_MSG_CHECKING(type of readline support)
  READLINE_DIR=""
  READLINE_INC=""
  GNU_READLINE=0
  READLINE_LIB="# -lreadline"
  case "x$with_readline_arg" in
    xslang|xno )
      AC_MSG_RESULT(slang);
      ;;
    * )
      termcap_lib=""
      for lib in termcap ncurses curses ; do
          AC_CHECK_LIB(${lib},tgetent,[termcap_lib="-l${lib}"; break])
      done
      case "x$with_readline_arg" in
        x|xyes|xgnu|xGNU )
           AC_CHECK_LIB(readline,main,
              [GNU_READLINE=1],[GNU_READLINE=0],[$termcap_lib])
           if test "$GNU_READLINE" == 1 ; then
              READLINE_LIB="-lreadline ${termcap_lib}"
           else
              READLINE_LIB="# -lreadline ${termcap_lib}"
           fi
           ;;
         * )
           for ext in a so dylib ; do
               if test -r "${with_readline_arg}/lib/libreadline.${ext}" ; then
                   GNU_READLINE=1
                   JD_SET_RPATH($with_readline_arg/lib)
                   READLINE_DIR="-L${with_readline_arg}/lib"
                   READLINE_INC="-I${with_readline_arg}/include"
                   READLINE_LIB="${READLINE_DIR} -lreadline ${termcap_lib}"
                   break
               fi
           done
           ;;
      esac
      if test "$GNU_READLINE" == 1 ; then
         AC_MSG_RESULT([Looks like you have GNU readline.  Can we use it?])
         tmp_libs=$LIBS
         LIBS="$LIBS ${READLINE_LIB}"
         tmp_cflags=$CFLAGS
         CFLAGS="$CFLAGS ${READLINE_INC}"
         AC_LINK_IFELSE(
            AC_LANG_PROGRAM([[#include <stdio.h>
                              #include <readline/readline.h>
                              #include <readline/history.h>
                              ]],
                            [[char *prompt=NULL, *line=NULL;
                              char **matches=NULL;
                              const char *text=NULL;
                              line = readline (prompt);
                              add_history (line);
                              rl_delete_text (0, rl_end);
                              rl_point = rl_end = 0;
                              rl_forced_update_display();
                              rl_refresh_line (0,0);
                              matches = rl_completion_matches (text,NULL);]]),
                             AC_MSG_RESULT([test program linked successfully -- readline looks ok.]),
                             [GNU_READLINE=0])
         CFLAGS=$tmp_cflags
         LIBS=$tmp_libs
      fi
      if test "$GNU_READLINE" == 1 ; then
         AC_DEFINE(HAVE_GNU_READLINE)
         AC_MSG_RESULT(gnu);
      else
         AC_MSG_RESULT([Failed configuring GNU readline.  Will use slang readline instead.]);
      fi
      ;;
  esac
  AC_SUBST(READLINE_LIB)
  AC_SUBST(READLINE_DIR)
  AC_SUBST(READLINE_INC)
])
#}}}

AC_DEFUN(JH_WITH_SUBDIR, dnl#{{{
[
AC_ARG_WITH(subdir,
  [  --with-subdir=DIR       Install in $prefix/isis/$DIR ],
  [ jh_subdir=$withval ], [jh_subdir=no])
if test "x$jh_subdir" != "xno"
then
  subdir="$jh_subdir"
else
  jh_isis_version=`grep " ISIS_VERSION_PREFIX " $CONFIG_DIR/src/isis.h | awk '{print [$]3}'`
  subdir="$jh_isis_version"
fi
if test "x$prefix" = "xNONE" ; then
   prefix_input="$ac_default_prefix"
   prefix="$ac_default_prefix/isis/$subdir"
else
   prefix_input="$prefix"
   prefix="$prefix/isis/$subdir"
fi
AC_ARG_WITH(dep,
  [  --with-dep=DIR          slang/cfitsio/pgplot in DIR/lib, DIR/include, DIR/pgplot ],
  [ jh_depdir=$withval ], [jh_depdir="$prefix_input"])
AC_SUBST(subdir)
AC_SUBST(prefix_input)
])

dnl#}}}

AC_DEFUN(SAFE_AC_FC_MANGLE, dnl #{{{
[
   AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
   pushdef([AC_MSG_ERROR], [FC_MANGLE_UPCASE=0; FC_MANGLE_SUFFIX= ;])
   AC_FC_FUNC(a,ma)
   case $ma in
      A)	FC_MANGLE_UPCASE=1; FC_MANGLE_SUFFIX="";;
      A_)	FC_MANGLE_UPCASE=1; FC_MANGLE_SUFFIX="_";;
      A__)	FC_MANGLE_UPCASE=1; FC_MANGLE_SUFFIX="__";;
      a)	FC_MANGLE_UPCASE=0; FC_MANGLE_SUFFIX="";;
      a_)	FC_MANGLE_UPCASE=0; FC_MANGLE_SUFFIX="_";;
      a__)	FC_MANGLE_UPCASE=0; FC_MANGLE_SUFFIX="__";;
      *)	AC_MSG_WARN([unknown Fortran mangling convention])
                FC_MANGLE_UPCASE=0; FC_MANGLE_SUFFIX="_"
                ;;
   esac

   popdef([AC_MSG_ERROR])
   AC_SUBST(FC_MANGLE_UPCASE)
   AC_SUBST(FC_MANGLE_SUFFIX)
])dnl #}}}
