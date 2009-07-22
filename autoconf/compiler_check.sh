#! /bin/sh

if test $# -lt 2 ; then
  echo "Usage compiler_check.sh \$CC \$FC \$HEADAS"
  exit 0
fi

user_cc="$1"
user_fc="$2"

user_cc_version=`$user_cc -dumpversion 2>&1`
if test $? -eq 0 ; then
  cc_dumpversion_worked=1
  user_cc_version_string=`$user_cc --version 2>&1 | head -1`
else
  cc_dumpversion_worked=0
  user_cc_version=`$user_cc -V 2>&1`
  if test $? -eq 0 ; then
    user_cc_version=`$user_cc -V 2>&1 | head -1`
    user_cc_version_string="$user_cc_version"
  else
    user_cc_version=UNKNOWN
    user_cc_version_string=UNKNOWN
  fi
fi

user_fc_version=`$user_fc --version 2>&1`
if test $? -eq 0 ; then
  user_fc_version=`$user_fc --version 2>&1 | head -1`
else
  user_fc_version=`$user_fc -V 2>&1`
  if test $? -eq 0 ; then
    user_fc_version=`$user_fc -V 2>&1 | head -1`
  else
    user_fc_version=UNKNOWN
  fi
fi

# Check for obvious compiler mismatches

printed_warning=0
line='----------------------------------------------------------'
explain=`cat <<EOF

Explanation:
  Mismatched compilers can cause the isis build to fail.
  Although this script attempts to detect a compiler mismatch,
  it is not infallible. If you think your compilers are compatible,
  then continue with the installation. If the build fails,
  it might be because your compilers are mismatched.
  For suggestions on how to solve some common problems, see:
      http://space.mit.edu/cxc/isis/faq.html

EOF`

if test $cc_dumpversion_worked -eq 1 ; then
   user_cc_fc_match=`echo $user_fc_version | sed s,"$user_cc_version",XXXXX,g | grep XXXXX`
   if test x"$user_cc_fc_match" = x"" ; then
      echo ""
      echo "$line"
      echo "WARNING: mismatched C and Fortran compilers?"
      echo "        CC version = $user_cc_version_string"
      echo "        FC version = $user_fc_version"
      printed_warning=1
   fi
fi

if test $# -lt 3 ; then
  if test $printed_warning -ne 0 ; then
      echo "$explain"
      echo "$line"
      echo ""
  fi
  exit 0
fi

headas="$3"

headas_config_status="$headas/BUILD_DIR/config.status"
if test ! -r "$headas_config_status" ; then
   echo '*** Cannot check compatibility of HEADAS compilers ($HEADAS/BUILD_DIR/config.status is unreadable)'
   exit 0
fi

headas_cc=`grep @CC@ $headas_config_status | cut -d, -f3`
headas_fc=`grep @FC@ $headas_config_status | cut -d, -f3`

headas_cc_version=`$headas_cc -dumpversion 2>&1`
if test $? -eq 0 ; then
  cc_dumpversion_worked=1
  headas_cc_version_string=`$headas_cc --version 2>&1 | head -1`
else
  cc_dumpversion_worked=0
  headas_cc_version=`$headas_cc -V 2>&1`
  if test $? -eq 0 ; then
    headas_cc_version=`$headas_cc -V 2>&1 | head -1`
    headas_cc_version_string="$headas_cc_version"
  else
    headas_cc_version=UNKNOWN
    headas_cc_version_string=UNKNOWN
  fi
fi

headas_fc_version=`$headas_fc --version 2>&1`
if test $? -eq 0 ; then
  headas_fc_version=`$headas_fc --version 2>&1 | head -1`
else
  headas_fc_version=`$headas_fc -V 2>&1`
  if test $? -eq 0 ; then
    headas_fc_version=`$headas_fc -V 2>&1 | head -1`
  else
    headas_fc_version=UNKNOWN
  fi
fi

if test ! x"$user_cc_version" = x"$headas_cc_version" ; then
   if test $printed_warning -eq 0 ; then
      echo ""
      echo "$line"
   else
      echo ""
   fi
   echo "WARNING: C compiler differs from HEADAS C compiler?"
   echo "           CC version = $user_cc_version"
   echo "    HEADAS CC version = $headas_cc_version"
   printed_warning=1
fi

if test ! x"$user_fc_version" = x"$headas_fc_version" ; then
   if test $printed_warning -eq 0 ; then
      echo ""
      echo "$line"
   else
      echo ""
   fi
   echo "WARNING: Fortran compiler differs from HEADAS Fortran compiler?"
   echo "           FC version = $user_fc_version"
   echo "    HEADAS FC version = $headas_fc_version"
   printed_warning=1
fi

if test $printed_warning -ne 0 ; then
      echo "$explain"
      echo "$line"
      echo ""
fi
