#! /bin/sh

if test $# -lt 2 ; then
  echo "Usage compiler_check.sh \$CC \$FC \$HEADAS"
  exit 0
fi

user_cc="$1"
user_fc="$2"

user_cc_version=`$user_cc -dumpversion 2>&1`
if test $? -eq 0 ; then
  user_cc_version_string=`$user_cc --version 2>&1 | head -1`
else
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

 If the build fails, this may be the cause.
 For help resolving common installation problems, see:
    http://space.mit.edu/cxc/isis/faq.html

EOF
`

user_cc_fc_match=`echo $user_fc_version | sed s,"$user_cc_version",XXXXX,g | grep XXXXX`
if test x"$user_cc_fc_match" = x"" ; then
   echo ""
   echo "$line"
   echo "WARNING: possibly mismatched C and Fortran compilers:"
   echo "        CC version = $user_cc_version_string"
   echo "        FC version = $user_fc_version"
   printed_warning=1
fi

if test $# -lt 3 ; then
  if test $printed_warning -ne 0 ; then
      echo "$explain"
      echo "$line"
      echo ""
  fi
  exit 0
fi

###############################################################
# 4/2010:  As of heasoft-6.9 compiler information is no longer
# stored in $HEADAS/BUILD_DIR/config.status so we can no longer
# check the compatibility of HEADAS compilers.
exit 0
###############################################################

headas="$3"

headas_config_status="$headas/BUILD_DIR/config.status"
if test ! -r "$headas_config_status" ; then
   echo '*** Cannot check compatibility of HEADAS compilers ($HEADAS/BUILD_DIR/config.status is unreadable)'
   exit 0
fi

# heasoft-6.7 introduced a new format for config.status
symbol=`grep @CC@ $headas_config_status`
if test "x$symbol" != "x" ; then
   headas_cc=`grep @CC@ $headas_config_status | cut -d, -f3`
   headas_fc=`grep @FC@ $headas_config_status | cut -d, -f3`
else
   headas_cc=`grep '"CC"' $headas_config_status | cut -d'"' -f4`
   headas_fc=`grep '"FC"' $headas_config_status | cut -d'"' -f4`
fi

headas_cc_version=`$headas_cc -dumpversion 2>&1`
if test $? -eq 0 ; then
  headas_cc_version_string=`$headas_cc --version 2>&1 | head -1`
else
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
   echo "WARNING: C compiler does not match HEADAS C compiler:"
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
   echo "WARNING: Fortran compiler does not match HEADAS Fortran compiler:"
   echo "           FC version = $user_fc_version"
   echo "    HEADAS FC version = $headas_fc_version"
   printed_warning=1
fi

if test $printed_warning -ne 0 ; then
      echo "$explain"
      echo "$line"
      echo ""
fi
