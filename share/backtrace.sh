#! /bin/sh

#   Author:  John C. Houck <houck@space.mit.edu>
#  Written:  March 2007
#
#  This script filters backtrace output to map hex addresses
#  to the associated function name, file name, and line number.
#
#  Usage:  backtrace.sh trace.txt $bindir/isis
#
#  Requirements:  glibc, binutils (addr2line)
#

ADDR2LINE="addr2line"

if test $# -eq 0 ; then
  echo "Usage:  $0 <backtrace.txt> [path-to-executable]"
  exit 0
fi

FILE="$1"

if test $# -eq 1 ; then
  BINARY=`which isis`
else
  BINARY="$2"
  if test "`basename $BINARY`" = "$BINARY" ; then
     BINARY=`which $BINARY`
  fi
fi

if ! test -x "$BINARY" ; then
   echo "*** Error: wrong path to executable:  $2"
   exit 1
fi

echo "#"
echo "# --- isis backtrace ---"
echo "#"
echo "# `date`"
echo "# `uname -a`"
echo "# $BINARY"
echo "#"
exec <"$FILE"
while read LINE
do
  case "$LINE" in

     *"("*")["*"]" )
        address=`echo "$LINE" | tr '[]' , | cut -d ',' -f2`
        echo "     $LINE" | sed "s%\[$address\]%%"
        ;;

     *"["*"]" )
        address=`echo "$LINE" | tr '[]' , | cut -d ',' -f2`
        info=`$ADDR2LINE --functions --basenames -e "$BINARY" "$address" | tr '\n' ' '`
        if ! test "$info" = "?? ??:0 " ; then
           echo "$info"
        fi
         ;;

     "*"* )
       ;;

     "" )
       ;;

     * )
        echo "$LINE"
        ;;
  esac
done

echo "#"
echo "#-end-backtrace"
