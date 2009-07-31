#! /bin/sh

if test $# -eq 1 ; then
  if test -x "$1" ; then
     ISIS="$1"
  else
     ISIS=`which isis`
  fi
fi

if test -x "$ISIS" ; then
   # solaris /bin/sh fails with 'if ! test' and seems
   # to require an executable statement between
   # 'then' and 'else'.
   notused=0
else
   exit 0
fi

string=`$ISIS -n - <<EOT
() = fprintf(stdout,"just testing");
EOT
`

if test "$string" = "just testing" ; then
  echo "testing here documents.... ok"
  exit 0
else
  echo "testing here documents.... FAILED"
  exit 1
fi
