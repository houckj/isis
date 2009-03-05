#! /bin/sh

# This script copies a list of files into a specified
# directory tree, then makes a tar file.
#
# Usage:
#   ./mkdist.sh  /tmp/foo-dir
#
# Here's one way to make a draft input file-list:
#  du -a | grep -v .svn | tac | awk '{print $2}'
#

if ! test $# -eq 1 ; then
   echo "Usage:" `basename $0` "<target-dir>"
   exit 0
fi

target_dir="$1"
source_dir="."
file_list="files.lis"

if ! test -f "$file_list" ; then
  echo '*** Wrong source directory?' "-- I don't see $file_list"
  exit 1
fi

if test -d "$target_dir" ; then
  echo '*** Error:' "$target_dir already exists."
  exit 1
else
  mkdir -p "$target_dir" || exit 1
fi

while read item
do
  src="$source_dir/$item"
  dest="$target_dir/$item"

  # If a file or directory is listed before the directory
  # that contains it, we have to create the destination
  # directory before copying the file.

  dir=`dirname $dest`
  if ! test -d "$dir" ; then
      mkdir -p $dir || exit 1
  fi

  if test -d "$src" ; then
     mkdir $dest || exit 1
  elif test -f "$src" ; then
     /bin/cp -a $src $dest || exit 1
  else
     echo '*** Error:' "$src is missing."
     exit 1
  fi
done < $file_list

rootdir=`dirname $target_dir`
dirname=`basename $target_dir`
tarfile="${dirname}.tar.gz"
cd $rootdir ; tar czf ${tarfile} ${dirname}

echo "Created $rootdir/$tarfile"
