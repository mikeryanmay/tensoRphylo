#!/bin/bash

# make a temporary file directory
echo "Creating temporary directory."
mkdir -p tmp
cd tmp

# first, clone an up-to-date version of NCL
echo "Cloning NCL."
git clone https://github.com/mtholder/ncl.git

# check if inst/include/ncl exists
if [ -d "../../inst/include/ncl" ]
then
  echo "Found inst/include/ncl"
else
  echo "Did not find inst/include/ncl/"
  mkdir -p ../../inst/include/ncl/ncl/
fi

# move into ncl/ncl
cd ncl/ncl

# transfer the header files to inst/include/ncl
echo "Updating header files."
for header in *.h; do
  # make sure the file exists
  [ -f "$header" ] || break
  # move the file
  echo "Updating $header".
  cp $header ../../../../inst/include/ncl/ncl/
done;

# check if the whitelist exists
whitelist="../../../ncl_sources.txt"
if [ -f "$whitelist" ]
then
  echo -n > $whitelist
else
  touch $whitelist
fi

# transfer the implementation files to inst/include/ncl
echo "Updating header files."
for impl in *.cpp; do
  # make sure the file exists
  [ -f "$impl" ] || break
  # move the file
  echo "Updating $impl".
  cp $impl ../../../../inst/include/ncl/ncl/
  # keep track of the file in the whitelist
  echo "../inst/include/ncl/ncl/$impl" >> $whitelist
done;

# cleanup
echo "Cleaning up."
cd ../../..
rm -rf tmp

# update makevars
bash regenerateMakevars.sh
