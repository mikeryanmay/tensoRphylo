#!/bin/bash

# specify the makefile
source="MakevarsTemplate"
target="../src/Makevars"

# copy the source file
cp $source $target

# find all inst/include cpp files
instFiles=($(find "../inst/include" -type f -name "*.cpp"))
sed -i '' -e "s|SRCPLACEHOLDER|${instFiles[*]}|" $target

# find all source cpp files
srcFiles=($(find "../src" -type f -name "*.cpp"))
sed -i '' -e "s|INCPLACEHOLDER|${srcFiles[*]}|" $target
