#!/bin/bash

# make a temporary file directory
echo "Creating temporary directory."
mkdir -p tmp
cd tmp

# first, clone an up-to-date version of NCL
echo "Cloning TensorPhylo."
git clone https://mrmay@bitbucket.org/mrmay/tensorphylo.git

# check if inst/include/TensorPhylo exists
if [ -d "../../inst/include/TensorPhylo" ]
then
  echo "Found inst/include/TensorPhylo"
else
  echo "Did not find inst/include/TensorPhylo"
  mkdir -p ../../inst/include/TensorPhylo
fi
target="../../../../inst/include/TensorPhylo"

# move into tensorphylo/src
cd tensorphylo/src

# check if the whitelist exists
whitelist="../../../tp_sources.txt"
if [ -f "$whitelist" ]
then
  echo -n > $whitelist
else
  touch $whitelist
fi

########
# DATA #
########

# transfer data files
echo "Updating Data."
mkdir -p ../../../../inst/include/TensorPhylo/Data

# find matching files
IFS=$'\n'
headerFiles=($(find "Data" -type f -name "*.h*"))
sourceFiles=($(find "Data" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;


#############
# INTERFACE #
#############

# transfer data files
echo "Updating Interface."
mkdir -p ../../../../inst/include/TensorPhylo/Interface

# find matching files
IFS=$'\n'
headerFiles=($(find "Interface" -type f -name "*.h*"))
sourceFiles=($(find "Interface" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

##############
# LIKELIHOOD #
##############

# transfer data files
echo "Updating Likelihood."
mkdir -p ../../../../inst/include/TensorPhylo/Likelihood

# find matching files
IFS=$'\n'
headerFiles=($(find "Likelihood" -type f -name "*.h*"))
sourceFiles=($(find "Likelihood" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

##############
# PARAMETERS #
##############

# transfer data files
echo "Updating Parameters."
mkdir -p ../../../../inst/include/TensorPhylo/Parameters

# find matching files
IFS=$'\n'
headerFiles=($(find "Parameters" -type f -name "*.h*"))
sourceFiles=($(find "Parameters" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

################
# SYNCH EVENTS #
################

# transfer data files
echo "Updating SynchronousEvents."
mkdir -p ../../../../inst/include/TensorPhylo/SynchronousEvents

# find matching files
IFS=$'\n'
headerFiles=($(find "SynchronousEvents" -type f -name "*.h*"))
sourceFiles=($(find "SynchronousEvents" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

##########
# TENSOR #
##########

# transfer data files
echo "Updating Tensor."
mkdir -p ../../../../inst/include/TensorPhylo/Tensor

# find matching files
IFS=$'\n'
headerFiles=($(find "Tensor" -type f -name "*.h*"))
sourceFiles=($(find "Tensor" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

########
# TEST #
########

# transfer data files
echo "Updating Test."
mkdir -p ../../../../inst/include/TensorPhylo/Test

# find matching files
IFS=$'\n'
headerFiles=("Test/Utils.h")
sourceFiles=("Test/Utils.cpp")
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

#########
# Utils #
#########

# transfer data files
echo "Updating Utils."
mkdir -p ../../../../inst/include/TensorPhylo/Utils

# find matching files
IFS=$'\n'
headerFiles=($(find "Utils" -type f -name "*.h*"))
sourceFiles=($(find "Utils" -type f -name "*.cpp"))
unset IFS

# loop over header files
for file in "${headerFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
done;

# loop over source files
for file in "${sourceFiles[@]}"; do
  echo "Updating $file"
  # move the file
  rsync -R $file $target
  # also store it in the whitelist
  echo "../inst/include/TensorPhylo/$file" >> $whitelist
done;

# cleanup
echo "Cleaning up."
cd ../../..
rm -rf tmp

# update makevars
bash regenerateMakevars.sh
