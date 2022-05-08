#!/bin/bash

# specify the makefile
target="../src/Makevars"

# PKG_CXXFLAGS
echo "# PACKAGE FLAGS"                                             >  $target
echo "PKG_CXXFLAGS += -I../inst/include/ncl/"                      >> $target
echo "PKG_CXXFLAGS += -I../inst/include/TensorPhylo/"              >> $target
echo "PKG_CXXFLAGS += -I../inst/include/package/"                  >> $target
echo "PKG_CXXFLAGS += -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS" >> $target

# cxx
echo -e ""             >> $target
echo "# CXX STD"       >> $target
echo "CXX_STD = CXX17" >> $target

# package source
echo -e ""                              >> $target
echo "# SOURCES"                        >> $target
echo 'SOURCES = tensoRphyloModules.cpp' >> $target

# NCL sources
while IFS= read -r line || [[ -n "$line" ]]; do
    # echo "Text read from file: $line"
    echo "SOURCES += $line" >> $target
done < "ncl_sources.txt"

# tensorphylo sources
while IFS= read -r line || [[ -n "$line" ]]; do
    # echo "Text read from file: $line"
    echo "SOURCES += $line" >> $target
done < "tp_sources.txt"

# finish up
echo -e ""                                        >> $target
echo "# OBJECTS etc."                             >> $target
echo 'OBJECTS = RcppExports.o $(SOURCES:.cpp=.o)' >> $target
