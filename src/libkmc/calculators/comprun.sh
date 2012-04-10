#!/bin/sh
echo "---- Compiling ----"
echo ""
g++ kmcstandalone.cc -I /people/thnfs/homes/kordt/votca/include/ -L /people/thnfs/homes/kordt/votca/lib -lvotca_tools -o kmcstandalone.exe -Wall -g
echo ""
echo "---- Compiling done. ----"

echo ""
echo "Press enter to continue executing the file."
read nix
echo ""

echo "---- Executing ----"
echo ""
chmod 755 kmcstandalone.exe
LD_LIBRARY_PATH=/people/thnfs/homes/kordt/votca/lib
./kmcstandalone.exe
echo ""
echo "---- Executing done. ----"



