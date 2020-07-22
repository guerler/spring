#!/bin/bash
#### help message ####
CURDIR=`dirname "$(readlink -e $0)"`
LIBDIR="$CURDIR/database"
if [ ! -s "$CURDIR/runSPRING.pl" ];then
    echo "ERROR! Please put this script inside SPRING suite directory"
    exit
fi
mkdir -p $LIBDIR

#### download amino-library database ####
echo "Downloading amino-library database to $LIBDIR/SPRING"
cd $LIBDIR
echo "wget http://zhanglab.ccmb.med.umich.edu/spring/package/amino-library.tar.bz2"
wget -o log -c http://zhanglab.ccmb.med.umich.edu/spring/package/amino-library.tar.bz2
tar -xvf amino-library.tar.bz2
rm -f amino-library.tar.bz2 log

#### download pdb database ####
echo "Downloading SPRING database to $LIBDIR/pdb"
cd $LIBDIR
echo "wget http://zhanglab.ccmb.med.umich.edu/spring/package/pdb.tar.bz2"
wget -o log -c http://zhanglab.ccmb.med.umich.edu/spring/package/pdb.tar.bz2
tar -xvf pdb.tar.bz2
rm -f pdb.tar.bz2 log

#### download pdb database ####
echo "Downloading SPRING database to $LIBDIR/pdb"
cd $LIBDIR
echo "wget http://zhanglab.ccmb.med.umich.edu/spring/package/SPRING.tar.bz2"
wget -o log -c http://zhanglab.ccmb.med.umich.edu/spring/package/SPRING.tar.bz2
tar -xvf SPRING.tar.bz2
rm -f SPRING.tar.bz2 log
