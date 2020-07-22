#PBS -A zhng_flux 
#PBS -q flux
#PBS -d /tmp/
#PBS -l nodes=1:ppn=1,pmem=200mb,walltime=40:00:00
#PBS -o /nfs/amino-home/bgovi/projects/StructureSearch/output/0/structSearch.out
#PBS -e /nfs/amino-home/bgovi/projects/StructureSearch/output/0/structSearch.err
#PBS -N structSearch_0
#PBS -M bgovi@umich.edu
#PBS -V
python /nfs/amino-home/bgovi/projects/StructureSearch/structureCompare.py 0
