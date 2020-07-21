------------------------------------------------------------------------------------------
 INSTRUCTIONS TO SET UP THE SPRING-DATABASE (i.e. PDB-REPOSITORY) with scripts
------------------------------------------------------------------------------------------
Part 1.A - Setting up the indexing
 1. Generate a text file containing all proteins (pdbcode+chain) separated by lines
 1. Download original PDB files with spring_install/pdbload.sh (if not already done)
 2. Execute spring_install/pdbchains.sh (it will show which parameters are required)
 3. Execute spring_install/pdbfasta.sh
 4. Execute spring_install/make_fastadb.sh
 5. Execute spring_install/pdbcomplex.sh
 6. Execute spring_install/pdbsplit.sh
 7. Execute spring_install/pdbindex.sh (PSI-BLAST required)

Part 1.B - Running HHsearch (read HHsearch manual in hhsearch/hhsuite if neccessary)
 8. Update local pathes in hhsearch/hhsuite/lib/hh/scripts/HHPaths.pm
 9. Update hhsearch/hhmake.sh
10. Test HHsearch by running hhsearch/hhmake.sh
11. Execute hhsearch/hhmakeall.sh
12. Test HHsearch by running hhsearch/hhsearch.sh
13. Execute hhsearch/hhsearchall.sh
14. Execute spring_install/make_hhmdb.sh

------------------------------------------------------------------------------------------
 INSTRUCTIONS TO SET UP A TARGET-LIBRARY (i.e. TARGET-REPOSITORY) with scripts
------------------------------------------------------------------------------------------
Part 2.A - Setting up the target repository
 1. Download a sequence database in fasta format
 2. Execute spring_other/fastadb_sequences.sh
 3. Execute spring_other/fastadb_list.sh
 4. Execute HHsearch (see Part 1.B above)

Part 2.B - For benchmarking
 5. Put crystal structure chains into your target-repository in the following format
    {your-target-repository}/chains/{first to chain name characters}/{chain name}
 6. Execute blast/psiblast_mpi.sh
 7. Execute blast/psiblast_exclusion.sh

------------------------------------------------------------------------------------------
 INSTRUCTIONS TO RUN SPRING
------------------------------------------------------------------------------------------
 1. Get familiar with springconfig.h 
 2. Update config.cnf
 3. Type 'make' in SPRING main directory
 4. Execute ./spring or spring_mpi.sh for larger sets (it will required parameters)

------------------------------------------------------------------------------------------
 HOW TO COMPARE RESULTS
------------------------------------------------------------------------------------------
 1. Get familiar with spring_compare.sh, which compares spring_mpi.sh outputs
 2. A MS-Excel template for the compare_x.txt file is in the performance directory
