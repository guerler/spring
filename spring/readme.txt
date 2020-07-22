SPRING: protein dimer structure prediction by mapping monomeric threading
to protein-protein interaction

1. How to install the SPRING Suite?

   a) unpack 'source.tar.bz2' by
      $ tar -xvf source.tar.bz2

      This will create the "source" directory, which hosts the SPRING program.

   b) download SPING libraries from the SPRING template library website
      and unpack them to "source/database"
      $ source/download_lib.sh

   c) (optional) recompile shared libraries used by SPRING
      $ cd source
      $ make

2. How to run SPRING?
   
   a) Make a folder to store input and output files. Put your target sequence
      in a fasta format sequence file "seq.fasta". Put this fasta file under
      the folder you have just created. The name you gave to your target
      sequence inside "seq.fasta" will be completely ignored by SPRING. Make
      sure there are two sequences in the "seq.fasta" file. An examples is
      given at "source/example"
   
   b) Main script for running SPRING is source/runSPRING.pl
      (it is a wrapper script around source/SPRING/spring.py).
      Run it directly without arguments will output the help information.

      The following arguments are mandatory. One example is:

      "sourced/runSPRING.pl -seqname sequence_name -datadir data_dir"

      -datadir means the directory which contains your sequence "seq.fasta"
      -seqname means the unique name of your query sequence. This name is
               used to name temporary folders and to tag PBS jobs. It does 
               not have to be the same as -datadir or the name you gave to
               target sequence inside "seq.fasta". If you want to run more
               than one targets at the same time, make sure different
               targets has different -seqname and different -datadir values.

3. How to interpretate the output files?

   a) Threading output:
      *.hhr               - HHsearch output for target sequences
      TemplateSummary.txt - list of top templates
      init.dat            - target-template alignment in pseudo-PDB format

4. How to report bugs?

   SPRING standalone suite is experimental. If you find bugs (and we
   expect plenty), please contact Chengxin Zhang (zcx@umich.edu) or Yang
   Zhang (zhng@umich.edu). Please make sure you include the following:
   a) input files 
   b) command you used to run the software and its screen output
   c) The operating system you used (e.g. Ubuntu Linux Trusty 14.04 amd6)
      and the output of:
      uname -a; lsb_release -a; head -n 25 /proc/cpuinfo; ulimit -a; free
   d) output files

5. System requirement:

   a) x86_64 Linux (kenerl version >= 2.6.18).
   b) Perl and Python (version 2.7) interpreters should be installed.
      cython and numpy module for Python is requred. We strongly recommend
      using Anaconda Python 2.7, which comes with cython and numpy
   c) Basic compress and decompress package should be installed to support: 
      tar and bunzip2

6. How to cite SPRING?

   SPRING standalone suite is packaged by Chengxin Zhang, Brandon Govindarajoo,
   and Yang Zhang. If you are using the SPRING package, you can cite:
   
   Aysam Guerler, Brandon Govindarajoo, Yang Zhang. Mapping Monomeric
   Threading to Protein-Protein Structure Prediction. Journal of Chemical
   Information and Modeling, 53, 717-725 (2013)

7. License:

   SPRING: Accurate Template Recognition for Protein Complex Structure
   Copyright (C) 2016 Brandon Govindarajoo, Chengxin Zhang, Yang Zhang
   Copyright (C) 2013 Aysam Guerler, Brandon Govindarajoo, Yang Zhang

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
