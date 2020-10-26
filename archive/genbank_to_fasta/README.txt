Installation Instructions:

1) Make sure you have a working version of Biopython installed. This software was
tested on version 1.70, though other versions may work. It can be found
at: http://biopython.org/.

2) The script will automatically use the version of python set in your environment. If your 
Biopython version uses a different version of python, adjust the first line of the 
genbank_to_fasta.py file to point to your preferred python binary.

3) Move the genbank_to_fasta.py file to an executable directory. On many systems,
it is /usr/local/bin/ or /usr/bin/. You may have to make the file executable like so:
$ chmod ug+rx /usr/local/bin/genbank_to_fasta.py

4) From the command line type 'genbank_to_fasta.py -h' for usage information.

