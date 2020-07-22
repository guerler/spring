# HHPaths.pm 
# (c) J. Soeding, A. Hauser 2012

# HHsuite version 2.0

# PLEASE INSERT CORRECT PATHS AT POSITIONS INDICATED BY ... BELOW
# THE ENVIRONMENT VARIABLE HHLIB NEEDS TO BE SET TO YOUR LOCAL HH-SUITE DIRECTORY, 
# AS DESCRIBED IN THE HH-SUITE USER GUIDE AND README FILE

package HHPaths;

# This block can stay unmodified
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION = "version 2.0.12 (Feb 2012)";
our @ISA     = qw(Exporter);
our @EXPORT  = qw($VERSION $hhlib $hhdata $hhbin $hhscripts $execdir $datadir $ncbidir $dummydb $pdbdir $dsspdir $dssp $cs_lib $context_lib);

##############################################################################################
# PLEASE COMPLETE THE PATHS ... TO PSIPRED AND OLD-STYLE BLAST (NOT BLAST+) (NEEDED FOR PSIPRED) 
our $execdir = "~/package/spectrum.next/tool/spring/scripts/psipred/bin";  # Where the PSIPRED V2 programs have been installed
our $datadir = "~/package/spectrum.next/tool/spring/scripts/psipred/data"; # Where the PSIPRED V2 data files have been installed
our $ncbidir = "~/package/spectrum.next/tool/spring/scripts/blast/bin";    # Where the NCBI programs have been installed (for PSIPRED in addss.pl)

##############################################################################################
# PLEASE COMPLETE THE PATHS ... TO YOUR LOCAL PDB FILES, DSSP FILES ETC.
our $pdbdir  =  "";  # where are the pdb files? Used in hhmakemodel.pl
our $dsspdir =  "";  # where are the dssp files? Used in addss.pl
our $dssp    =  "";  # where is the dssp binary? Used in addss.pl
##############################################################################################

# The lines below probably do not need to be changed

# Setting paths for hh-suite perl scripts
our $hhlib    = $ENV{"HHLIB"};     # main hh-suite directory
our $hhdata   = $hhlib."/data";    # path to data directory for hhblits, example files
our $hhbin    = $hhlib."/bin";     # path to cstranslate (path to hhsearch, hhblits etc. should be in $PATH)
our $hhscripts= $hhlib."/scripts"; # path to hh perl scripts (addss.pl, reformat.pl, hhblitsdb.pl etc.)
our $dummydb  = $hhdata."/do_not_delete"; # Name of dummy blast db for PSIPRED (single sequence formatted with NCBI formatdb)

# HHblits data files
our $cs_lib = "$hhdata/cs219.lib";
our $context_lib = "$hhdata/context_data.lib";

# Add hh-suite scripts directory to search path
$ENV{"PATH"} = $hhscripts.":".$ENV{"PATH"}; # Add hh scripts directory to PATH

return 1;
