// product details
#define PRODUCT_NAME 	"SPRING-Release"
#define PRODUCT_COPY 	"University of Michigan (Zhang Lab.)"
#define PRODUCT_AUTHORS	"Aysam Guerler (aysam.guerler@gmail.com)"
#define PRODUCT_VERSION	"3.12"

// const
#define TOLERANCE       0.01
#define EPSILON         0.00000000001
#define MINUSEPSILON    -0.00000000001
#define INT_MAX         999999999
#define LARGE           INT_MAX
#define NOCOORD         999
#define BSIZE           1048576

// natural constants
#define COULOMB         332.0636
#define PI              M_PI

// generic string buffer
char strbuf[2048];

// std includes
using namespace std;
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <cctype>
#include <map>
#include <string.h>
#include <errno.h>
#include <string>
#include <list>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <cstdarg>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <stddef.h>
#include <dirent.h>

//
// external
//

// kabsch
#include "components/kabsch.h"

//
// internal
//

// general library
#include "components/clock.h"
#include "components/vec.h"
#include "components/vec3.h"
#include "components/msg.h"
#include "components/convert.h"
#include "components/fname.h"			 				 
#include "components/lib.h"
#include "components/file.h"
#include "components/molecule/specmolecule.h"
#include "components/format.h"
#include "components/sequencealignment.h"
#include "components/trans.h"
#include "components/transmatrix.h"
#include "components/storage.h"
#include "components/hash.h"
#include "components/listmin.h"
#include "components/exclude.h"

// load configuration
#include "config.h"

// include plugins
#include "components/tmalign.h"
#include "components/tmscore.h"
#include "components/pulchra.h"

// molecular interaction quality measures
#include "components/interfacequality/fnat.h"
#include "components/interfacequality/irmsd.h"
#include "components/interfacequality/irmsdinterface.h"
#include "components/interfacequality/mintm.h"
#include "components/interfacequality/interfacequality.h"

// pymol
#include "components/pymol/pymol.h"

// threading
#include "components/threading/hhsearch.h"

