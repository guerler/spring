// define release
#define SPRING_RELEASE

// parameters
struct SpringConfig
{
    // main thresholds
    static const double MINSPRING = -LARGE;
    static const double MINZ      = 20.0;
    static const bool   WITHREF   = true;
	
    // benchmark threshold
    static const double SEQCUTOFF = 1.0;

    // search depth
    static const int    NTMPL     = 100;
    static const int    NTMPL_PDB = 10;
    
    // threshold for main features
    static const double MINTM     = 0.1;
    static const double MINENERGY = 0.1;

    // clash configuration
    static const double RADIUS2   = 9.0;
    static const double MAXCLASH  = 0.1;
        
    // score configuration
    static const double W1        = 0.1;	// 0.1 energy weighting factor
    static const double NORMALIZE = 0.32;
	    		    
    #ifdef SPRING_RELEASE
    static const int    NTOP      = 1;
    static const bool   DETAILS   = false;
    static const bool   WRITE     = true; // writes PDB-FILE
    #else
    static const int    NTOP      = 100;
    static const bool   DETAILS   = false;
    static const bool   WRITE     = false;
    #endif  
};
