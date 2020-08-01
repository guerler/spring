// define release
#define SPRING_RELEASE

// parameters
struct SpringConfig
{
    // main thresholds
    static const double MINSPRING     = -LARGE;
    static const double MINZ          = -LARGE;
    static const bool   WITHREF       = true;

    // benchmark threshold
    static const double SEQCUTOFF     = LARGE;

    // search depth
    static const int NTMPL            = 200;
    static const int NTMPL_PER_TARGET = 1000;
    static const int NTMPL_PER_PDB    = 10;

    // clash configuration
    static const double RADIUS2       = 9.0;

    // minimum rmsd of resulting models
    static const double MINRMSD       = 4.0;

    // score configuration
    static const double W0            = 12.0;  // 14.5, 12.0 no refinement (CHECK SCORING FUNCTION)
    static const double W1            = 1.4;   // 0.95, 1.4 no refinement
    static const double W2            = 0.0;   // -8.8, -1800 no refinement
    static const double NORMALIZE     = 6.5;   // 6.7, 6.5 no refinement

    // refinement
    static const bool   REFINEMENT    = true;
    static const bool   SUGGESTEDMONO = true;

    #ifdef SPRING_RELEASE
    static const int    NTOP          = 5;
    static const bool   DETAILS       = false;
    static const bool   WRITE         = false; // writes PDB-FILE
    #else
    static const int    NTOP          = 100;
    static const bool   DETAILS       = false;
    static const bool   WRITE         = false;
    #endif  
};
