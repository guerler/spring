// atom
struct SpecAtom
{   
    // const
    static const char SSECOIL   = '?';
	static const int  SSEEMPTY  = -1;
	
	// data
	Vec3       pos;
	string     name;
	string     res;
	int        rescode;	
	int        resno;
	int        respdb;
	int        sseno;	
	char       sse;

    // additional pdb
	double     occupancy;
    double     bfactor;
    
	// construct
	SpecAtom(Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, char sse = SSECOIL, int sseno = SSEEMPTY)
	{
	        set (pos, name, res, resno, respdb, occupancy, bfactor, sse, sseno);      
	}
    
	// setup atom
	inline void set (Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, char sse = SSECOIL, int sseno = SSEEMPTY)
	{
		// data
		this->pos        = pos;
		this->name       = Lib::strclean(name);
		this->res        = res;
		this->resno      = resno;
		this->occupancy  = occupancy;
		this->bfactor    = bfactor;

        // original residue identifier from PDB
        if (respdb == LARGE)
            this->respdb = resno;
		else
            this->respdb = respdb;
	
        // secondary information
		this->sse        = sse;
		this->sseno      = sseno;

        // get residue code
		this->rescode    = SpecDetails::getResCode(res);
	}          
};

