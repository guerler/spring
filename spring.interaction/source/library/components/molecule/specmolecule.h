// secondary structure element (SSE)
#include "sse.h"
// parameters
#include "specdetails.h"
// atoms
#include "specatom.h"

// molecule object
class SpecMolecule
{
	// parameters
	int 		n;		
	public:
	
	// molecule name or pdbcode
	string             name, info;

	// lists
	Vec < SpecAtom >   latom;
	Vec < SpecAtom >   lhetatom;
	Vec < SSE  >       lsse;
	Vec < int >        lcalpha;
	Vec < char >       lssebyres;

	// tables
	Vec < Vec < SpecAtom* > > tresidue;

	// type tags
	static const char SSECOIL   = '?';
	static const char SSEALPHA  = 'H';
	static const char SSEBETA   = 'S';
	static const char SSETURN   = 'T';	
	static const int  SSEEMPTY  = -1;
	string ATOMCA;

	// construct
	SpecMolecule()
	{
		construct();
	}

	// construction
	void construct (string name = "", string info = "")
	{
		// initialize
		initialize (info, name);
	}

	// open initialize
	void construct(SpecMolecule* mol, string name = "")
	{
		// change name
		if (name == "")
			name = mol->name;
			
		// initialize
		initialize (mol->info, name);

		// setup size and center
		n = mol->size();
	
		// copy lists
		latom.insert(latom.begin(), mol->latom.begin(), mol->latom.end());
		lhetatom.insert(lhetatom.begin(), mol->lhetatom.begin(), mol->lhetatom.end());
		lsse.insert(lsse.begin(), mol->lsse.begin(), mol->lsse.end() );
		lcalpha.insert(lcalpha.begin(), mol->lcalpha.begin(), mol->lcalpha.end() );
		lssebyres.insert(lssebyres.begin(), mol->lssebyres.begin(), mol->lssebyres.end() );
	}

    // merge two molecules
	void construct (SpecMolecule* a, SpecMolecule* b)
	{
       	// get size
       	int na = a->lcalpha.size();

		// append as atom
		construct();
		for (int i = 0; i < a->latom.size(); i++)
		{
      		SpecAtom* atom = &a->latom[i];
       		append (atom->pos, atom->name, atom->res, atom->resno, atom->respdb);
      	}

		for (int i = 0; i < b->latom.size(); i++)
		{
       		SpecAtom* atom = &b->latom[i];
       		append (atom->pos, atom->name, atom->res, na + atom->resno, atom->respdb);
      	}
       	finalize();
   	}
    	
	// copy coordinates
	void readCoordinates(SpecMolecule* mol)
	{
		// assert
		assert (n == mol->n);
		
		// copy
		for (int i = 0; i < n; i++)
		    latom[i].pos = mol->latom[i].pos;
	}

	// copy coordinates
	void readCAlpha(SpecMolecule* mol)
	{
		// assert
		assert (lcalpha.size() == mol->lcalpha.size());
		
		// copy
		for (int i = 0; i < lcalpha.size(); i++)
		    latom[lcalpha[i]].pos = mol->latom[mol->lcalpha[i]].pos;
	}

	// load molecule
	void readStructure(SpecMolecule* mol)
	{	
		// initialize
		int cres = 0;
		if (latom.size() > 0)
			cres = latom[latom.size() - 1].resno;

		// add atoms		
		int nres = 0;        
		for (int i = 0; i < mol->latom.size(); i++)
		{
		    // atom
		    SpecAtom* a = &mol->latom[i];

		    // check
		    if(nres != a->resno)
		    {
                nres = a->resno;
                cres++;
		    }

		    // append
		    append(a->pos, a->name, a->res, cres);
		}

		// prepare
		lssebyres.insert(lssebyres.end(), mol->lssebyres.begin(), mol->lssebyres.end() );
	}

	// read
	bool readSequence (string sequence, string structure = "", string name = "Molecule")
	{        
		// initialize
		initialize (info, name);

		// setup molecule
		for (int i = 0; i < sequence.size(); i++)
		    append(Vec3(NOCOORD), "CA", SpecDetails::getResName(sequence[i]), i);

		// prepare
		lssebyres.resize(sequence.size());
		lssebyres.fill(SSECOIL);

		// finalize
		if (finalize())
			return true;   
		else
			return false; 
	}

	// copy calphas
	void copyCAlpha(SpecMolecule* selection)
	{
		// setup molecule
		selection->construct(name);
	
		// setup atoms
		SpecAtom* a;
		for (int i = 0; i < lcalpha.size(); i++)
		{
		    a = &latom[lcalpha[i]];
		    selection->append (a->pos, a->name, a->res, a->resno, a->respdb, a->occupancy, a->bfactor, a->sse);
		}

		// update
		selection->finalize();
	}

	// copy calphas
	void copyBone(SpecMolecule* selection)
	{
		// setup molecule
		selection->construct(name);
	
		// setup atoms
		SpecAtom* a;
		for (int i = 0; i < latom.size(); i++)
		{
		    a = &latom[i];
		    if (a->name == "CA" || a->name == "CB" || a->name == "N" || a->name == "O" || a->name == "C")
    		    selection->append (a->pos, a->name, a->res, a->resno, a->respdb, a->occupancy, a->bfactor, a->sse);
		}

		// update
		selection->finalize();
	}

	// copy calphas
	void copyHeavy(SpecMolecule* selection)
	{
		// setup molecule
		selection->construct(name);
	
		// setup atoms
		SpecAtom* a;
		for (int i = 0; i < latom.size(); i++)
		{
		    a = &latom[i];
            if (a->name.size() > 0)
		    if (a->name[0] != 'H')
    		    selection->append (a->pos, a->name, a->res, a->resno, a->respdb, a->occupancy, a->bfactor, a->sse);
		}

		// update
		selection->finalize();
	}
    	
	// copy surface
	void copySequence(string& sequence)
	{
		// setup atoms
		sequence = "";
		for (int j = 0; j < lcalpha.size(); j++)
			sequence += SpecDetails::getCharCode(latom[lcalpha[j]].rescode);
	}

	// number of elements
	inline int size()
	{
		return n;
	}

	// residues
	inline int sizeRes()
	{
		return (int) lcalpha.size();
	}

	// append atom
	inline void append (SpecAtom& a)
	{
		append (a.pos, a.name, a.res, a.resno, a.respdb, a.occupancy, a.bfactor, a.sse);
	}

	// append atom
	inline void append (Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, int sse = SSECOIL)
	{
		// increase counter
		n++;
		    	
		// resize
		latom.append (SpecAtom(pos, name, res, resno, respdb, occupancy, bfactor, sse));
	}

	// append hetatom
	inline void appendhet (Vec3 pos, string name, string res, int resno = 0, int respdb = LARGE)
	{
		lhetatom.append (SpecAtom(pos, name, res, resno, respdb));
	}

	// set
	inline void set (int i, Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, int sse = SSECOIL)
	{
		latom[i].set(pos, name, res, resno, respdb, occupancy, bfactor, sse);
	}

	// center of gravity
	inline Vec3 cog()
	{
		// mean vector
		Vec3 tmp (0, 0, 0);
		int cx = 0;
		for (int i = 0; i < n; i++)
		{
            if (latom[i].pos != NOCOORD)
            {
    			tmp += latom[i].pos;
    	   		cx++;
            }
        }
		tmp /= cx;
			
		// return 
		return tmp;
	}

	// center of gravity
	inline void setcenter()
	{
		// mean vector
		Vec3 tmp (0, 0, 0);
		for (int i = 0; i < n; i++)
			tmp += latom[i].pos;
		tmp /= n;
		
        for (int i = 0; i < n; i++)
            latom[i].pos += tmp;
	}
	
	// log
	void lg(bool details = false)
	{
		// sse seq
		string seq;        
		for (int i = 0; i < lsse.size(); i++)
		    seq += lsse[i].getType();

		// title
		Msg::write ("Molecule %s %s [%i, %i]", name.c_str(), seq.c_str(), lsse.size(), lcalpha.size());

		if (details)
		{
		    for (int i = 0; i < n; i++)
			Msg::write ("ATOM %i \t %f,%f,%f,%s,%s,%i,%c", i, latom[i].pos.x, latom[i].pos.y, latom[i].pos.z, latom[i].name.c_str(), latom[i].res.c_str(), latom[i].resno, latom[i].sse);	
		    
		    for (int i = 0; i < lsse.size(); i++)
			Msg::write ("SSE  %i \t %c \t %i \t %i", i, lsse[i].getType(), lsse[i].getStart(), lsse[i].getEnd());
		}
	}

	// finalize
	inline bool finalize()
	{   
		// verify
		if (n == 0)
		{
			Msg::write("No atoms or secondary structure were defined in molecule %s.", name.c_str());
			return false;
		}

		// verify consistency
		if (latom[0].resno != 0)
		{
			Msg::write("Residue numbers inconsistent.");
			return false;
		}
		for (int i = 0; i < latom.size() - 1; i++)
		{
			int diff = latom[i+1].resno - latom[i].resno;
			if (!(diff == 0 || diff == 1))
			{
				Msg::write("Residue numbers inconsistent.");
				return false;
			}
		}

		// load structure
		makeSecondary();
		
		// update secondary structure information if missing
		if (lssebyres.size() == 0)
		{
			lssebyres.resize(latom[n - 1].resno + 1);
			lssebyres.fill(SSECOIL);
		}

		// verify consistency
		if (lssebyres.size() != latom[n - 1].resno + 1)
		{
			Msg::write("Residue numbers inconsistent with secondary information per residue (%i, %i).", lssebyres.size(), latom[n - 1].resno + 1);
			return false;
		}

		// setup calpha	
		finalizeAlpha();

		// setup secondary structure details
		finalizeSecondary();
			
		// return
		return true;
	}

	// check for c-alpha atom
	bool isCA(SpecAtom& a)
	{
		// check
		if(a.name == ATOMCA)
			return true;

		// return
		return false;
	}

	// n-defined
	int ndefined()
	{
		// calculate
		int ndefined = 0;
		for(int i = 0; i < lcalpha.size(); i++)
		{
            if (latom[lcalpha[i]].pos == NOCOORD)
                continue;
            ndefined++;			
        }

		// return
		return ndefined;
	}
	
	// rmsd upon forcetopology
	double rmsd(SpecMolecule* mol)
	{
		// verify
		if (lcalpha.size() != mol->lcalpha.size())
			Msg::error ("SpecMolecule::rmsd", "Number of atoms not equivalent");
		
		// calculate
		double dist = 0;
		int naligned = 0;
		for(int i = 0; i < lcalpha.size(); i++)
		{
            if (latom[lcalpha[i]].pos == NOCOORD || mol->latom[mol->lcalpha[i]].pos == NOCOORD)
                continue;
			dist += pow (latom[lcalpha[i]].pos.dist(mol->latom[mol->lcalpha[i]].pos), 2);
            naligned++;			
        }

		// return
		return sqrt ( dist / naligned );
	}
	
    // template model score
    double tms(SpecMolecule* mol)
	{
		// verify
		if (lcalpha.size() != mol->lcalpha.size() && lcalpha.size() != 0)
			Msg::error ("SpecMolecule::rmsd", "Number of atoms not equivalent");

        // calculate tms
        double Ln = lcalpha.size();
        double d0 = 1.24 * pow(Ln - 15, (double) (1.0 / 3.0)) - 1.8;
        double tms = 0;
        for(int i = 0; i < Ln; i++)
        {
            if (latom[lcalpha[i]].pos == NOCOORD || mol->latom[mol->lcalpha[i]].pos == NOCOORD)
                continue;            
            tms += 1 / (1 + pow(latom[lcalpha[i]].pos.dist(mol->latom[mol->lcalpha[i]].pos) / d0, 2));
        }
        return tms / Ln;
    }	

    // this function overwrites the original pdb residue numbers
    void enumerateResidues()
	{
        for(int i = 0; i < latom.size(); i++)
            latom[i].respdb = latom[i].resno;
    }    
private:	
	// initialize class.x
	void initialize(string& info, string& name)
	{
		// definitions
		ATOMCA = "CA";
		 
		// initialize
		this->info  	= info;
		this->name  	= name;
		this->n     	= 0;
	
		// clear all lists
		latom.clear();
		lhetatom.clear();		
		lsse.clear();
		lssebyres.clear();
		lcalpha.clear();

		// clear all tables
		tresidue.clear();
	}

	// setup calpha list
	void finalizeAlpha()
	{
		// tresidue
		tresidue.resize(lssebyres.size());

		// lcalphas
		lcalpha.resize(lssebyres.size());
		lcalpha.fill(0);

		int cres = 0; 
		if (n > 0 && lssebyres.size() > 0)
		{
			// last residue
			int last = latom[0].resno;
	    
			// loop
			for (int i = 0; i < n; i++)
			{								
				// residue counter
				if (last != latom[i].resno)
				{
					last = latom[i].resno;
					cres++;
				}

				// is c-alpha atom
				if (isCA(latom[i]))
					lcalpha[cres] = i;
			
				// append residue table
				tresidue[cres].push_back(&latom[i]);
			
				// sse
				latom[i].sse = lssebyres[cres];
			}
		}
	}

	// setup secondary structure details
	void finalizeSecondary()
	{
		char last = 0;
		int nsse = 0;
		bool valid = false;

		// over all atoms
		for (int i = 0; i < n; i++)
		{            
			// append sses
			if (last != latom[i].sse)
			{                    
				if (valid)
					lsse[nsse-1].setEnd(latom[i-1].resno);

				if (latom[i].sse == SSEALPHA || latom[i].sse == SSEBETA)
				{
					lsse.push_back(SSE (latom[i].sse, latom[i].resno));
					nsse++;
				}

				// update
				last    = latom[i].sse;
				valid   = (last == SSEALPHA || last == SSEBETA);
			}

			// backup sse index
			if (latom[i].sse == SSEALPHA || latom[i].sse == SSEBETA)
				latom[i].sseno = nsse - 1;
		}

		// finalize
		if (valid)
			lsse[nsse-1].setEnd(latom[n-1].resno);
	}
	
	// get secondary structure
	void makeSecondary()
	{
    	// variables
    	int j1, j2, j3, j4, j5;
    	double dis13, dis14, dis15, dis24, dis25, dis35;
        
        // copy calpha positions to array
        Vec <Vec3> lpos;
        for (int i = 0; i < latom.size(); i++)
            if(isCA(latom[i]))
                lpos.push_back(latom[i].pos);        
        
		// prepare sse by residue list
		int n = lpos.size();

		// prepare
		lssebyres.resize(n);
		lssebyres.fill(SSECOIL);	
		
		// return
		if (n == 0)
			return;	

		// 1->coil, 2->helix, 3->turn, 4->strand
		for (int i = 2; i < n - 2; i++)
		{
     		j1=i-2;
     		j2=i-1;
     		j3=i;
     		j4=i+1;
     		j5=i+2;
     		if(j1 > 0 && j5 < n)
            {
                // get distances
			    dis13 = lpos[j1].dist(lpos[j3]);
			    dis14 = lpos[j1].dist(lpos[j4]);
			    dis15 = lpos[j1].dist(lpos[j5]);
			    dis24 = lpos[j2].dist(lpos[j4]);
			    dis25 = lpos[j2].dist(lpos[j5]);
			    dis35 = lpos[j3].dist(lpos[j5]);
			    
                // helix
        		double delta = 2.5; //2.1 ->2.5
        		if(fabs(dis15-6.37) < delta)
        		if(fabs(dis14-5.18) < delta)
        		if(fabs(dis25-5.18) < delta)
        		if(fabs(dis13-5.45) < delta)
        		if(fabs(dis24-5.45) < delta)
        		if(fabs(dis35-5.45) < delta)
        		{
        			lssebyres[i] = SSEALPHA;
        			continue;
                }

                // strand
        		delta = 2.0; //1.42 ->2.0
        		if(fabs(dis15-13.0) < delta)
        		if(fabs(dis14-10.4) < delta)
        		if(fabs(dis25-10.4) < delta)
        		if(fabs(dis13-6.1) < delta)
        		if(fabs(dis24-6.1) < delta)
        		if(fabs(dis35-6.1) < delta)
        		{
        			lssebyres[i] = SSEBETA;
        			continue;
                }        			
        
                // turn
                if (dis15 < 8.0)
                {
                    lssebyres[i] = SSETURN;
                    continue;
                }
        
                // coil
           		lssebyres[i] = SSECOIL;				    
            }
        }

    	// x%x >> xxx
    	for (int i = 0; i < lssebyres.size() - 2; i++)
    	{
         	if(lssebyres[i] != SSECOIL)
    		{
            	int j = lssebyres[i];
    			if(lssebyres[i+2] == j)
    				lssebyres[i+1] = j;
    		}
        }
    
    	// -x- >> -----
    	for (int i = 0; i < lssebyres.size() - 2; i++)
    	{
         	if(lssebyres[i] == SSECOIL)
    		if(lssebyres[i+1] != SSECOIL)
    		if(lssebyres[i+2] == SSECOIL)
    			lssebyres[i+1] = SSECOIL;
        }
    }  	
};

