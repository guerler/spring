// native contacts
class EnergyContacts
{
    // dimensions
    static const int    NTYPE  = 21;
    static const int    NDIST  = 20;
	static const double NSCALE = 2.0;	 

    // matrix
   	Vec < double > mat;
   	
   	// size
   	int ntarget, nmol;
   	
   	// template
   	SpecMolecule* targettmp;
    SpecMolecule* moltmp;
    
    // alignment
    Vec < int > restarget, resmol;
public:

	// obtain interface
	bool construct (SpecMolecule* target, SpecMolecule* mol)
	{
		// replace
		getSCM(target);
		getSCM(mol);

        // backup
        this->targettmp = target;
        this->moltmp    = mol;
        
		// backup size
		ntarget = target->lcalpha.size();
		nmol = mol->lcalpha.size();
        
        // set
        restarget.resize(ntarget);    
        resmol.resize(nmol);

        // setup scoring table
        Format::readList("data/dfire.txt", mat);
                        		
        // return
        return true;
	}
    
	// get energy
	double get(SpecMolecule* target, SpecMolecule* mol, Vec < int >& alntarget, Vec < int >& alnmol)
	{
        // verify
        if (alntarget.size() != ntarget || alnmol.size() != nmol)
        {
            Msg::write("EnergyContacts::construct() : Alignment does not fit molecule size.");
            return 0.0;
        }

        // transfer
        for (int i = 0; i < alntarget.size(); i++)
        {
            if (alntarget[i] != LARGE)
                restarget[i] = target->latom[target->lcalpha[alntarget[i]]].rescode;
            else
                restarget[i] = LARGE;
        }

        // transfer
        for (int i = 0; i < alnmol.size(); i++)
        {
            if (alnmol[i] != LARGE)
                resmol[i] = mol->latom[mol->lcalpha[alnmol[i]]].rescode;
            else
                resmol[i] = LARGE;
        }
                        
        // determine contact energy
        double e = 0.0;
  	    for (int i = 0; i < ntarget; i++)
            for (int j = 0; j < nmol; j++)
                if (restarget[i] != LARGE && resmol[j] != LARGE)
                {
					// index
					int resa = restarget[i];
					int resb = resmol[j];
					int dist = int (targettmp->latom[targettmp->lcalpha[i]].pos.dist(moltmp->latom[moltmp->lcalpha[j]].pos) * NSCALE);
					
					// index
					if (dist < NDIST)
					{
						int index = getIndex(resa, resb, dist);
                    	e += mat[index];
					}
				}

        // return
		return e;
	}
	
    // get index
    int getIndex(int i, int j, int k)
    {
		return i * NTYPE * NDIST + j * NDIST + k;
	}
    
    // get center of mass for side chains
    void getSCM (SpecMolecule* mol)
    {
        // get size
        int n = mol->lcalpha.size();
        
        // reset
		Vec < Vec3 > scm;
        scm.resize(n);
        scm.fill(0.0);

        // center of side chains
        Vec < int > scn;
		scn.resize(n);
		
        // loop
        for (int i = 0; i < mol->latom.size(); i++)
        {
            // get atom
            SpecAtom* atom = &mol->latom[i];
             
            // skip backbone atoms
            if (atom->name == "C" || atom->name == "O" || atom->name == "N")
				continue;
             
            // get residue number
            int nres = atom->resno;
            
            // add coordinates
            scm[nres]+= atom->pos;
            scn[nres]++;          
        }
        
        // calculate average of side chain mass
        for (int i = 0; i < n; i++)
        {
            if (scn[i] > 0)
                scm[i] /= scn[i];
        }
        
        // copy
        SpecMolecule tmp;
        mol->copyCAlpha(&tmp);
        mol->construct(&tmp);
        for (int i = 0; i < n; i++)
        	mol->latom[i].pos = scm[i];
    }
};
