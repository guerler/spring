// gobal interaction rmsd
class IRmsd
{
	// interface
	Vec < Vec3 > inter, decoy;	
		
	// parameters
	Vec < Vec3 > r;
	Vec3 t;
	
	// reference
    SpecMolecule *reftarget, *refmol;	
	
	// interface residues
	int ninterface;
public:
    		
	// obtain interface
	bool construct (SpecMolecule* target, SpecMolecule* mol)
	{	
        // reference
        this->reftarget = target;
        this->refmol = mol;
        
        // compare
        ninterface = target->lcalpha.size() + mol->lcalpha.size();
        
        // verify
        if (ninterface == 0)
            return false;
        
        // return
        return true;
	}
    
    // initialize orientation
    bool initialize(SpecMolecule* target, SpecMolecule* mol, bool translate = true)
	{
        // verify size
		if (target->lcalpha.size() != reftarget->lcalpha.size() || mol->lcalpha.size() != refmol->lcalpha.size())
        {
            Msg::error ("IRmsd::initialize()", "Size of reference and target differs.");
            return false;
        }
				
		// reset
		inter.clear();
		decoy.clear();
		
		// collect assigned interface residues
		for (int i = 0; i < target->lcalpha.size(); i++)
		{
			if (target->latom[target->lcalpha[i]].pos != NOCOORD)
			{
				inter.push_back(reftarget->latom[reftarget->lcalpha[i]].pos);
				decoy.push_back(target->latom[target->lcalpha[i]].pos);				
			}
		}
		for (int i = 0; i < mol->lcalpha.size(); i++)
		{
			if (mol->latom[mol->lcalpha[i]].pos != NOCOORD)
			{
                inter.push_back(refmol->latom[refmol->lcalpha[i]].pos);
				decoy.push_back(mol->latom[mol->lcalpha[i]].pos);
			}
		}

		// verify
		if (decoy.size() < 5)
			return false;
		
		// align interface
		if (translate)
		{
            Trans::align (decoy, inter, r, t);
            Trans::apply (decoy, r, t);
        }
        
        // return
        return true;
    }

    // coverage    
	double getcoverage()
	{
        // no interface found
        if (ninterface == 0)
            return 0.0;

		// return
        return (double) decoy.size() / (double) ninterface;
	}
        
	// root mean square deviation
	double getrmsd()
	{
        // no interface found
        if (ninterface == 0 || decoy.size() == 0)
            return 0.0;

		// calculate
		double d = 0;
		for (int i = 0; i < decoy.size(); i++)
   			d += pow (decoy[i].dist(inter[i]), 2);
		
		// return
        return sqrt ( d / decoy.size() );
	}

    // template model score
    double gettms()
	{
        // no interface found
        if (ninterface == 0 || decoy.size() == 0)
            return 0.0;

        // calculate tms            
        double Ln = ninterface;
        double d0 = 1.24 * pow(Ln - 15, (double) (1.0 / 3.0)) - 1.8;
        double tms = 0;
        for(int i = 0; i < decoy.size(); i++)
            tms += 1 / (1 + pow(decoy[i].dist(inter[i]) / d0, 2));
        return tms / Ln;
    }
    
    // apply transformation
    void align (SpecMolecule* target, SpecMolecule* mol, bool reverse = false)
	{
        // check
        if (r.size() == 0 || t.size() == 0)
            return;
                    
        // reverse transformation
        if (reverse)
            Trans::reverse(r, t);
            
        // apply transformation
        Trans::apply (target, r, t);
        Trans::apply (mol, r, t);        
    }      
};

