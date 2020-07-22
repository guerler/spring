// interface rmsd
class IRmsdInterface
{
	// reference
    SpecMolecule *reftarget, *refmol;	

	// interface
	Vec < Vec3 > inter, decoy;
public:   	
	// interface residues
	int ninterface;
	int nmol, ntarget;

	// interface flags
	Vec < bool > molalpha, targetalpha;
	
	// parameters
	Vec < Vec3 > r;
	Vec3 t;
    		
	// obtain interface
	bool construct (SpecMolecule* target, SpecMolecule* mol, double cut = 10.0)
	{	
        // reference
        this->reftarget = target;
        this->refmol = mol;
        
		// prepare flags
		targetalpha.resize(target->lcalpha.size());
		targetalpha.fill(false);
		molalpha.resize(mol->lcalpha.size());
		molalpha.fill(false);
				
		// determine interface residues
		for (int i = 0; i < target->latom.size(); i++)
		{
			SpecAtom* a = &target->latom[i];
			for (int j = 0; j < mol->latom.size(); j++)
			{
				SpecAtom* b = &mol->latom[j];				
				if (a->pos.dist(b->pos) < cut && a->pos != NOCOORD && b->pos != NOCOORD)
					targetalpha[a->resno] = molalpha[b->resno] = true;
			}
		}		
		
		// count interface residues
		ntarget = nmol = 0;
		for (int i = 0; i < targetalpha.size(); i++)
            if (targetalpha[i])
                ntarget++;
		for (int i = 0; i < molalpha.size(); i++)
            if (molalpha[i])
                nmol++;

        // target
        ninterface = ntarget + nmol;

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
            Msg::error ("IRmsdInterface::initialize()", "Size of reference and target differs.");
            return false;
        }
				
		// reset
		inter.clear();
		decoy.clear();
		
		// collect assigned interface residues
		for (int i = 0; i < target->lcalpha.size(); i++)
		{
			if (targetalpha[i] && target->latom[target->lcalpha[i]].pos != NOCOORD)
			{
				inter.push_back(reftarget->latom[reftarget->lcalpha[i]].pos);
				decoy.push_back(target->latom[target->lcalpha[i]].pos);				
			}
		}
		for (int i = 0; i < mol->lcalpha.size(); i++)
		{
			if (molalpha[i] && mol->latom[mol->lcalpha[i]].pos != NOCOORD)
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
        if (decoy.size() == 0)
            return LARGE;

		// calculate
		double d = 0;
		for (int i = 0; i < decoy.size(); i++)
   			d += pow (decoy[i].dist(inter[i]), 2);
		
		// return
        return floor (sqrt ( d / decoy.size() ) * 100.0) / 100.0;
	}
	
    // template model score
    double gettms(bool ind = false)
	{
        // no interface found
        int minsize = 20;
        if (ninterface < minsize || decoy.size() == 0)
            return 0.0;

        // calculate tms            
        double Ln = ninterface;
        double d0 = 1.24 * pow(Ln - 15, (double) (1.0 / 3.0)) - 1.8;
        double tms = 0;
        for(int i = 0; i < decoy.size(); i++)
            tms += 1 / (1 + pow(decoy[i].dist(inter[i]) / d0, 2));
            
        // correct length bias and return
        if (ind)
            return (tms / Ln) - (0.0875 - 0.3899 / sqrt(Ln)) + 0.10;
        return (tms / Ln);
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
    
	// interface residues
	static bool reduce (SpecMolecule* target, SpecMolecule* mol, double cutoff = 10.0, bool both = true)
	{
        // interface residues
    	Vec < int > molalpha, targetalpha;

		// prepare flags
		targetalpha.resize(target->lcalpha.size());
		targetalpha.fill(-1);
		molalpha.resize(mol->lcalpha.size());
		molalpha.fill(-1);

		// determine interface residues
		int ncontact = 0;
		for (int i = 0; i < target->latom.size(); i++)
		{
			SpecAtom* a = &target->latom[i];
			for (int j = 0; j < mol->latom.size(); j++)
			{
				SpecAtom* b = &mol->latom[j];
                if (a->pos != NOCOORD && b->pos != NOCOORD)
    				if (a->pos.dist(b->pos) < cutoff)
    				{
	       				targetalpha[a->resno] = molalpha[b->resno] = 0;
	       				ncontact++;
                    }
			}
		}
		if (ncontact == 0)
    		return false;

		// construct interface molecules
		if (both)
		{
            // index residues
            int nres = 0;
    		for (int i = 0; i < targetalpha.size(); i++)
                if (targetalpha[i] != -1)
                    targetalpha[i] = nres++;

    		// construct interface molecules
            SpecMolecule targetint;
            targetint.construct(target->name);

    		// construct interface molecules
            for (int i = 0; i < target->latom.size(); i++)
            {
    			SpecAtom* a = &target->latom[i];
                if (targetalpha[a->resno] != -1 && a->pos != NOCOORD)
                    targetint.append (a->pos, a->name, a->res, targetalpha[a->resno], a->respdb, a->occupancy, a->bfactor);
            }
            if (!targetint.finalize())
                return false;
            target->construct(&targetint);
        }
    
        // index residues
        int nres = 0;
        for (int i = 0; i < molalpha.size(); i++)
            if (molalpha[i] != -1)
                molalpha[i] = nres++;

        // molecule interface 
        SpecMolecule molint;
        molint.construct(mol->name);
        for (int i = 0; i < mol->latom.size(); i++)
        {
            SpecAtom* a = &mol->latom[i];
            if (molalpha[a->resno] != -1 && a->pos != NOCOORD)
                molint.append (a->pos, a->name, a->res, molalpha[a->resno], a->respdb, a->occupancy, a->bfactor);
        }
        if (!molint.finalize())
            return false;

        // change
        mol->construct(&molint);

        // validate
        return true;
    }
};

