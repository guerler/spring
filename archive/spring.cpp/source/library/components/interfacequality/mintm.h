// interface rmsd
class MinTM
{
	// reference
    SpecMolecule mref, mtmp;	

	// interface residues
	int ninterface;
	int nmol, ntarget;

	// interface flags
	Vec < bool > molalpha, targetalpha;

public:   	
	// obtain interface
	bool construct (SpecMolecule* target, SpecMolecule* mol, double cut = 10.0)
	{	
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
        if (ntarget == 0 || nmol == 0)
            return false;
        
        // generate merged molecule
        merge(target, targetalpha, mol, molalpha, &mref);
        
        // return
        return true;
	}
	
	// get tm-score
	double get (SpecMolecule* target, SpecMolecule* mol)
	{
		// check if an interface is available
        if (ntarget == 0 || nmol == 0)
			return 0.0;
					
	    // generate merged molecule
        merge(target, targetalpha, mol, molalpha, &mtmp);
        
        // perform alignment
        TMScore::align(&mtmp, &mref, true);        
        
        // get score
        return min(tmscore(0, ntarget), tmscore(ntarget, ntarget + nmol));
	}	

private:
	// get tm-score
	double tmscore(int nstart, int nend)
	{
        // get diff
        double Ln = nend - nstart;
        
        // length check
        if (Ln == 0)
        	return 0.0;
        	
        // calculate score
        double d0 = 1.24 * pow(Ln - 15, (double) (1.0 / 3.0)) - 1.8;
        double score = 0.0;
		for(int i = nstart; i < nend; i++)
		{
			// get coordinates
			Vec3 ref = mref.latom[mref.lcalpha[i]].pos;
			Vec3 tmp = mtmp.latom[mtmp.lcalpha[i]].pos;			
			
			// check if both are valid
            if (ref != NOCOORD && tmp != NOCOORD)
                score += 1 / (1 + pow(ref.dist(tmp) / d0, 2));
        }
		return score / Ln;
	}        
    
    // copy residues
	void merge (SpecMolecule* a, Vec <bool>& aalpha, SpecMolecule* b, Vec <bool>& balpha, SpecMolecule* m)
	{
        // get size
        int na = a->lcalpha.size();
        int nb = b->lcalpha.size();        
		int nc = 0;

		// append as atom
		m->construct();
		for (int i = 0; i < na; i++)
		{
			if (aalpha[i])
			{
	            SpecAtom* atom = &a->latom[a->lcalpha[i]];
    	        m->append (atom->pos, atom->name, atom->res, nc++);
			}
        }
		for (int i = 0; i < nb; i++)
		{
			if (balpha[i])
			{			
	            SpecAtom* atom = &b->latom[b->lcalpha[i]];
    	        m->append (atom->pos, atom->name, atom->res, nc++);
			}
        }
        m->finalize();       
    }
};

