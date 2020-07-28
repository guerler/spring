// interface rmsd
struct FNat
{
    // cut
    double cut;
    
	// interface residues
	double ncontact;

	// interface flags
	Vec < int >  molcontact, targetcontact;
	
    // reference
    SpecMolecule* reftarget;
    SpecMolecule* refmol;
	
	// obtain interface
	bool construct (SpecMolecule* target, SpecMolecule* mol, double cut = 10.0)
	{	
        // cut off
        this->cut = cut;
        reftarget = target;
        refmol = mol;        
            
        // calculate contacts
        targetcontact.clear();        
        molcontact.clear();
        
		// determine interface residues
		for (int i = 0; i < target->lcalpha.size(); i++)
		{
			SpecAtom* a = &target->latom[target->lcalpha[i]];
			for (int j = 0; j < mol->lcalpha.size(); j++)
			{
				SpecAtom* b = &mol->latom[mol->lcalpha[j]];				
				if (a->pos.dist(b->pos) < cut && a->pos != NOCOORD && b->pos != NOCOORD)
				{
                    targetcontact.push_back(i);
                    molcontact.push_back(j);                    
                }
			}
		}

        // number of contacts
        ncontact = (double) targetcontact.size();

        // verify        
        if (ncontact == 0)
            return false;
        
        // return
        return true;
	}

    // fraction of native interfacial contacts
    double get(SpecMolecule* target, SpecMolecule* mol)
	{
        // no contacts found
        if (ncontact == 0.0)
            return 0.0;
            
        // verify size
        if (reftarget->lcalpha.size() != target->lcalpha.size() || refmol->lcalpha.size() != mol->lcalpha.size())
        {
            Msg::error ("FNat::get()", "Size of reference and target differs.");
            return 0.0;
        }
        
		// determine interface residues
		double cx = 0;
		for (int i = 0; i < ncontact; i++)
		{
            // count contacts
			SpecAtom* a = &target->latom[target->lcalpha[targetcontact[i]]];
			SpecAtom* b = &mol->latom[mol->lcalpha[molcontact[i]]];			
			if (a->pos.dist(b->pos) < cut && a->pos != NOCOORD && b->pos != NOCOORD)
                cx++;
		}

        // return fraction of conserved contacts
        return cx / ncontact;
    }       
};

