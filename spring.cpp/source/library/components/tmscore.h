// interface isscore
struct TMScore
{
	// obtain alignment and score
	static double align (SpecMolecule* mol, SpecMolecule* tmpl, bool transform = false)
	{
		// get score
        TransMatrix tmat;
        double score = align (mol, tmpl, tmat);
        
        // orient
        if (transform)
			tmat.apply(mol);
		
		// return
        return score;
    }

	// obtain alignment and score
	static double align (SpecMolecule* mol, SpecMolecule* tmpl, TransMatrix& tmat)
	{
        // verify
        if (min(mol->lcalpha.size(), tmpl->lcalpha.size()) <= 5)
            return 0.0;	
                           
        // pars
        Storage store;
        string uuid;
        string workpath;
        string unique;
    
        // obtain
        Config cnf;
        workpath = cnf.ptTemp;        
        uuid = cnf.get_uuid();

        // unqiue filename combo        
        unique = workpath + uuid + ".";        

        // p-value
        double score = 0.0;
        
        // overwrite possible enumeration of pdb file
        // this is important particularly for tmscore,
        // since it relies on the right enumeration
        mol->enumerateResidues();
        tmpl->enumerateResidues();
        
        // save
        store.save (mol, unique + "mol.pdb");
        store.save (tmpl, unique + "tmpl.pdb");

        // calc        
        string execute = "./plugins/tmscore/tmscore ";
        execute += unique + "mol.pdb " + unique + "tmpl.pdb > " + unique + "out.txt";
        system (execute.c_str());
        
        // get
        File f;
        f.open(unique + "out.txt");

        // aligned segment
        while (f.good())
        {
            f.readLine();
            
            // get score
            if (f.get(0, 13) == "TM-score    =")
            {
                score = Convert::toDbl(f.get(14, 6));
                break;
            }
        }

        // verify
        if (score == 0)
            return 0.0;

        // while
        while (f.good())
        {
            f.readLine();

            // read out transformation matrix
            if (f.get(0, 10) == " -------- ")
            {
                // read next line
                f.readLine();
                break;
            }
        }
    
        // get rotation matrix
        int i = 0;
        while (i < 3 && f.good())
        {
            // readlin
            f.readLine();

            // matrix
            tmat.t[i] = Convert::toDbl(f.get(3, 17));
            tmat.r[i][0] = Convert::toDbl(f.get(21, 14));
            tmat.r[i][1] = Convert::toDbl(f.get(36, 14));
            tmat.r[i][2] = Convert::toDbl(f.get(51, 14));

            // increase counter
            i++;
        }
                    
        // clean
        system (string("rm " + unique + "*").c_str());
 
        // return
        return score;
    }
};

