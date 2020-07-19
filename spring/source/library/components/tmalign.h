// interface isscore
struct TMAlign
{
	// obtain alignment and score
	static double align (SpecMolecule* mol, SpecMolecule* tmpl, int l = 0)
	{
        Vec <int> aln;
        TransMatrix tmat;
        return TMAlign::align (mol, tmpl, aln, tmat, l);
    }

	// obtain alignment and score
	static double align (SpecMolecule* mol, SpecMolecule* tmpl, Vec <int>& aln, int l = 0)
	{
        TransMatrix tmat;
        return TMAlign::align (mol, tmpl, aln, tmat, l);
    }
	
    // obtain alignment and score
	static double align (SpecMolecule* mol, SpecMolecule* tmpl, TransMatrix& tmat, int l = 0)
	{
        Vec <int> aln;
        return TMAlign::align (mol, tmpl, aln, tmat, l);
    }
            
	// obtain alignment and score
	static double align (SpecMolecule* mol, SpecMolecule* tmpl, Vec <int>& aln, TransMatrix& tmat, int l = 0)
	{
        // verify
        if (min(mol->lcalpha.size(), tmpl->lcalpha.size()) <= 5)
            return 0.0;	
                    
        // pars
        Storage store;
        string uuid;
        string workpath;
        string unique;
                    
        // prepare mapping
        aln.resize(tmpl->lcalpha.size());
        aln.fill(LARGE);
    
        // obtain
        Config cnf;
        workpath = cnf.ptTemp;        
        uuid = cnf.get_uuid();

        // unqiue filename combo        
        unique = workpath + uuid + ".";        

        // p-value
        double score = 0.0;
        
        // save
        store.save (mol, unique + "mol.pdb");
        store.save (tmpl, unique + "tmpl.pdb");
        
        // calc        
        string execute = "./plugins/tmalign/tmalign ";
        if (l != 0)
            execute += " -L " + Convert::toString(l) + " ";
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
            if (f.get(0, 7) == "Aligned")
            {
                score = Convert::toDbl(f.get(43, 7));
                break;
            }
        }

        // verify
        if (score == 0)
        {                            
            system (string("rm " + unique + "*").c_str());                    
            return 0.0;
        }

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

        // apply transformation
        if (i == 3)
            tmat.apply(mol);
        else {
            //Msg::write("TMAlign: Transformation matrix invalid.");
                    
            // clean
            system (string("rm " + unique + "*").c_str());                    
            return 0.0;            
        }

        // get alignment sequence
        while (f.good())
        {
            f.readLine();

            // read out transformation matrix
            if (f.get(2, 1) == ":")
            {
                // alignment data
                string targetaln = "";
                string templatealn = "";
                                
                // read next line
                if (f.good())
                    if (f.read())
                        targetaln = f.get(0);
                
                // skip
                f.readLine();
                
                // read next line
                if (f.good())
                    if (f.read())
                        templatealn = f.get(0);

                // verify
                if (templatealn.size() != targetaln.size())
                {
                    // message
                    Msg::write ("TMalign: Misdefined alignment.");
                    
                    // clean
                    system (string("rm " + unique + "*").c_str());
                    return 0.0;
                }
    
                // count residues
                int nres = 0;
                for (int i = 0; i < templatealn.size(); i++)
                    if (templatealn[i] != '-')
                        nres++;
                
                // verify
                if (nres != (int) tmpl->lcalpha.size())
                {
                    // message
                    Msg::write ("TMalign: alignment size and template size do not match. %i, %i.", nres, tmpl->lcalpha.size());
                    
                    // clean
                    system (string("rm " + unique + "*").c_str());
                    return 0.0;
                }

                // get original index
                Vec < int > index;
                index.resize(mol->lcalpha.size());
                index.fill(LARGE);
                int nmol = 0;
                for (int i = 0; i < mol->lcalpha.size(); i++)
                    if (mol->latom[mol->lcalpha[i]].pos != NOCOORD)
                        index[nmol++] = i;
                
                // load alignment
                int targetindex = 0;
                int templateindex = 0;
                for (int i = 0; i < templatealn.size(); i++)
                {
                    if (targetaln[i] != '-' && templatealn[i] != '-')
                        aln[templateindex] = index[targetindex];
                    if (targetaln[i] != '-')
                        targetindex++;
                    if (templatealn[i] != '-')
                        templateindex++;
                }
            }
        }
        
        // clean
        system (string("rm " + unique + "*").c_str());
                        
        // return
        return score;
    }
};

