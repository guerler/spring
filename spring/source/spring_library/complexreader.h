// generates complex models
struct ComplexReader
{   
    // list of all molecules
    Vec < Complex > lcomplex;
    
    // cross reference
    CrossReference reference;

    // molecules
    SpecMolecule *target;
    SpecMolecule *mol;

    // path
    string pdbpath;
    
    // storage
    Storage store;
    
    // size
    int n_mol, n_target;
    
    // score
	Vec < Vec < int > > subst;
    
    // construct
    void construct (string pdbpath, string matrix = "data/matblosum.txt")
    {
        // format
		Format::readMatrix(matrix, subst);        

        // link
        this->pdbpath = pdbpath;
        
        // load cross
        reference.construct(pdbpath + "index.txt");
    }

    // reset    
    void initialize(SpecMolecule* target, SpecMolecule* mol)
    {
        // link
        this->target = target;
        this->mol    = mol;

        // calphas
        n_target     = target->lcalpha.size();        
        n_mol        = mol->lcalpha.size();
        
        // reset        
        lcomplex.clear();
    }

    // read model from coordinates file
    bool add (string corename, double corescore, Vec <string>& lpartnername, Vec <double>& lpartnerscore, bool swap = false)
    {
        // log
        Msg::rewrite ("Template %s.", corename.c_str());

        // read directory
        Vec < CrossTemplate > ltemplate;
        if (!reference.get (corename, ltemplate))
            return false;

        // top complex
        Complex c;

        // make a hash
        Hash < string, double > hindex;

        // generate hash list for partner names
        for (int i = 0; i < lpartnername.size(); i++)
            hindex.insert(lpartnername[i], lpartnerscore[i]);

        // align and generate complexes
        int n = SpringConfig::NTMPL_PER_PDB;
        for (int i = 0; i < ltemplate.size(); i++)
        for (int j = 0; j < ltemplate[i].lcore.size(); j++)
        for (int k = 0; k < ltemplate[i].lpartner.size(); k++)
        {
            // skip
            if (ltemplate[i].lcore[j] == ltemplate[i].lpartner[k])
                continue;

            // partner score
            double partnerscore = hindex.get(ltemplate[i].lmatch[k]);

            // check
            if (partnerscore == -1)
                continue;

            // counter
            if (n-- < 0)
            	goto stop;

            // set score
            double zscore = min(corescore, partnerscore);

            // make name
            string core = ltemplate[i].lcore[j];
            string partner = ltemplate[i].lpartner[k];

            // update best complex
            c.corename = corename;
            c.core     = core;
            c.partner  = partner;
            c.zscore   = zscore;
            c.swap     = swap;

            // add to list
            lcomplex.push_back(c);
        }

        // stop loop
        stop:

        // return score
        return true;
    }

    // compare by score
	static bool sortbyzscore(const Complex& a, const Complex& b)
	{
		return a.zscore > b.zscore;
	}

    // compare by score
	static bool sortbyspringscore(const Complex& a, const Complex& b)
	{
		return a.score > b.score;
	}
	
    // read model from coordinates file
    void makemodels()
    {
        // sort
        std::sort (lcomplex.begin(), lcomplex.end(), sortbyzscore);

        // maximum number of templates
        int n = 0;
        for (int i = 0; i < lcomplex.size(); i++)
        {
            // top complex
            Complex* c = &lcomplex[i];

            // load templates
            SpecMolecule targettmp, moltmp;
        
            // generate file and directory names
            string splitpath 	= pdbpath + "splits/" + c->corename.substr(0, 2) + "/";
            string corepdb 		= splitpath + c->corename + "/" + c->core + ".pdb";
            string partnerpdb 	= splitpath + c->corename + "/" + c->partner + ".pdb";

            // switch between target = core and mol = core
            if (!c->swap)
            {
                store.read(&targettmp, corepdb, c->core);
                store.read(&moltmp, partnerpdb, c->partner);
            } else {
                store.read(&targettmp, partnerpdb, c->partner);
                store.read(&moltmp, corepdb, c->core);
            }
            	
			// check
			if (targettmp.lcalpha.size() == 0 || moltmp.lcalpha.size() == 0)
			{
				//Msg::write ("Requested split of %s not found.", c->core.c_str());
				continue;
			}

            // reduce
            if(!IRmsdInterface::reduce(&targettmp, &moltmp))
			{
				//Msg::write ("Requested split has not interface.", c->core.c_str());
				continue;
			}

            // alignment
            Vec < int > targetaln, molaln;
            double tms = min(TMAlign::align(target, &targettmp, targetaln), TMAlign::align(mol, &moltmp, molaln));
            if (tms == 0.0)
                continue;

            // get global score
            EnergyContacts e;
            e.construct(&targettmp, &moltmp);
            double energy = e.get(target, mol, targetaln, molaln);
            if (energy == 0.0)
                continue;

            // count
            if (n++ >= SpringConfig::NTMPL)
                break;

            // check for clashes
            double clashes = getclashes(target, mol);
                    
            // update score
            double score = (c->zscore + SpringConfig::W0 * tms + SpringConfig::W1 * energy) / SpringConfig::NORMALIZE;

            // sum aligned residues
            int aln_total = 0;
            int aln_same  = 0;
            
            // copy molecules
            c->target.construct(target);
            c->mol.construct(mol);
                    
            // get details
            c->info =  getalignment ("query_a", &targettmp, target, targetaln, aln_total, aln_same);
            c->info += getalignment ("query_b",&moltmp, mol, molaln, aln_total, aln_same);
                    
            // get energy / tms
            c->aln_total    = aln_total;
            c->aln_same     = aln_same;
            c->energy       = energy;
            c->tms          = tms;
            c->clashes      = clashes;
            c->score        = score;
            
            // log
            Msg::rewrite ("Template %i: %s with %s, %s (score=%4.2f, minz=%4.2f).", n, c->corename.c_str(), c->core.c_str(), c->partner.c_str(), c->score, c->zscore);
        }
        
        // log
        Msg::rewrite("");
        Msg::write("Model generation completed.");
        
        // sort
        std::sort (lcomplex.begin(), lcomplex.end(), sortbyspringscore);
    }

    // reweight sasa
    void updatevdw (SpecMolecule* mol, double sasa_core_treshold = 0.1, double sasa_low = 0.3, double sasa_high = 1.3)
    {
        // make pulchra
        Pulchra::make(mol);
        
        // get solvent
        NAccess::make(mol);
        
        // reduce to backbone at solvent accessible residues
        int nchange = 0;
        for (int i = 0; i < mol->latom.size(); i++)
        {
            // core
            if(mol->latom[i].occupancy > sasa_core_treshold)
            {
                mol->latom[i].bfactor *= sasa_low;
                nchange++;
            } else
            {
                mol->latom[i].bfactor *= sasa_high;
                nchange++;
            }
        }
        Msg::write ("Updated vdw-radii for %i atoms.", nchange);
        
        // reduce to backbone
        SpecMolecule tempmol;
        mol->copyBone(&tempmol);
        mol->construct(&tempmol);
    }

    // refine models
    // HEY BRANDON - THE DEBUG STATEMENT PUTS OUT THE TRAJECTORIES...ALWAYS REMOVE WHEN SEND TO CLUSTER
    void refinemodels(bool debug = false)
    {
        // refining models
        Msg::write("Refining models...");
        
        // copy/make side
        SpecMolecule cptarget, cpmol;
        cptarget.construct(target);
        cpmol.construct(mol);
        
        // update vdw-radii
        updatevdw(&cpmol);
        updatevdw(&cptarget);

        // refine/unclash
        int nmodels = 0;
        for (int i = 0; i < lcomplex.size(); i++)
        {
            // top complex
            Complex* c = &lcomplex[i];
            
            // check
            if (c->score == -LARGE)
                continue;
            
            // counter
            nmodels++;
            
            // align full-chain copies to target/mol
            Trans::align(&cptarget, &c->target);
            Trans::align(&cpmol, &c->mol);
            
            // backup initial orientation
            SpecMolecule cpstart;
            cpstart.construct(&cpmol);
            
            // check for clashes
            if (!debug)
            {
                // optimize binding mode
                Refinement::minimize(&cptarget, &cpmol);

                // update deviation
                c->deviation = cpstart.rmsd(&cpmol);

                // log
                Msg::rewrite ("Refined model %i with a deviation of %4.2f...", nmodels, c->deviation);
            } else {
                double before_clash = getclashes(&c->target, &c->mol);
                
                // filename
                string fname = c->target.name + "_" +  c->mol.name + "_" + Convert::toString(nmodels) + ".pdb";
                Vec < SpecMolecule > pdb;
                pdb.push_back(cptarget);
                pdb.push_back(cpmol);

                // refine
                Refinement::minimize(&cptarget, &cpmol, &pdb);

                // update deviation
                c->deviation = cpstart.rmsd(&cpmol);
                
                // save trajectory
                store.save (pdb, "results/refinement" + fname, 2);
                
                // get after values
                double after_clash = getclashes(&c->target, &cpmol);
                
                // log
                Msg::rewrite ("Template %i: Clashes before %4.2f and after %4.2f.", nmodels, before_clash, after_clash);
            }

            // update spring score
            c->score += SpringConfig::W2 * c->deviation / SpringConfig::NORMALIZE;

            // copy refined result back
            c->mol.construct(&cpmol);
        }
        
        // sort
        std::sort (lcomplex.begin(), lcomplex.end(), sortbyspringscore);
        
        // model refinement completed
        Msg::rewrite("");
        Msg::write("Model refinment completed.");
    }

    // calculate clashes
    static double getclashes (SpecMolecule* target, SpecMolecule* mol)
    {
        double clashes = 0.0;
        for (int m = 0; m < target->lcalpha.size(); m++)
        for (int n = 0; n < mol->lcalpha.size(); n++)
            if (target->latom[target->lcalpha[m]].pos != NOCOORD && mol->latom[mol->lcalpha[n]].pos != NOCOORD)
            {
                if (target->latom[target->lcalpha[m]].pos.dist2(mol->latom[mol->lcalpha[n]].pos) < SpringConfig::RADIUS2)
                    clashes++;
            }
        clashes /= min (target->lcalpha.size(), mol->lcalpha.size());
        return clashes;
    }
    
    // get alignment
    string getalignment (string name, SpecMolecule* tmpl, SpecMolecule* mol, Vec<int>& aln, int& aln_total, int& aln_same)    
    {   
        // result 
        string summary = "";
            
        // alignment
        int start     = 0;
        int block     = 70;
        int n_aligned = 0;
        int n_same    = 0;        
        while (start < aln.size())
        {
            // alignment
            string tmpl_sequence      = "";
            string target_sequence    = "";
            string alignment_sequence = "";   
                    
            // get end
            int end = min(start + block, aln.size());
            
            // alignment
            for (int i = start; i < end; i++)
                tmpl_sequence += SpecDetails::getCharCode(tmpl->latom[tmpl->lcalpha[i]].rescode);
        
            // alignment        
            int minindex = LARGE, maxindex = 0;
            for (int i = start; i < end; i++)
            {
                if (aln[i] != LARGE)
                {
                    // get alignment
                    int tmpl_code = tmpl->latom[tmpl->lcalpha[i]].rescode;
                    int mol_code  = mol->latom[mol->lcalpha[aln[i]]].rescode;
    
                    // get target sequence
                    target_sequence += SpecDetails::getCharCode(mol_code);
    
                    // aligned
                    n_aligned++;
    
                    // get alignment                
                    if (tmpl_code == mol_code)
                    {
                        alignment_sequence += ":";
                        n_same++;
                    } else {
                        if (subst[tmpl_code][mol_code] > 0)
                            alignment_sequence += ".";
                        else
                            alignment_sequence += "_";
                    }
                    
                    // get index
                    minindex = min (minindex, aln[i]);
                    maxindex = max (maxindex, aln[i]);
                } else {
                    target_sequence    += "-";
                    alignment_sequence += "_";
                }
            }
            
            // format
            if (minindex != LARGE && maxindex >= minindex)
            {
                sprintf (strbuf, "Q %-15s %-5i %-70s %-5i", name.c_str(), minindex + 1, target_sequence.c_str(), maxindex + 1);
                summary += string(strbuf) + "\n";
                sprintf (strbuf, "I %-15s %-5i %-70s %-5i", "similarity", start + 1, alignment_sequence.c_str(), end);
                summary += string(strbuf) + "\n";
                sprintf (strbuf, "T %-15s %-5i %-70s %-5i", tmpl->name.c_str(), start + 1, tmpl_sequence.c_str(), end);
                summary += string(strbuf) + "\n\n";
            }
                    
            // update
            start += end;            
        }
        
        // update
        aln_total += n_aligned;
        aln_same  += n_same;
                
        // return
        string info = " Aligned=(" + Convert::toString(n_aligned) + "/" + Convert::toString((int)aln.size()) + ")";
        info += " Identities=" + Convert::toString(Lib::percent(n_same, n_aligned)) + "%";        
        return ">" + tmpl->name + info + "\n" + summary;
    }	
};
