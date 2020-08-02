// read hhsearch result files
struct HHsearch
{
    // index
    Vec < double > lscore;
    
    // index
    Vec < string > ltemplate;
        
    // filename
    string fname, chain, modelpath;
    
    // hhsearch internal chain name
    string query;

    // storage
    Storage store;
        
    // construct on the basis of tmalign table
    bool construct(string chain, string modelpath, double minscore = -LARGE, double seqcut = LARGE)
    {
        // names
        this->modelpath = modelpath;        
        this->fname     = chain + ".hhr";
        this->chain     = chain;
    
        // return    
        return construct(chain, modelpath, fname, minscore, seqcut);
    }                

    // construct on the basis of tmalign table
    bool construct(string chain, string modelpath, string fname, double minscore = -LARGE, double seqcut = LARGE)
    {
        // fname
        fname           = modelpath + fname;
        
        // set
        this->modelpath = modelpath;
        this->fname     = fname;
        this->chain     = chain;
        
        // file
        File f;
            
        // open file
        if (!f.open (fname))
        {
            Msg::write ("WARNING: HHsearch file %s not found!", fname.c_str());
            return false;
        }

        // hash
        Hash < string, int > lhash;
        
        // reset
        query = "";
        lscore.clear();
        ltemplate.clear();
        int n = 0;
        while (f.good())
        {
            if (!f.read())
                continue;

            // get query name
            if (f.get(0) == "Query")
                query = f.get(1);
            
            // skip comments
    		if (f.get(0)[0] == '>')
    		{
                // get name         
                string templatename = Lib::strclean(f.get(0).substr(1));

                // verify name
                if (templatename == "")
                    continue;
                
                // next line
                double score = 0.0;
                double identity = 0.0;
                if (f.read())
                {
                    // get score
                    if (f.get(1).substr(0, 8) == "E-value=")
                    {
                        score = Convert::toDbl(f.get(1).substr(8));
                        if (score > 0.0)
                            score = -log10(score);
                    }
                    
                    // get identity
                    if (f.get(4).substr(0, 11) == "Identities=")
                    {
                        string tmpstr = f.get(4).substr(0, f.get(4).length() - 1);
                        identity = Convert::toDbl(tmpstr.substr(11)) / 100.0;
                    } else {
                        Msg::error ("HHsearch::construct()", "Identity not found.");
                    }
                }

                // verify
                if (score > minscore)
                {
                    // score
                    if (lhash.get(templatename) == -1)
                    {
                        // insert
                        lhash.insert(templatename, 1);
                        ltemplate.push_back(templatename);
                        lscore.push_back(score);
                        
                        // count
                        n++;
                    }
                }
            }
        }

        // close file
        f.close();
        
        // verify
        if (ltemplate.size() == 0 || query == "")
        {
            Msg::write ("No template or query definition found in hhr (%s).", fname.c_str());
            return false;
        }

        // sort
        Lib::sort_greater(lscore, ltemplate);
        
        // return
        return true;	
    }       

    // exclude
    bool getmodel (SpecMolecule* mol, string templatename, string templatepath)
    {       
        // get template
        SpecMolecule moltmpl;
        string fname = templatepath + "chains/" + templatename.substr(0, 2) + "/" + templatename;
        if(!store.read(&moltmpl, fname, templatename))
            Msg::write ("Template file %s not found.", fname.c_str());
        else
            moltmpl.lg();
            
        // getmodel
        return getmodel (mol, templatename, &moltmpl);
    }
        
    // exclude
    bool getmodel (SpecMolecule* mol, string templatename, SpecMolecule* moltmpl)
    {
        // temporary
        SpecMolecule model;

        // generate empty molecule
        model.construct(templatename);
        Msg::write("Loading template %s.", templatename.c_str());
        // open sequence alignment file
        File f;
        if (!f.open(fname))
            Msg::write("HHsearch::getmodel()","HHsearch file not found %s.", fname.c_str());

        // steps
        int step = 0;
        while (f.good())
        {
            // read line
            if (!f.read())
                continue;

            if (f.get(0) == ">" + templatename)
            {
                step = 1;
                break;
            }
        }
        
        // verify
        if (step == 0)
        {
            // close file
            f.close();
            
            // message
            Msg::write ("Sequence alignment of %s not found!", templatename.c_str());
            return false;
        }
        
        // read sequence alignment
        int querystart = 0, tmplstart = 0;
        string queryseq, tmplseq;        
        while (f.good())
        {
            // read
            if (!f.read())
                continue;
            
            // next code
            string key = f.get(0);
            char k = key[0];
            if (k == '>')
                break;

            // get target offset
            if (step == 1 && k == 'Q' && f.get(1) == query)
            {
                querystart = f.getInt(2);
                queryseq = f.get(3);
                step = 2;
            }
            
            // copy from template
            if (step == 2 && k == 'T' && f.get(1) == templatename)
            {
                tmplstart = f.getInt(2);
                tmplseq = f.get(3);
                step = 1;

                // verify
                if (tmplseq.size() != queryseq.size())
                {
                    Msg::write ("Sizes differ! %i %i.", (int) tmplseq.size(), (int) queryseq.size());
                    Msg::error ("HHsearch::thread()","Sizes differ.");
                }

                // copy assignment
                int tmplindex = 0;
                for (int i = 0; i < tmplseq.size(); i++)
                {
                    if (queryseq[i] != '-' && tmplseq[i] != '-' && queryseq[i] != '.' && tmplseq[i] != '.')
                    {
                        // convert indices                    
                        int tmplres = tmplstart + tmplindex - 1;

                        // template residues
                        if (tmplres < moltmpl->lcalpha.size())
                        {
                            // get atom
                            SpecAtom* atom = &moltmpl->latom[moltmpl->lcalpha[tmplres]];
    
                            // check
                            if (atom->res != SpecDetails::getResName(tmplseq[i]))
								Msg::write("HHsearch::getmodel() : Template sequence inconsistent : %s, %s, %i", atom->res.c_str(), SpecDetails::getResName(tmplseq[i]).c_str(), i);
                        
                            // name
                            string resname = SpecDetails::getResName(queryseq[i]);
                           
                            // setup atom
                            model.append (atom->pos, "CA ", resname, model.size());                        
                        } else {
                            Msg::write("HHsearch::getmodel() : Invalid template index.");
                        }
                    }
                    if (tmplseq[i] != '-' && tmplseq[i] != '.')
                        tmplindex++;                        
                }
            }
        }
                
        // close file      
        f.close();

        // finalize
        model.finalize();

        // copy molecules back
        mol->construct(&model);	

        // return
        return true;        
    }
};


