// complex def
struct CrossTemplate
{
    // list of all molecules
    Vec < string > lcore, lpartner, lmatch;
    
    // construct
    CrossTemplate()
    {
        // reset
        lcore.clear();
        lpartner.clear();
        lmatch.clear();
    }
};

// cross referencing
struct CrossReference
{
    // list of all lists
    Vec < Vec < CrossTemplate > > list;
    
    // hash index of available lists
    Hash < string, int > index;
    
    // construct
    void construct (string filename)
    {
        // reset
        list.clear();
        index.clear();

		// library size
		int nlib = 0;

        // cross reference
        File f;
        f.open (filename);
        Vec < string > lcompound, lmatch;
        string newbio = "", oldbio = "";
        string newcore = "", oldcore = "";        
        Vec < CrossTemplate > ltemplate;      
        while (f.read())
        {
            // get core name
            newcore = f.get(0);

            // first core
            if (oldcore == "")
                oldcore = newcore;

            // check for new core
            if (oldcore != newcore)
            {
                // add final biomolecule of core to list of templates
                addtolist (lcompound, lmatch, ltemplate);
								
                // add new core
                if (ltemplate.size() > 0)
                {
					// check if empty
					if (ltemplate[0].lcore.size() == 0)
					{
						cout << newcore << endl;
						Msg::error("CrossReference::construct()", "No core found.");
					}

					// count
					for (int i = 0; i < ltemplate.size(); i++)
						nlib += ltemplate[i].lcore.size() * ltemplate[i].lpartner.size();

                    // add
                    list.push_back(ltemplate);
                    index.insert(oldcore, (int) index.size());
                    
                    // reset
                    ltemplate.clear();
                    lcompound.clear();
                    lmatch.clear();             
                }
                    
                // reset name
                oldcore = newcore; 
                oldbio  = "";                   
            }
                            
            // get biomol id
            Vec < string > tagarr; 
            string tagstr = f.get(1);
            Lib::split(tagarr, &tagstr[0], '_');
            if (tagarr.size() != 3)
                continue;
            
            // get biomolecule identifier
            newbio = newcore + tagarr[0];
            
            // first biomolecule
            if (oldbio == "")
                oldbio = newbio;
            
            // check for change
            if (oldbio != newbio)
            {
                // add to list of templates
                addtolist (lcompound, lmatch, ltemplate);
                
                // reset
                lcompound.clear();
                lmatch.clear();
                oldbio = newbio;
            }

            // backup
            lcompound.push_back(f.get(1));           
            lmatch.push_back(f.get(2));           
        }
        f.close();
        
        // check
        if (list.size() == 0)
            Msg::error ("CrossReference::construct()", "Index file not available or empty.");
        else
        	Msg::write ("Estimated library size %i.", nlib);
    }
    
    // add
    void addtolist(Vec < string > lcompound, Vec < string > lmatch, Vec < CrossTemplate >& ltemplate)
    {       
        // new cross template definition
        CrossTemplate tmpl;

        // loop
        Vec < string > tag;
        bool valid = false;
        for (int i = 0; i < lcompound.size(); i++)
        {
            // split name
            Lib::split(tag, &lcompound[i][0], '_');
            
            // collect cores
            if (tag[1] == "0")
            {
                tmpl.lcore.push_back(lcompound[i]);
				valid = true;                
			}

            // collect everything
            tmpl.lpartner.push_back(lcompound[i]);
            tmpl.lmatch.push_back(lmatch[i]);         
        }
        
        // add to list
        if (valid)
	        ltemplate.push_back(tmpl);
    }
    
    // read model from coordinates file
    bool get (string corename, Vec < CrossTemplate >& tmpl)
    {
        // make name
        int i = index.get(corename);
        if (i != -1)
        {
            tmpl = list[i];
            return true;
        }
        // core not found
        return false;
    }	
};
