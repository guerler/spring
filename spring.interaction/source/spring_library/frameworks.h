// cross referencing
struct Act_CrossReference
{
    // list of all lists
    Vec < Vec < string > > list;
    
    // hash index of available lists
    Hash < string, int > index;
    
    // construct
    void construct (string filename)
    {
        // reset
        list.clear();
        index.clear();

        // cross reference
        File f;
        f.open (filename);
        Vec < string > lmatch;
        string newcore = "", oldcore = "";
        string newbio = "", oldbio = "";
        bool skip_core = false;
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
                // add new core
                if (lmatch.size() > 0)
                {
                    // add
                    list.push_back(lmatch);
                    index.insert(oldcore, (int) index.size());
                    
                    // reset
                    lmatch.clear();                    
                }
                
                // reset name
                oldcore = newcore;
                oldbio  = "";
            }
            
            // get match
            string match = f.get(2);
            
            // check
            if (match == "")
                continue;
            
            // is not first core
            Vec < string > tagarr;
            string tagstr = f.get(1);
            Lib::split(tagarr, &tagstr[0], '_');
            if (tagarr.size() != 3)
                continue;

            // get biomolecule identifier
            newbio = newcore + tagarr[0];

            // check for change
            if (oldbio != newbio)
            {
                oldbio = newbio;
                skip_core = true;
            }
            
            // skip first core in every biomolecule
            if (skip_core && tagarr[1] == "0")
            {
                skip_core = false;
                continue;
            }
            
            // is unique
            bool found = false;
            for (int i = 0; i < lmatch.size(); i++)
            {
                if (lmatch[i] == match)
                {
                    found = true;
                    break;
                }
            }
            
            // backup
            if (!found)            
                lmatch.push_back(match);           
        }
        f.close();
        
        // check
        if (list.size() == 0)
            Msg::error ("Act_CrossReference::construct()", "Index file not available or empty.");
    }
    
    // read model from coordinates file
    bool get (string corename, Vec < string >& tmpl)
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


// generates complex models
struct Frameworks
{
    // get top ranks
    double maxscore;
    
    // cross reference
    Act_CrossReference reference;
        
    // construct
    void construct (string filename)
    {
        // load cross
        reference.construct(filename);

        // reset
        initialize();
    }

    // reset    
    void initialize()
    {
        maxscore = -LARGE;
    }        


    // identify frameworks
    void getframeworks (Vec < string > ltarget, Vec < ltarget_score)
    {
    }
    
    // read model from coordinates file
    bool add (ThreadingInfo& core, Vec < ThreadingInfo >& partner)
    {
        // read directory
        Vec < string > ltemplate;
        if (!reference.get (core.tmpl, ltemplate))
            return false;

        // make a hash
        Hash < string, double > hscore;
        
        // generate hash list for identified templates
        for (int i = 0; i < partner.size(); i++)
            hscore.insert(partner[i].tmpl, partner[i].score);

        // align and generate complexes
        for (int i = 0; i < ltemplate.size(); i++)
        {
            // start with the most similar template
            double partnerscore = hscore.get(ltemplate[i]);

            // check
            if (partnerscore == -1)
                continue;

            // score and energy
            double score = min(partnerscore, core.score);

            // backup maximal score
            maxscore = max(maxscore, score);
        }
                
        // return score
        return true;
    }	
};
