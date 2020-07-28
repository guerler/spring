// threading info
struct ThreadingInfo
{
    string tmpl;
    double score;
};

// cross referencing
struct ThreadingReference
{
    // list of all lists
    Vec < Vec < ThreadingInfo > > list;
    
    // hash index of available lists
    Hash < string, int > index;
    
    // construct
    void construct (string filename)
    {
        // reset
        list.clear();
        index.clear();

        // reference
        Vec < ThreadingInfo > linfo;
        string newcore = "", oldcore = "";  
        
        // cross reference
        File f;
        f.open (filename);
        while (f.read())
        {
            // get core name
            newcore = f.get(0);

            // first core
            if (oldcore == "")
                oldcore == newcore;

            // check for new core
            if (oldcore != newcore)
            {
                // add new core
                if (linfo.size() > 0)
                {
                    // add
                    list.push_back(linfo);
                    index.insert(oldcore, (int) index.size());
                    
                    // reset
                    linfo.clear();                    
                }
                
                // reset name
                oldcore = newcore; 
            }
            
            // get match
            string tmpl = f.get(1);
            double score = Convert::toDbl(f.get(2));
                        
            // check
            if (tmpl == "")
                continue;
            
            // backup
            ThreadingInfo ti;
            ti.tmpl = tmpl;
            ti.score = score;   
            linfo.push_back(ti);
        }
        f.close();
        
        // check
        if (list.size() == 0)
            Msg::error ("ThreadingReference::construct()", "Index file not available or empty.");
    }
    
    // read model from coordinates file
    bool get (string corename, Vec < ThreadingInfo >& linfo)
    {
        // make name
        int i = index.get(corename);
        if (i != -1)
        {
            linfo = list[i];
            return true;
        }
        // core not found
        return false;
    }	
};
