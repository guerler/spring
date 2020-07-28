// exclusion for tmscore similarity
struct Exclude
{
    // exlcusion hash
    Hash < string, double > hlist;

    // file
    File f;

    // load file
    bool construct(string fname)
    {
        // reset
        hlist.clear();

        // open file
        if (!f.open (fname))
        {
            //Msg::write ("WARNING: Exclusion file %s not found!", fname.c_str());
            return false;
        }
        
        // reset
        while (f.read())
        {
            // skip comments
    		if (f.get(0) == "#")
	       		continue;

            // string
            string key = f.get(0);
            double score = f.getDbl(1);
                    
            // list
            if (key != "")
                hlist.insert(key, score);
        }
        
        // close file
        f.close();
        return true;	
    }
    
    // exclude
    double get (string search)
    {
        return hlist.get(search);
    }
};


