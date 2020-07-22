// generates complex models
struct StructuralNeighbor
{   
    // hash of all cluster centers
    Hash < string, int > hindex;

    // clusters
    Vec < Vec < string > > lcluster;

    // path
    string pdbpath, modelpath;

    // construct
    void construct (string modelpath, string pdbpath)
    {
        // set up
        this->pdbpath = pdbpath;
        this->modelpath = modelpath;
        
        // reset
        hindex.clear();
        lcluster.clear();

        // cluster counter
        int ncluster = 0;

        // read in clusters
        File fc;
        fc.open(pdbpath + "/neighbors/cluster");
        while (fc.read())
        {
            // get number of columns
            int ncol = fc.size();
            if (ncol == 0)
                continue;

            // collect cluster entries
            Vec < string > lcolumns;
            
            // loop through columns
            for (int i = 0; i < fc.size(); i++)
            {
                string id = fc.get(i);
                hindex.insert(id, ncluster);
                lcolumns.push_back(id);
            }

            // add to cluster
            lcluster.push_back(lcolumns);

            // next cluster
            ncluster++;
        }
        fc.close();
        
        // log
        Msg::write ("Number of recognized structural clusters is %i for %i proteins.", ncluster, hindex.size());
    }

    // reset    
    bool get (string chain, string toptemplate, Vec < string >& lchains, double seqcut = LARGE, double thresh = 0.40)
    {
        // exclude
        Exclude exclude;
        exclude.construct(modelpath + "exclusion/identity/" + chain.substr(0, 2) + "/" + chain);

        // make score list for sorting
        Vec < double > lscore;

        // prepare list of neighbors
        lchains.clear();
        
        // get cluster center
        int index = hindex.get(toptemplate);
        if (index != -1)
        {
            // the first entry is the cluster center
            string clustername = lcluster[index][0];
            
            // open file
            string filename = pdbpath + "/neighbors/tmalign/" + clustername.substr(0, 2) + "/" + clustername + ".tmalign";
            File fn;
            fn.open(filename);
            while (fn.read())
            {
                // check
                if (fn.size() < 2)
                    continue;
                    
                // get content
                string nchain = fn.get(0);
                double nscore = fn.getDbl(1);
                
                // check threshold
                if (nscore < thresh)
                    continue;
                    
                // get cluster center for nchain
                int nindex = hindex.get(nchain);
                if (nindex == -1)
                    continue;
                    
                // add to neighbor list
                for (int i = 0; i < lcluster[nindex].size(); i++)
                {
                    // get template
                    string templatename = lcluster[nindex][i];
                    
                    // verify
	       			if (seqcut < 1.0)
		      			if (templatename.substr(0, 4) == chain.substr(0, 4) || exclude.get(templatename) > seqcut)
		  	       			continue;

                    // check
                    lchains.push_back(templatename);
                    lscore.push_back(nscore);
                }
            }
        }
        
        // check
        if (lchains.size() == 0)
            return false;

        // sort
        Lib::sort_greater(lscore, lchains);

        // return
        return true;
    }        
};
