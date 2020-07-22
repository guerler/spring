// local library
struct Lib
{
    // make directory
	static void makedir(string d)
	{                
        mkdir(d.c_str(), 0777);
	}
	
    // read quality
    template < typename t_num >
    static int percent(t_num a, t_num b)
    {        
        return int(double(a) * 100.0 / double(b));
    }

    // read quality
    static void getquality(Vec < bool >& lref, Vec < bool >& lpred, double& cov, double& acc, double& mcc)
    {
        // get values
        double tp = 0, fp = 0, tn = 0, fn = 0;
        
        // get interface size
        for (int i = 0; i < lref.size(); i++)
        {
            if (lref[i])
            {
                if (lpred[i])
                    tp++;
                else
                    fn++;
            } else {
                if (lpred[i])
                    fp++;
                else
                    tn++;
            }
        }
        
        // get coverage
        if (tp + fn > 0)
            cov = tp / (tp + fn);
        else
            cov = 0.0;

        // get accuracy
        if (tp + fp > 0)
            acc = tp / (tp + fp);
        else
            acc = 0.0;

        // get mathew correlation coefficiant
        mcc = (tp+fn)*(tp+fp)*(tn+fp)*(tn+fn);
        if (mcc > 0)
            mcc = (tp*tn - fp*fn) / sqrt(mcc);
        else
            mcc = 0.0;
    }

	// normalize deviation
	static bool normalizemax(Vec <double>& list)
	{
		// normalize
		if (list.size() == 0)
			return false;

		// minimum
		double mi = LARGE;
		for (int i = 0; i < list.size(); i++)
            mi = min(mi, list[i]);
        
        // apply
		for (int i = 0; i < list.size(); i++)
			list[i] -= mi;

		// maximum
		double ma = -LARGE;
		for (int i = 0; i < list.size(); i++)
			ma = max(ma, list[i]);

		// apply
		if (ma > 0.0)
            for (int i = 0; i < list.size(); i++)
                list[i] /= ma;

		// return true
		return true;
	}

	// normalize deviation
	static bool normalize(Vec <double>& list, int nsize = 0)
	{
		// normalize
		if (nsize == 0)
            nsize = list.size();
            
        // verify
		if (list.size() == 0)
			return false;

		// average
		double key = 0;
		for (int i = 0; i < nsize; i++)
			key += list[i];
		key /= nsize;

		// apply
		for (int i = 0; i < nsize; i++)
			list[i] -= key;		

		// standard deviation
		key = 0;
		for (int i = 0; i < nsize; i++)
			key += pow(list[i], 2.0);
		key /= nsize - 1.0;
		key = sqrt(key);

		// apply
		if (key != 0)
            for (int i = 0; i < nsize; i++)
                list[i] /= key;

		// return true
		return true;
	}

	// str find and replace
	static string strreplace(string search, string replace, string str)
	{
        int found = str.find(search);
        while (found != string::npos)
        {
            str.replace(found, search.size(), replace);
            found = str.find(search);
        }
	    return str;
	}

	// str clean
	static string strclean(string str)
	{
        string output = "";
		for (int i = 0; i < str.length(); i++)
		{
			if (str[i] != ' ' && str[i] != '\n' && str[i] != '\r' && str[i] != '\t' && str[i] != 0)
				output += str[i];
		} 
		return output;
	}
	
	// str resize
	static string strresize(string str, int size, char c = ' ')
	{
		// enlarge
		while (str.length() < size)
			str += c;
		if (str.length() > size)
			str.resize(size);		
		return str;
	}
    	    
	// open file to get list
	template < typename TValue >
	static TValue maximum (TValue a, TValue b, TValue c)
	{
		if (a >= b && a >= c)
			return a;
		if (b > c)
			return b;
		return c;
	}

	// open file to get list
	template < typename TValue >
	static int maxindex (TValue a, TValue b, TValue c)
	{
		if (a >= b && a >= c)
			return 0;
		else if (b >= a && b >= c)
			return 1;
		return 2;
	}

    /**
        structural score
    */    
    static double stscore (double rmsd, int pairs, int sifactor = 100)
    {
        // scaled structural alignment score (SAS)
		return  (rmsd * (double) sifactor) / (double) pairs;
    }

   	// tm score
    template < typename TMol >
	static double tmscore(TMol* mola, TMol* molb, double& rmsd, int& nres, double limit = 5.0)
	{
        // loop
		int n = mola->lcalpha.size();
		int m = molb->lcalpha.size();

		// prepare contact matrix
		Vec < Vec < int > > d;
		
		// setup sizes
		d.resize(n);
		for (int i = 0; i < n; i++)
		{
			d[i].resize(m);
			d[i].fill(0);
		}			

		// generate matrix
		for (int i = 1; i < n; i++)
		for (int j = 1; j < m; j++)
		{   
			// check distance
			int weight = 0;
			if (mola->latom[mola->lcalpha[i]].pos.dist(molb->latom[molb->lcalpha[j]].pos) < limit)
			    weight = 1;
			
            // set left gap, upper gap, match/mismatch    
			d[i][j] = maximum (d[i-1][j-1] + weight, d[i-1][j], d[i][j-1]);
		}

        // reproduce alignment
        int ia = n - 1;
        int ib = m - 1;
        Vec < int > lista, listb;
        while (ia > 0 && ib > 0)
        {
            // value
            int maxvalue = maximum (d[ia-1][ib-1], d[ia-1][ib], d[ia][ib-1]);
            if (d[ia-1][ib-1] == maxvalue)
            {
                // check for match
                if (d[ia][ib] > maxvalue)
                {
                    lista.push_back(mola->lcalpha[ia]);
                    listb.push_back(molb->lcalpha[ib]);
                }
                
                // update
                ia--;
                ib--;
            } else if (d[ia-1][ib] == maxvalue)
                ia--;
            else
                ib--;
        }

        // calculate rmsd
		double dist = 0;
		for(int i = 0; i < lista.size(); i++)
			dist += pow (mola->latom[lista[i]].pos.dist(molb->latom[listb[i]].pos), 2);
		rmsd = sqrt ( dist / lista.size() );
		
		// number of residues
		nres = lista.size();
		
        // calculate score
        double Ln = min (n, m);
        double d0 = 1.24 * pow(Ln - 15, (double) (1.0 / 3.0)) - 1.8;
        double score = 0;
		for(int i = 0; i < lista.size(); i++)
		{
            dist = mola->latom[lista[i]].pos.dist(molb->latom[listb[i]].pos);
			score += 1 / (1 + pow(dist / d0, 2));
        }
		return score / Ln;
	}
	
	/**
		read directory
	*/     
	static void readDirectory(string path, Vec < string >& files, string type = "")
	{
		// clear index
		files.clear();
		
		// read directory
		DIR *dp;
		struct dirent *ep;
		dp = opendir (path.c_str());
		if (dp != NULL)
		{
			while ((ep = readdir (dp)))
			{
                // get name
				string name = ep->d_name;
				
				// check type
				if (type != "")
				    if (type != FName::gettype(ep->d_name))
				        continue;
				        
				// verify and insert
				if (name != "." && name != "..")
					files.push_back(name);
			}
			(void) closedir (dp);
		}
	}
    	
	/**
		tokenizer
	*/
	template < typename TChar >
	static void split(Vec < string >& data, TChar* ptr, char key = ' ')
	{
       	// reset
		data.clear();
		string str = "";
        
        // loop
		bool valid = false;
		while (*ptr != '\0')
		{
			// add data
			if (*ptr == key || *ptr == '\n' || *ptr == '\r' || *ptr == '\t')
			{
				if (valid)
				{
					data.push_back(str);
					str = "";
				}
				
				valid = false;
			} else {
				str += *ptr;
				valid = true;
			}

            // next character
			ptr++;
		}
		
		// add final item
		if (str != "")
			data.push_back(str);	
	}

	// sort by
	template < typename TA, typename TB >
	static void sort (Vec < TA >& data, Vec < TB >& index, int size = 0)
	{
		// get size
		int n = size;
        if (n == 0)
            n = data.size();
        if (index.size() < n)
            Msg::error("Lib::sort()", "Size of index array is too small!");
		
		// sort data
		TA	temp;
		TB	tint;
		for (int i = 0; i < n-1; i++)
		for (int j = i+1; j < n; j++)
		{
			if (data[i] > data[j])
			{
				temp = data[i];
				data[i] = data[j];
				data[j] = temp;

				tint = index[i];
				index[i] = index[j];
				index[j] = tint;
			}
		}
	}
	
	// sort by
	template < typename TA, typename TB >
	static void sort_greater (Vec < TA >& data, Vec < TB >& index, int size = 0)
	{
		// get size
		int n = size;
        if (n == 0)
            n = data.size();
        if (index.size() < n)
            Msg::error("Lib::sort_greater()", "Size of index array is too small!");
		
		// sort data
		TA  temp;
		TB	tint;
		for (int i = 0; i < n-1; i++)
		for (int j = i+1; j < n; j++)
		{
			if (data[i] < data[j])
			{
				temp = data[i];
				data[i] = data[j];
				data[j] = temp;

				tint = index[i];
				index[i] = index[j];
				index[j] = tint;
			}
		}
	}
    	
	/**
		Compare
	*/
	inline static bool compare (Vec <int> d, Vec <int> d0)
	{
		// assert
		assert (d.size() == d0.size());
		
		// compare both arrays
		int cx = 0;
		for (int i = 0; i < d.size(); i++)
		{
			if (d[i] > 0 && d[i] == d0[i])
                cx++;
			
			if (cx > 3)
                return true;
		}
		
		// return true
		return (cx == d.size());
	}
	
    // sign
    template < typename TItem >
    inline static int sign (TItem f)
    {
        // setup sign
	    if (f > 0)
	        return 1;
    	else
    	    if (f < 0)
	        	return -1;
	        else
	        	return 0;
    }
	
    // square
    inline static double sqr(double val)
    {
	    return pow(val, 2);
    }
		
    // square
    inline static double floor(double val, int index = 2)
    {
        double factor = pow (10.0, index);
	    return double (int (val * factor)) / factor;
    }	
};
