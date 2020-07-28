// local library
struct FName
{
	// returns path/path/name.typ.typ -> name.typ.typ
    static string getfullname(string path)
	{
		return path.substr(path.rfind("/") + 1);
	}
	
	// returns path/path/name.typ.typ -> name
	static string getname(string path)
	{
		// find
		int start = path.rfind("/") + 1;
		if (start != -1)
            path = path.substr(start);

        // find
	    int end = path.find(".");
        if (end != -1)
            path = path.substr(0, end);
		
		// return		
		return path;
	}

	static string reduceprefix(string path)
	{
		return path.substr(path.find("_") + 1);
	}

	// returns name.typ.typ -> name.typ
	static string reduce(string fname)
	{
        int pos = fname.rfind(".");
        if (pos != -1)
    		return fname.substr(0, fname.rfind("."));
    	else
            return fname;
	}

	// returns name.typ.typ -> .typ	
	static string gettype(string fname)
	{
        int pos = fname.rfind(".");
        if (pos != -1)
    		return fname.substr(pos);
    	else
            return "";
	}
};
