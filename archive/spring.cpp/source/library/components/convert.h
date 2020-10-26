// local library
struct Convert
{
	// anything to string
    template < typename TAny >
    inline static string toString(TAny t)
    {
        std::ostringstream sstrm;
        sstrm << t;
        return sstrm.str();
    }

	// format
	inline static string toString(double t)
	{
		// check
		if (t == LARGE)
			return "undefined";
		
		// return
        std::ostringstream sstrm;
        sstrm << t;
        return sstrm.str();
	}
	
	// format
	inline static float toFloat(string str)
	{
		return atof(str.c_str());
	}
	
	inline static double toDbl(string str)
	{
		return atof(str.c_str());
	}
	
	inline static char toChar (string str)
	{
		return str.c_str()[0];
	}
	
	inline static int toInt (string str)
	{
		return atoi(str.c_str());
	}	
};
