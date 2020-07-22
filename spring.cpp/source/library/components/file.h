// parameters
char buffer[BSIZE];

// file handle
class File
{
	// file stream
	ifstream* ifs;
	
    // maxsize
    int ndata;	    
public:
	// data
	Vec < string > data;
	
	// real construct
	File()
	{
		ifs = 0;
		ndata = 0;
	}

	// construct
	bool construct(string filename)
	{
		return open (filename);
	}
	
	// open file	
	bool open (string filename)
	{
		// file stream
		if (ifs != 0)
			delete ifs;

		// open stream
		ifs = new ifstream;

		// open filestream
		ifs->open (filename.c_str(), ios::in);
		if (!ifs->good())		
		    return false;
		else
		    return true;        
	}
	
	// destructor
	~File()
	{
        close();
    }
	
	// close
	void close()
	{
		if (ifs != 0)
		{
			delete ifs;
			ifs = 0;
		}
	}
    
	// open file
	static bool exists (string filename)
	{
		// open file
		ifstream ifs;
		ifs.open (filename.c_str(), ios::in);
		if (ifs.good())
			return true;
		else
			return false;        
	}
		
	// good
	inline bool good()
	{
		return ifs->good();
	}

	// default get
	inline string get(int i = 0)
	{
		if (i < ndata)
			return data[i];
		else
			return "";
	}
    
	inline string get(int i, int j)
	{
        if (i >= 0 && i < ndata)
    		return data[0].substr(i, j);
        else
            return "";
	}
    
	// get
	inline char getChar (int i = 0)
	{
		return Convert::toChar(get(i));
	}
	
	inline int getInt (int i = 0)
	{
		return Convert::toInt(get(i));
	}
	
	inline double getDbl(int i = 0)
	{
		return Convert::toDbl(get(i));
	}

	inline double getDouble(int i = 0)
	{
		return getDbl(i);
	}
        
	// reader
	inline bool readLine()
	{
		// read
		ifs->getline(buffer, BSIZE, '\n');
		data.resize(1);
		data[0] = string (buffer);
		ndata = data[0].length();
		
		// return
		if (ndata > 0)
			return true;
		else
			return false;
	}
    
	inline bool read()
	{
        // load
        ifs->getline(buffer, BSIZE, '\n');
    
        // split
        Lib::split (data, buffer);
        ndata = data.size();
        
		// return
		if (ndata > 0)
			return true;
		else
			return false;		
	}

	// return size of entries
	inline int size()
	{
		return ndata;
	}
};
