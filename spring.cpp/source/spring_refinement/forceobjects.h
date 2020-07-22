// force set
struct ForceSet
{    
	// parameter set
	Vec3 	force, amom, cog;
	double 	angle, energy;

	// load
	ForceSet()
	{
        // reset
		reset();
	}
	
	// reset parameters
	void reset()
	{
        // results
        cog     = Vec3(0.0, 0.0, 0.0);
		force	= Vec3(0.0, 0.0, 0.0);
		amom	= Vec3(0.0, 0.0, 0.0);
		energy  = 0.0;
		angle   = 1.0;		
	}
};

// force topology
struct ForceTopo
{
	// index arrays
	Vec < Vec < int > > a;
	Vec < Vec < int > > b;

	// resize
	template < typename TTopo >
	void construct (TTopo* ft)
	{
		// clear
		a.clear();
		b.clear();
		
		// resize
		a.insert(a.begin(), ft->a.begin(), ft->a.end());
		b.insert(b.begin(), ft->b.begin(), ft->b.end());
	}
    
	// clear
	void clear()
	{
		a.clear();
		b.clear();
	}

	// resize
	void resize(int n)
	{
		// clear
		a.clear();
		b.clear();
		
		// resize
		a.resize(n);
		b.resize(n);
	}

	// reset
	void reset(int n)
	{
        resize(n);
	}
    	
	// number of groups
	int nGroups()
	{
		assert (a.size() == b.size());
		return a.size();
	}
	
	// count interaction pairs
	int npairs()
	{
		int n = 0;
		for (int k = 0; k < a.size(); k++)
			n += a[k].size() * b[k].size();
		return n;
	}
};


