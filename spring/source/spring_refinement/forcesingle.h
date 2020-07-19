// force
template < typename TPotential >
class ForceSingle
{
	// molecule
	SpecMolecule* mol;
	SpecMolecule* target;

  	// force set
   	ForceSet fs;
    
	// potential function
	TPotential potential;
	
	// minimizer
	StepBFGS	adAngle;
	StepRPROP	adForce;
	
	// parameters
	double 	epsilon, weight;
	int 	maxSteps;
public:		
	// topologies
	ForceTopo 	ft;
	
    // constructor
	void construct (SpecMolecule* target, SpecMolecule* mol, double epsilon = 0.001, double weight = 0.5, int maxSteps = 50)
	{
		// parameters		
		this->target	= target;
		this->mol		= mol;

		// potential
		potential.construct(&fs, target, mol);
		
		// load parameters
		this->epsilon 	= epsilon;
		this->weight	= weight;
		this->maxSteps	= maxSteps;
		
    	// construct adjuster
   		adAngle.construct(weight);
		adForce.construct(weight, 3);
	}
	
	// improve positioning
	void minimize()
	{
        // initialize force
		initialize();
		
    	// apply force
		int steps = maxSteps;
		while (steps > 0 && !next())
		{
        	Trans::apply (mol, fs);
			steps--;
    	}
    }

	// improve positioning
	void minimize_pairwise()
	{
        // initialize force
		initialize();
		
    	// apply force
		int steps = maxSteps;
		while (steps > 0 && !next_pairwise())
		{
        	Trans::apply (mol, fs);
			steps--;
    	}
    }
        
    // return energy
    double energy()
    {
        return fs.energy;
    }
private:
	// get force
	bool next_pairwise()
	{	
        // invalid
        if (target->lcalpha.size() != mol->lcalpha.size())
            return false;
        
     	// reset force params
		fs.reset();
    							
    	// calc center of gravity
		fs.cog = mol->cog();
		
		// evaluate
		for(int i = 0; i < target->lcalpha.size(); i++)
            potential.evaluate(target->lcalpha[i], mol->lcalpha[i]);

		// setup
        fs.angle = adAngle.adjust(fs.amom);
		
        // setup force
		for (int i = 0; i < 3; i++)
            fs.force[i] = adForce.adjust(fs.force[i], i);        
        
        // criteria
        if (fs.force.length() < epsilon && fs.amom.length() < epsilon)
            return true;
        else
            return false;
	}
        
	// get force
	int next()
	{	
     	// reset force params
		fs.reset();
    							
    	// calc center of gravity
		fs.cog = mol->cog();
		
		// evaluate
		for (int k = 0; k < ft.a.size(); k++)
			for (int j = 0; j < ft.a[k].size(); j++)
				for(int i = 0; i < ft.b[k].size(); i++)
					potential.evaluate(ft.a[k][j], ft.b[k][i]);

		// setup
        fs.angle = adAngle.adjust(fs.amom);
		
        // setup force
		for (int i = 0; i < 3; i++)
            fs.force[i] = adForce.adjust(fs.force[i], i);        
        
        // criteria
        if (fs.force.length() < epsilon && fs.amom.length() < epsilon)
            return true;
        else
            return false;
	}

	void initialize()
	{
		// initialize step manager
		adForce.initialize();
		adAngle.initialize();
	}
};

