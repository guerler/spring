// job
struct TransMatrix
{
	// values
	Vec < Vec3 >   r, from, to;
	Vec3           t, center;

    // static
    static const int nalpha = 5;

	// construct
	TransMatrix()
	{
		t = 0.0;
		r.resize(3);
		r.fill(0.0);
	}
    
	// set reference
	void construct (SpecMolecule* mol)
	{
		// read from
		from.resize(nalpha);
		for (int i = 0; i < nalpha; i++)
			from[i] = mol->latom[mol->lcalpha[i]].pos;
	}
	
	// copy constructor
	void construct (TransMatrix& matrix)
	{
		from.assign(matrix.from.begin(), matrix.from.end());
		to.assign(matrix.to.begin(), matrix.to.end());
		r.assign(matrix.r.begin(), matrix.r.end());
		t = matrix.t;
	}

	// set orientations
	void prepare (SpecMolecule* mol)
	{
		// read to
		to.resize(nalpha);
		for (int i = 0; i < nalpha; i++)
			to[i] = mol->latom[mol->lcalpha[i]].pos;
              
		// get matrix
		Trans::align(from, to, r, t);
		
		// center
		center = mol->cog();
	} 

	// set orientations
	void prepare (Vec <Vec3 > pos)
	{        
		// read to
		to.resize(nalpha);
		for (int i = 0; i < nalpha; i++)
			to[i] = pos[i];
              
		// get matrix
		Trans::align(from, to, r, t);      
	}

	// set orientations
	void apply (SpecMolecule* mol, Vec <Vec3 >& pos)
	{
		// read from
		from.resize(nalpha);
		for (int i = 0; i < nalpha; i++)
			from[i] = mol->latom[mol->lcalpha[i]].pos;

		// from
		to.resize(nalpha);
		for (int i = 0; i < nalpha; i++)
			to[i] = pos[i];
                          
		// get matrix
		Trans::align(from, to, r, t);    
		Trans::apply(mol, r, t);
	}
            
	// apply
	void apply (SpecMolecule* mol)
	{
		Trans::apply(mol, r, t);
	} 

	// apply
	void apply (Vec <Vec3>& pos)
	{
		Trans::apply(pos, r, t);
	} 
		
	// print
	void print(double epsilon = 0.000001)
	{
		// clean up zeros
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				if (fabs(r[i][j]) < epsilon)
					r[i][j] = 0.0;
			if (fabs(t[i]) < epsilon)
				t[i] = 0.0;
		}

		// write
		printf (" i          t(i)         u(i,1)         u(i,2)         u(i,3)\n");
		printf (" 1    %14.10f %14.10f %14.10f %14.10f\n", t[0], r[0][0], r[0][1], r[0][2]);
		printf (" 2    %14.10f %14.10f %14.10f %14.10f\n", t[1], r[1][0], r[1][1], r[1][2]);
		printf (" 3    %14.10f %14.10f %14.10f %14.10f\n", t[2], r[2][0], r[2][1], r[2][2]);
	}
    
    // reverse
    void reverse()
    {        
        Trans::reverse(r, t);
    }
            	    
	// get string
	string toString()
	{
		return	Convert::toString(r[0].x) + " " +
				Convert::toString(r[0].y) + " " +
				Convert::toString(r[0].z) + " " +
				Convert::toString(r[1].x) + " " +
				Convert::toString(r[1].y) + " " +
				Convert::toString(r[1].z) + " " +
				Convert::toString(r[2].x) + " " +
				Convert::toString(r[2].y) + " " +
				Convert::toString(r[2].z) + " " +				
				Convert::toString(t.x) + " " +
				Convert::toString(t.y) + " " +
				Convert::toString(t.z);
	}
		
	// prepare matrix
   	static double load (SpecMolecule* mol, string input, double& rmsd, bool invert = false)
	{
		// split data
		Vec < string > data;
		Lib::split (data, input.c_str());
        
		// verify
		if (data.size() <= 12)
			return LARGE;
            
		// load data
		Vec < Vec3 > r;
		r.resize(3);
		r.fill(0.0);
		Vec3 t = 0.0;
        
		// read		
		r[0].x    = Convert::toDbl(data[0]);
		r[0].y    = Convert::toDbl(data[1]);
		r[0].z    = Convert::toDbl(data[2]);
		r[1].x    = Convert::toDbl(data[3]);
		r[1].y    = Convert::toDbl(data[4]);
		r[1].z    = Convert::toDbl(data[5]);
		r[2].x    = Convert::toDbl(data[6]);
		r[2].y    = Convert::toDbl(data[7]);
		r[2].z    = Convert::toDbl(data[8]);
		t.x       = Convert::toDbl(data[9]);
		t.y       = Convert::toDbl(data[10]);
		t.z       = Convert::toDbl(data[11]);
        
		// apply
		Trans::apply(mol, r, t);

		// read last column rmsd
		if (data.size() <= 13)
			rmsd = (double) LARGE;
		else
			rmsd = Convert::toDbl(data[13]);
                
		// return
		return Convert::toDbl(data[12]);
	}	
};
