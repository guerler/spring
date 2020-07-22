// transformation
struct Trans
{
	/**
		apply transformation
	**/

	// apply translation on whole molecule
	template < typename TVector >
	static void apply (SpecMolecule* mol, TVector& shift)
	{
		// compute new positions
		for(int i = 0; i < mol->size(); i++)
		{
			mol->latom[i].pos.x += shift.x;
			mol->latom[i].pos.y += shift.y;			
			mol->latom[i].pos.z += shift.z;			
		}			
	}

	// transformation matrix
	inline static void apply (Vec < Vec3>& pos, Vec < Vec3 >& m, Vec3& t)
	{
		// apply translation on every atom
		double x, y, z;
		for (int i = 0; i < pos.size(); i++)
		{
			// verify
			if (pos[i] == NOCOORD)
				continue;
                            
			// calculate
			x = pos[i].x * m[0][0] + pos[i].y * m[0][1] + pos[i].z * m[0][2] + t[0];
			y = pos[i].x * m[1][0] + pos[i].y * m[1][1] + pos[i].z * m[1][2] + t[1];
			z = pos[i].x * m[2][0] + pos[i].y * m[2][1] + pos[i].z * m[2][2] + t[2];
						
			// update
			pos[i].x = x;
			pos[i].y = y;
			pos[i].z = z;
		}
	}
	        
	// transformation matrix
	inline static void apply (SpecMolecule* mol, Vec < Vec3 >& m, Vec3& t)
	{
		// 
		// apply translation on every atom
		double x, y, z;
		for (int i = 0; i < mol->size(); i++)
		{
			// verify
			if (mol->latom[i].pos == NOCOORD)
				continue;

			// calculate
			x = mol->latom[i].pos.x * m[0][0] + mol->latom[i].pos.y * m[0][1] + mol->latom[i].pos.z * m[0][2] + t[0];
			y = mol->latom[i].pos.x * m[1][0] + mol->latom[i].pos.y * m[1][1] + mol->latom[i].pos.z * m[1][2] + t[1];
			z = mol->latom[i].pos.x * m[2][0] + mol->latom[i].pos.y * m[2][1] + mol->latom[i].pos.z * m[2][2] + t[2];
						
			// update
			mol->latom[i].pos.x = x;
			mol->latom[i].pos.y = y;
			mol->latom[i].pos.z = z;
		}

		// all other atoms		
		for (int i = 0; i < mol->lhetatom.size(); i++)
		{            
			// calculate
			x = mol->lhetatom[i].pos.x * m[0][0] + mol->lhetatom[i].pos.y * m[0][1] + mol->lhetatom[i].pos.z * m[0][2] + t[0];
			y = mol->lhetatom[i].pos.x * m[1][0] + mol->lhetatom[i].pos.y * m[1][1] + mol->lhetatom[i].pos.z * m[1][2] + t[1];
			z = mol->lhetatom[i].pos.x * m[2][0] + mol->lhetatom[i].pos.y * m[2][1] + mol->lhetatom[i].pos.z * m[2][2] + t[2];
						
			// update
			mol->lhetatom[i].pos.x = x;
			mol->lhetatom[i].pos.y = y;
			mol->lhetatom[i].pos.z = z;
		}
	}

	        
	// transformation matrix
	inline static void reverse (Vec < Vec3 >& r, Vec3& t)
	{
        // swap 0
        double s = 0;
        s = r[1][0];
        r[1][0] = r[0][1];
        r[0][1] = s;

        // swap 0
        s = r[2][0];
        r[2][0] = r[0][2];
        r[0][2] = s;

        // swap 0
        s = r[1][2];
        r[1][2] = r[2][1];
        r[2][1] = s;
        
        // Invert the translation part by negating it; and then rotate it by the new rotation part
        Vec3 v = Vec3(-t.x, -t.y, -t.z);
        t.x = v.x * r[0][0] + v.y * r[0][1] + v.z * r[0][2];
        t.y = v.x * r[1][0] + v.y * r[1][1] + v.z * r[1][2];
        t.z = v.x * r[2][0] + v.y * r[2][1] + v.z * r[2][2];
    }
    		
	/**
		align
	**/

	// kabsch alignment
	inline static bool align (SpecMolecule* mola, SpecMolecule* molb)
	{
        // perform sequence alignment
        Vec < int > lmola, lmolb;
        SequenceAlignment sa;
        sa.construct();
        sa.getSharedResidues(mola, molb, lmola, lmolb);

		// list
		Vec < Vec3 > lista, listb;
        		
		// copy target atoms to lista
		for (int i = 0; i < lmola.size(); i++)
		{
            // get assigned residues
            SpecAtom* atoma = &mola->latom[mola->lcalpha[lmola[i]]];
            SpecAtom* atomb = &molb->latom[molb->lcalpha[lmolb[i]]];

            // check if they have valid coordinates
            if (atoma->pos != NOCOORD && atomb->pos != NOCOORD)
            {
    			lista.push_back(atoma->pos);
                listb.push_back(atomb->pos);    			
            }
        }

		// apply alignment
		return Trans::align(lista, listb, mola);
	}

	// kabsch alignment
	inline static bool align (SpecMolecule* mola, SpecMolecule* molb, Vec < Vec3 >& rot, Vec3& trans)
	{
		// list
		Vec < Vec3 > lista, listb;
		
		// copy target atoms to lista
		int minlength = min (mola->lcalpha.size(), molb->lcalpha.size());
		for (int i = 0; i < minlength; i++)				
		{
            if (mola->latom[mola->lcalpha[i]].pos != NOCOORD && molb->latom[molb->lcalpha[i]].pos != NOCOORD)
            {
    			lista.push_back(mola->latom[mola->lcalpha[i]].pos);
                listb.push_back(molb->latom[molb->lcalpha[i]].pos);    			
            }
        }

		// apply alignment
		return align(lista, listb, rot, trans);
	}
	
	// kabsch alignment
	inline static bool align (Vec < Vec3 >& a, Vec < Vec3 >& b)
	{
		// parameters
		Vec < Vec3 > r;
		Vec3 t;

		// apply
		if (!align(a, b, r, t))
            return false;
        
        // apply
		apply(a, r, t);
		
		// return
		return true;
	}

	// kabsch alignment
	inline static bool align (Vec < Vec3 >& a, Vec < Vec3 >& b, SpecMolecule* mol)
	{
		// parameters
		Vec < Vec3 > r;
		Vec3 t;

		// apply
		if (!align(a, b, r, t))
            return false;
        
        // apply
		apply(mol, r, t);
		
		// return
		return true;
	}		

	// kabsch alignment
	inline static bool align (Vec < Vec3 >& a, Vec < Vec3 >& b, Vec < Vec3 >& rot, Vec3& trans)
	{
		// amount
		int n = a.size();

		// results
        double rms     = LARGE;
        double t[3]    = {0,0,0};
        double r[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
            
		// reset
		if (n == b.size() && n >= 3)
		{
			// setup matrices a and b			
			double **ref;
			double **mov;
			
			// initialize			
			KabschArray (&ref, n, 3);
			KabschArray (&mov, n, 3);
						
			for (int i = 0; i < n; i++)
				for (int j = 0; j < 3; j++)
				{
					mov[i][j] = a[i][j];
					ref[i][j] = b[i][j];
				}
  
			// get matrix
			Kabsch (mov, ref, n, 1, &rms, t, r);
			
			// delete array
			KabschDeleteArray (&ref, n);
			KabschDeleteArray (&mov, n);
		} else
            return false;

		// copy transformation
		rot.resize(3);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
   				rot[i][j] = r[i][j];
		for (int i = 0; i < 3; i++)	    	
			trans[i] = t[i];
			
		// return
		return true;
	}
};

