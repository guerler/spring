// merged molecules
struct MergeMol
{
    // get data
    int na, nb;
    SpecMolecule mab, mba;
    SpecMolecule refa, refb;
    SpecMolecule ma, mb;
    
    // storage
    Storage store;    

    // construct
    void construct(SpecMolecule* mola, SpecMolecule* molb)
    {
        // full
        refa.construct(mola);
        refb.construct(molb);

        // interface
        ma.construct(mola);
        mb.construct(molb);

        // get size
        na = mola->lcalpha.size();
        nb = molb->lcalpha.size();

        // generate merged molecules
        merge (mola, molb, &mab);
        merge (molb, mola, &mba);

        // set center
        mab.setcenter();
        mba.setcenter();
    }
    
    // score
    double tmscore_symmetry (string& log)
    {
        // alignment
        Vec < int > aln;
        TransMatrix tmat;
        
        // first mode
        TMAlign::align(&mab, &mba, aln, tmat);
        double score = min(tmscore(&mab, &mba, aln, 0, na), tmscore(&mab, &mba, aln, na, na + nb));
        double ident = identity (&mab, &mba, aln);
        
        // convert rotation to axis-angle
        Vec3 rotAxis;
        double rotAngle;
        toAxisAngle (tmat.r, rotAxis, rotAngle);
        
        // log
        log = Convert::toString(score) + " " + Convert::toString(ident) + " " + Convert::toString(rotAxis.x) + " " + Convert::toString(rotAxis.y) + " " + Convert::toString(rotAxis.z) + " " + Convert::toString(rotAngle);
        
        // return score
        return score;
    }
    
private:

    void toAxisAngle(Vec < Vec3 >& m, Vec3& v, double& angle)
    {
        // reset
        v.x = 0.0;
        v.y = 1.0;
        v.z = 0.0;
        angle = 0.0;
        
        // margins
        double epsilon = 0.01;

        // check for singularity
        if ((fabs(m[0][1]-m[1][0])< epsilon) && (fabs(m[0][2]-m[2][0])< epsilon)
        	  && (fabs(m[1][2]-m[2][1])< epsilon))
        {
            // first check for identity matrix which must have +1 for all terms
            // in leading diagonaland zero in other terms
            if ((fabs(m[0][1]+m[1][0]) < epsilon)
                && (fabs(m[0][2]+m[2][0]) < epsilon)
                && (fabs(m[1][2]+m[2][1]) < epsilon)
                && (fabs(m[0][0]+m[1][1]+m[2][2]-3) < epsilon))
                // this singularity is identity matrix so angle = 0
                return;

            // otherwise this singularity is angle = 180
            angle = PI;
            double xx = (m[0][0]+1)/2;
            double yy = (m[1][1]+1)/2;
            double zz = (m[2][2]+1)/2;
            double xy = (m[0][1]+m[1][0])/4;
            double xz = (m[0][2]+m[2][0])/4;
            double yz = (m[1][2]+m[2][1])/4;
            if ((xx > yy) && (xx > zz)) { // m[0][0] is the largest diagonal term
                if (xx< epsilon)
                {
				    v.x = 0;
				    v.y = 0.7071;
				    v.z = 0.7071;
                } else {
                    v.x = sqrt(xx);
                    v.y = xy/v.x;
                    v.z = xz/v.x;
                }
            } else if (yy > zz) { // m[1][1] is the largest diagonal term
                if (yy< epsilon) {
                    v.x = 0.7071;
                    v.y = 0;
                    v.z = 0.7071;
                } else {
                    v.y = sqrt(yy);
                    v.x = xy/v.y;
                    v.z = yz/v.y;
                }
            } else { // m[2][2] is the largest diagonal term so base result on this
                if (zz< epsilon) {
                    v.x = 0.7071;
				    v.y = 0.7071;
				    v.z = 0;
                    } else {
				        v.z = sqrt(zz);
				        v.x = xz/v.z;
				        v.y = yz/v.z;
                    }
                }
                return; // return 180 deg rotation
    	}

        // as we have reached here there are no singularities so we can handle normally
    	double s = sqrt((m[2][1] - m[1][2])*(m[2][1] - m[1][2])
                   + (m[0][2] - m[2][0])*(m[0][2] - m[2][0])
                   + (m[1][0] - m[0][1])*(m[1][0] - m[0][1])); // used to normalise
                
    	if (fabs(s) > 0)
        {
            angle = acos(( m[0][0] + m[1][1] + m[2][2] - 1)/2);
            v.x = (m[2][1] - m[1][2])/s;
            v.y = (m[0][2] - m[2][0])/s;
            v.z = (m[1][0] - m[0][1])/s;
        }
    }

	void merge (SpecMolecule* a, SpecMolecule* b, SpecMolecule* m)
	{
        // get size
        int na = a->lcalpha.size();
        int nb = b->lcalpha.size();

		// append as atom
		m->construct();
		for (int i = 0; i < na; i++)
		{
            SpecAtom* atom = &a->latom[a->lcalpha[i]];
            m->append (atom->pos, atom->name, atom->res, i);
        }
		for (int i = 0; i < nb; i++)
		{
            SpecAtom* atom = &b->latom[b->lcalpha[i]];
            m->append (atom->pos, atom->name, atom->res, na + i);
        }
        m->finalize();
    }
    
   	// tm score
	double tmscore(SpecMolecule* merge, SpecMolecule* tmpl, Vec <int>& aln, int nstart, int nend)
	{
        // loop
		int n = tmpl->lcalpha.size();
		int m = nend - nstart;
		
        // calculate score
        double Ln = min (n, m);
        double d0 = 1.24 * pow(Ln - 15, (double) (1.0 / 3.0)) - 1.8;
        double score = 0.0;
		for(int i = nstart; i < nend; i++)
		{
            if (aln[i] != LARGE)
            {
                double dist = tmpl->latom[tmpl->lcalpha[aln[i]]].pos.dist(merge->latom[merge->lcalpha[i]].pos);
                score += 1 / (1 + pow(dist / d0, 2));
            }
        }
		return score / Ln;
	}        
	
   	// tm score
	double identity(SpecMolecule* tmpl, SpecMolecule* merge, Vec <int>& aln)
	{
        // loop
		int n = tmpl->lcalpha.size();
		int m = merge->lcalpha.size();
		
        // calculate score
        double ln = min (n, m);
        double score = 0.0;
		for(int i = 0; i < (int) aln.size(); i++)
		{
            if (aln[i] != LARGE)
                if (tmpl->latom[tmpl->lcalpha[aln[i]]].res == merge->latom[merge->lcalpha[i]].res)
                    score++;
        }
		return score / ln;
	}        	
};
