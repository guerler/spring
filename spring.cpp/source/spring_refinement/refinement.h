#include "vec4.h"
#include "transmath.h"
#include "rotation.h"
#include "forceobjects.h"
#include "minimizer.h"
#include "forcesingle.h"

// force
struct Refinement
{
	// constructor
	static void minimize (SpecMolecule* target, SpecMolecule* mol, Vec < SpecMolecule >* traj = 0, double weight = 0.01, double epsilon = 0.001, int maxSteps = 100)
	{
        // attraction
        SpecMolecule molstart;
        molstart.construct(mol);
        
        // force set
        ForceSet fs;

        // minimizer
        StepBFGS   adAngle;
        StepRPROP  adForce;

		// construct adjuster
   		adAngle.construct(weight);
		adForce.construct(weight, 3);		

		// initialize minimization schemes
		adForce.initialize();
		adAngle.initialize();

    	// add start
        if (traj != 0)
            traj->push_back(*mol);

		// apply force
		for (int steps = 0; steps < maxSteps; steps++)
		{
			// reset force params
			fs.reset();
								
			// calc center of gravity
			fs.cog = mol->cog();
			
            // calculate
            Vec3 v;
            double r;
            double vdw;
            
            // target repulsion
            for (int i = 0; i < target->latom.size(); i++)
            for (int j = 0; j < mol->latom.size(); j++)
            {
                // get distance
                if (target->latom[i].pos == NOCOORD || mol->latom[j].pos == NOCOORD)
                    continue;

                // get distance
				v = target->latom[i].pos - mol->latom[j].pos;
				r = v.length();

                // get clash
                vdw = target->latom[i].bfactor + mol->latom[j].bfactor;
                if (r < vdw)
                {
                    // potential function
                    v *= -1.0 / (1.0 + r*r);
                
                    // sum rigid body forces
                    fs.force	+= v;
                    fs.amom	    += v.cross(mol->latom[j].pos - fs.cog);
                    fs.energy   += -0.1 / (0.1 + r*r);
                }
            }

            // template attraction
            for (int i = 0; i < molstart.latom.size(); i++)
            {
                if (molstart.latom[i].pos == NOCOORD)
                    continue;

                // get distance
                Vec3 p = mol->latom[i].pos;
				v = molstart.latom[i].pos - p;
				r = v.length();

                // potential function
                v *= 0.001 / (1.0 + r*r);

                // sum rigid body forces
                fs.force	+= v;
                fs.amom	    += v.cross(p - fs.cog);
            }

			// setup
			fs.angle = adAngle.adjust(fs.amom);
			
            // setup force
			for (int i = 0; i < 3; i++)
                fs.force[i] = adForce.adjust(fs.force[i], i);                
					
			// criteria
            if (fs.amom.length() == 0 || fs.force.length() <= epsilon)
                break;

			// apply transformation
            apply (mol, fs);
            
        	// add
            if (traj != 0)
                traj->push_back(*mol);
		}
	}
	

	// apply on whole molecule
	static void apply (SpecMolecule* mol, ForceSet& fs)
	{
		// parameter
		Vec3 tmp;
		Mat4 trans;
		Rotation rot (fs.amom, -fs.angle);

		// generate transform matrix (rigid transformation)
		trans.setTransform(fs.force, rot.q, fs.cog);

		// compute new positions
		for(int i = 0; i < mol->size(); i++)
		{
			trans.multVecMatrix(mol->latom[i].pos, tmp);
			mol->latom[i].pos = tmp;
		}
	}
};


