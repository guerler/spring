// top model information
struct InterfaceQuality_Info
{
    // properties
    double interface_contacts, interface_isscore, interface_min_tm, global_max_rmsd, global_tm, interface_rmsd, global_rmsd, interface_coverage, global_coverage, global_min_tm;
    
    // backup transformation matrix
    TransMatrix tmat;
    
    // info object
    InterfaceQuality_Info()
    {
        reset();
    }
    
    // info object
    void reset()
    {
        interface_rmsd = global_rmsd = interface_coverage = global_coverage = LARGE;
        interface_contacts = interface_isscore = interface_min_tm = global_min_tm = global_tm = global_max_rmsd = -LARGE;
    }
        
    // make string
    string toString()
    {
        return Convert::toString(global_tm) + " " + Convert::toString(interface_contacts) + " " +
		       Convert::toString(global_min_tm) + " " +
               Convert::toString(interface_rmsd) + " " + Convert::toString(interface_coverage) + " " +
               Convert::toString(global_rmsd) + " " + Convert::toString(global_coverage);
    }
    
    // get title
    static string getHeader()
    {
        return "global_tm interface_contacts global_min_tm interface_rmsd interface_coverage global_rmsd global_coverage";
    }  
};

// gobal interaction rmsd
class InterfaceQuality
{
	// reference
    SpecMolecule *targetref, *molref;	
public:

	// obtain interface
	bool construct (SpecMolecule* target, SpecMolecule* mol)
	{	
        // reference
        this->targetref = target;
        this->molref = mol;
        
        // return
        return true;
	}
    
    // initialize orientation
    InterfaceQuality_Info get(SpecMolecule* targetpredicted, SpecMolecule* molpredicted)
	{
        /** COMPLEX QUALITY **/
				
		// make a copy
		SpecMolecule target;
		target.construct(targetpredicted);

		// make a copy
		SpecMolecule mol;
		mol.construct(molpredicted);
				
        // info object		
        InterfaceQuality_Info info;

        // interface contacts
        FNat fnat;
        fnat.construct(targetref, molref);
        info.interface_contacts = fnat.get(&target, &mol);
        
        // construct
        IRmsdInterface qint;
        qint.construct(targetref, molref);
        if(qint.initialize(&target, &mol))
        {
        	info.interface_rmsd = qint.getrmsd();
        	info.interface_coverage = qint.getcoverage();
		}

        // construct
        IRmsd qglobal;
        qglobal.construct(targetref, molref);
        if(qglobal.initialize(&target, &mol))
        {
        	info.global_rmsd = qglobal.getrmsd();
	    	info.global_coverage = qglobal.getcoverage();
		}

        // construct tmscore
        SpecMolecule mreference, mtarget, mmol;
        mreference.construct(targetref, molref);
        mtarget.construct(&target, &mol);
        info.global_tm = TMScore::align(&mtarget, &mreference, info.tmat); 

        // minimum interfacial tm-score
		MinTM qmtm;
        qmtm.construct(targetref, molref);
        info.interface_min_tm = qmtm.get(&target, &mol);

        /** MONOMERIC QUALITY **/
        
        // minimum monomeric tm to native with no transformation
        info.global_min_tm = min(TMScore::align(&target, targetref), TMScore::align(&mol, molref));       

        // maximum with transformation
        Trans::align (&target, targetref);
        Trans::align (&mol, molref);
        info.global_max_rmsd = max(targetref->rmsd(&target), molref->rmsd(&mol));       
                                                                       
        // return
        return info;
    }
};

