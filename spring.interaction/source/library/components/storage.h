class Storage
{
	// file
	File f;	

    // write all atoms also for models
    bool nocoordmode;
public:
	// molecule
	string path, pathalternative;
	
	// construct
	Storage()
	{
        construct();
    }

	// construct
	void construct (string path = "", string pathalternative = "")
	{
		// initialize
		this->path              = path;
		this->pathalternative   = pathalternative;
		this->nocoordmode       = false;
	}   
	
	// set nocoordmode
	void setnocoord(bool mode)
	{
        nocoordmode = mode;
    }
	
	// read
	bool read (SpecMolecule* mol, string pdbname, string name = "Molecule")
	{        
		// setup molecule
		mol->construct(name);

		// load pdb with/without hydrogens
		pdb(mol, pdbname);

		// finalize
		if (mol->finalize())
		{
			if (mol->lcalpha.size() == 0)
			{
				Msg::write("Molecule has no residues.");
				return false;
			} else
				return true;   
		} else
			return false; 
	}

	// pdb reader
	void pdb(SpecMolecule* mol, string name, bool verbose = true)
	{   
		// read file     
		load (name);

        // key
		string key = "";
		
		// loop
        int nres = 0;
		int resno = 0;
		string chain = "";
		string resname = "";
		
		// residue information
        Vec3 respos;
		
		// validate residue
		SpecMolecule residue;
		residue.construct();

		// read file
		while (f.good())
		{
			// read
			f.readLine();

			// get pdb key
			key = Lib::strclean(f.get(0, 6));

			// end model
			if (key == "END" || key == "ENDMDL" || key == "TER")
				break;

			// atom
			if (key == "ATOM" || (key == "HETATM" && f.get(21, 1) == chain))
			{
				// generate index for sse information				
				if (resname == "")
				{
					// backup residue details
                    resname = f.get(17, 3);
					chain   = f.get(21, 1);
	       			resno   = Convert::toInt(f.get(22, 4));
					
                    // reset
					respos  = 0.0;
				}
		        
				// verify new residue
				if (resname != f.get(17, 3) || chain != f.get(21, 1) || resno != Convert::toInt(f.get(22, 4)))
				{
					// verify
					if (respos == 0.0)
					{
                        if (verbose)
                            Msg::write("Residue %i incomplete in %s.", nres, name.c_str());
				    } else {
    					// copy residue into molecule
                        for (int i = 0; i < residue.size(); i++)
                            mol->append(residue.latom[i]);
                        
                        // count
                        nres++;
                    }
                
					// reset residue
					residue.construct();
	
					// backup residue details
                    resname = f.get(17, 3);
					chain   = f.get(21, 1);
	       			resno   = Convert::toInt(f.get(22, 4));
					
                    // reset
					respos     = 0.0;					
				}

				// get atom name and coordinate
				string atomname = f.get(12, 4);
				Vec3 atompos = Vec3( Convert::toDbl(f.get(30, 8)) , Convert::toDbl(f.get(38, 8)), Convert::toDbl(f.get(46, 8)) );

				// check calpha coordinate
				if (Lib::strclean(atomname) == "CA")
				{
                    // calpha already set
                    if (respos == 0.0)
                    {  
                        // get first calpha coordinate  
                        respos = atompos;
                    
                        // append calpha
                        residue.append(atompos, atomname, resname, nres, resno, Convert::toDbl(f.get(54, 6)), Convert::toDbl(f.get(60, 6)));
                    }
                } else {
                    // append			
                    residue.append(atompos, atomname, resname, nres, resno, Convert::toDbl(f.get(54, 6)), Convert::toDbl(f.get(60, 6)));
                }
			} else if (key == "HETATM")
				mol->appendhet(Vec3( Convert::toDbl(f.get(30, 8)) , Convert::toDbl(f.get(38, 8)), Convert::toDbl(f.get(46, 8)) ), f.get(12, 4), f.get(17, 3));
		}
        		
		// append last residue
        if (respos != 0.0)
        {
            // add atoms
            for (int i = 0; i < residue.size(); i++)
                mol->append(residue.latom[i]);
        }

		// return
		f.close();
	}

	// save
	void save (Vec < SpecMolecule >& molb, string name)
	{
		FILE *fp;
		open (fp, name);
	
		for (int i = 0; i < molb.size(); i++)
			print (fp, &molb[i], Convert::toString(i+1));
	
		// finalize
		close (fp);
	}
	
	// save
	void save (Vec < SpecMolecule* >& molb, string name)
	{
		FILE *fp;
		open (fp, name);
	
		for (int i = 0; i < molb.size(); i++)
			print (fp, molb[i], Convert::toString(i+1));
	
		// finalize
		close (fp);
	}
	
	// save
	void save (SpecMolecule* mola, Vec < SpecMolecule* >& molb, string name)
	{
		FILE *fp;
		open(fp, name);    
	
		// print
		print (fp, mola, "A", "receptor");
		for (int i = 0; i < molb.size(); i++)
		{
			if (i == 0)
				print (fp, molb[i], "B", "ligand reference");
			else
				printMinimal (fp, molb[i], "B", "ligand " + Convert::toString(i-1));                
		}
	
		// finalize
		close (fp);
	}
		   
	void save (SpecMolecule* mola, SpecMolecule* molb, string name)
	{
		FILE *fp;
		open(fp, name);                
		print (fp, mola);
		print (fp, molb, "B");
		close (fp);
	}

	void save (SpecMolecule* mol, string name, bool minimal = false)
	{
		FILE *fp;
		open (fp, name);
		
		// print
		if (minimal)
    		printMinimal (fp, mol);		
		else
    		print (fp, mol);

		// finalize
		close (fp);
	}

private:
	// open file
	void open(FILE*& fp, string name)
	{
		string nm (name);
	    if ((fp = fopen(nm.c_str(), "w")) == NULL)
			Msg::error("open()", "Can not open file " + string(nm.c_str()));
	}
	
	void load(string name)
	{
		// setup path		
		if (!f.open(path + name))
			if (!f.open(pathalternative + name))
				Msg::write ("Unable to open file %s.", string (path + name).c_str());
	}
            
	// close file
	void close(FILE* fp)
	{
		fflush(fp);
		fclose(fp);
	}
	
	// write molecule to file
	void print (FILE *fp, SpecMolecule* mol, string chain = "A", bool details = false)
	{                  
        // info
        if (mol->info != "")
            fprintf (fp, "%s", string (mol->info + "\n").c_str());            
            
		// print more details
		if (details)
		{
			// are there sses
			if (mol->lsse.size() > 0)
			{
				fprintf (fp, "REM Secondary structure information\n");
				fprintf (fp, "REM No. Type \t From \t To\n");
				for (int i = 0; i < mol->lsse.size(); i++)
					fprintf (fp, "REM %i \t %c \t %i \t %i\n", i, mol->lsse[i].getType(), mol->lsse[i].getStart(), mol->lsse[i].getEnd());
			}
			
			// calphas
			fprintf (fp, "REM C-Alpha table\n");
			for (int i = 0; i < mol->lcalpha.size(); i++)
				fprintf (fp, "REM %i \t %i \n", i, mol->lcalpha[i]);
	
			// sse per residues
			fprintf (fp, "REM SSE-Residue table\n");
			for (int i = 0; i < mol->lssebyres.size(); i++)
				fprintf (fp, "REM %i \t %c \n", i, mol->lssebyres[i]);
		}
		
		// print atoms
		for (int i = 0; i < mol->size(); i++)
		{
            if (mol->latom[i].pos.x != NOCOORD || (mol->latom[i].pos.x == NOCOORD && nocoordmode))
            {
                fprintf (fp, "ATOM  %5d  %-4s%-3s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12c    %12c\n",
                        i,
				        mol->latom[i].name.c_str(),
						mol->latom[i].res.c_str(),
						chain.c_str(),
						min (9999, mol->latom[i].respdb),
						mol->latom[i].pos.x, mol->latom[i].pos.y, mol->latom[i].pos.z,
						min(99.0, mol->latom[i].occupancy), max (-99.0, min(99.0, mol->latom[i].bfactor)), mol->latom[i].name[0], mol->latom[i].sse);
            }
        }

		// print hetatoms
		for (int i = 0; i < mol->lhetatom.size(); i++)
		{
			fprintf (fp, "HETATM%5d  %-4s%-4s%1s%4d    %8.3f%8.3f%8.3f  1.00  1.00%12c    %12c\n",
						i,
						mol->lhetatom[i].name.c_str(),
						mol->lhetatom[i].res.c_str(),
						"H",
						min (9999, mol->lhetatom[i].respdb),
						mol->lhetatom[i].pos.x, mol->lhetatom[i].pos.y, mol->lhetatom[i].pos.z,
						mol->lhetatom[i].name[0], mol->lhetatom[i].sse);
        }
        
        // terminate
        fprintf(fp, "TER\n");
        
        // make model tag
        if (Convert::toInt(chain) > 0)
            fprintf (fp, "END\n");        
	}
	
	// write molecule in minimal mode
	void printMinimal (FILE *fp, SpecMolecule* mol, string chain = "A", string tag = "")
	{
		// log
		if (tag != "")
		{
			tag = "REM " + tag + "\n";			
			fprintf (fp, "%s", tag.c_str());
		}
			
		// loop
		for (int index = 0; index < mol->lcalpha.size(); index++)
		{
			int i = mol->lcalpha[index];
			fprintf(fp,"ATOM  %5d  %-4s%-4s%1s%4d    %8.3f%8.3f%8.3f\n",
						i, 	mol->latom[i].name.c_str(),
						mol->latom[i].res.c_str(),
						chain.c_str(),
						min (9999, mol->latom[i].respdb),
						mol->latom[i].pos.x, mol->latom[i].pos.y, mol->latom[i].pos.z);
		}
        
        // terminate
        fprintf(fp, "TER\n");
	}
};


