// interface isscore
struct NAccess
{
	// obtain alignment and score
	static bool make (SpecMolecule* mol)
	{
        // verify
        if (mol->lcalpha.size() <= 5)
            return false;

        // pars
        Storage store;
        string uuid;
        string workpath;
        string unique;

        // obtain
        Config cnf;
        workpath = cnf.ptTemp;
        uuid = cnf.get_uuid();

        // unqiue filename combo
        unique = workpath + uuid;

        // save
        store.save (mol, unique + ".pdb");

        // calc
        chdir("plugins/naccess");
        string execute = "./naccess.sh " + unique + ".pdb";
        system (execute.c_str());
        chdir("../../");
        
        // read
        SpecMolecule molrebuilt;
        bool success = store.read(&molrebuilt, unique + ".asa");

        // copy pulchra model
        if (success)
        {
            // copy molecule
            mol->construct(&molrebuilt);

            if (molrebuilt.lcalpha.size() != mol->lcalpha.size())
                Msg::write ("naccess ignored residues and was canceled.");
        } else
            Msg::write ("naccess failed.");

        // clean
        system (string("rm " + unique + "*").c_str());

        // return
        return success;
    }
};

