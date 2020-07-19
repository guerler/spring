// interface isscore
struct Pulchra
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
        unique = workpath + uuid + ".";

        // save
        store.save (mol, unique + "mol.pdb");

        // calc
        string execute = "./plugins/pulchra/pulchra " + unique + "mol.pdb";
        system (execute.c_str());

        // read
        SpecMolecule molrebuilt;
        bool success = store.read(&molrebuilt, unique + "mol.rebuilt.pdb");

        // copy pulchra model
        if (success)
        {
            // copy molecule
            mol->construct(&molrebuilt);

            // check size for PULCHRA may fail for some residues
            if (molrebuilt.lcalpha.size() != mol->lcalpha.size())
                Msg::write ("Pulchra ignored residues and was canceled.");
        } else
            Msg::write ("Pulchra failed.");

        // clean
        system (string("rm " + unique + "*").c_str());

        // return
        return success;
    }
};

