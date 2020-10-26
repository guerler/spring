// check for symmetry
struct Symmetry
{
    // link
    SpecMolecule target, mol;

    // construct
    void construct (SpecMolecule* reftarget, SpecMolecule* refmol, string outputname)
    {
        // link
        target.construct(reftarget);
        mol.construct(refmol);

        // molecules
        MergeMol merge;
        merge.construct(&target, &mol);

        /*/ align
        string log = "";
        merge.tmscore_symmetry(log);
        cout << log << endl;*/
        makemultimeric(outputname);
    }

    // write higher order
    void makemultimeric(string outputname, int maximum_symmetry = 50)
    {
        // number of c-alpha atoms
        double ncalpha = (double) mol.lcalpha.size();
            
        // list of molecules
        Vec < SpecMolecule > lcomplex;

        // set first two monomers into list
        lcomplex.push_back(target);
        lcomplex.push_back(mol);

        // generation loop
        int i = 0;
        while (i < maximum_symmetry)
        {
            // increase counter
            i++;

            // re-orient
            TransMatrix tmat;
            TMAlign::align(&target, &lcomplex[(int) lcomplex.size() - 1], tmat);
            tmat.apply(&mol);

            // check for clashes
            double clash_fract  = 0.0;
            double clash_radius = sqrt(SpringConfig::RADIUS2);
            for (int j = 0; j < lcomplex.size(); j++)
            for (int m = 0; m < lcomplex[j].lcalpha.size(); m++)
            for (int n = 0; n < mol.lcalpha.size(); n++)
                if (lcomplex[j].latom[lcomplex[j].lcalpha[m]].pos.dist(mol.latom[mol.lcalpha[n]].pos) < clash_radius)
                    clash_fract++;

            // normalize
            if (clash_fract > 0.0)
                clash_fract /= ncalpha;

            // check
            if (clash_fract > 0.10)
                break;

            // add to list
            lcomplex.push_back(mol);
        }

        // write multi-complex
        Msg::write("Multimer with %i units generated.", (int) lcomplex.size());
        Storage store;
        store.save(lcomplex, outputname);
    }
};
