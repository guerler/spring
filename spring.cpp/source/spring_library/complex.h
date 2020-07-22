// complex def
struct Complex
{
    // complex definition
    string corename, core, partner, match;
    double zscore, score, energy, tms, clashes, deviation;
    int    aln_total, aln_same;
    bool   swap;
    string info;
        
    // molecules
    SpecMolecule target, mol;
    
    // construct
    Complex()
    {
        corename = info = core = partner = match = "";
        aln_total = aln_same = 0;
        tms = score = zscore = -LARGE;
        clashes = energy = deviation = LARGE;
        swap = false;    
    }
};
