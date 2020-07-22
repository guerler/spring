// fasta module
struct SequenceAlignment
{
	// scoring matrix
	Vec < Vec < int > > subst;
	
	// dynamical programming matrix
	Vec < Vec < int > > d;

	// automatic
	SequenceAlignment()
	{
		construct();
	}

	// construct
	void construct (string matrix = "data/matblosum.txt")
	{ 
		Format::readMatrix(matrix, subst);
	}

	// get score for a single aminoacid pair
	int getScore (int rescodea, int rescodeb)
	{
		return subst[rescodea][rescodeb];
	}

	// get alignment score for two molecules
	int getScore (SpecMolecule* mola, SpecMolecule* molb)
	{
		// read sequence a
		Vec < int > vseqa;
		vseqa.resize(mola->lcalpha.size());
		for (int i = 0; i < vseqa.size(); i++)
			vseqa[i] = mola->latom[mola->lcalpha[i]].rescode;
            
		// read sequence b
		Vec < int > vseqb;
		vseqb.resize(molb->lcalpha.size());
		for (int i = 0; i < vseqb.size(); i++)
			vseqb[i] = molb->latom[molb->lcalpha[i]].rescode;
            
		// return score
		return getScore(vseqa, vseqb);
    }

	// sequence alignment score
	int getScore (string seqa, string seqb, int gap = 5)
	{
		// read sequence a
		Vec < int > vseqa;
		vseqa.resize(seqa.size());
		for (int i = 0; i < seqa.size(); i++)
			vseqa[i] = SpecDetails::getResCode(seqa[i]);
            
		// read sequence b
		Vec < int > vseqb;
		vseqb.resize(seqb.size());
		for (int i = 0; i < seqb.size(); i++)
			vseqb[i] = SpecDetails::getResCode(seqb[i]);
            
		// return score
		return getScore(vseqa, vseqb, gap);
    }
   	
	// sequence alignment score
	int getScore (Vec <int>& vseqa, Vec<int>& vseqb, int gap = 5)
	{	
		// loop
		int n = vseqa.size();
		int m = vseqb.size();
                
		// verify
		if (m == 0 || n == 0)
			return 0;

		// resize
		n++;
		m++;
		
		// setup sizes
		d.resize(n);
		for (int i = 0; i < n; i++)
			d[i].resize(m);

		// initialize values
		for (int i = 0; i < n; i++)
			d[i][0] = 0;		
		for (int j = 0; j < m; j++)
			d[0][j] = - j * gap;		
		
		// generate matrix
		// set left gap, upper gap, match/mismatch
		for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)		
			d[i][j] = Lib::maximum ( d[i-1][j-1] + subst[vseqa[i-1]][vseqb[j-1]], d[i-1][j], d[i][j-1] - gap);

		// return score
		return d[n-1][m-1];
	}
	
	// reduce to shared residues
	void reduceSharedResidues(SpecMolecule* mola, SpecMolecule* molb)
	{	
		// determine residues
		Vec < int > lmola, lmolb;
		getSharedResidues(mola, molb, lmola, lmolb);

		// log
		Msg::write("%i shared residues out of %i and %i atoms detected.", lmola.size(), mola->size(), molb->size());

		// recopy mola
		SpecMolecule m0;
		m0.construct(mola->name);
		for (int i = lmola.size() - 1; i >= 0; i--)
			m0.append (mola->latom[lmola[i]]);
		m0.finalize();
		mola->construct(&m0);

		// recopy molb
		m0.construct(molb->name);
		for (int i = lmolb.size() - 1; i >= 0; i--)
			m0.append (molb->latom[lmolb[i]]);
		m0.finalize();
		molb->construct(&m0);
	}

	// sequence alignment score
	void getSharedResidues (SpecMolecule* mola, SpecMolecule* molb, Vec < int >& lmola, Vec < int >& lmolb)
	{
		// read sequence a
		Vec < int > vseqa;
		vseqa.resize(mola->lcalpha.size());
		for (int i = 0; i < vseqa.size(); i++)
			vseqa[i] = mola->latom[mola->lcalpha[i]].rescode;
        
		// read sequence b
		Vec < int > vseqb;
		vseqb.resize(molb->lcalpha.size());
		for (int i = 0; i < vseqb.size(); i++)
			vseqb[i] = molb->latom[molb->lcalpha[i]].rescode;
            
		// return score
		getSharedResidues(vseqa, vseqb, lmola, lmolb);
	}
    
	// sequence alignment score
	void getSharedResidues (Vec <int>& vseqa, Vec<int>& vseqb, Vec < int >& lmola, Vec < int >& lmolb)
	{
		// loop
		int n = vseqa.size();
		int m = vseqb.size();
                
		// verify
		if (m == 0 || n == 0)
			return;

		// resize
		n++;
		m++;
		
		// setup sizes
		d.resize(n);
		for (int i = 0; i < n; i++)
			d[i].resize(m);

		// initialize values
		for (int i = 0; i < n; i++)
			d[i][0] = 0;		
		for (int j = 0; j < m; j++)
			d[0][j] = 0;		
		
		// generate matrix
		// set left gap, upper gap, match/mismatch
		for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)	
			d[i][j] = Lib::maximum ( d[i-1][j-1] + ((vseqa[i-1] == vseqb[j-1])? 1 : 0),  d[i-1][j], d[i][j-1]);

		// index
		int i = n - 1;
		int j = m - 1;

		// retrieve shared residues
		int a, b, c, max;	
		while (i > 0 && j > 0)
		{	
			// read values		
			a = d[i-1][j-1];
			b = d[i-1][j];
			c = d[i][j-1];
			
			// read left gap, upper gap, match/mismatch
			max = Lib::maximum (a, b, c);			
			if (max == a)
			{
				lmola.push_back(--i);
				lmolb.push_back(--j);
			} else {
				if (max == b)
					i--;
				else
					j--;
			}				
		}
	}

	// copy model
	bool makemodel(SpecMolecule* model, string sequencefile)
	{
        // read sequence
        string sequence = "";
        Format::readSequence(sequencefile, sequence);

        // read sequence
        return makemodel_sequence (model, sequence);
    }

	// copy model
	bool makemodel_sequence (SpecMolecule* model, string sequence)
	{
        // sequence empty
        if (sequence == "")
            Msg::error("SequenceAlignment::makemodel_sequence()", "Sequence not found or empty.");
        
        // read sequence
        SpecMolecule mol;
        mol.readSequence(sequence);

        // perform sequence alignment
        Vec < int > lmodel, lmol;
        getSharedResidues(&mol, model, lmol, lmodel);
      
        // transfer residues
        for (int i = 0; i < lmol.size(); i++)
            mol.latom[mol.lcalpha[lmol[i]]].pos = model->latom[model->lcalpha[lmodel[i]]].pos;
        
        // construct
        model->construct(&mol);

        /*/ check
        if (model->lcalpha.size() != lmodel.size())
            Msg::write ("SequenceAlignment::makemodel() : Model/Template alignment inconsistency! (%i, %i)", model->lcalpha.size(), lmodel.size());*/

        // check
        if (model->size() != lmodel.size())
            return false;
        
        // renew model
        return true;
    }	
};
