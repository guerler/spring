// pymol script writer
class PyMol
{
	// pymol
	fstream pms;
	string path;
	string fname;
public:
	// colors: lightblue, density, lightorange, orange
	// constructor
	void construct(string fname)
	{
		// write file
		string fn = string(FName::reduce(fname) + ".pml");
		pms.open(fn.c_str(),ios_base::out);		
		pms << "bg_color white" << endl;
	}

    void load (string chain_file)
    {
        pms << "load " << chain_file << endl;
    }
	
    // write all
	void add (SpecMolecule* mol, char chain_id, string chain_color)
	{
		// add
		pms << "color " << chain_color << ", chain " << chain_id << endl;
        		
		// sse
		for (int i = 0; i < mol->lsse.size(); i++)
		{
			int s = mol->lsse[i].getStart();
			int e = mol->lsse[i].getEnd();			
			sse(s, e-1, mol->lsse[i].getType(), chain_id);
		}
	}

	// close file
	void close()
	{
		pms << "hide everything" << endl;
		pms << "show cartoon" << endl; 
		pms.close();
    }
    
	// colors
	inline static string getcolor(int i)
	{
        // switch
		static string colors[] = {"yellow","marine","green","lightmagenta","cyan","density","red",
					"forest","deepolive","lime","sand","deeppurple","palecyan","orange"};
		i = i % 13;
		return colors[i];
	}    
private:	
	
	// write residue
	void residue(int start, int end, char chain, string item)
	{
		pms << "show cartoon, resi " << start << "-" << end << " AND chain " << chain << endl;
		pms << item << ", resi " << start << "-" << end << " AND chain " << chain << endl;
	}

	// write residue
	void residue(int residue, char chain, string item)
	{
		pms << "show cartoon, resi " << residue << " AND chain " << chain << endl;
		pms << item << ", resi " << residue << " AND chain " << chain << endl;
	}

	// write sse	
	void sse(int start, int end, char type, char chain)
	{
		pms << "alter " << chain << "/" << min (start, end) << "-" << max (start, end) << "/, ss='" << type << "'" << endl;
	}	
};

