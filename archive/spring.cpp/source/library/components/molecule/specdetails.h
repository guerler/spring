struct SpecDetails
{    	    	
	// types
	static const int AMINO = 20;
	static const int ATOMS = 37;

	// setup atom
	inline static int getResCode (string res)
	{
		// set amino acid code sorted by single code
		if (res == "ALA")			return 0; 	// 'A'
		else if (res == "CYS")		return 1; 	// 'C'
		else if (res == "ASP")		return 2; 	// 'D'
		else if (res == "GLU")		return 3; 	// 'E'
		else if (res == "PHE")		return 4; 	// 'F'
		else if (res == "GLY")		return 5; 	// 'G'
		else if (res == "HIS")		return 6; 	// 'H'
		else if (res == "ILE")		return 7; 	// 'I'
		else if (res == "LYS")		return 8; 	// 'K'
		else if (res == "LEU")		return 9; 	// 'L'
		else if (res == "MET")		return 10; 	// 'M'
		else if (res == "ASN")		return 11; 	// 'N'
		else if (res == "PRO")		return 12; 	// 'P'
		else if (res == "GLN")		return 13; 	// 'Q'
		else if (res == "ARG")		return 14; 	// 'R'
		else if (res == "SER")		return 15; 	// 'S'
		else if (res == "THR")		return 16; 	// 'T'
		else if (res == "VAL")		return 17; 	// 'V'
		else if (res == "TRP")		return 18; 	// 'W'
		else if (res == "TYR")		return 19; 	// 'Y'
		return 20;
	}
    
	// setup atom
	inline static int getResCode (char res)
	{
		// set amino acid code sorted by single code
		if (res == 'A')				return 0; 	// 'A'
		else if (res == 'C')		return 1; 	// 'C'
		else if (res == 'D')		return 2; 	// 'D'
		else if (res == 'E')		return 3; 	// 'E'
		else if (res == 'F')		return 4; 	// 'F'
		else if (res == 'G')		return 5; 	// 'G'
		else if (res == 'H')		return 6; 	// 'H'
		else if (res == 'I')		return 7; 	// 'I'
		else if (res == 'K')		return 8; 	// 'K'
		else if (res == 'L')		return 9; 	// 'L'
		else if (res == 'M')		return 10; 	// 'M'
		else if (res == 'N')		return 11; 	// 'N'
		else if (res == 'P')		return 12; 	// 'P'
		else if (res == 'Q')		return 13; 	// 'Q'
		else if (res == 'R')		return 14; 	// 'R'
		else if (res == 'S')		return 15; 	// 'S'
		else if (res == 'T')		return 16; 	// 'T'
		else if (res == 'V')		return 17; 	// 'V'
		else if (res == 'W')		return 18; 	// 'W'
		else if (res == 'Y')		return 19; 	// 'Y'
		return 20;
	}
    
	// setup atom
	inline static string getResName (char res)
	{
		// set amino acid code sorted by single code
		if (res == 'A')				return "ALA"; 	// 'A'
		else if (res == 'C')		return "CYS"; 	// 'C'
		else if (res == 'D')		return "ASP"; 	// 'D'
		else if (res == 'E')		return "GLU"; 	// 'E'
		else if (res == 'F')		return "PHE"; 	// 'F'
		else if (res == 'G')		return "GLY"; 	// 'G'
		else if (res == 'H')		return "HIS"; 	// 'H'
		else if (res == 'I')		return "ILE"; 	// 'I'
		else if (res == 'K')		return "LYS"; 	// 'K'
		else if (res == 'L')		return "LEU"; 	// 'L'
		else if (res == 'M')		return "MET"; 	// 'M'
		else if (res == 'N')		return "ASN"; 	// 'N'
		else if (res == 'P')		return "PRO"; 	// 'P'
		else if (res == 'Q')		return "GLN"; 	// 'Q'
		else if (res == 'R')		return "ARG"; 	// 'R'
		else if (res == 'S')		return "SER"; 	// 'S'
		else if (res == 'T')		return "THR"; 	// 'T'
		else if (res == 'V')		return "VAL"; 	// 'V'
		else if (res == 'W')		return "TRP"; 	// 'W'
		else if (res == 'Y')		return "TYR"; 	// 'Y'
		return "   ";
	}
        
	// setup atom
	inline static char getCharCode(int rescode)
	{
		// set amino acid code sorted by single code
		if (rescode == 0)			return 'A';
		else if (rescode == 1)		return 'C';
		else if (rescode == 2)		return 'D';
		else if (rescode == 3)		return 'E';
		else if (rescode == 4)		return 'F';
		else if (rescode == 5)		return 'G';
		else if (rescode == 6)		return 'H';
		else if (rescode == 7)		return 'I';
		else if (rescode == 8)		return 'K';
		else if (rescode == 9)		return 'L';
		else if (rescode == 10)		return 'M';
		else if (rescode == 11)		return 'N';
		else if (rescode == 12)		return 'P';
		else if (rescode == 13)		return 'Q';
		else if (rescode == 14)		return 'R';
		else if (rescode == 15)		return 'S';
		else if (rescode == 16)		return 'T';
		else if (rescode == 17)		return 'V';
		else if (rescode == 18)		return 'W';
		else if (rescode == 19)		return 'Y';
        else                        return 'X';
	}              
	  
};

