// file formats
class Format
{

public:
	// print		
	template < typename TList >
	static void writeList(TList& list, string fname)
	{
		// write file
		fstream pms;
		pms.open(fname.c_str(),ios_base::out);
		for (int x = 0; x < list.size(); x++)
			pms << list[x] << endl;
		pms.close();
	}	
	
	// print		
	template < typename TMatrix >
	static void writeMatrix(TMatrix& matrix, string fname)
	{
		// write file
		fstream pms;
		pms.open(fname.c_str(),ios_base::out);
		for (int x = 0; x < matrix.size(); x++)
		{
    			for (int y = 0; y < matrix[x].size(); y++)
				pms << setw(15) << matrix[x][y];
			pms << endl;
		}
		pms.close();
	}	
		
	// print		
	template < typename TCube >
	static void writeCube(TCube& matrix, string fname)
	{
		// verify
		if (matrix.size() == 0)
			Msg::error("Format::writeCube", "Cube is empty.");
		if (matrix[0].size() == 0)
			Msg::error("Format::writeCube", "Cube is empty.");
		
		// write file
		fstream pms;
		pms.open(fname.c_str(),ios_base::out);
		
		// write header
		pms << "DIMX " << matrix.size() << endl;
		pms << "DIMY " << matrix[0].size() << endl;
		pms << "DIMZ " << matrix[0][0].size() << endl;		
				
		// write data
		for (int x = 0; x < matrix.size(); x++)
    		for (int y = 0; y < matrix[x].size(); y++)
			{
				for (int z = 0; z < matrix[x][y].size(); z++)
			        pms << setw(15) << matrix[x][y][z];
                pms << endl;
            }
		pms.close();
	}	

	// print		
	static void writeString(string& str, string fname)
	{
		// write file
		fstream pms;
		pms.open(fname.c_str(),ios_base::out);
        pms << str;
        pms.close();
	}	

	// open list
	static void readList (string fname, Vec < string >& list)
	{
        // open stride
        File f;
		f.open(fname);
        
        // loop
		list.clear();
		while (f.good())
        {
            // read
			if (!f.read())
                continue;
            
            // merge
            for (int i = 0; i < min((int) f.size(), 2); i++)
                list.push_back(f.get(i));
        }			
        f.close();
	}

	// open list
	static void readList (string fname, Vec < int >& list)
	{
        // open stride
        File f;
		f.open(fname);
        
        // loop
		list.clear();
		while (f.good())
        {
            // read
			if (!f.read())
                continue;
            
            // merge
            for (int i = 0; i < f.size(); i++)
                list.push_back(Convert::toInt(f.get(i)));            
        }			
        f.close();
	}

	// open list
	static void readSequence (string fname, string& sequence)
	{
        // open stride
        File f;
		f.open(fname);
        
        // loop
		sequence = "";
		while (f.good())
        {
            // read
			if (!f.read())
                continue;
            
            // skip
            if (f.get(0).substr(0, 1) == ">")
                continue;
            
            // merge
            sequence += f.get(0);   
        }			
        f.close();
	}
	
	// open list
	static void readList (string fname, Vec < double >& list)
	{
        // open stride
        File f;
		f.open(fname);
        
        // loop
		list.clear();
		while (f.good())
        {
            // read
			if (!f.read())
                continue;
            
            // merge
            for (int i = 0; i < f.size(); i++)
                list.push_back(Convert::toDbl(f.get(i)));   
        }			
        f.close();
	}
	
	// open file to get list
	static void readPairs(string fname, Vec < string >& lista, Vec < string >& listb)
	{
		// locals
		string ida, idb;
		
		// reset
		lista.clear();
		listb.clear();		
		
		// open file to get list
        File f;
		f.open(fname);
		while (f.good())
        {
            // read
			f.read();
			
			// get
			ida = f.get(0);
			idb = f.get(1);
			
			// verify
			if (ida != "" && idb != "")
			{
            	lista.push_back(ida);
            	listb.push_back(idb);
			}
        }			
        f.close();
	}

	// read int matrix
	static void readMatrix (string fname, Vec < Vec <int> >& data)
	{						
        // string content
        Vec < Vec <string> > dstr;
        readTable (fname, dstr);

		// reset
		data.clear();

		// loop
		for (int i = 0; i < dstr.size(); i++)
		{
			// read
			Vec < int > line;
			for (int j = 0; j < dstr[i].size(); j++)
				line.push_back(Convert::toInt(dstr[i][j]));

			// write
			data.push_back(line);
		}
		
		// verify
		int n = data.size();
		for (int i = 0; i < n; i++)
			if (n != data[i].size())
			{
                Msg::write("File %s inconsistency.", fname.c_str());
				Msg::write("Length difference in row %i, %i instead of %i", i, data[i].size(), n);
				Msg::error("Format::readMatrix()", "Inconsistent matrix definition. Dimensions differ.");
			}
	}
	
	// read a matrix
	static void readMatrix (string fname, Vec < Vec <double> > & data)
    {
        // string content
        Vec < Vec <string> > dstr;
        readTable (fname, dstr);
        
		// reset
		data.clear();

		// loop
		for (int i = 0; i < dstr.size(); i++)
		{
			// read
			Vec < double > line;
			for (int j = 0; j < dstr[i].size(); j++)
				line.push_back(Convert::toDbl(dstr[i][j]));

			// write
			data.push_back(line);
		}
		
		// verify
		int n = data.size();
		for (int i = 0; i < n; i++)
			if (n != data[i].size())
			{
                Msg::write("File %s inconsistency.", fname.c_str());
				Msg::write("Length difference in row %i, %i instead of %i", i, data[i].size(), n);
				Msg::error("Format::readMatrix()", "Inconsistent matrix definition. Dimensions differ.");
			}
	}

	// read int matrix
	static string readTable (string fname, Vec < Vec <string> >& data)
	{
        // string tablename
        string name = "not_available";

		// makup data matrix
		File file;
		file.construct (fname);

		// reset
		data.clear();

		// loop
		while (file.good())
		{
			// read
			if (file.read())
			{
				// skip comments
				if (file.get(0) == "#")
				{
                    name = file.get(1);
					continue;
                }

				// read
				Vec < string > line;
				for (int i = 0; i < file.data.size(); i++)
					line.push_back(file.get(i));

				// write
				data.push_back(line);
			}
		}

		// close file
		file.close();

		// log
        if (data.size() == 0)
        {
            Msg::write("File %s not found.", fname.c_str());
            Msg::error("Format::readMatrix()", "File not found.");
        }
        
		// verify
		int n = data[0].size();
		for (int i = 0; i < data.size(); i++)
			if (n != data[i].size())
			{
                Msg::write("File %s inconsistency.", fname.c_str());
				Msg::write("Length difference in row %i, %i instead of %i", i, data[i].size(), n);
				Msg::error("Format::readMatrix()", "Inconsistent matrix definition. Dimensions differ.");
			}

        // retur name
        return name;
	}

	// read a cube
	template < typename TCube >
	static void readCube (string fname, TCube& data)
	{						
		// makup data matrix
		File file;
		file.construct (fname);

		// dimensions
		int dimx = 0, dimy = 0, dimz = 0;
		
		// read header
		while (file.good())
		{
			// read
			file.read();
			
			// skip comments
			if (file.get(0) == "#")
				continue;
						
			// dimensions
			if (file.get(0) == "DIMX")
				dimx = file.getInt(1);			
			if (file.get(0) == "DIMY")
				dimy = file.getInt(1);
			if (file.get(0) == "DIMZ")
				dimz = file.getInt(1);
						
			// verify
			if (dimx != 0 && dimy != 0 && dimz != 0)
				break;
		}
								
		// prepare matrix
		data.clear();
		data.resize(dimx);
		for (int x = 0; x < dimx; x++)
		{
			data[x].resize(dimy);
			for (int y = 0; y < dimy; y++)
			{
				data[x][y].resize(dimz);
				data[x][y].fill(0.0);
			}
		}
					
		// read header	
		int x = 0, y = 0, z = 0;
		while (file.good())
		{
			// read
			file.read();
			
			// skip comments
			if (file.get(0) == "#")
				continue;
						
			// read
			for (z = 0; z < dimz; z++)
			{
				// verify
				if (file.get(z) == "")
				{
					Msg::write("Invalid matrix entry at %i, %i, %i", x, y, z);
					Msg::error("Format::readCube()", "Invalid matrix entry"); 
				}
				
				// read data
				data[x][y][z] = file.getDbl(z);
			}
			
			// update counter
			y++;
			if (!(y < dimy))
			{
				y = 0;
				x++;				
				if (!(x < dimx))
					break;
			}
		}
				
		// close file
		file.close();
		
		// log
		Msg::write("Cube loaded from %s [%i]", fname.c_str(), dimx*dimy*dimz);		
	}	
};

