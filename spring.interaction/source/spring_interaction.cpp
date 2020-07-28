// log on screen
#define LG_ON

// include
#include "library/interface.h"

// prefilter interacting pairs
struct ModSpring_Interaction
{
    // output file
    fstream ftrain;

    // complex reader
    Act_ComplexReader cr;

    // threading results
    ThreadingReference tr;

    // construct
    void construct (string target_index_str, string ftarget, string targetpath, string pdbpath)
    {
        // construct molecules
        Vec < string > lchains;

        // read list
        Format::readList (ftarget, lchains);

        // get index
        int target_index = Convert::toInt (target_index_str);

        // check
        string target_name = "";
        if (target_index >= 0 && target_index < lchains.size())
            target_name = lchains[target_index];
        else
            Msg::error("ModSpring_Interaction::construct()", "Invalid target index.");

        // report
        string opath = targetpath + "/interactions/" + target_name.substr(0, 2) + "/";
        string fname = opath + target_name;

        // check if output exists already
        if (File::exists(fname))
            Msg::error ("ModSpring_Interaction::construct()", "Output file already exists.");

        // make directory
        Lib::makedir(targetpath + "/interactions/");
        Lib::makedir(opath);

        // threading reference
        string threading_reference = targetpath + "/summary.txt";

        // validate that threading reference exists
        if (!File::exists(threading_reference))
        {
            Msg::write ("Generating threading reference...please wait.");

            // open file
            fstream fref;
            fref.open(threading_reference.c_str(), ios_base::out);

            // open report file
            HHsearch pptarget;
            for (int i = 0; i < lchains.size(); i++)
            {
                // target
                if(!pptarget.construct(lchains[i], targetpath, SpringConfig::MINZ, SpringConfig::SEQCUTOFF))
                    continue;

                // get scores and template names
                Vec < string > ltarget = pptarget.ltemplate;
                Vec < double > ltargetscore = pptarget.lscore;

                // write
                for (int j = 0; j < ltarget.size(); j++)
                    fref << lchains[i] << " " << ltarget[j] << " " << ltargetscore[j] << endl;

                // rewrite
                Msg::rewrite ("Processing %i of %i.", i, (int) lchains.size());
            }

            // close file
            fref.close();

            // done
            Msg::write ("Threading reference written to %s.", threading_reference.c_str());
        }

        // threading reference
        tr.construct(threading_reference);

        // complex reader
        cr.construct(pdbpath + "index.txt");

        // open file
        ftrain.open(fname.c_str(), ios_base::out);

        // loop
        for (int j = target_index; j < lchains.size(); j++)
            construct_pair (target_name, lchains[j]);

        // close
        ftrain << "DONE";
        ftrain.close();

        // log
        Msg::write ("Done.");
    }

    // construct
    bool construct_pair(string ida, string idb)
    {
        // read
        Vec < ThreadingInfo > linfoa;
        if (!tr.get (ida, linfoa))
            return false;

        // read
        Vec < ThreadingInfo > linfob;
        if (!tr.get (idb, linfob))
            return false;

        // reset complex list
        cr.initialize();

        // generate complexes using description from pdb files
        for (int i = 0; i < linfoa.size(); i++)
            cr.add (linfoa[i], linfob);

        // generate complexes using description from pdb files
        for (int i = 0; i < linfob.size(); i++)
        {
            // search
            bool found = false;
            for (int j = 0; j < linfoa.size(); j++)
                if (linfoa[j].tmpl == linfob[i].tmpl)
                {
                    found = true;
                    break;
                }

            // add
            if (!found)
                cr.add (linfob[i], linfoa);
        }

        // show top
        if (cr.maxscore != -LARGE)
            ftrain << ida + " " + idb + " " + Convert::toString(cr.maxscore) << endl;

        // return
        return true;
    }
};

// main 
int main(int n, char* arg[]) 
{ 
	// parameters
	Vec < string > input;

	// load parameters
	for (int i = 0; i < n; i++)
		input.push_back((string) arg[i]);
		
	// result 
	cout << " **************************************************************************" << endl; 
	cout << " *                          SPRING_INTERACTION                            *" << endl;
	cout << " *                                                                        *" << endl; 
	cout << " * Reference: A. Guerler et al.                                           *" << endl;
	cout << " * Comments on the program? Please contact: aysam.guerler@gmail.com       *" << endl;
	cout << " **************************************************************************" << endl; 
	cout << endl; 

	// run 
    if (n < 5)
        Msg::error ("Missing information", "[target index] [file list] [target path] [pdb path]");

    // check args
    ModSpring_Interaction spec;
    spec.construct(input[1], input[2], input[3], input[4]);

	// return 
	return 0; 
} 
