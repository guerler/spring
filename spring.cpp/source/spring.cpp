// log on screen
#define LG_ON

// include
#include "library/interface.h"

// spring specific
#include "springconfig.h"

// energy function
#include "spring_library/energycontacts.h"

// refinement strategy
#include "spring_refinement/refinement.h"

// spring modeling includes
#include "spring_library/complex.h"
#include "spring_library/structuralneighbor.h"
#include "spring_library/crossreference.h"
#include "spring_library/complexreader.h"
#include "spring_library/mergemol.h"

// multimeric threading with spring
struct ModSpring_Model
{
    // storage
    Storage store;

    // output file
    fstream ftrain;

    // complex reader
    ComplexReader cr;
    
    // construct
    void construct (string ftarget, string modelpath, string pdbpath, string output)
    {
        // make sure output is a directory
        output += "/";

        // complex reader
        cr.construct(pdbpath);

        // make directory
        Lib::makedir(output);

        // report
        string fname = FName::getname(ftarget);

		// open file
		Msg::write ("Evaluating %s.", string (output + fname).c_str());

        // open file
        ftrain.open(string (output + fname + ".results").c_str(), ios_base::out);

        // if details are not available then print header
        if (!SpringConfig::DETAILS)
            ftrain << "#query_a query_b template_a template_b score " << InterfaceQuality_Info::getHeader() << endl;

        // construct molecules
        Vec < string > lchains;

        // read list
        Format::readList (ftarget, lchains);

        // open report file
        for (int m = 0; m < lchains.size(); m+=2)
            construct_pair (lchains[m], lchains[m+1], modelpath, pdbpath, output);

        // close
        ftrain << "DONE";
        ftrain.close();
    }

    // make full-atom model
    void getsidechains (SpecMolecule* target, SpecMolecule* mol, SpecMolecule* targetfull = 0, SpecMolecule* molfull = 0)
    {
        // replace models with full models
        bool done = false;
        if (targetfull != 0 && molfull != 0)
        {
            if(targetfull->latom.size() > 0 && molfull->latom.size() > 0)
            {
                // get side chains from full chains
                TMAlign::align(targetfull, target);
                TMAlign::align(molfull, mol);

                // replace
                target->construct(targetfull);
                mol->construct(molfull);
                
                // set flag
                done = true;
            }
        }

        if (!done)
        {
            // build side chains
            Pulchra::make (target);
   	        Pulchra::make (mol);
        }
    }

    // construct
    bool construct_pair(string ida, string idb, string modelpath, string pdbpath, string output)
    {
        // new line
        Msg::next();

		// log
        Msg::write ("Targets %s, %s.", ida.c_str(), idb.c_str());

        /**
            GET THREADING FILES
        **/

        //
        // threading results
        //
        HHsearch pptarget, ppmol;

        // target
        if(!pptarget.construct(ida, modelpath, SpringConfig::MINZ, SpringConfig::SEQCUTOFF))
            return false;

        // molecule
        if(!ppmol.construct(idb, modelpath, SpringConfig::MINZ, SpringConfig::SEQCUTOFF))
            return false;

        // double check
        if (pptarget.ltemplate.size() == 0 || ppmol.ltemplate.size() == 0)
        {
            Msg::write("Target threading results incomplete.");
            return false;
        }

        // get names of monomeric templates
        string targetmononame = pptarget.ltemplate[0];
        string molmononame    = ppmol.ltemplate[0];

        // get template names
        Vec < string > ltarget = pptarget.ltemplate;
        Vec < double > ltargetscore = pptarget.lscore;
        Vec < string > lmol = ppmol.ltemplate;
        Vec < double > lmolscore = ppmol.lscore;

        /**
            CONSTRUCT REFERENCE AND MODEL
        **/
        // make model
        Msg::write ("Accessing template repository at %s.", string(pdbpath).c_str());

        // sequence alignment tool
        SequenceAlignment sa;

        // check for suggested monomers
        SpecMolecule targetmono, molmono;
        SpecMolecule targetfullmono, molfullmono;
        bool withmono = false;
        if (SpringConfig::SUGGESTEDMONO && store.read(&targetfullmono, modelpath + "monomers/" + ida, ida) &&
            store.read(&molfullmono, modelpath + "monomers/"  + idb, idb))
        {
            // both loaded
            withmono = true;

            // align model to sequence
            sa.makemodel(&targetfullmono, modelpath + "fasta/" + ida);
            sa.makemodel(&molfullmono, modelpath + "fasta/" + idb);

            // copy calphas
            targetfullmono.copyCAlpha(&targetmono);
            molfullmono.copyCAlpha(&molmono);
            
            // log
            Msg::write ("Suggested monomeric structures available.");
        } else {
        	// check
            if (!pptarget.getmodel (&targetmono, targetmononame, pdbpath))
            {
                Msg::write ("Requested monomeric template for %s not found.", targetmononame.c_str());
                return false;
            }

            // check
            if (!ppmol.getmodel (&molmono, molmononame, pdbpath))
            {
                Msg::write ("Requested monomeric template for %s not found.", molmononame.c_str());
                return false;
            }

            // align model to sequence
            //sa.makemodel(&targetmono, modelpath + "fasta/" + ida);
            //sa.makemodel(&molmono, modelpath + "fasta/" + idb);

            // log
            Msg::write ("Suggested monomeric structures not available.");
        }

        // update names
        targetmono.name = ida;
        molmono.name    = idb;
        
        // construct reference molecules
        SpecMolecule targetref, molref;
        bool withref = SpringConfig::WITHREF;
        if (withref)
        {
	        // construct reference molecules
	        if (!store.read(&targetref, modelpath + "chains/" + ida, ida))
    	        withref = false;

	        // construct reference molecules
            if (!store.read(&molref, modelpath + "chains/" + idb, idb))
        	    withref = false;

            // proceeding without reference
            if (withref)
            {
                // align model to sequence
                sa.makemodel(&targetref, modelpath + "fasta/" + ida);
                sa.makemodel(&molref, modelpath + "fasta/" + idb);

                // check
                if (targetmono.lcalpha.size() != targetref.lcalpha.size() || molmono.lcalpha.size() != molref.lcalpha.size())
                {
            	    Msg::error ("ModSpring_Model::construct()", "Size of reference and target differs.");
                	return false;
                }
            } else
                Msg::write ("Proceeding without reference complex.");
		}

        /**
            EVALUATE INTERACTION
        **/

        // log
        Msg::write ("Monomeric templates %s, %s.", targetmononame.c_str(), molmononame.c_str());
        Msg::write ("Template names loaded. %i, %i.", ltarget.size(), lmol.size());

        // reset complex list
        cr.initialize(&targetmono, &molmono);

        // identify templates by min z-score
        for (int i = 0; i < min((int) SpringConfig::NTMPL_PER_TARGET, ltarget.size()); i++)
            cr.add (ltarget[i], ltargetscore[i], lmol, lmolscore);

        // generate complexes using description from pdb files
        for (int i = 0; i < min((int) SpringConfig::NTMPL_PER_TARGET, lmol.size()); i++)
        {
            // search
            bool found = false;
            for (int j = 0; j < ltarget.size(); j++)
                if (ltarget[j] == lmol[i])
                {
                    found = true;
                    break;
                }

            // add
            if (!found) {
                cr.add (lmol[i], lmolscore[i], ltarget, ltargetscore, true);
            }
        }

        // check
        Msg::write("Total number of identified templates is %i.", cr.lcomplex.size());

        // return
        if (cr.lcomplex.size() == 0)
            return false;

        // sort templates by score
        cr.makemodels();

        /**
            SUMMARIZE AND EVALUATE RESULTS
        **/

        // get top ranks
        Vec < Complex > lrank;

        // check
        int nredundant = 0;
        for (int i = 0; i < cr.lcomplex.size(); i++)
        {
            // read complex
            Complex* clx = &cr.lcomplex[i];

            // done in sorted list
            if (clx->score == -LARGE || lrank.size() >= SpringConfig::NTOP)
                break;

            // get details
            string name  = clx->corename.substr(0, 4);
            double score = clx->score;

            // fix modelsizes
            sa.makemodel(&clx->mol, modelpath + "fasta/" + idb);
            
            // placed
            bool done = false;

            // first placement option
            for (int j = 0; j < lrank.size(); j++)
            {
                // search lower scoring structure of same template
                if (lrank[j].corename.substr(0, 4) == name)
                {
                    // replace if score is higher
                    if (score > lrank[j].score)
                        lrank[j] = cr.lcomplex[i];

                    // placement completed
                    done = true;
                    break;
                }
            }

            // second placement option
            if (!done)
            {
                // load rmsd calculator
                IRmsd irmsd;
                irmsd.construct(&clx->target, &clx->mol);

                // loop
                for (int j = 0; j < lrank.size(); j++)
                {
                    // search similar structure
                    irmsd.initialize(&lrank[j].target, &lrank[j].mol);
                    if (irmsd.getrmsd() < SpringConfig::MINRMSD)
                    {
                        nredundant++;
                        done = true;
                        break;
                    }
                }
            }

            // add to list
            if(!done)
                lrank.push_back(cr.lcomplex[i]);
        }

        // sort list
        std::sort (lrank.begin(), lrank.end(), ComplexReader::sortbyspringscore);

        // models clustered
        Msg::write("Model clustering and sorting completed. Found %i redundant models at %2.2fA RMSD.", nredundant, SpringConfig::MINRMSD);

        // info
        InterfaceQuality_Info top_quality;

        // table
        string outtable   = "";
        string outdetails = "";

        // loop models
        int nmodel = 0;
        for (int i = 0; i < lrank.size(); i++)
        {
            if (lrank[i].score != -LARGE)
            {
                // read complex
                Complex* clx = &lrank[i];

                // count model
                nmodel++;

                // filename
                string fname = ida + "_" +  idb + "_" + Convert::toString(nmodel) + ".pdb";

                // template identifier
                string tmplnamecore = clx->corename + "/" + clx->core;
                string tmplnamepartner = clx->corename + "/" + clx->partner;
                
                // check minimum
                //if (clx->score < SpringConfig::MINSPRING)
                //	continue;

                // write models
                Msg::write ("Score %f.", clx->score);

                // trajectory
                if (SpringConfig::WRITE)
                {
                    // make side chains
                    getsidechains(&clx->target, &clx->mol, &targetfullmono, &molfullmono);

                    // check for clashes
                    double clash_radius = sqrt(SpringConfig::RADIUS2);
                    for (int m = 0; m < clx->target.lcalpha.size(); m++)
                    for (int n = 0; n < clx->mol.lcalpha.size(); n++)
                    {
                        if (clx->target.latom[clx->target.lcalpha[m]].pos.dist(clx->mol.latom[clx->mol.lcalpha[n]].pos) < clash_radius)
                        {
                            if (clx->target.lcalpha.size() > clx->mol.lcalpha.size())
                                clx->target.latom[clx->target.lcalpha[m]].pos = NOCOORD;
                            else
                                clx->mol.latom[clx->mol.lcalpha[n]].pos = NOCOORD;
                        }
                    }

                    // alignment details
                    outdetails += "\nNo. " + Convert::toString(nmodel) + "\n" + clx->info;
                    outdetails += "----------------------------------------------------------------------------------------------------\n";

                    // finalize
                    clx->target.info = "REM spring " + Convert::toString(clx->score) + "\n";
                    clx->target.info += "REM templates " + tmplnamecore + " " + tmplnamepartner + "\n";
                    store.save (&clx->target, &clx->mol, output + fname);
                    Msg::write ("Predicted model %i stored.", nmodel);
                }

                // write table
                if(SpringConfig::DETAILS)
                    sprintf (strbuf, "%-5i %-15s %-15s %10.1f %10.1f %10.2f %10i %9i%%\n", nmodel, tmplnamecore.c_str(), tmplnamepartner.c_str(), clx->score, clx->energy, clx->tms, clx->aln_total, Lib::percent(clx->aln_same, clx->aln_total));
                else
                    sprintf (strbuf, "%-5i %-15s %-15s %-15s %-15s %10.1f %10.1f %10.2f %10i %9i%%\n", nmodel, ida.c_str(), idb.c_str(),tmplnamecore.c_str(), tmplnamepartner.c_str(), clx->score, clx->energy, clx->tms, clx->aln_total, Lib::percent(clx->aln_same, clx->aln_total));
                outtable += strbuf;
            }
        }

        // verify
        if (outtable == "")
        {
            Msg::write ("No acceptable prediction made.");
            return false;
        }

        // write header
        if (SpringConfig::DETAILS && outtable != "" && !withref)
        {
            string outheader = "";
            outheader += "----------------------------------------------------------------------------------------------------\n";
            outheader += "SPRING - [S]ingle-chain based [PR]ediction of [IN]teractions and [G]eometries\n";
            outheader += "\n";
            outheader += "Reference : Mapping monomeric threading to protein-protein structure prediction\n";
            outheader += "            Guerler, Govindarajoo, Zhang, Journal of Chemical Information and Modeling\n";
            outheader += "\n";
            outheader += "Please contact aysam.guerler@gmail.com for questions and comments.";
            outheader += "\n----------------------------------------------------------------------------------------------------\n";
            outheader += "Query (with monomeric template) " + ida + " (" + targetmononame + ") / " + idb + " (" + molmononame + ")";
            outheader += "\n----------------------------------------------------------------------------------------------------\n";
            sprintf (strbuf, "%-5s %-15s %-15s %10s %10s %10s %10s %10s\n", "No.", "Template_a", "Template_b", "Score", "Energy", "TM-score", "Aligned", "Ident.");
            outheader += strbuf;
            ftrain << outheader;
        }

        // write features
        ftrain << outtable;

        // check details
        if (SpringConfig::DETAILS && outtable != "" && !withref)
        {
            // separate table from details
            string outlegend = "";
            outlegend += "----------------------------------------------------------------------------------------------------\n";
            outlegend += "Score    = SPRING-score (confident predictions usually have a score above 10)\n";
            outlegend += "Energy   = Contact-based interaction potential\n";
            outlegend += "TM-score = TM-score of target model to the interface region of the dimeric template\n";
            outlegend += "Aligned  = Number of aligned residues\n";
            outlegend += "Identity = Fraction of identical residues\n";
            outlegend += "----------------------------------------------------------------------------------------------------\n";
            outlegend += "Below, the alignments of the individual monomer models (of query a and b) to the interfacial regions\n";
            outlegend += "of the dimer templates are shown. The interfacial regions are named by the chain, which was used to \n";
            outlegend += "identify the template, followed by a SPRING specific interface identifier.\n";
            outlegend += "----------------------------------------------------------------------------------------------------\n";

            // write
            ftrain << outlegend;
            ftrain << outdetails;

        }
        ftrain << endl;

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
	cout << " *                                SPRING                                  *" << endl;
	cout << " *                                                                        *" << endl; 
	cout << " * Reference: A. Guerler et al.                                           *" << endl;
	cout << " * Comments on the program? Please contact: aysam.guerler@gmail.com       *" << endl;
	cout << " **************************************************************************" << endl; 
	cout << endl; 
	
	// run 
    if (n < 5)
        Msg::error ("Missing information", "[target list] [target repository] [pdb repository] [output directory]");

    // check args
    ModSpring_Model spec;
    spec.construct(input[1], input[2], input[3], input[4]);

	// return 
	return 0; 
} 
