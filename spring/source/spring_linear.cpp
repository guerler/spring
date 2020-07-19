// log on screen
#define LG_ON

// include
#include "library/interface.h"

// item
struct LinearItem
{
    // double
    double score, ref; 
    string tmpl, line;    
    
    // linear item
    LinearItem()
    {
        score = -LARGE;
        ref = 0.0;
        tmpl = "";
        line = "";
    }
};

// clean module
struct ModLinear
{
    static const double q = 0.3;
			
	// construct
	void construct(string flist, string fname)
	{
        // accepted ids
        Vec < string> lselect;
        Format::readList(flist, lselect);
        Hash < string, int  > hselect;
        for (int i = 0; i < lselect.size(); i++)
        {
            string key = lselect[i].substr(0, 4);
            if (hselect.get(key) == -1)
                hselect.insert(key, 1);
        }

        // features
        Vec < string > lline;
        Vec < string > ltemplate;
        Vec < double > lzscore;
        Vec < double > lenergy;
        Vec < double > ltmscore;      
        Vec < double > ldeviation;
        Vec < double > lquality;
        Vec < string > lid;
        
        // lines
        int n = 0;

        // feature file        
        File fl;
        fl.construct (fname.c_str());
        while (fl.good())
			if (fl.readLine())
				lline.push_back(fl.get());
		fl.close();
		
		// read details        
        fl.construct (fname.c_str());
        while (fl.good())
        {
            // read line
            if (!fl.read())
                continue;
            
            // feature
            double zscore = Convert::toDbl(fl.get(4));

            // energy
            double energy = Convert::toDbl(fl.get(5));

            // tmscore
            double tmscore = Convert::toDbl(fl.get(6));

            // clash
            double deviation = Convert::toDbl(fl.get(7));

			// quality
            double quality = Convert::toDbl(fl.get(9));
				
			/*/ check for redundancy
			bool found = false;
			for (int i = 0; i < n; i++)
			{
				if (fabs(lquality[i]-quality) > 0.05)
					continue;

				if (fabs(lfeature[i]-feature) > 5.0)
					continue;

				if (fabs(lenergy[i]-energy) > 5.0)
					continue;

				if (fabs(ltmscore[i]-tmscore) > 0.05)
					continue;

				if (fabs(lclash[i]-clash) > 0.05)
					continue;
					
				found = true;
				break;
			}
			if (found)
				continue;*/

			// add to list
            lzscore.push_back(zscore);
            lenergy.push_back(energy);
            ltmscore.push_back(tmscore);
            ldeviation.push_back(deviation);
            lquality.push_back(quality);
            
            // read pair identifier
            lid.push_back(fl.get(0).substr(0, 4) + " " + fl.get(1).substr(0, 4));

            // add
            ltemplate.push_back(fl.get(2) + " " + fl.get(3));
                   
            // increase counter
            n++;
        }
        fl.close();
        
        // verify
        if (n == 0)
            Msg::error("ModLinear::construct()", " Feature file not found or empty.");
        
		/*/ normalize
		Lib::normalize(lzscore);
		Lib::normalize(lenergy);
		Lib::normalize(ltmscore);		
		Lib::normalize(ldeviation);	*/
                
        /*/ file
        for (int i = 0; i < lzscore.size(); i++)
        {
    	    string qstr = "-1";
			if (lquality[i] >= q)
				qstr = "+1";
            cout << qstr << " 1:" << lzscore[i] << " 2:" << lenergy[i] << " 3:" << ltmscore[i] << " 4:" << ldeviation[i]  << endl;
		}
        exit(0);*/
   
        // hit count
        int nhits = 0;
        
        // counter
        int nentry = 0;
    
        // top selection
        int ntop = 5;
        
        // min quality
        double opt = 0.0;
        
        // data string
        string sout;
                
        // output
        int NMAX = 100000;
        Vec < double > lscore;
        Vec < double > lref;
        lscore.resize(NMAX);
        lref.resize(NMAX);
        
        // weighting
        double w0  = 9;// without deviation 12.0, 1.4, 1800;
        double w1  = 1.4;
        double w2  = 12.1;
        
//#define TRAIN
#ifdef TRAIN
        for (double w0 = 0.0; w0 < 15.0; w0 += 0.5)
        for (double w1 = 0.0; w1 < 1.5;  w1 += 0.05)
        for (double w2 = 0.0; w2 < 15.0; w2 += 0.1)
        //for (double w2 = 500.0; w2 < 2000.0; w2 += 20.0)
        {
#endif
            // reset
            sout    = "";
            nentry  = 0;
            nhits   = 0;
        
            // min list
            ListMin < LinearItem > models;
            models.construct(ntop);
            
            // loop
            string id = lid[0];   
            for (int i = 0; i < n; i++)
            {         
                // next target
                if (lid[i] != id)
                {
                    // get best model
                    double ref = 0.0;
                    double score = 0.0;
                    string tmpl = "";
                    string line = "";
                    for (int j = 0; j < models.size(); j++)
                        if (models.lscore[j] != LARGE)
                        {
                            if (models.litem[j].ref >= ref)
                            {
                                score = min(models.lscore[j], score);
                                ref = models.litem[j].ref;
                                tmpl = models.litem[j].tmpl;
                                line = models.litem[j].line;
                            }
                        }
                    
                    // add to results
                    lscore[nentry] = score;
                    lref[nentry] = ref;
                    nentry++;
     
                    // backup
                    //sout += id + " " + Convert::toString(score) + " " + Convert::toString(ref) + " " + tmpl + "\n";
    				sout += line + "\n";
    
                    // reset
                    models.initialize();
                        
                    // reset
                    id = lid[i];
                }
    
                // make item   
                LinearItem tmpitem;
                tmpitem.score = lzscore[i] + w0 * ltmscore[i] + w1 * lenergy[i] - w2  * ldeviation[i];
                tmpitem.ref = lquality[i];
                tmpitem.tmpl = ltemplate[i];
                tmpitem.line = lline[i];
                models.update(-tmpitem.score, tmpitem);
            }

#ifdef TRAIN
            // roc
            Lib::sort(lscore, lref, nentry);
                    
            // area under curve
            double average = 0.0;
            for (int i = 0; i < nentry; i++)            
            {
                // count
                if (lref[i] > q)
                    nhits++;
    
                // sum                    
                average += lref[i];
            }
                    
            // normalize
            average /= nentry;
    
            // get mcc
            double tt = nhits;
            double tf = nentry - nhits;
    
            // get mcc
            double tp = 0;
            double fp = 0;
            double maxmcc = 0.0;
            double maxmcc_score = 0.0;
            double maxmcc_tp = 0.0;
            for (int i = 0; i < nentry; i++)
            {
                // prediction
                if (lref[i] > q)
                    tp++;
                else
                    fp++;
    
                // mcc
                double fn = tt - tp;
                double tn = tf - fp;
                double mcc = 0.0;
                if ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0.0)
                    mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));                
                if (mcc > maxmcc)
                {
                    maxmcc = mcc;
                    maxmcc_score = lscore[i];
                    maxmcc_tp = tp;
                }
            }
    
            // maximize
            if (opt < maxmcc)
            {
                opt = maxmcc;
                cout << w0 << " " << w1 << " " << w2 << " " << nentry << " " << nhits << " " << nhits/double(nentry) << " " << average << " " << maxmcc_tp << " " << maxmcc_score << " " << maxmcc << endl;
            }
        }
#else
        cout << sout.c_str();
#endif
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
	cout << " *                            SPRING_LINEAR                               *" << endl;
	cout << " *                                                                        *" << endl;
	cout << " * Reference: A. Guerler et al.                                           *" << endl;
	cout << " * Comments on the program? Please contact: aysam.guerler@gmail.com       *" << endl;
	cout << " **************************************************************************" << endl;
	cout << endl;

	// run
    if (n < 3)
        Msg::error ("Missing information", "[target list] [feature file]");

    // check args
    ModLinear spec;
    spec.construct(input[1], input[2]);

	// return
	return 0;
}
