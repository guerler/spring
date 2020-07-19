/**
	Minimization Interface
 */
template < typename TInput >
class IStep
{
	public :
		// weight or step size
		double weight;

		// interval
		double wMin, wMax;

		// weights
		Vec <double> w;

		// pre forces
		Vec <TInput> pre;

		// rho
		double rhoplus, rhominus;

		// step size in percentage of weight
		void construct (double weight, int size = 1, double max = 2.0, double min = EPSILON, double rhoplus = 1.5, double rhominus = 0.5)
		{
			// init params
			this->weight	= weight;
			this->wMin		= min;
			this->wMax		= max;

			// adaption parameters
			this->rhoplus	= rhoplus;
			this->rhominus	= rhominus;

			// resize arrays
			w.resize(size);
			pre.resize(size);
		}

		// reset values
		void initialize ()
		{
			w.fill(weight);
			pre.fill(TInput (0.0));
		}
};

//	broyden-fletcher-goldfarb-shanno (bfgs) method
class StepBFGS
{
	// parameters
	double grad, fLength;
	
	// parent
	IStep < Vec3 > p;
public:
	
	void construct(double weight)
	{
		p.construct(weight);
	}
	
	void initialize()
	{
		p.initialize();
	}
	
	inline double adjust(Vec3& f)
	{
		// calc. gradient
		grad = p.pre[0].scalar(f);

		// switch
		if (grad > 0.0)
			p.w[0] = min(p.w[0] * p.rhoplus, p.wMax);
		else 
		if (grad < 0.0)
			p.w[0] = max(p.w[0] * p.rhominus, p.wMin);

		// backup f
		p.pre[0] = f;
		
		// return
		return p.w[0];
	}
};

// rprop
class StepRPROP
{
	// parameters
	double grad;

	// direction
	int f1;
	
	// interface
	IStep < double > p;
	
	public :		
		void construct(double weight, int size = 1)
		{
			p.construct(weight, size);
		}
	
		void initialize()
		{
			p.initialize();
		}
	
		inline double adjust(double& f, int i)
		{
			// setup sign
			f1 = Lib::sign(f);
		
			// calc. gradient
			grad = f1 * p.pre[i];

			// rprop
			if (grad > 0)
			{
				// increase learning rate
				p.w[i]		= min(p.w[i] * p.rhoplus, p.wMax);
				p.pre[i]	= f1;
			} else if (grad < 0)
			{
				// actual shift in direction
				p.w[i]		= max(p.w[i] * p.rhominus, p.wMin);
				p.pre[i]	= 0;
			} else
			{
				// no direction or lastly shift in direction
				p.pre[i]	= f1;
			}

			// result
			return f1 * p.w[i];
		}
};
