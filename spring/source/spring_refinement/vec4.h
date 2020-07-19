// 4d vector
struct Vec4
{
	// values
	double w, x, y, z;
	
	// construct
	Vec4()
	{
		initialize(0, 0, 0, 0);
	}
			
	Vec4(double w, double x, double y, double z)
	{
		initialize(w, x, y, z);
	}
	
	// operator
	void operator*=(const double &d)
	{
		w *= d;
		x *= d;
		y *= d;
		z *= d;
	}
	
	// operator
	void operator/=(const double &d)
	{
		w /= d;
		x /= d;
		y /= d;
		z /= d;
	}
	
	// get length
    double length() const
    {
        double factor = w*w + x*x + y*y + z*z;
        if (factor > 0)
		  return sqrt(factor);
		else
		  return 0;
    }
	
    // inefficient for compatibility only
    inline double& operator[](int i)
	{
       	assert(i < 4);
		assert(i >= 0);
		return (&w)[i];
	}
private:
	void initialize(double w, double x, double y, double z)
	{
		this->w = w;
		this->x = x;
		this->y = y;
		this->z = z;						
	}		
};
