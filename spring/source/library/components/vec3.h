// default vector class
struct Vec3
{
	//
	// data
	//
	double x, y, z;

	//
	// construction
	//
	Vec3(double x, double y, double z)
	{
		initialize (x, y, z);
	}
	
	Vec3(double x = 0.0)
	{
		initialize (x, x, x);
	}
	
	int size()
	{
		return 3;
	}
	
	//
	// operators
	//
	Vec3 operator+(const Vec3 &v) const
	{
		Vec3 v0;
		v0.x = x + v.x;
		v0.y = y + v.y;
		v0.z = z + v.z;
		return v0;
	}
	
	Vec3 operator-(const Vec3& v) const
	{
		Vec3 v0;
		v0.x = x - v.x;
		v0.y = y - v.y;
		v0.z = z - v.z;
		return v0;
	}
	
	void operator+=(const Vec3 &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}
	
	void operator-=(const Vec3 &v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	void operator/=(const double &d)
	{
		x /= d;
		y /= d;
		z /= d;
	}
	
	void operator*=(const double &d)
	{
		x *= d;
		y *= d;
		z *= d;
	}

	void operator-=(const double &d)
	{
		x -= d;
		y -= d;
		z -= d;
	}

	void operator+=(const double &d)
	{
		x += d;
		y += d;
		z += d;
	}

	void operator*=(const Vec3 &v)
	{
		x = y * v.z - z * v.y;
        y = z * v.x - x * v.z;
        z = x * v.y - y * v.x;
	}
	

	double operator*(const Vec3& v) const
	{
		return dot(v);
	}
	
	Vec3 operator*(const double& d) const
	{
		return Vec3 (x*d, y*d, z*d);
	}
	
	bool operator==(const Vec3 &v) const
	{
		return (length() == v.length());
	}

	bool operator==(const double &d) const
	{
		return (x == d && y == d && z == d);
	}
    
	bool operator!=(const double &d) const
	{
		return !operator==(d);
	}
    
	// unary minus
	friend Vec3 operator-(const Vec3& v)
	{
		return Vec3(-v.x, -v.y, -v.z);
	}
	
    // inefficient for compatibility only
    inline double& operator[](int i)
	{
       	assert(i < 3);
		assert(i >= 0);
		return (&x)[i];
	}
	
    //
	// methods
	//

	// distance
	double dist2(const Vec3& v) const
	{
		Vec3 d = (*this) - v;		
		return d.x*d.x + d.y*d.y + d.z*d.z;
	}
    	
	// distance
	double dist(const Vec3& v) const
	{
		return (Vec3 ((*this) - v).length());
	}
	
	// dot product
	double dot(const Vec3& v) const
	{
		return (x*v.x + y*v.y + z*v.z);
	}

    // cross product
	Vec3 cross(const Vec3& v) const
	{
		return Vec3 (y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}
	
	// one cross product  self = cross(v1, v2)
	template < typename TA, typename TB >
	void cross(TA& v1, TB& v2)
	{
		x= v1.y*v2.z-v2.y*v1.z;
		y= v2.x*v1.z-v1.x*v2.z;
		z= v1.x*v2.y-v2.x*v1.y;
	}

	// multiplying a cross product by a scalar is very common
	// one cross product  v3 = k*cross(v1, v2)
	template < typename TA, typename TB >
	static Vec3 cross (double k, TA& v1, TB& v2)
	{
		return Vec3( k*(v1.y*v2.z-v2.y*v1.z),
		      k*(v2.x*v1.z-v1.x*v2.z),
		      k*(v1.x*v2.y-v2.x*v1.y) );
	}

    // cross
	void neg()
	{
		x = -x;
		y = -y;
		z = -z;
	}
		  
	int equals(Vec3& v)
	{
		return ( (x == v.x) && (y == v.y) && (z == v.z) );
	}

	// set all
	inline void set(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	// get length
	double length() const
	{
		double factor = x*x + y*y + z*z;
        if (factor > 0)
			return sqrt(factor);
		else
			return 0.0;
	}

	// is normalized
	bool isNormalized() const
	{
		return (length() == 1);
	}
    
	// normalize
	void normalize()
	{
		(*this) /= length();
	}
	
	// scalar
	double scalar(const Vec3 &v) const
	{
		return (x * v.x) + (y * v.y) + (z * v.z);
	}

	// get angle
	inline double angle(Vec3 v)
	{
		double lng1 = this->length();
		double lng2 = v.length();
		if ( lng1 == 0 || lng2 == 0 )
			return 0;
	
		double a = acos( this->dot(v) / ( lng1 * lng2 ) );
		return ( a / M_PI ) * 180.0;
	}

private:
	// initialize
	void initialize (double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
};

// output
ostream& operator<<(ostream& stream, const Vec3& v)
{
	stream << setw(15) << v.x << setw(15) << v.y << setw(15) << v.z;
	return stream;
}

