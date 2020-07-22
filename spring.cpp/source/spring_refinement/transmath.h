// matrix
class Mat4
{
	// matrix
	double m[4][4];
public:	
	void identity()
	{
		m[0][0] = 1.0;  m[0][1] = 0.0;  m[0][2] = 0.0;  m[0][3] = 0.0;
		m[1][0] = 0.0;  m[1][1] = 1.0;  m[1][2] = 0.0;  m[1][3] = 0.0;
		m[2][0] = 0.0;  m[2][1] = 0.0;  m[2][2] = 1.0;  m[2][3] = 0.0;
		m[3][0] = 0.0;  m[3][1] = 0.0;  m[3][2] = 0.0;  m[3][3] = 1.0;
	}
	
	void setRotate (double* q)
	{
		m[0][0] = (double)(1 - 2.0 * (q[1] * q[1] + q[2] * q[2]));
		m[0][1] = (double)(    2.0 * (q[0] * q[1] + q[2] * q[3]));
		m[0][2] = (double)(    2.0 * (q[2] * q[0] - q[1] * q[3]));
		m[0][3] = 0.0;
		m[1][0] = (double)(    2.0 * (q[0] * q[1] - q[2] * q[3]));
		m[1][1] = (double)(1 - 2.0 * (q[2] * q[2] + q[0] * q[0]));
		m[1][2] = (double)(    2.0 * (q[1] * q[2] + q[0] * q[3]));
		m[1][3] = 0.0;
		m[2][0] = (double)(    2.0 * (q[2] * q[0] + q[1] * q[3]));
		m[2][1] = (double)(    2.0 * (q[1] * q[2] - q[0] * q[3]));
		m[2][2] = (double)(1 - 2.0 * (q[1] * q[1] + q[0] * q[0]));
		m[2][3] = 0.0;
		m[3][0] = 0.0;
		m[3][1] = 0.0;
		m[3][2] = 0.0;
		m[3][3] = 1.0;
	}
	
	void setTranslate (Vec3& t)
	{
		m[0][0] = 1.0;
		m[0][1] = 0.0;
		m[0][2] = 0.0;
		m[0][3] = 0.0;
		m[1][0] = 0.0;
		m[1][1] = 1.0;
		m[1][2] = 0.0;
		m[1][3] = 0.0;
		m[2][0] = 0.0;
		m[2][1] = 0.0;
		m[2][2] = 1.0;
		m[2][3] = 0.0;
		m[3][0] = t[0];
		m[3][1] = t[1];
		m[3][2] = t[2];
		m[3][3] = 1.0;
	}
	
	void setTransform(Vec3 &t, double* q, Vec3 &c)
	{
		Mat4 matrix;
		
		// reset
		identity();
		
		// encode translation
		if (t != 0)
		{
    		matrix.setTranslate(t);
            set(matrix * (*this));
        }
		
		// shift to center
		if (c != 0)
		{
            matrix.setTranslate(c);
	       	set(matrix * (*this));
        }
		
		// rotate
		matrix.setRotate(q);
		set(matrix * (*this));
		
		// shift back
		if (c != 0)
		{
            c.neg();
    		matrix.setTranslate(c);
            set(matrix * (*this));
        }
	}
	
	friend Mat4 operator*(const Mat4& m1, const Mat4& m2)
	{
		Mat4 out;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				out.m[i][j] = 0.0;
				for (int k = 0; k < 4; k++)
					out.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
		return out;
	}
	
	void reset()
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				m[i][j] = 0.0;
	}
	
	void set(const Mat4& matrix )
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				m[i][j] = matrix.m[i][j];
	}
	
	void multVecMatrix(Vec3 &src, Vec3 &dst)
	{
		double f = 1.0 / (src[0]*m[0][3]+src[1]*m[1][3]+src[2]*m[2][3]+m[3][3]);		
		dst[0]  =    f*(src[0]*m[0][0]+src[1]*m[1][0]+src[2]*m[2][0]+m[3][0]);
		dst[1]  =    f*(src[0]*m[0][1]+src[1]*m[1][1]+src[2]*m[2][1]+m[3][1]);
		dst[2]  =    f*(src[0]*m[0][2]+src[1]*m[1][2]+src[2]*m[2][2]+m[3][2]);
	}
};

