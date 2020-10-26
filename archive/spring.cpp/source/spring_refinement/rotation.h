// rotation
struct Rotation
{
	public:
		double q[4];
		
		// axis and angle to quaternion
		Rotation (Vec3& axis, double angle)
		{
			// prepare
			axis.normalize();
			axis *= sin(angle / 2.0);

			// set
			q[0] = axis[0];
			q[1] = axis[1];
			q[2] = axis[2];
			q[3] = cos(angle / 2.0);
		}
		
		// axis and angle to quaternion
		Rotation (Vec4& v)
		{
			q[0] = v.w;
			q[1] = v.x;
			q[2] = v.y;
			q[3] = v.z;
		}
};

// matrix
class RotationMatrix
{

	// matrix
	double m[4][4];
public:
	// setup
	void setup(Vec4& q)
	{		
		// reset
		identity();
		
		// rotate
		RotationMatrix matrix;
		matrix.setRotate(q);
		set(matrix * (*this));
	}
	
	void identity()
	{
		m[0][0] = 1.0;  m[0][1] = 0.0;  m[0][2] = 0.0;  m[0][3] = 0.0;
		m[1][0] = 0.0;  m[1][1] = 1.0;  m[1][2] = 0.0;  m[1][3] = 0.0;
		m[2][0] = 0.0;  m[2][1] = 0.0;  m[2][2] = 1.0;  m[2][3] = 0.0;
		m[3][0] = 0.0;  m[3][1] = 0.0;  m[3][2] = 0.0;  m[3][3] = 1.0;
	}
	
	void setRotate (Vec4& q)
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
	
	friend RotationMatrix operator*(const RotationMatrix& m1, const RotationMatrix& m2)
	{
		RotationMatrix out;
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
	
	void set(const RotationMatrix& matrix )
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

