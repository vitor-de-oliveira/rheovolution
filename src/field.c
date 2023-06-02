#include "field.h"

int
func (double t, const double y[], double f[],
      void *params)
{
	(void)(t);
	
	double G = 1.0;
	double m1 = 1.0;
	double m2 = 1.0 - m1;
	
  	double total_mass = m1 + m2;
  	
  	double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  	double r_cube = r * r * r;
  
  	/**
  	 * y[0]  = x
  	 * y[1]  = y
  	 * y[2]  = z
  	 * y[3]  = x_dot
  	 * y[4]  = y_dot
  	 * y[5]  = z_dot
  	 * y[6]  = l_x
  	 * y[7]  = l_y
  	 * y[8]  = l_z
  	 * y[9]  = b_0_11
  	 * y[10] = b_0_12
  	 * y[11] = b_0_13
  	 * y[12] = b_0_21
  	 * y[13] = b_0_22
  	 * y[14] = b_0_23
  	 * y[15] = b_0_31
  	 * y[16] = b_0_32
  	 * y[17] = b_0_33
  	 * y[18] = B_11
  	 * y[19] = B_12
  	 * y[20] = B_13
  	 * y[21] = B_21
  	 * y[22] = B_22
  	 * y[23] = B_23
  	 * y[24] = B_31
  	 * y[25] = B_32
  	 * y[26] = B_33
	*/
  
	f[0] = y[3];
	f[1] = y[4];
	f[2] = y[5];
	f[3] = -1.0 * G * total_mass * y[0] / r_cube;
	f[4] = -1.0 * G * total_mass * y[1] / r_cube;
	f[5] = -1.0 * G * total_mass * y[2] / r_cube;
	
	return GSL_SUCCESS;
}
