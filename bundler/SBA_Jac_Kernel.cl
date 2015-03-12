#define JAC_DEBUG 0

//Structure copied out of sfm.h, no need to worry
typedef struct {
    double R[9];     /* Rotation */
    double t[3];     /* Translation */
    double f;        /* Focal length */
    double k[2];     /* Undistortion parameters */
    double k_inv[6]; /* Inverse undistortion parameters */
    char constrained[9];
    double constraints[9];  /* Constraints (if used) */
    double weights[9];      /* Weights on the constraints */
    double K_known[9];  /* Intrinsics (if known) */
    double k_known[5];  /* Distortion params (if known) */

    char fisheye;            /* Is this a fisheye image? */
    char known_intrinsics;   /* Are the intrinsics known? */
    double f_cx, f_cy;       /* Fisheye center */
    double f_rad, f_angle;   /* Other fisheye parameters */
    double f_focal;          /* Fisheye focal length */

    double f_scale, k_scale; /* Scale on focal length, distortion params */
} camera_params_t;

void matrix_product33(double *A, double *B, double *R)
{
    R[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    R[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    R[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];    
    
    R[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    R[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    R[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];    

    R[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    R[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    R[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];            
} 


void matrix_product33_g(double *A, __global double *B, __global double *R)
{
    R[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    R[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    R[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];    
    
    R[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    R[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    R[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];    

    R[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    R[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    R[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];            
} 

void matrix_product331_gg(__global double *A, __global double *b, double *r)
{
    r[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
    r[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
    r[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}

void matrix_product331_g(__global double *A, double *b, double *r)
{
    r[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
    r[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
    r[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}


/* Scale a matrix by a scalar */
void matrix_scale(int m, int n, double *A, double s, double *R) {
    int i;
    int entries = m * n;
    
    for (i = 0; i < entries; i++) {
	R[i] = A[i] * s;
    }
}


/* Compute the matrix sum R = A + B */
void matrix_sum(int Am, int An, int Bm, int Bn, 
                double *A,  double *B, double *R) {
    int r = Am;
    int c = An;
    int n = r * c, i;
    
    if (Am != Bm || An != Bn) {
	//printf("[matrix_sum] Error: mismatched dimensions\n");
	return;
    }

    for (i = 0; i < n; i++) {
	R[i] = A[i] + B[i];
    }
}


/* Compute an updated rotation matrix given the initial rotation (R)
 * and the correction (w) */
void rot_update(__global double *R, __global double *w, __global double *Rnew) 
{
	
    if (JAC_DEBUG) printf("[rot_update] %s\n", "Calling Rotate");

    double theta, sinth, costh, n[3];
    double nx[9], nxsq[9];
    double term2[9], term3[9];
    double tmp[9], dR[9];

    double ident[9] = 
	{ 1.0, 0.0, 0.0,
	  0.0, 1.0, 0.0,
	  0.0, 0.0, 1.0 };

    theta = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);

    if (theta == 0.0) {
    	int i;
    	for (i=0;i<9;i++) {
	    Rnew[i]=R[i];
	}
	return;
    }

    n[0] = w[0] / theta;
    n[1] = w[1] / theta;
    n[2] = w[2] / theta;

    nx[0] = 0.0;   nx[1] = -n[2];  nx[2] = n[1];
    nx[3] = n[2];  nx[4] = 0.0;    nx[5] = -n[0];
    nx[6] = -n[1]; nx[7] = n[0];   nx[8] = 0.0;

    matrix_product33(nx, nx, nxsq);

    sinth = sin(theta);
    costh = cos(theta);

   if (JAC_DEBUG) printf("[rot_update] sinth=%f, costh=%f\n", sinth, costh);

    matrix_scale(3, 3, nx, sinth, term2);
    matrix_scale(3, 3, nxsq, 1.0 - costh, term3);

    matrix_sum(3, 3, 3, 3, ident, term2, tmp);
    matrix_sum(3, 3, 3, 3, tmp, term3, dR);

    matrix_product33_g(dR, R, Rnew);
}


int sfm_project_rd(__global camera_params_t *init, double *K, __global double *k,	//FIXME: init will cause issues, maybe not
                   __global double *R, __global double *dt, __global double *b, double *p)
{
    double tnew[3];
    double b_cam[3];


    tnew[0] = dt[0];
    tnew[1] = dt[1];
    tnew[2] = dt[2];

    int in_front = 1;

    /* Project! */
    double b2[3];
    b2[0] = b[0] - tnew[0];
    b2[1] = b[1] - tnew[1];
    b2[2] = b[2] - tnew[2];
    matrix_product331_g(R, b2, b_cam);
    
    if (b_cam[2] >= 0.0)
        in_front = 0; // cheirality violation...

    if (!init->known_intrinsics) {
        p[0] = -b_cam[0] * K[0] / b_cam[2];
        p[1] = -b_cam[1] * K[0] / b_cam[2];
    } else {
        /* Apply intrinsics */
        double x_n = -b_cam[0] / b_cam[2];
        double y_n = -b_cam[1] / b_cam[2];

        __global double *k = init->k_known;
	double rsq = x_n * x_n + y_n * y_n;
	double factor = 1.0 + k[0] * rsq + 
            k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

        double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
        double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

	double x_d = x_n * factor + dx_x;
	double y_d = y_n * factor + dx_y;

        __global double *K = init->K_known;
        p[0] = K[0] * x_d + K[1] * y_d + K[2];
        p[1] = K[4] * y_d + K[5];
    }

    /* Apply radial distortion */
        double k1 = k[0] / init->k_scale;
        double k2 = k[1] / init->k_scale;

	double rsq = (p[0] * p[0] + p[1] * p[1]) / (K[0] * K[0]);
	double factor = 1.0 + k1 * rsq + k2 * rsq * rsq;

	p[0] *= factor;
	p[1] *= factor;

    return in_front;
}

void sfm_project_point3(int j, int i, __global double *aj, __global double *bi, 
			       double *xij, __global camera_params_t *init_params, __global double* last_ws, __global double* last_Rs)	{

    double K[9] = 
	{ 1.0, 0.0, 0.0,
	  0.0, 1.0, 0.0,
	  0.0, 0.0, 1.0 };

    __global double *w, *dt, *k;

    /* Compute intrinsics */
    K[0] = K[4] = aj[6] / init_params[j].f_scale;
    if (JAC_DEBUG) printf("K[0] is %f\n", K[0]);
    
    /* Compute translation, rotation update */
    dt = aj + 0;
    w = aj + 3;
    k = aj + 7;

    if (JAC_DEBUG) printf("W[0]=%f, W[1]=%f, W[2]=%f, last_ws[%d]=%f, last_ws[%d]=%f, last_ws[%d]=%f\n", w[0], w[1], w[2], 3*j, last_ws[3 * j + 0], 3*j+1, last_ws[3 * j + 1], 3*j+2, last_ws[3 * j + 2]);

    if (w[0] != last_ws[3 * j + 0] ||
	w[1] != last_ws[3 * j + 1] ||
	w[2] != last_ws[3 * j + 2]) {

	rot_update(init_params[j].R, w, last_Rs + 9 * j);
	last_ws[3 * j + 0] = w[0];
	last_ws[3 * j + 1] = w[1];
	last_ws[3 * j + 2] = w[2];
    }

    if (JAC_DEBUG) printf("Global_last_Rs {%f %f %f  %f %f %f %f %f %f}\n", *(last_Rs + 9 * j), *(last_Rs + 9 * j + 1), *(last_Rs + 9 * j +2), *(last_Rs + 9 * j + 3), *(last_Rs + 9 * j + 4), *(last_Rs + 9 * j + 5), *(last_Rs + 9 * j + 6), *(last_Rs + 9 * j + 7), *(last_Rs + 9 * j + 8));
    
    sfm_project_rd(&init_params[j], K, k, last_Rs + 9 * j, 
                   dt, bi, xij);
}
 

__kernel void SBA_Jac(
    __global double *p,          	/* I: current parameter estimate, (m*cnp+n*pnp)x1 */

    //idxij parameters
    uint nr, uint nc, uint nnz,	 	/*I*/
    __global int *val,			/*I*/
    __global int *colidx,		/*I*/
    __global int *rowptr,		/*I*/

    //Program shared global variables
    __global double *last_Rs,	/*IO*/	
    __global double *last_ws,	/*IO*/	

    __global double *jac,               /* O: array for storing the approximated jacobian */

    //Dat parameters	
    __global void *init_params		/*I*/	
) {
  int i, j, ii, jj;
  __global double *pa, *pb, *paj, *pbi;
  __global double *pAB;
  int n, m, Asz, Bsz, ABsz;

  double tmp;
  double d, d1;

  //Get the work item or thread id
  uint thread_ID = get_global_id(0);

#define cnp 9
#define pnp 3
#define mnp 2

  /* allocate memory for hxij, hxxij */
  double hxij[2*mnp];
  double *hxxij=hxij+mnp;

/*
  if (thread_ID==0) {
      int iiii;
      for (iiii=0; iiii<nr*pnp + nc*cnp; iiii++)  
  	printf("[SBA_Jac] thread %d, p[%d] = %f\n", thread_ID, iiii, p[iiii]);	
  }
*/

  /* retrieve problem-specific information passed in *dat */
  if (JAC_DEBUG) printf("[SBA_Jac] thread %d, m = %d\n", thread_ID, m);	

  n=nr; m=nc;
  pa=p; pb=p+m*cnp;
  Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;

    /* compute A_ij */
    j=thread_ID; 
      /*for (j=0; j<m; ++j)*/ {
      paj=pa+j*cnp; // j-th camera parameters

      for(jj=0; jj<cnp; ++jj){
        d=(double)(0.0001)*paj[jj]; // force evaluation
        d=(((d)>=0)? (d) : -(d));
        if(d< 0.000001) d = 0.000001;
        d1=1.0/d; /* invert so that divisions can be carried out faster as multiplications */

	int low, high, diff, mid;
	for (i=0; i<nr; ++i){
	    //Do binary search tree to get the columns we want
	    low = rowptr[i];
	    high = rowptr[i+1] - 1;
   	    while(low<=high){
		mid=(low+high)>>1; //(low+high)/2;
		diff=j-colidx[mid];
      		if(diff<0)
        	    high=mid-1;
      		else if(diff>0)
        	    low=mid+1;
      		else{ /* found */
		    pbi=pb + i*pnp; // i-th point parameters
	  	    if (JAC_DEBUG) printf("[SBA_Jac] thread %d, will be call sfm_project_point3 for calculating A_ij, with %d, %d, %f, %f\n", thread_ID, j, i, *paj, *pbi);	
          	    sfm_project_point3(j, i, paj, pbi, hxij, init_params, last_ws, last_Rs); // evaluate supplied function on current solution

          	    tmp=paj[jj];
          	    paj[jj]+=d;
          	    sfm_project_point3(j, i, paj, pbi, hxxij, init_params, last_ws, last_Rs);
          	    paj[jj]=tmp; /* restore */

          	    pAB=jac + val[mid]*ABsz; // set pAB to point to A_ij
          	    for(ii=0; ii<mnp; ++ii) {
            	         pAB[ii*cnp+jj]=(hxxij[ii]-hxij[ii])*d1;
			 if (JAC_DEBUG) printf("[SBA_Jac] thread %d wrote Jac[%d] = %f, d1 = %f\n", thread_ID, val[mid]*ABsz + ii*cnp+jj, (hxxij[ii]-hxij[ii])*d1, d1);
		    }	
        	    break;
      		} 
	    }
	}
      }
    }

    /* compute B_ij */
    for(i=0; i<n; ++i){
      pbi=pb+i*pnp; // i-th point parameters

      for(jj=0; jj<pnp; ++jj){
        d=(double)(0.0001)*pbi[jj]; // force evaluation
        d=(((d)>=0)? (d) : -(d));
        if(d< 0.000001) d = 0.000001;
        d1=1.0/d; /* invert so that divisions can be carried out faster as multiplications */

	for(j=rowptr[i]; j < rowptr[i+1]; ++j) {
	    if (colidx[j]==thread_ID) {
		if (JAC_DEBUG) printf("[SBA_Jac] thread %d, will be call sfm_project_point3 for calculating B_ij, with %d, %d, %f, %f\n", thread_ID, colidx[j], i, *paj, *pbi);	
          	paj=pa + colidx[j]*cnp; // j-th camera parameters
          	sfm_project_point3(colidx[j], i, paj, pbi, hxij, init_params, last_ws, last_Rs); // evaluate supplied function on current solution

          	tmp=pbi[jj];
          	pbi[jj]+=d;
          	sfm_project_point3(colidx[j], i, paj, pbi, hxxij, init_params, last_ws, last_Rs);
          	pbi[jj]=tmp; /* restore */

          	pAB=jac + val[j]*ABsz + Asz; // set pAB to point to B_ij
          	for(ii=0; ii<mnp; ++ii) {
            	    pAB[ii*pnp+jj]=(hxxij[ii]-hxij[ii])*d1;
		    if (JAC_DEBUG) printf("[SBA_Jac] thread %d wrote Jac[%d] = %f\n", thread_ID, val[j]*ABsz + Asz + ii*pnp+jj, (hxxij[ii]-hxij[ii])*d1);
		}
	    }
	}
      }
    }

    if (thread_ID==0) {
    	int iiii;
    	for (iiii=0; iiii<nr*pnp; iiii++)  
  	   if (JAC_DEBUG) printf("[SBA_Jac] thread %d, Jac[%d] = %f\n", thread_ID, iiii, jac[iiii]);	
    }
}
