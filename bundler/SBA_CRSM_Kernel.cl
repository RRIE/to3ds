__kernel void SBA_CRSM(__global int *rcsubs_array, __global int *rcidxs_array,
			__global int *val, __global int *nnz_array, 
			__global int *l_array, uint nvis, uint maxCPvis, uint nr,
			uint m,
			__global double* V, __global double* Yj, __global double* W,
			uint cnp, uint pnp, uint Vsz, uint Wsz, uint Ysz)
{
	uint idj = get_global_id(0);
	uint idk = get_global_id(1);
	uint idi = get_global_id(2);

	__global int *colidx = val+nvis;
	__global int *rowptr = colidx+nvis;
	__global int *temp_ptr1, *temp_ptr2;

	__global double *ptr1, *ptr2, *ptr3, *ptr4;
	double sum = 0.0;

	uint i, l, k, ii, jj;
	uint nnz;

/* Computing nnz_array values */

	temp_ptr1 = rcsubs_array+idj*maxCPvis;
	temp_ptr2 = rcidxs_array+idj*maxCPvis;
		
	for(i = l = 0; i < nr; ++i)
	{
		for(k = rowptr[i]; k < rowptr[i+1]; ++k)
		{
			if(colidx[k] == idj)
			{
				temp_ptr2[l] = k;
				temp_ptr1[l] = i;
				l++;
			}
		}
	}
	nnz = l;	
	



/* Computing l_array values */
	if(idi < nnz)
	{
		i = rcsubs_array[idj*maxCPvis+idi];
		for(k = rowptr[i]; k < rowptr[i+1]; ++k)
		{
			if(idj == colidx[k])
				l_array[idj*maxCPvis*m+idk*maxCPvis+idi] = k;
		}
		l_array[idj*maxCPvis*m+idk*maxCPvis+idi] = -1;
	}

	
/* Compute Yj */

	if(idk == idj)
	{
		if(idi < nnz)
		{
			ptr3 = V+rcsubs_array[idj*maxCPvis+idi]*Vsz;

			ptr1 = Yj+idj*maxCPvis*Ysz+idi*Ysz;
			ptr2 = W+val[rcidxs_array[idj*maxCPvis+idi]]*Wsz;
			for(ii = 0; ii < cnp; ++ii)
			{
				ptr4 = ptr2+ii*pnp;
				for(jj = 0; jj < pnp; ++jj)
				{
					for(l = 0, sum = 0.0; l <= jj; ++l)
						sum+=ptr4[l]*ptr3[jj*pnp+l];
					for( ;l <pnp; ++l)
						sum+=ptr4[l]*ptr3[l*pnp+jj];
					ptr1[ii*pnp+jj] = sum;
				}
			}
		}
	}

	if(idk == 0 && idi == 0)
		nnz_array[idj] = nnz;

}
