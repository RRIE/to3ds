__kernel void SBA_CRSM(__global int *rcsubs_array, __global int *rcidxs_array,
			__global int *val, uint nvis, uint maxCPvis, uint nr,
			uint m,
			__global double* V, __global double* Yj, __global double* W,
			 uint Vsz, uint Wsz, uint Ysz,
			__global double *YWt, uint YWtsz,
			__global double *S, uint mcon, uint mmconxUsz,
			__global double *U, uint Usz, uint Sdim,
			__global double *E, __global double *eab)
{
	uint idj = get_global_id(0);
	uint idk = get_global_id(1);
	uint idi = get_global_id(2);

	__global int *colidx = val+nvis;
	__global int *rowptr = colidx+nvis;
	__global int *temp_ptr1, *temp_ptr2;

	__global double *ptr1, *ptr2, *ptr3, *ptr4, *pYWt, *ea, *eb;
	double sum = 0.0;

	uint i, l, k, ii, jj;
	uint nnz;
	int l_val;

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
				l_val = k; //l_array[idj*maxCPvis*m+idk*maxCPvis+idi] = k;
		}
		l_val = -1;//l_array[idj*maxCPvis*m+idk*maxCPvis+idi] = -1;

	
/* Compute Yj */

		ptr3 = V+rcsubs_array[idj*maxCPvis+idi]*Vsz;

		ptr1 = Yj+idj*maxCPvis*Ysz+idi*Ysz;
		ptr2 = W+val[rcidxs_array[idj*maxCPvis+idi]]*Wsz;
		#pragma unroll 9
		for(ii = 0; ii < 9; ++ii) // cnp = 9
		{
			ptr4 = ptr2+ii*3; // pnp = 3 
			#pragma unroll 3
			for(jj = 0; jj < 3; ++jj) // pnp = 3
			{
				for(l = 0, sum = 0.0; l <= jj; ++l)
					sum+=ptr4[l]*ptr3[jj*3+l]; //pnp = 3
				for( ;l < 3; ++l) // pnp = 3
					sum+=ptr4[l]*ptr3[l*3+jj]; //pnp = 3
				ptr1[ii*3+jj] = sum; // pnp = 3
			}
		}
		if(l_val != -1)
		{
			ptr2 = W+val[l_val]*Wsz;
			for(i = 0; i < nnz; i++)
			{
				ptr1 = Yj+idj*maxCPvis*Ysz+i*Ysz;
			
				for(ii = 0; ii < 9; ++ii) // cnp = 9
				{
					pYWt = YWt+idj*m*YWtsz*maxCPvis + idk*YWtsz*maxCPvis + idi*YWtsz + ii*9; // cnp = 9
					ptr3 = ptr1+ii*3; // pnp = 3
	
					ptr4 = ptr2;
					#pragma unroll 9
					for(jj = 0; jj < 9; ++jj) // cnp = 9
					{
						#pragma unroll 3
						for(l = 0, sum = 0.0; l < 3; ++l) // pnp = 3
							sum+=ptr3[l]*ptr4[l];
						if(i == 0)
							pYWt[jj] = sum;//YWt[idj*m*YWtsz*maxCPvis+idk*YWtsz*maxCPvis+idi*YWtsz+ii*cnp+jj] = sum;
						else 
							pYWt[jj] += sum;
						ptr4 += 3; // pnp = 3
					}
				}
			}
		}
	}

/* Compute S */

	ptr2 = S + (idk-mcon)*mmconxUsz + (idj - mcon)*9; // cnp = 9
	ptr4 = YWt;//+idj*m*YWtsz*maxCPvis + idk*YWtsz*maxCPvis + idi*YWtsz;
	
	if(idj!=idk)
	{
		#pragma unroll 9
		for(ii = 0; ii < 9; ++ii, ptr2+=Sdim) // cnp = 9
		{
			#pragma unroll 9
			for(jj = 0; jj < 9; ++jj) // cnp = 9
			{
				ptr2[jj] = -ptr4[jj*9+ii]; // cnp = 9
			}
		}
	}
	else
	{
		ptr1 = U + idj*Usz;
		#pragma unroll 9
		for(ii = 0; ii < 9; ++ii, ptr2+=Sdim) // cnp = 9
		{
			#pragma unroll 9
			for(jj = 0; jj < 9; ++jj) // cnp = 9
			{
				ptr2[jj] = ptr1[jj*9+ii]-ptr4[jj*9+ii]; // cnp = 9
			}
		}
	}

	/* Compute e_j=ea_j */
	ptr1 = E + idj*9*maxCPvis+idi*9; // easz = 9
	ea=eab; 
	eb=eab+m*9; // cnp = 9
	
	if((idk == 0) && (idi < nnz))
	{
		ptr2 = Yj+idj*maxCPvis*Ysz+idi*Ysz;
		ptr3 = eb+rcsubs_array[idj*maxCPvis+idi]*3; //ebsz = 3
		#pragma unroll 9
		for(ii = 0; ii < 9; ++ii) // cnp = 9
		{
			ptr4 = ptr2+ii*3;//pnp = 3
			#pragma unroll 3
			for(jj=0, sum = 0.0; jj < 3; ++jj) //  pnp = 3
				sum+=ptr4[jj]*ptr3[jj];
			ptr1[ii] = sum;
		}
	}

}
