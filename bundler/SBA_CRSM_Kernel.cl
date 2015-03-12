__kernel void SBA_CRSM(__local int *rcsubs_array, __local int *rcidxs_array,
			__global int *val, uint nvis, uint maxCPvis, uint nr,
			uint m,
			__global double* V, __global double* Yj, __global double* W,
			__global double *S, uint mmconxUsz,
			__global double *U, uint Sdim,
			__global double *E, __global double *eab)
{
	uint idj = get_global_id(0);
	uint idk = get_global_id(1);
	uint idi = get_global_id(2);

	uint localidj = get_local_id(0);
	uint localidk = get_local_id(1);
	uint localidi = get_local_id(2);

	__global int *colidx = val+nvis;
	__global int *rowptr = colidx+nvis;
	__local int *temp_ptr1, *temp_ptr2;

	__global double *ptr1, *ptr2, *ptr3, *ptr4, *eb;
	__local double *pYWt;
	__local double YWt_test[81];
	double sum = 0.0;

	uint i, l, k, ii, jj;
	uint nnz;
	int l_val;

/* Computing nnz_array values */

	temp_ptr1 = rcsubs_array;
	temp_ptr2 = rcidxs_array;
	
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
	
	#pragma unroll 81
	for(ii = 0; ii < 81; ++ii)
		YWt_test[ii] = 0.0;

/* Computing l_array values */
	if(idi < nnz)
	{
		i = rcsubs_array[idj];
		for(k = rowptr[i]; k < rowptr[i+1]; ++k)
		{
			if(idj == colidx[k])
				l_val = k; //l_array[idj*maxCPvis*m+idk*maxCPvis+idi] = k;
		}
		l_val = -1;//l_array[idj*maxCPvis*m+idk*maxCPvis+idi] = -1;

	
/* Compute Yj */

		ptr3 = V+rcsubs_array[idi]*9; //pnp*pnp = 3*3 = 9

		ptr1 = Yj+idj*maxCPvis*27+idi*27; // cnp*pnp = 9*3 = 27
		ptr2 = W+val[rcidxs_array[idi]]*27; // cnp*pnp = 9*3 = 27
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
			ptr2 = W+val[l_val]*27; // cnp*pnp = 9*3 = 27
			if(localidi == 0 && localidj == 0 && localidk == 0)
			{
				for(i = 0; i < nnz; i++)
				{
					ptr1 = Yj+idj*maxCPvis*27+i*27; // cnp*pnp = 9*3 = 27
					pYWt = YWt_test;//+idj*m*81*maxCPvis + idk*81*maxCPvis + idi*81;
					for(ii = 0; ii < 9; ++ii) // cnp = 9
					{
						 // cnp*cnp = 9*9 = 81, cnp = 9
						ptr3 = ptr1+ii*3; // pnp = 3
		
						ptr4 = ptr2;
						#pragma unroll 9
						for(jj = 0; jj < 9; ++jj) // cnp = 9
						{
							#pragma unroll 3
							for(l = 0, sum = 0.0; l < 3; ++l) // pnp = 3
								sum+=ptr3[l]*ptr4[l];
							if(i == 0)
								pYWt[ii*9+jj] = sum;//YWt[idj*m*YWtsz*maxCPvis+idk*YWtsz*maxCPvis+idi*YWtsz+ii*cnp+jj] = sum;
							else 
								pYWt[ii*9+jj] += sum;
							ptr4 += 3; // pnp = 3
						}
					}
				}
			}
		}
	/*	else
		{
			if(localidi == 0 && localidj == 0 && localidk == 0)
			{
				
			}
		}*/
	}
	barrier(CLK_LOCAL_MEM_FENCE);

/* Compute S */

	ptr2 = S + (idk)*mmconxUsz + (idj)*9; // cnp = 9
	ptr4 = YWt_test;//+idj*m*YWtsz*maxCPvis + idk*YWtsz*maxCPvis + idi*YWtsz;
	
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
		ptr1 = U + idj*81; // Usz = cnp*cnp = 9*9 = 81
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
	eb=eab+m*9; // cnp = 9
	
	if((idk == 0) && (idi < nnz))
	{
		ptr2 = Yj+idj*maxCPvis*27+idi*27; // Ysz = cnp*pnp = 9*3 = 27
		ptr3 = eb+rcsubs_array[idi]*3; //ebsz = 3
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
