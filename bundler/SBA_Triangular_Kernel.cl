__kernel void SBA_Triangular(__global double *src_Yj, __global double *src_W, 
			__global double *src_YWt, __global double *src_S,
			__global double *src_U, __global int* crsm_array,
			uint mcon, uint mmconxUsz, uint j,uint nnz, uint cnp, 
			uint pnp, uint Wsz, uint Ysz, uint Usz, uint Sdim,
			__global int *idxij_val, uint YWtsz, uint m)
{
	uint idk = get_global_id(0);
	uint idi = get_global_id(1);
	uint idii = get_global_id(2);
	uint jj = 0;	
	uint h = 0;
	double sum = 0.0;

	
	if((idk < m) && (idi < nnz) && (idii < cnp))
	{
		__global double *ptr1, *ptr2, *ptr3, *ptr4, *pYWt;

		int l = crsm_array[idk*nnz+idi];

		if(l != -1)
		{
			ptr2 = src_W + idxij_val[l]*Wsz;
			ptr1 = src_Yj + idi*Ysz;
			
			ptr3 = ptr1 + idii*pnp;
			pYWt = src_YWt + idk*YWtsz*nnz + idi*YWtsz + idii*cnp;

			ptr4 = ptr2;

			for(jj = 0; jj < cnp; ++jj)
			{
				for(h = 0, sum = 0.0; h < pnp; ++h)
					sum+=ptr3[h]*ptr4[h];
				pYWt[jj] = sum; 
				ptr4 += pnp;
				//printf("SBA_Triangular sum is %f\n", sum);
			}	
	
		//	printf("SBA_Triangular_Kernel read back valid l value %d\n", l);
		}
		//else
		//{
		//	printf("Wrong result %d\n", l);
		//}
	/*	barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

		ptr2 = src_S + (idk-mcon)*mmconxUsz + (j-mcon)*cnp;

		if(idk != j)
		{
			ptr2 = ptr2 + idii*Sdim;
			for(jj = 0; jj < cnp; ++jj)
				ptr2[jj] = -src_YWt[idk*YWtsz+jj*cnp+idii];
		}
		else
		{
			ptr1 = src_U + j*Usz;
			ptr2 = ptr2 + idii*Sdim;
			for(jj = 0; jj < cnp; ++jj)
				ptr2[jj] = ptr1[jj*cnp+idii]-src_YWt[idk*YWtsz+jj*cnp+idii];
		}*/
	}
}
