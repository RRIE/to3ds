__kernel void SBAOpenCL(__global double* src_V, __global double* src_Yj,
			__global double* src_W, __global int* src_rcidxs,
			__global int* src_rcsubs, __global int* src_idxij_val,
			uint nnz, uint cnp, uint pnp, uint Vsz, 
			uint Wsz, uint Ysz)
{
	uint idi = get_global_id(0);
	uint idii = get_global_id(1);
	uint idjj = get_global_id(2);

	__global double *ptr3, *ptr1, *ptr2, *ptr4;
	double sum = 0.0;
	uint k;
	if(idi < nnz && idii < cnp && idjj < pnp)
	{
		ptr3 = src_V+src_rcsubs[idi]*Vsz;
		ptr1 = src_Yj+idi*Ysz;
		ptr2 = src_W + src_idxij_val[src_rcidxs[idi]]*Wsz;
		ptr4 = ptr2+idii*pnp;
				
		for(k = 0, sum = 0.0; k <= idjj; ++k)
			sum+=ptr4[k]*ptr3[idjj*pnp+k];
		for( ;k<pnp; ++k)
			sum+=ptr4[k]*ptr3[k*pnp+idjj];
		ptr1[idii*pnp+idjj] = sum;
	}
}
