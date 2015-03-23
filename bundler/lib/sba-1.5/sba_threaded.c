#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include "sba.h"
#include "sba_levmar.h"

#define NR_THREADS 8
#define BUFFER_SIZE 8

#define cnp 9
#define pnp 3
#define Vsz 9
#define Ysz 27
#define Wsz 27
#define Usz 81
#define Sblsz 81
#define easz 9
#define ebsz 3
#define YWtsz 81

#define FABS(x)           (((x)>=0)? (x) : -(x))

pthread_mutex_t *global_lock = NULL;

struct threadParams
{
	int mmconxUsz;
	int Sdim;
	int m; 
	int maxPvis;
	int maxCPvis;
	int beg;
	int end;
	double *V;
	double *W;
	double *E;
	double *ea;
	double *eb;
	double *S;
	double *U;
	struct sba_crsm idxij;
};

void *stub(void *ptr)
{
	double *V = ((struct threadParams*)ptr)->V;
	double *W = ((struct threadParams*)ptr)->W;
	double *E = ((struct threadParams*)ptr)->E;
	double *ea = ((struct threadParams*)ptr)->ea;
	double *eb = ((struct threadParams*)ptr)->eb;
	double *S = ((struct threadParams*)ptr)->S;
	double *U = ((struct threadParams*)ptr)->U;
	struct sba_crsm idxij = ((struct threadParams*)ptr)->idxij;	

	int mmconxUsz = ((struct threadParams*)ptr)->mmconxUsz;
	int Sdim = ((struct threadParams*)ptr)->Sdim;
	int m = ((struct threadParams*)ptr)->m;
	int maxPvis = ((struct threadParams*)ptr)->maxPvis;
	int maxCPvis = ((struct threadParams*)ptr)->maxCPvis;

	double *Yj = (double *)malloc(sizeof(double)*maxPvis*Ysz);
	double YWt[81];
	int *rcidxs = (int *)malloc(sizeof(int)*maxCPvis);
	int *rcsubs = (int *)malloc(sizeof(int)*maxCPvis);
	
	double *ptr1, *ptr2, *ptr3, *ptr4;
	int j;
	
	register int k, i, ii, jj, l;
	register double *pYWt;

	double sum;
	int nnz;

	int beg = ((struct threadParams *)ptr)->beg;
	int end = ((struct threadParams *)ptr)->end;

	for(j = beg; j < end; j++)
	{

		// need V, W, E, eab, S, mmconxUsz, Sdim
		memset(Yj, 0.0, maxPvis*Ysz*sizeof(double));
		
		nnz = sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs);
		for(i = 0; i < nnz; ++i)
		{
			ptr3 = V + rcsubs[i]*Vsz;
			ptr1 = Yj + i*Ysz;
			ptr2 = W + idxij.val[rcidxs[i]]*Wsz;
			
			for(ii = 0; ii < cnp; ++ii)
			{
				ptr4 = ptr2 + ii*pnp;
				for(jj = 0; jj < pnp; ++jj)
				{
					for(k = 0, sum = 0.0; k <= jj; ++k)
						sum += ptr4[k]*ptr3[jj*pnp+k];
					for( ; k < pnp; ++k)
						sum += ptr4[k]*ptr3[k*pnp+jj];
					ptr1[ii*pnp+jj] = sum;
				}
			}
		}
		for(k=j; k<m; ++k)
		{ 		  
			/* Recall that Y_ij is cnp x pnp and W_ik is 
			* cnp x pnp */ 
			memset(YWt, 0, YWtsz*sizeof(double));

			for(i=0; i<nnz; ++i)
			{
				ii=idxij.colidx[idxij.rowptr[rcsubs[i]]];
				jj=idxij.colidx[idxij.rowptr[rcsubs[i]+1]-1];
				if(k<ii || k>jj) continue; /* W_ik == 0 */

				/* set ptr2 to point to W_ik */
				l=sba_crsm_elmidxp(&idxij, rcsubs[i], k, j, rcidxs[i]);
		
				if(l==-1) continue; /* W_ik == 0 */

				ptr2=W + idxij.val[l]*Wsz;
				/* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
				ptr1=Yj + i*Ysz;

				for(ii=0; ii<cnp; ++ii)
				{
					ptr3=ptr1+ii*pnp;
					pYWt=YWt+ii*cnp;

					ptr4=ptr2;
					for(jj=0; jj<cnp; ++jj)
					{
						for(l=0, sum=0.0; l<pnp; ++l)
							sum+=ptr3[l]*ptr4[l];
						pYWt[jj]+=sum; 
						ptr4+=pnp;
					}
				}
			}
			ptr2=S + (k)*mmconxUsz + (j)*cnp; // set ptr2 to point to the beginning of block j,k in S
		 
			if(j!=k)
			{ /* Kronecker */
				for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
					for(jj=0; jj<cnp; ++jj)
						ptr2[jj]=-YWt[jj*cnp+ii];
			}
			else
			{
				ptr1=U + j*Usz; // set ptr1 to point to U_j

				for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
					for(jj=0; jj<cnp; ++jj)
						ptr2[jj] = ptr1[jj*cnp+ii] - YWt[jj*cnp+ii];
			}
		}

		/* compute e_j=ea_j - \sum_i Y_ij eb_i */
		/* Recall that Y_ij is cnp x pnp and eb_i is pnp x 1 */
		ptr1=E + j*easz; // set ptr1 to point to e_j

		for(i=0; i<nnz; ++i)
		{
			ptr2=Yj + i*Ysz;

			/* set ptr3 to point to eb_i */
			ptr3=eb + rcsubs[i]*ebsz;
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=ptr2+ii*pnp;
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr3[jj];
				ptr1[ii]+=sum;
			}
		}

		ptr2=ea + j*easz; // set ptr2 to point to ea_j
		for(i=0; i<easz; ++i)
			ptr1[i]=ptr2[i] - ptr1[i];
	}
	free(Yj);
	free(rcidxs);
	free(rcsubs);
	return;
}

int SEJ(int mmconxUsz, int Sdim, int m, int maxPvis, int maxCPvis, double *V, double *W, double *E, double *ea, double *eb, double *S, double *U, struct sba_crsm *idxij)
{
	
	int i;
	int num_t = NR_THREADS;
	if(m > NR_THREADS)
		num_t = NR_THREADS;
	else
		num_t = m;
	
	pthread_t thread[NR_THREADS];
	struct threadParams val[NR_THREADS];

	for(i = 0; i < num_t; i++)
	{
		val[i].mmconxUsz = mmconxUsz;
		val[i].Sdim = Sdim;
		val[i].m = m; 
		val[i].maxPvis = maxPvis;
		val[i].maxCPvis = maxCPvis;
		val[i].V = V;
		val[i].W = W;
		val[i].E = E;
		val[i].ea = ea;
		val[i].eb = eb;
		val[i].S = S;
		val[i].U = U;
		val[i].idxij = *idxij;
		val[i].beg = i*m/num_t;
		val[i].end = (i+1)*m/num_t;
		pthread_create(&thread[i], NULL, &stub, &(val[i]));

	}
	for(i = 0; i < num_t; i++)
		pthread_join(thread[i], NULL);
	return 0;
}

struct jac_threadParams
{
	int beg;
	int end;
	int maxCPvis;
	int mnp;
	int Asz;
	int ABsz;
	double *hxij;
	double *pb;
	double *pa;
	double *jac;
	int *rcsubs;
	int *rcidxs;
	struct sba_crsm *idxij;
	void *proj;
	void *adata;
};

void *stub_jac_AIJ(void *ptr)
{
	register int i, jj, j, ii;
	int nnz;
	double d, d1, tmp;
	int maxCPvis = ((struct jac_threadParams*)ptr)->maxCPvis;

	int *rcsubs = (int *)malloc(sizeof(int)*maxCPvis);
	int *rcidxs = (int *)malloc(sizeof(int)*maxCPvis);

	double *pbi, *paj;

	double *hxij;
	double *hxxij;
	double *pb = ((struct jac_threadParams*)ptr)->pb;
	double *pa = ((struct jac_threadParams*)ptr)->pa;
	double *jac = ((struct jac_threadParams*)ptr)->jac;
	double *pAB;
	
	int beg = ((struct jac_threadParams*)ptr)->beg;
	int end = ((struct jac_threadParams*)ptr)->end;

	void (*proj)(int j, int i, double *aj, double *bi, double *xij, void *adata);

	proj = ((struct jac_threadParams*)ptr)->proj;

	void *adata = ((struct jac_threadParams*)ptr)->adata;

	struct sba_crsm *idxij = ((struct jac_threadParams*)ptr)->idxij;
	
	int mnp = ((struct jac_threadParams*)ptr)->mnp;

	int ABsz = ((struct jac_threadParams*)ptr)->ABsz;
	int Asz = ((struct jac_threadParams*)ptr)->Asz;
	
	if((hxij=malloc(2*mnp*sizeof(double)))==NULL){
		fprintf(stderr, "memory allocation request failed in sba_motstr_Qs_fdjac()!\n");
		exit(1);
	}
	hxxij=hxij+mnp;

	/* compute A_ij */
	for(j=beg; j<end; ++j){
		paj=pa+j*cnp; // j-th camera parameters
		printf("j is %d\n", j);
		nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
		for(jj=0; jj<cnp; ++jj){
			d=(double)(SBA_DELTA_SCALE)*paj[jj]; // force evaluation
			d=FABS(d);
			if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
			d1=1.0/d; /* invert so that divisions can be carried out faster as multiplications */

			for(i=0; i<nnz; ++i){
				pbi=pb + rcsubs[i]*pnp; // i-th point parameters
				(*proj)(j, rcsubs[i], paj, pbi, hxij, adata); // evaluate supplied function on current solution

				tmp=paj[jj];
				paj[jj]+=d;
				(*proj)(j, rcsubs[i], paj, pbi, hxxij, adata);
				paj[jj]=tmp; /* restore */

				pAB=jac + idxij->val[rcidxs[i]]*ABsz; // set pAB to point to A_ij
				for(ii=0; ii<mnp; ++ii) {
					pAB[ii*cnp+jj]=(hxxij[ii]-hxij[ii])*d1;
				}
			}
		}
	}
	free(rcsubs);
	free(rcidxs);
	free(hxij);
}

void *stub_jac_BIJ(void *ptr)
{
	register int i, jj, j, ii;
	int nnz;
	double d, d1, tmp;
	int maxCPvis = ((struct jac_threadParams*)ptr)->maxCPvis;

	int *rcsubs = (int *)malloc(sizeof(int)*maxCPvis);
	int *rcidxs = (int *)malloc(sizeof(int)*maxCPvis);

	double *pbi, *paj;

	double *hxij;
	double *hxxij;
	double *pb = ((struct jac_threadParams*)ptr)->pb;
	double *pa = ((struct jac_threadParams*)ptr)->pa;
	double *jac = ((struct jac_threadParams*)ptr)->jac;
	double *pAB;
	
	int beg = ((struct jac_threadParams*)ptr)->beg;
	int end = ((struct jac_threadParams*)ptr)->end;

	void (*proj)(int j, int i, double *aj, double *bi, double *xij, void *adata);

	proj = ((struct jac_threadParams*)ptr)->proj;

	void *adata = ((struct jac_threadParams*)ptr)->adata;

	struct sba_crsm *idxij = ((struct jac_threadParams*)ptr)->idxij;
	
	int mnp = ((struct jac_threadParams*)ptr)->mnp;

	int ABsz = ((struct jac_threadParams*)ptr)->ABsz;
	int Asz = ((struct jac_threadParams*)ptr)->Asz;
	
	if((hxij=malloc(2*mnp*sizeof(double)))==NULL){
		fprintf(stderr, "memory allocation request failed in sba_motstr_Qs_fdjac()!\n");
		exit(1);
	}
	hxxij=hxij+mnp;

	/* compute B_ij */
	for(i=beg; i<end; ++i)
	{
		pbi=pb+i*pnp; // i-th point parameters

		nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
		for(jj=0; jj<pnp; ++jj)
		{
			/* determine d=max(SBA_DELTA_SCALE*|pbi[jj]|, SBA_MIN_DELTA), see HZ */
			d=(double)(SBA_DELTA_SCALE)*pbi[jj]; // force evaluation
			d=FABS(d);
			if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
			d1=1.0/d; /* invert so that divisions can be carried out faster as multiplications */

			for(j=0; j<nnz; ++j)
			{
				paj=pa + rcsubs[j]*cnp; // j-th camera parameters
				pthread_mutex_lock(&(global_lock[rcsubs[j]]));
				(*proj)(rcsubs[j], i, paj, pbi, hxij, adata); // evaluate supplied function on current solution
				tmp=pbi[jj];
				pbi[jj]+=d;
				(*proj)(rcsubs[j], i, paj, pbi, hxxij, adata);
				pthread_mutex_unlock(&(global_lock[rcsubs[j]]));
				pbi[jj]=tmp; /* restore */
				
				pAB=jac + idxij->val[rcidxs[j]]*ABsz + Asz; // set pAB to point to B_ij
				for(ii=0; ii<mnp; ++ii)
					pAB[ii*pnp+jj]=(hxxij[ii]-hxij[ii])*d1;
				}
		}
	}
	free(rcsubs);
	free(rcidxs);
	free(hxij);
}

int jac_threaded(double *p, struct sba_crsm *idxij, int maxCPvis, double *jac, void *dat)
{
	register int i, j, ii, jj;
	double *pa, *pb, *paj, *pbi;
	register double *pAB;
	int n, m, nnz, Asz, Bsz, ABsz;

	double tmp;
	register double d, d1;

	struct wrap_motstr_data_ *fdjd;
	void (*proj)(int j, int i, double *aj, double *bi, double *xij, void *adata);
	int mnp;
	void *adata;	

	global_lock = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t)*idxij->nc);

	for(i = 0; i < idxij->nc; i++)
		pthread_mutex_init(&(global_lock[i]),NULL);
	
	/* retrieve problem-specific information passed in *dat */
	fdjd=(struct wrap_motstr_data_ *)dat;
	proj=fdjd->proj;
	mnp=fdjd->mnp;
	adata=fdjd->adata;

	n=idxij->nr; m=idxij->nc;
	pa=p; pb=p+m*cnp;
	Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;

	int num_t = NR_THREADS;
	if(m > NR_THREADS)
		num_t = NR_THREADS;
	else
		num_t = m;
	
	pthread_t thread[NR_THREADS];
	struct jac_threadParams val[NR_THREADS];
	//num_t = 1;
	for(i = 0; i < num_t; i++)
	{
		val[i].beg = i*m/num_t;
		val[i].end = (i+1)*m/num_t;
		val[i].idxij = idxij;
		val[i].pa = pa;
		val[i].pb = pb;
		val[i].proj = proj;
		val[i].maxCPvis = maxCPvis;
		val[i].jac = jac; 
		val[i].adata = adata;
		val[i].mnp = mnp;
		val[i].ABsz = ABsz;
		val[i].Asz = Asz;
		pthread_create(&thread[i], NULL, &stub_jac_AIJ, &(val[i]));
	}
	for(i = 0; i < num_t; i++)
		pthread_join(thread[i], NULL);

	if(n > NR_THREADS)
		num_t = NR_THREADS;
	else
		num_t = n;
	for(i = 0; i < num_t; i++)
	{
		val[i].beg = i*n/num_t;
		val[i].end = (i+1)*n/num_t;
		val[i].idxij = idxij;
		val[i].pa = pa;
		val[i].pb = pb;
		val[i].proj = proj;
		val[i].maxCPvis = maxCPvis;
		val[i].jac = jac; 
		val[i].adata = adata;
		val[i].mnp = mnp;
		val[i].ABsz = ABsz;
		val[i].Asz = Asz;
		pthread_create(&thread[i], NULL, &stub_jac_BIJ, &(val[i]));
	}
	for(i = 0; i < num_t; i++)
		pthread_join(thread[i], NULL);
	
	for(i = 0; i < idxij->nc; i++)
		pthread_mutex_destroy(&(global_lock[i]));
	free(global_lock);
	return 0;
}
