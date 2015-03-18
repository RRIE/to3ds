#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include "sba.h"


#define NR_THREADS 8
#define BUFFER_SIZE 8

pthread_cond_t notempty, notfull, done;
pthread_mutex_t lock; 

int num_requests = 0;
int *buf = NULL;
int count = 0;
int bcast = NR_THREADS;

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


struct threadParams
{
	int mmconxUsz;
	int Sdim;
	int m; 
	int maxPvis;
	int maxCPvis;
	double *V;
	double *W;
	double *E;
	double *ea;
	double *eb;
	double *S;
	double *U;
	struct sba_crsm idxij;
};

int add_queue(int* queue, int max_length, int val)
{
	int i = 0;
	for(i = 0; i < max_length; i++)
	{
		if(queue[i] == -1)
		{
			queue[i] = val;
			return 0;
		}
	}
	return -1;
}

int remove_queue(int * queue, int max_length, int val)
{
	int i = 0;
	int j = 0;
	for(i = 0; i < max_length; i++)
	{
		if(queue[i] == val)
		{
			for(j = i+1; j < max_length; j++)
				queue[j-1] = queue [j];
			queue[max_length-1] = -1;
			return 0;
		}
	}
	return -1;
}

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
	double *YWt = (double *)malloc(sizeof(double)*YWtsz);
	int *rcidxs = (int *)malloc(sizeof(int)*maxCPvis);
	int *rcsubs = (int *)malloc(sizeof(int)*maxCPvis);
	
	double *ptr1, *ptr2, *ptr3, *ptr4;
	int j;
	
	register int k, i, ii, jj, l;
	register double *pYWt;

	double sum;
	int nnz;
	printf("[sba_threaded] in stub\n");

	
	while(1)
	{
		pthread_mutex_lock(&lock);
		while(num_requests == 0 && count < m)
			pthread_cond_wait(&notempty, &lock);
		if(count >= m && num_requests == 0)
		{
			pthread_mutex_unlock(&lock);
			goto DONEWORK;
		}
		j = buf[0];
		remove_queue(buf, BUFFER_SIZE, j);
		num_requests--;
		if(num_requests < BUFFER_SIZE)
			pthread_cond_signal(&notfull);
		pthread_mutex_unlock(&lock);

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
	DONEWORK: 
	pthread_mutex_lock(&lock);
	bcast--;
	pthread_cond_signal(&done);
	pthread_mutex_unlock(&lock);
	free(Yj);
	free(YWt);
	free(rcidxs);
	free(rcsubs);
	return;
}

int SEJ(int mmconxUsz, int Sdim, int m, int maxPvis, int maxCPvis, double *V, double *W, double *E, double *ea, double *eb, double *S, double *U, struct sba_crsm *idxij)
{
	buf = (int *)malloc(sizeof(int)*BUFFER_SIZE);
	pthread_t *thread = (pthread_t *)malloc(sizeof(pthread_t)*NR_THREADS);
	pthread_cond_init(&notempty, NULL);
	pthread_cond_init(&notfull, NULL); 
	pthread_cond_init(&done, NULL); 
	pthread_mutex_init(&lock, NULL);

	struct threadParams val; 
	val.mmconxUsz = mmconxUsz;
	val.Sdim = Sdim;
	val.m = m; 
	val.maxPvis = maxPvis;
	val.maxCPvis = maxCPvis;
	val.V = V;
	val.W = W;
	val.E = E;
	val.ea = ea;
	val.eb = eb;
	val.S = S;
	val.U = U;
	val.idxij = *idxij;
	int i;
	for(i = 0; i < BUFFER_SIZE; i++)
		buf[i] = -1;
	int num_t = NR_THREADS;
	if(m > NR_THREADS)
	{
		for(i = 0; i < NR_THREADS; i++)
			pthread_create(&thread[i],NULL, &stub, &val);
	}
	else
	{
		for(i = 0; i < m; i++)
			pthread_create(&thread[i],NULL, &stub, &val);
		bcast = m;
		num_t = m;
	}
	while(count < m)
	{
		pthread_mutex_lock(&lock);		
		while(num_requests == BUFFER_SIZE)
			pthread_cond_wait(&notfull, &lock);
		add_queue(buf, BUFFER_SIZE, count);
		count++;
		num_requests++;
		if(num_requests > 0)
			pthread_cond_signal(&notempty);
		pthread_mutex_unlock(&lock);
	}
	pthread_mutex_lock(&lock);
	pthread_cond_broadcast(&notempty);
	while(bcast != 0)
		pthread_cond_wait(&done, &lock);
	pthread_mutex_unlock(&lock);
	printf("[sba_threaded] joining threads\n");
	for(i = 0; i < num_t; i++)
		pthread_join(thread[i], NULL);
	free(buf);
	pthread_mutex_destroy(&lock);
	pthread_cond_destroy(&done);
	pthread_cond_destroy(&notempty);
	pthread_cond_destroy(&notfull);
	buf = NULL;
	num_requests = 0;
	bcast = NR_THREADS;
	count = 0;
	return 0;
}
