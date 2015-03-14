/*
 * key_match.cpp
 *
 *  Created on: 2014-12-30
 *      Author: edward
 */

#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>


#include "key_util.h"
#include "keys2a.h"
#include "key_match.h"

#include "tbb/parallel_for.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/task_scheduler_init.h"




typedef tbb::enumerable_thread_specific< std::vector<KEY_MATCH_RESULT> > KEY_MATCH_RESULT_TLS;

/********************************************************************************************************
 *  Match Keys algorithm
 *  	Run Parallel For on outerloop of original match keys algorithm
 *  	Outerloop will make a FLANN Structure
 *  	Innerloop will do matching
 *  	Good for small sets and / or limited number of cores
 ********************************************************************************************************/
std::vector<KEY_MATCH_RESULT> KEY_MATCHER::match_keys_v1(
		std::vector<unsigned char*> & keys,
		std::vector<int> & num_keys,
		KEY_MATCH_PARAMS match_params)
{
	tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());  // Explicit number of threads
	std::vector<KEY_MATCH_RESULT> matches;
	KEY_MATCH_RESULT_TLS key_match_result_tls (matches);

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0,num_keys.size()),
		[&num_keys, &keys, &match_params, &key_match_result_tls](const tbb::blocked_range<size_t>& r) {
		for (int i = r.begin(); i < r.end(); i++){
			if(num_keys[i] == 0) continue;
			// set up the index that will be queried

			MATRIX_CHAR key_points( keys[i], num_keys[i], KEY_DIM);
			KEY_INDEX_L2_CHAR key_index(key_points,  flann::KDTreeIndexParams(1));

			key_index.buildIndex();

			// Compute the start index for inner loop
			int start_idx = 0;
			//Commenting out Window Radius
			//if (match_params.window_radius() > 0) start_idx = std::max( (int) i - match_params.window_radius(), 0);

			//set up the query
			for(int j = start_idx; j < i; j++)
			{
				if(num_keys[j] == 0) continue;

				std::vector<KeypointMatch> matches = KEY_MATCHER::search(keys[j], num_keys[j], KEY_DIM,  key_index, match_params, j, i);

				// push the matches
				int num_matches = (int) matches.size();

				if (num_matches >= match_params.min_matches()) {
					KEY_MATCH_RESULT result(j, i, matches);
					key_match_result_tls.local().push_back(result);
				}
			}
		}
	});


	std::vector<KEY_MATCH_RESULT> match_results = key_match_result_tls.combine (
			[](std::vector<KEY_MATCH_RESULT>x, std::vector<KEY_MATCH_RESULT> y)
			{ x.insert( x.end(), y.begin(), y.end() ); return x; }
	);

	return match_results;

}


/********************************************************************************************************
 *  Match Keys algorithm
 *  	- Run Parallel For  in fine grain
 *  	- Overhead step of allocating FLANN Indexes for all images (time / memory)
 *  	- Parallel For will do matching. Search tree with larger number of keys for marginal speedup
 *  	- Good for large sets and / or many cores / threads available
 ********************************************************************************************************/
std::vector<KEY_MATCH_RESULT> KEY_MATCHER::match_keys_v2(
		std::vector<unsigned char*> & keys,
		std::vector<int> & num_keys,
		KEY_MATCH_PARAMS match_params)
{
		tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());  // Explicit number of threads
		std::vector<KEY_MATCH_RESULT> matches;
		KEY_MATCH_RESULT_TLS key_match_result_tls (matches);

		std::vector<KEY_INDEX_L2_CHAR> key_index_db;
		printf("[KeyMatchFull] Building Index...\n");
		// Build All the Indexes
		for( int iter = 0; iter < num_keys.size(); iter++)
		{
			MATRIX_CHAR key_points( keys[iter], num_keys[iter], KEY_DIM);
			KEY_INDEX_L2_CHAR key_index(key_points,  flann::KDTreeIndexParams(1));
			key_index.buildIndex();
			key_index_db.push_back(key_index);
		}

		printf("[KeyMatchFull] Searching...\n");
		int num_checks = (num_keys.size()) * (num_keys.size()-1) / 2;

		tbb::parallel_for(
			tbb::blocked_range<size_t>(0,num_checks),
			[&num_keys, &keys, &key_index_db, &match_params, &key_match_result_tls](const tbb::blocked_range<size_t>& r) {

			for (int n = r.begin(); n < r.end(); n++){

				// we need to determine the pair to match (i,j)
				// map single number n to (i,j) (consecutive sums / combinatorics)
				int m = 0;
				int i = 0; int j = 0;
				while( m < (n+1))
				{
					m+= (i+1);
					i++;
				}

				j = n - ( (i * ( i - 1) ) / 2 );

				// now we will pick the index to query and the query index
				int n_index = (num_keys[j] < num_keys[i]) ? i : j;
				int n_query = (num_keys[j] < num_keys[i]) ? j : i;

				// skip in nothing to search
				if(num_keys[n_query] == 0) continue;

				std::vector<KeypointMatch> matches = KEY_MATCHER::search(keys[n_query], num_keys[n_query], KEY_DIM, key_index_db[n_index], match_params, n_query, n_index);

				int num_matches = (int) matches.size();

				if (num_matches >= match_params.min_matches()) {
					if(n_query < n_index) key_match_result_tls.local().push_back(KEY_MATCH_RESULT(n_query, n_index, matches));
					else key_match_result_tls.local().push_back(KEY_MATCH_RESULT(n_index, n_query, matches));
				}
			}
		});

		std::vector<KEY_MATCH_RESULT> match_results = key_match_result_tls.combine (
				[](std::vector<KEY_MATCH_RESULT>x, std::vector<KEY_MATCH_RESULT> y)
				{ x.insert( x.end(), y.begin(), y.end() ); return x; }
		);

		std::sort( match_results.begin(), match_results.end(),
				[](KEY_MATCH_RESULT a, KEY_MATCH_RESULT b)
				{
					if(a.img_id2() < b.img_id2()) return 1;
					if(a.img_id2() == b.img_id2() && a.img_id1() < b.img_id1()) return 1;
					return 0;
				});


		return match_results;

}


/********************************************************************************************************
 *  Match Keys Search
 * 		Will run search FLANN Key Index.
 * 		Runs two searches.
 * 			1. First search  is used to test the sets of keys to see if they will match in probability
 * 			2. Second search performed if first search returns true
 ********************************************************************************************************/
std::vector<KeypointMatch> KEY_MATCHER::search(
		unsigned char * keys, int num_keys,
		int key_dim, KEY_INDEX_L2_CHAR & key_index,
		KEY_MATCH_PARAMS match_params,
		int n_query,
		int n_index)
{
	assert(match_params.search_ratio() > 0 && match_params.search_ratio() <= 1.0);

	std::vector<KeypointMatch> matches;

	// general search params
	int nn = 2; // number of nearest neighbors to search for
	size_t * index_ptr;
	float * dist_ptr;

	// params for searches search
	int num_keys1 = (int) num_keys * match_params.search_ratio();


	MATRIX_CHAR query1( keys, num_keys1, key_dim); // search 1

	index_ptr= new size_t[query1.rows * nn];
	dist_ptr = new float[query1.rows * nn];
	flann::Matrix<size_t> indices(index_ptr, query1.rows * nn, nn);
	flann::Matrix<float> dists(dist_ptr, query1.rows * nn, nn);

	// search
	key_index.knnSearch(query1, indices, dists, nn, flann::SearchParams(match_params.flann_checks()));

	for( int k = 0; k < query1.rows; k++)
	{
			if( dists[k][0] < match_params.dist_ratio()*match_params.dist_ratio()* dists[k][1])
			{
				if(n_query < n_index)	matches.push_back(KeypointMatch(k, (int) indices[k][0]));
				else matches.push_back(KeypointMatch((int) indices[k][0], k));
			}
	}

	delete[] index_ptr;
	delete[] dist_ptr;

	// push the matches
	int num_matches = (int) matches.size();

	// do the second search if first search meets the expected minimum and has less keys than the total
	// we don't want to search twice if we have search_ratio = 1.0
	if (num_matches >= match_params.min_matches() * match_params.search_ratio() && num_keys1 < num_keys) {

		int num_keys2 = num_keys - num_keys1;
		MATRIX_CHAR query2( keys + sizeof(unsigned char) * key_dim * num_keys1, num_keys2, key_dim); // search 2

		index_ptr= new size_t[query2.rows * nn];
		dist_ptr = new float[query2.rows * nn];

		flann::Matrix<size_t> indices(index_ptr, query2.rows * nn, nn);
		flann::Matrix<float> dists(dist_ptr, query2.rows * nn, nn);

		// search
		key_index.knnSearch(query2, indices, dists, nn, flann::SearchParams(match_params.flann_checks()));

		for( int k = 0; k < query2.rows; k++)
		{
				if( dists[k][0] < match_params.dist_ratio()*match_params.dist_ratio()* dists[k][1])
				{
					// we add an offset to the index based on first search
					if(n_query < n_index)	matches.push_back(KeypointMatch(num_keys1 + k, (int) indices[k][0]));
					else matches.push_back(KeypointMatch( (int) indices[k][0], num_keys1 + k));
				}
		}

		delete[] index_ptr;
		delete[] dist_ptr;
	}
	else if(num_matches < match_params.min_matches() * match_params.search_ratio())
	{
		matches.clear();
	}

	return matches;
}
