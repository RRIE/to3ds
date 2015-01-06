/*
 * key_match.cpp
 *
 *  Created on: 2014-12-30
 *      Author: edward
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>

#include "key_util.h"
#include "keys2a.h"
#include "key_match.h"

#include "tbb/parallel_for.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/task_scheduler_init.h"

#include "flann/flann.hpp"


typedef tbb::enumerable_thread_specific< std::vector<KEY_MATCH_RESULT> > KEY_MATCH_RESULT_TLS;

std::vector<KEY_MATCH_RESULT> KEY_MATCHER::match_keys(std::vector<unsigned char*> & keys, std::vector<int> & num_keys, int window_radius, double ratio)
{
	tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());  // Explicit number of threads
	std::vector<KEY_MATCH_RESULT> matches;
	KEY_MATCH_RESULT_TLS key_match_result_tls (matches);

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0,num_keys.size()),
		[&num_keys, &keys, &window_radius, &ratio, &key_match_result_tls](const tbb::blocked_range<size_t>& r) {
		for (int i = r.begin(); i < r.end(); i++){
			if(num_keys[i] == 0) continue;
			// set up the index that will be queried
			flann::Matrix<unsigned char> key_points( keys[i], num_keys[i], 128);
			flann::Index < flann::L2<unsigned char>> key_index(key_points,  flann::KDTreeIndexParams(1));
			//flann::Index < flann::L2<unsigned char>> key_index(key_points,  flann::LshIndexParams(1,20,2));
			//flann::Index < flann::L2<unsigned char>> key_index(key_points,  flann::AutotunedIndexParams(0.6,0.01,0,0.1));
			key_index.buildIndex();

			// Compute the start index for inner loop
			int start_idx = 0;
			if (window_radius > 0) start_idx = std::max( (int) i - window_radius, 0);

			//set up the query
			for(int j = start_idx; j < i; j++)
			{
				if(num_keys[j] == 0) continue;


				int nn = 2; // number of nearest neighbors to search for

				flann::Matrix<unsigned char> query( keys[j], num_keys[j], 128);

				size_t * index_ptr = new size_t[query.rows * nn];
				float * dist_ptr = new float[query.rows * nn];
				flann::Matrix<size_t> indices(index_ptr, query.rows * nn, nn);
				flann::Matrix<float> dists(dist_ptr, query.rows * nn, nn);

				// search

				key_index.knnSearch(query, indices, dists, nn, flann::SearchParams(32));
				std::vector<KeypointMatch> matches;

				for( int k = 0; k < query.rows; k++)
				{
						if( dists[k][0] < ratio*ratio* dists[k][1])
						{
							matches.push_back(KeypointMatch(k, (int) indices[k][0]));
						}
				}

				// push the matches
				int num_matches = (int) matches.size();

				if (num_matches >= 16) {
					KEY_MATCH_RESULT result(j, i, matches);
					key_match_result_tls.local().push_back(result);
				}
				delete[] index_ptr;
				delete[] dist_ptr;

			}
		}
	});


	std::vector<KEY_MATCH_RESULT> match_results = key_match_result_tls.combine (
			[](std::vector<KEY_MATCH_RESULT>x, std::vector<KEY_MATCH_RESULT> y)
			{ x.insert( x.end(), y.begin(), y.end() ); return x; }
	);

	return match_results;

}
