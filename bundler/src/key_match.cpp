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


typedef tbb::enumerable_thread_specific< std::vector<KEY_MATCH_RESULT> > KEY_MATCH_RESULT_TLS;

std::vector<KEY_MATCH_RESULT> KEY_MATCHER::match_keys(std::vector<unsigned char*> & keys, std::vector<int> & num_keys, int window_radius, double ratio)
{
	tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());  // Explicit number of threads

	std::vector<KEY_MATCH_RESULT> matches;
	KEY_MATCH_RESULT_TLS key_match_result_tls (matches);

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0,num_keys.size()),
		[&num_keys, &keys, &window_radius, &ratio, &key_match_result_tls](const tbb::blocked_range<size_t>& r) {
			for (size_t i=r.begin();i<r.end();++i){
				if (num_keys[i] == 0) continue;
				/* Create a tree from the keys */
				ANNkd_tree *tree = CreateSearchTree(num_keys[i], keys[i]);

				/* Compute the start index */
				int start_idx = 0;
				if (window_radius > 0) start_idx = std::max( (int) i - window_radius, 0);

				for (int j = start_idx; j < i; j++) {
					if (num_keys[j] == 0) continue;
					/* Compute likely matches between two sets of keypoints */
					std::vector<KeypointMatch> matches = MatchKeys(num_keys[j], keys[j], tree, ratio);

					int num_matches = (int) matches.size();

					if (num_matches >= 16) {
						KEY_MATCH_RESULT result(j, i, matches);
						key_match_result_tls.local().push_back(result);
					}
				}
				delete tree;
			}
		}
	);

	std::vector<KEY_MATCH_RESULT> match_results = key_match_result_tls.combine (
			[](std::vector<KEY_MATCH_RESULT>x, std::vector<KEY_MATCH_RESULT> y)
			{ x.insert( x.end(), y.begin(), y.end() ); return x; }
	);

	return match_results;

}
