/*
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* KeyMatchFull.cpp */
/* Read in keys, match, write results to a file */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "key_util.h"
#include "keys2a.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include <iostream>
#include <vector>


struct mytask {
  mytask(size_t n)
    :_n(n)
  {}
  void operator()() {
    for (int i=0;i<1000000;++i) {}  // Deliberately run slow
    std::cerr << "[" << _n << "]" << "\n";
  }
  size_t _n;
};

void parallel_for_test(){
	tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());  // Explicit number of threads

	std::vector<mytask> tasks;
	for (int i=0;i<16;++i)
	tasks.push_back(mytask(i));

	tbb::parallel_for(
			tbb::blocked_range<size_t>(0,tasks.size()),
		[&tasks](const tbb::blocked_range<size_t>& r) {
		for (size_t i=r.begin();i<r.end();++i) tasks[i]();
		}
	);
}

int main(int argc, char **argv) {
	//parallel_for_test();
    char *list_in;
    char *file_out;

    // Ratio Used For Matching (Replace with Macro?)
    double ratio;
    ratio = 0.6;

    if (argc != 3 && argc != 4) {
        printf("Usage: %s <list.txt> <outfile> [window_radius]\n", argv[0]);
        return EXIT_FAILURE;
    }

    list_in = argv[1];
    file_out = argv[2];
    FILE *f;
    if ((f = fopen(file_out, "w")) == NULL) {
        printf("Could not open %s for writing.\n", file_out);
        return EXIT_FAILURE;
    }

    int window_radius = -1;
    if (argc == 4) {
        window_radius = atoi(argv[3]);
    }

    clock_t start = clock();
    std::vector<unsigned char*> keys;
    std::vector<int> num_keys;
    if(KEY_UTIL::read_key_files_from_list(list_in, keys, num_keys) != 0) return EXIT_FAILURE;
    clock_t end = clock();
    printf("[KeyMatchFull] Reading keys took %0.3fs\n",
           (end - start) / ((double) CLOCKS_PER_SEC));

    for (int i = 0; i < num_keys.size(); i++) {
        if (num_keys[i] == 0)
            continue;

        printf("[KeyMatchFull] Matching to image %d\n", i);

        start = clock();

        /* Create a tree from the keys */
        ANNkd_tree *tree = CreateSearchTree(num_keys[i], keys[i]);

        /* Compute the start index */
        int start_idx = 0;
        if (window_radius > 0)
            start_idx = std::max(i - window_radius, 0);

        for (int j = start_idx; j < i; j++) {
            if (num_keys[j] == 0)
                continue;

            /* Compute likely matches between two sets of keypoints */
            std::vector<KeypointMatch> matches =
                MatchKeys(num_keys[j], keys[j], tree, ratio);

            int num_matches = (int) matches.size();

            if (num_matches >= 16) {
                /* Write the pair */
                fprintf(f, "%d %d\n", j, i);

                /* Write the number of matches */
                fprintf(f, "%d\n", (int) matches.size());

                for (int i = 0; i < num_matches; i++) {
                    fprintf(f, "%d %d\n",
                            matches[i].m_idx1, matches[i].m_idx2);
                }
            }
        }

        end = clock();
        printf("[KeyMatchFull] Matching took %0.3fs\n",
               (end - start) / ((double) CLOCKS_PER_SEC));
        fflush(stdout);

        delete tree;
    }

    /* Free keypoints */
    for (int i = 0; i < keys.size(); i++) {
        if (keys[i] != NULL)
            delete [] keys[i];
    }

    fclose(f);
    return EXIT_SUCCESS;
}

