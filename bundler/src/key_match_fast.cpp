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
#include "key_match.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
    // Ratio Used For Matching (Replace with Macro?)
    double ratio = 0.6;

    char *list_in;
    char *file_out;

    int window_radius = -1;

    // Process Arguments
    if (argc != 3 && argc != 4) {
        printf("Usage: %s <list.txt> <outfile> [window_radius]\n", argv[0]);
        return EXIT_FAILURE;
    }

    list_in = argv[1];
    file_out = argv[2];

    if (argc == 4) {
        window_radius = atoi(argv[3]);
    }

    // read in the keys
    clock_t start = clock();

    std::vector<unsigned char*> keys;
    std::vector<int> num_keys;
    if(KEY_UTIL::read_key_files_from_list(list_in, keys, num_keys) != 0) return EXIT_FAILURE;

    clock_t end = clock();
    printf("[KeyMatchFull] Reading keys took %0.3fs\n",  (end - start) / ((double) CLOCKS_PER_SEC));

    start = clock();
    std::vector<KEY_MATCH_RESULT> match_results =  KEY_MATCHER::match_keys_v2(keys, num_keys);
    end = clock();
    printf("[KeyMatchFull] Matching keys took %0.3fs\n",  (end - start) / ((double) CLOCKS_PER_SEC));
    // dump the results
    FILE *f;
    if ((f = fopen(file_out, "w")) == NULL) {
        printf("Could not open %s for writing.\n", file_out);
        return EXIT_FAILURE;
    }

    for(int i = 0; i < match_results.size(); i++)
    {
    	match_results[i].print(f);
    }
    fclose(f);

    /* Free keypoints */
    for (int i = 0; i < keys.size(); i++) {
        if (keys[i] != NULL)
            delete [] keys[i];
    }

    return EXIT_SUCCESS;
}

