/* 
*  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
*    and the University of Washington
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

/* keys2.cpp */
/* Class for SIFT keypoints */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <zlib.h>

#include "keys2a.h"




/* Create a search tree for the given set of keypoints */
ANNkd_tree *CreateSearchTree(int num_keys, unsigned char *keys)
{
    // clock_t start = clock();

    /* Create a new array of points */
    ANNpointArray pts = annAllocPts(num_keys, 128);

    for (int i = 0; i < num_keys; i++) {
        memcpy(pts[i], keys + 128 * i, sizeof(unsigned char) * 128);
    }

    /* Create a search tree for k2 */
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys, 128, 16);
    // clock_t end = clock();

    // printf("Building tree took %0.3fs\n", 
    //        (end - start) / ((double) CLOCKS_PER_SEC));    

    return tree;
}

std::vector<KeypointMatch> MatchKeys(int num_keys1, unsigned char *k1,
                                     ANNkd_tree *tree2,
                                     double ratio, int max_pts_visit)
{
    annMaxPtsVisit(max_pts_visit);
    std::vector<KeypointMatch> matches;

    /* Now do the search */
    // clock_t start = clock();
    for (int i = 0; i < num_keys1; i++) {
        ANNidx nn_idx[2];
        ANNdist dist[2];

        tree2->annkPriSearch(k1 + 128 * i, 2, nn_idx, dist, 0.0);

        if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) {
            matches.push_back(KeypointMatch(i, nn_idx[0]));
        }
    }
    // clock_t end = clock();

    // printf("Searching tree took %0.3fs\n", 
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    return matches;    
}

/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatch> MatchKeys(int num_keys1, unsigned char *k1,
                                     int num_keys2, unsigned char *k2,
                                     double ratio, int max_pts_visit) 
{
    annMaxPtsVisit(max_pts_visit);

    int num_pts = 0;
    std::vector<KeypointMatch> matches;

    num_pts = num_keys2;
    clock_t start = clock();

    /* Create a new array of points */
    ANNpointArray pts = annAllocPts(num_pts, 128);

    for (int i = 0; i < num_pts; i++) {
        memcpy(pts[i], k2 + 128 * i, sizeof(unsigned char) * 128);
    }

    /* Create a search tree for k2 */
    ANNkd_tree *tree = new ANNkd_tree(pts, num_pts, 128, 16);
    clock_t end = clock();

    // printf("Building tree took %0.3fs\n", 
    //	      (end - start) / ((double) CLOCKS_PER_SEC));

    /* Now do the search */
    start = clock();
    for (int i = 0; i < num_keys1; i++) {
        ANNidx nn_idx[2];
        ANNdist dist[2];

        tree->annkPriSearch(k1 + 128 * i, 2, nn_idx, dist, 0.0);

        if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) {
            matches.push_back(KeypointMatch(i, nn_idx[0]));
        }
    }
    end = clock();
    // printf("Searching tree took %0.3fs\n", 
    //        (end - start) / ((double) CLOCKS_PER_SEC));

    /* Cleanup */
    annDeallocPts(pts);
    // annDeallocPt(axis_weights);

    delete tree;

    return matches;
}
