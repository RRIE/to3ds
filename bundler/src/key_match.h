/*
 * key_match.h
 *
 *  Created on: 2014-12-30
 *      Author: edward
 */

#ifndef KEY_MATCH_H_
#define KEY_MATCH_H_

#include "flann/flann.hpp"

#define DEFAULT_MATCHES							16
#define DEFAULT_SEARCH_RATIO					1.0
#define DEFAULT_WINDOW_RADIUS				-1
#define DEFAULT_DIST_RATIO 						0.6

#define DEFAULT_FLANN_CHECKS				16

#define KEY_DIM											128

// Some useful typedefs for working with FLANN
typedef  flann::Index < flann::L2<unsigned char>> KEY_INDEX_L2_CHAR;
typedef  flann::Matrix<unsigned char> MATRIX_CHAR;

// PARAMS FOR MATCH
class KEY_MATCH_PARAMS{
private:
	int m_min_matches;
	double m_search_ratio;

	int m_window_radius;
	double m_dist_ratio;

	int m_flann_checks;
public:
	KEY_MATCH_PARAMS(
			int min_matches = DEFAULT_MATCHES,
			double search_ratio = DEFAULT_SEARCH_RATIO,
			int window_radius = DEFAULT_WINDOW_RADIUS,
			double dist_ratio = DEFAULT_DIST_RATIO,
			int flann_checks = DEFAULT_FLANN_CHECKS) :
	 m_min_matches(min_matches), m_search_ratio(search_ratio), m_window_radius(window_radius), m_dist_ratio(dist_ratio),
	 m_flann_checks(flann_checks) {}

	// Getters
	int 				min_matches() 		{return m_min_matches;}
	double 		search_ratio() 		{return m_search_ratio;}
	int 				window_radius() 	{return m_window_radius;}
	double 		dist_ratio() 			{return m_dist_ratio;}
	int				flann_checks()		{return m_flann_checks;}

	// Setters
	void 			min_matches(int min_matches) 		{ m_min_matches = min_matches;}
	void 			search_ratio(double search_ratio) 	{ m_search_ratio = search_ratio;}
	void 			window_radius(int window_radius) 	{ m_window_radius = window_radius;}
	void 			dist_ratio(double dist_ratio) 			{ m_dist_ratio = dist_ratio;}
	void				flann_checks(int flann_checks)			{ m_flann_checks = flann_checks;}
};

class KEY_MATCH_RESULT
{
private:
	int m_img_id1;
	int m_img_id2;
	std::vector<KeypointMatch> m_matches;

public:
	KEY_MATCH_RESULT(int img_id1, int img_id2, std::vector<KeypointMatch> matches) :
		m_img_id1(img_id1), m_img_id2(img_id2), m_matches(matches){}

	inline void print( FILE *f)
	{
        fprintf(f, "%d %d\n", m_img_id1, m_img_id2);

        /* Write the number of matches */
        fprintf(f, "%d\n", (int) m_matches.size());

        for (int i = 0; i < m_matches.size(); i++) {
            fprintf(f, "%d %d\n",
            		m_matches[i].m_idx1, m_matches[i].m_idx2);
        }
	}

	inline int img_id1() {return m_img_id1;}
	inline int img_id2() {return m_img_id2;}


};


class KEY_MATCHER
{
public:
	static std::vector<KEY_MATCH_RESULT> match_keys_v1( std::vector<unsigned char*> & keys, std::vector<int> & num_keys, KEY_MATCH_PARAMS match_params = KEY_MATCH_PARAMS());
	static std::vector<KEY_MATCH_RESULT> match_keys_v2(std::vector<unsigned char*> & keys, std::vector<int> & num_keys, KEY_MATCH_PARAMS match_params = KEY_MATCH_PARAMS());
	static std::vector<KeypointMatch> search( unsigned char * keys, int num_keys, int key_dim, KEY_INDEX_L2_CHAR & key_index, KEY_MATCH_PARAMS match_params, int n_query, int n_index );
};

#endif /* KEY_MATCH_H_ */
