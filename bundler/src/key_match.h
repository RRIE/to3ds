/*
 * key_match.h
 *
 *  Created on: 2014-12-30
 *      Author: edward
 */

#ifndef KEY_MATCH_H_
#define KEY_MATCH_H_

class KEY_MATCH_RESULT
{
public:
	KEY_MATCH_RESULT(int img_id1, int img_id2, std::vector<KeypointMatch> matches) :
		m_img_id1(img_id1), m_img_id2(img_id1), m_matches(matches){}

	inline void print( FILE *f)
	{
        fprintf(f, "%d %d\n", m_img_id2, m_img_id1);

        /* Write the number of matches */
        fprintf(f, "%d\n", (int) m_matches.size());

        for (int i = 0; i < m_matches.size(); i++) {
            fprintf(f, "%d %d\n",
            		m_matches[i].m_idx1, m_matches[i].m_idx2);
        }
	}

private:
	int m_img_id1;
	int m_img_id2;
	std::vector<KeypointMatch> m_matches;
};


class KEY_MATCHER
{
public:
	static std::vector<KEY_MATCH_RESULT> match_keys( std::vector<unsigned char*> & keys, std::vector<int> & num_keys, int window_radius = -1, double ratio = 0.6);
};;

#endif /* KEY_MATCH_H_ */
