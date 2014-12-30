/*
 * key_util.h
 *
 *  Created on: 2014-12-30
 *      Author: edward
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <zlib.h>
#include "keys2a.h"

#ifndef KEY_UTIL_H_
#define KEY_UTIL_H_

class KEY_UTIL
{
public:
	static int get_number_of_keys_normal(FILE *fp);
	static int get_number_of_keys_gzip(gzFile fp);
	static int get_number_of_keys(const char *filename);

	static int read_key_file_list(char* list_in, std::vector<std::string>& key_files);
	static int read_key_files_from_list(char* list_in,  std::vector<unsigned char*> & keys, std::vector<int> & num_keys );
	static int read_key_file(const char *filename, unsigned char **keys, keypt_t **info = NULL);
	static int read_keys(FILE *fp, unsigned char **keys, keypt_t **info = NULL);
	static int read_keys_gzip(gzFile fp, unsigned char **keys, keypt_t **info = NULL);
};


#endif /* KEY_UTIL_H_ */
