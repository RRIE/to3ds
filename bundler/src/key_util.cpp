#include "key_util.h"

/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::get_number_of_keys_normal(FILE *fp)
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    int num, len;

    if (fscanf(fp, "%d %d", &num, &len) != 2) {
        printf("Invalid keypoint file.\n");
        return 0;
    }

    return num;
}

/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::get_number_of_keys_gzip(gzFile fp)
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    int num, len;

    char header[256];
    gzgets(fp, header, 256);

    if (sscanf(header, "%d %d", &num, &len) != 2) {
        printf("Invalid keypoint file.\n");
        return 0;
    }

    return num;
}

/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::get_number_of_keys(const char *filename)
//			returns the number of keys in a file
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    FILE *file;

    file = fopen (filename, "r");
    if (! file) {
        /* Try to file a gzipped keyfile */
        char buf[1024];
        sprintf(buf, "%s.gz", filename);
        gzFile gzf = gzopen(buf, "rb");

        if (gzf == NULL) {
            printf("Could not open file: %s\n", filename);
            return 0;
        } else {
            int n = KEY_UTIL::get_number_of_keys_gzip(gzf);
            gzclose(gzf);
            return n;
        }
    }

    int n = KEY_UTIL::get_number_of_keys_normal(file);
    fclose(file);
    return n;
}


/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::read_key_file_list(char* list_in, std::vector<std::string>& key_files)
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    FILE* fp;

    if ((fp = fopen(list_in, "r")) == NULL) {
        printf("Error opening file %s for reading.\n", list_in);
        return 1;
    }

    char buf[512], *start;
    while (fgets(buf, 512, fp)) {
        // Remove trailing new-line
        if (buf[strlen(buf) - 1] == '\n') buf[strlen(buf) - 1] = '\0';

        // Find first non-space character
        start = buf;
        while(isspace(*start)) start++;

        // Skip empty lines
        if (strlen(start) == 0) continue;

        // Append file-name to key_files
        key_files.push_back(std::string(buf));
    }

    // Check we found input files
    if (key_files.size() == 0) {
        printf("No input files found in %s.\n", list_in);
        return 1;
    }

    return 0;
}

/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::read_key_files_from_list(char* list_in,  std::vector<unsigned char*> & keys, std::vector<int> & num_keys )
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{

    /* Read the list of files */
    std::vector<std::string> key_files;
    if (KEY_UTIL::read_key_file_list(list_in, key_files) != 0) return EXIT_FAILURE;
    int num_images = (int) key_files.size();


    keys.clear();
    num_keys.clear();

    keys.resize(num_images);
    num_keys.resize(num_images);

    /* Read all keys */
    for (int i = 0; i < num_images; i++) {
        keys[i] = NULL;
        num_keys[i] = KEY_UTIL::read_key_file(key_files[i].c_str(), &keys[i]);
    }

	return 0;
}


/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::read_key_file(const char *filename, unsigned char **keys, keypt_t **info)
//	This reads a keypoint file from a given filename and returns the list
//	of keypoints.
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    FILE *file;

    file = fopen (filename, "r");
    if (! file) {
        /* Try to file a gzipped keyfile */
        char buf[1024];
        sprintf(buf, "%s.gz", filename);
        gzFile gzf = gzopen(buf, "rb");

        if (gzf == NULL) {
            printf("Could not open file: %s\n", filename);
            return 0;
        } else {
            int n = KEY_UTIL::read_keys_gzip(gzf, keys, info);
            gzclose(gzf);
            return n;
        }
    }

    int n = KEY_UTIL::read_keys(file, keys, info);
    fclose(file);
    return n;
}

/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::read_keys(FILE *fp, unsigned char **keys, keypt_t **info)
// 	Read keypoints from the given file pointer and return the list of
// 	keypoints.  The file format starts with 2 integers giving the total
//	number of keypoints and the size of descriptor vector for each
//	keypoint (currently assumed to be 128). Then each keypoint is
//	specified by 4 floating point numbers giving subpixel row and
//	column location, scale, and orientation (in radians from -PI to
//	PI).  Then the descriptor vector for each keypoint is given as a
//	list of integers in range [0,255]. */
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    int i, num, len;

    std::vector<Keypoint *> kps;

    if (fscanf(fp, "%d %d", &num, &len) != 2) {
        printf("Invalid keypoint file\n");
        return 0;
    }

    if (len != 128) {
        printf("Keypoint descriptor length invalid (should be 128).");
        return 0;
    }

    *keys = new unsigned char[128 * num + 8];

    if (info != NULL)
        *info = new keypt_t[num];

    unsigned char *p = *keys;
    for (i = 0; i < num; i++) {
        /* Allocate memory for the keypoint. */
        // short int *d = new short int[128];
        float x, y, scale, ori;

        if (fscanf(fp, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
            printf("Invalid keypoint file format.");
            return 0;
        }

        if (info != NULL) {
            (*info)[i].x = x;
            (*info)[i].y = y;
            (*info)[i].scale = scale;
            (*info)[i].orient = ori;
        }

        char buf[1024];
        for (int line = 0; line < 7; line++) {
            fgets(buf, 1024, fp);

            if (line < 6) {
                sscanf(buf,
                    "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu "
                    "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
                    p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, p+8, p+9,
                    p+10, p+11, p+12, p+13, p+14,
                    p+15, p+16, p+17, p+18, p+19);

                p += 20;
            } else {
                sscanf(buf,
                    "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
                    p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7);
                p += 8;
            }
        }
    }

    return num; // kps;
}

/////////////////////////////////////FUNCTION HEADER START/////////////////////////////////////////////
int KEY_UTIL::read_keys_gzip(gzFile fp, unsigned char **keys, keypt_t **info)
/////////////////////////////////////FUNCTION HEADER END///////////////////////////////////////////////
{
    int i, num, len;

    std::vector<Keypoint *> kps;
    char header[256];
    gzgets(fp, header, 256);

    if (sscanf(header, "%d %d", &num, &len) != 2) {
        printf("Invalid keypoint file.\n");
        return 0;
    }

    if (len != 128) {
        printf("Keypoint descriptor length invalid (should be 128).");
        return 0;
    }

    *keys = new unsigned char[128 * num + 8];

    if (info != NULL)
        *info = new keypt_t[num];

    unsigned char *p = *keys;
    for (i = 0; i < num; i++) {
        /* Allocate memory for the keypoint. */
        // short int *d = new short int[128];
        float x, y, scale, ori;
        char buf[1024];
        gzgets(fp, buf, 1024);

        if (sscanf(buf, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
            printf("Invalid keypoint file format.");
            return 0;
        }

        if (info != NULL) {
            (*info)[i].x = x;
            (*info)[i].y = y;
            (*info)[i].scale = scale;
            (*info)[i].orient = ori;
        }

        for (int line = 0; line < 7; line++) {
            char *str = gzgets(fp, buf, 1024);
            assert(str != Z_NULL);

            if (line < 6) {
                sscanf(buf,
                    "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu "
                    "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
                    p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, p+8, p+9,
                    p+10, p+11, p+12, p+13, p+14,
                    p+15, p+16, p+17, p+18, p+19);

                p += 20;
            } else {
                sscanf(buf,
                    "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
                    p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7);
                p += 8;
            }
        }
    }

    assert(p == *keys + 128 * num);

    return num; // kps;
}


