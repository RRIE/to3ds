/*****************************************************************************/
/* SBA OpenCL include file						     */
/*****************************************************************************/

#include "CL/cl.h"

/* struct for storing OpenCL platform/device specific information 	     */

struct opencl_info{
	cl_platform_id platform;
	cl_context context;
	cl_command_queue queue;
	cl_device_id device;
	cl_kernel kernel[3];
};

/* Extracts the contents of the kernel and returns them in a string 	     */

char * extract_kernel(char *kernel_file);

/* OpenCL setup function 						     */
/* Returns the platform, context, command queue and device id in	     */
/* and opencl_info struct						     */

struct opencl_info sba_opencl_setup(void);

/* Creates an opencl kernel and stores it in info			     */
/* kernel_file: file name of the kernel					     */
/* kernel_name: name of the kernel					     */
/* kernel_index: index of the kernel, currently 0-2			     */

void sba_create_kernel(struct opencl_info *info, char *kernel_file, char *kernel_name, int kernel_index);

