/*****************************************************************************/
/* SBA OpenCL include file						     */
/*****************************************************************************/

#include "CL/cl.h"

/* struct for storing OpenCL platform/device specific information */

struct opencl_info{
	cl_platform_id platform;
	cl_context context;
	cl_command_queue queue;
	cl_device_id device;
};

/* OpenCL setup function */
/* Returns the platform, context, command queue and device id in	     */
/* and opencl_info struct						     */

struct opencl_info sba_opencl_setup(void);
