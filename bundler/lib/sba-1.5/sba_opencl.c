#include "stdio.h"
#include "string.h"
#include "sba_opencl.h"

#define SUCCESS 0
#define FAILURE 1

struct opencl_info sba_opencl_setup(void)
{
	/* OpenCL device variables */
	struct opencl_info info;
	cl_int error = 0;
	cl_context context;
	cl_command_queue queue;
	cl_device_id device;
	cl_uint numPlatforms;	
	cl_platform_id platform = NULL;
	
	/* Platform Information */
	char *platform_name = NULL;
	char *platform_vendor = NULL;
	size_t size;

	error = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (error != CL_SUCCESS)
	{
		printf("Error getting platforms\n");
		exit(error);
	}

	printf("[sba_opencl] number of platforms is: %d\n", numPlatforms);
	if(numPlatforms > 0)
	{
		cl_platform_id* platforms = (cl_platform_id* )malloc(numPlatforms* sizeof(cl_platform_id));
		error = clGetPlatformIDs(numPlatforms, platforms, NULL);
		platform = platforms[0];
		free(platforms);
	}

	error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
	if(error != CL_SUCCESS)
	{
		printf("[sba_opencl] Error getting device ids\n");
		exit(error);
	}

	context = clCreateContext(0, 1, &device, NULL, NULL, &error);
	if(error != CL_SUCCESS)
	{
		printf("[sba_opencl] Error creating context\n");
		exit(error);
	}
	
	queue = clCreateCommandQueue(context, device, 0, &error);
	if(error != CL_SUCCESS)
	{
		printf("[sba_opencl] Error creating command queue\n");
		exit(error);
	}
	
	error = clGetPlatformInfo(platform, CL_PLATFORM_NAME, NULL, platform_name, &size); // get size of platform name char array
	platform_name = (char*)malloc(size);
	error = clGetPlatformInfo(platform, CL_PLATFORM_NAME, size, platform_name, NULL); // get platform name char array

	//error = clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(char), platform_name, NULL);
	if(error != CL_SUCCESS)
	{
		printf("[sba_opencl] Error getting platform name\n");
		exit(error);
	}
	//error = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(char), platform_vendor, NULL);
	
	error = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, NULL, platform_vendor, &size); // get size of platform name char array
	platform_vendor = (char*)malloc(size);
	error = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, size, platform_vendor, NULL); // get platform name char array

	if(error != CL_SUCCESS)
	{
		printf("[sba_opencl] Error getting platform name\n");
		exit(error);
	}

	printf("[sba_opencl] Platform Name: %s, Platform Vendor: %s\n", platform_name, platform_vendor);
	/* Clean up */
	free(platform_name);
	free(platform_vendor);

	info.platform = platform;
	info.context = context;
	info.queue = queue;
	info.device = device;

	return info;
}
