#include "stdio.h"
#include "string.h"
#include "sba_opencl.h"
#include "assert.h"

char * extract_kernel(char *kernel_file)
{
	char * buffer = 0;
	long length;
	FILE * f = fopen (kernel_file, "r");

	if (f)
	{
  		fseek (f, 0, SEEK_END);
  		length = ftell (f);
  		fseek (f, 0, SEEK_SET);
  		buffer = malloc (length+1);
  		if (buffer)
  		{
    			fread (buffer, 1, length, f);
  		}
  		fclose (f);
		buffer[length] = '\0';
	}
	else
		printf("[sba_opencl] Invalid filename\n");
	return buffer;
}

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

void sba_create_kernel(struct opencl_info *info, char *kernel_file, char *kernel_name, int kernel_index)
{
	size_t src_size = 0;
	cl_int error;

	// Extract the kernel as a string
	char *kernel_string = extract_kernel(kernel_file);
	src_size = strlen(kernel_string);

	// Create opencl program
	cl_program program = clCreateProgramWithSource(info->context, 1, &kernel_string, &src_size, &error);

	assert(error == CL_SUCCESS);

	// Build the opencl program
	error = clBuildProgram(program, 1, &(info->device), NULL, NULL, NULL);

	// Get build log
	char *build_log;
	size_t log_size;

	clGetProgramBuildInfo(program, info->device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	build_log = (char *)malloc((log_size+1)*sizeof(char));

	clGetProgramBuildInfo(program, info->device, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
	build_log[log_size] = '\0';

	printf("[sba_opencl] Build log: %s\n", build_log);

	free(build_log);
	free(kernel_string);

	cl_kernel sba_kernel = clCreateKernel(program, kernel_name, &error);
	if(error != CL_SUCCESS)
		printf("[sba_opencl] Error is: %d\n", error);
	
	assert(error == CL_SUCCESS);

	info->kernel[kernel_index] = sba_kernel;
}
