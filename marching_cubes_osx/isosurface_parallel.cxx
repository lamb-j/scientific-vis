#include <iostream>
#include <vector>

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#include <OpenCL/cl.hpp>
#define MAX_SOURCE_SIZE (0x100000)

int GetNumberOfCells(const int *dims);
const char *getErrorString(cl_int error);

float *GenerateIsosurfaceParallel(float *X, float *Y, float *Z, float *F,
    const int *dims, float iso_value, int *num_tris_total, int device_num) {

  int Fsize = dims[0]*dims[1]*dims[2];
  int nCells = GetNumberOfCells(dims);
  int status;

  // Build OpenCL program
  std::vector<cl::Platform> all_platforms;
  cl::Platform::get(&all_platforms);
  if(all_platforms.size() == 0){
    std::cout << "No platforms found. Check OpenCL installation!" << std::endl;
    exit(1);
  }
  cl::Platform default_platform = all_platforms[0];
  std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

  std::vector<cl::Device> all_devices;
  default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
  if(all_devices.size() == 0){
    std::cout << "No devices found. Check OpenCL installation!" << std::endl;
    exit(1);
  }
  cl::Device default_device = all_devices[device_num];
  std::cout<< "Using device: "<< default_device.getInfo<CL_DEVICE_NAME>() << std::endl;

  cl::Context context(default_device);

  FILE *fp;
  char fileName[] = "./isosurface.cl";
  char *source_str;
  size_t source_size;

  fp = fopen(fileName, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel source.");
    exit(1);
  }
  source_str = (char*) malloc (MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  cl::Program::Sources sources;
  sources.push_back({source_str, source_size});

  cl::Program program(context, sources);

  if(program.build({default_device}) != CL_SUCCESS) {
    std::cout<< "Error building: "<< program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
    exit(1);
  }

  cl::CommandQueue queue(context, default_device);

  // Buffers 
  cl::Buffer bufferX(context, CL_MEM_READ_WRITE, dims[0] * sizeof(float));
  cl::Buffer bufferY(context, CL_MEM_READ_WRITE, dims[1] * sizeof(float));
  cl::Buffer bufferZ(context, CL_MEM_READ_WRITE, dims[2] * sizeof(float));
  cl::Buffer bufferF(context, CL_MEM_READ_WRITE, Fsize * sizeof(float));
  cl::Buffer bufferD(context, CL_MEM_READ_WRITE, 3 * sizeof(int));

  int *tri_occupied = (int *) malloc(sizeof(int) * nCells);
  cl::Buffer bufferOCC(context, CL_MEM_READ_WRITE, nCells * sizeof(int));

  queue.enqueueWriteBuffer(bufferX, CL_TRUE, 0, dims[0]*sizeof(float), X);
  queue.enqueueWriteBuffer(bufferY, CL_TRUE, 0, dims[1]*sizeof(float), Y);
  queue.enqueueWriteBuffer(bufferZ, CL_TRUE, 0, dims[2]*sizeof(float), Z);
  queue.enqueueWriteBuffer(bufferF, CL_TRUE, 0, Fsize*sizeof(float), F);
  queue.enqueueWriteBuffer(bufferD, CL_TRUE, 0, 3*sizeof(int), dims);

  // Create Prefix Kernel 
  cl::Kernel kernelPrefix(program, "CalculatePrefix");
  kernelPrefix.setArg(0, bufferX);
  kernelPrefix.setArg(1, bufferY);
  kernelPrefix.setArg(2, bufferZ);
  kernelPrefix.setArg(3, bufferF);
  kernelPrefix.setArg(4, bufferD);
  kernelPrefix.setArg(5, bufferOCC);
	kernelPrefix.setArg(6, iso_value);

  int global_work_size = nCells; 
  int local_work_size = dims[2] - 1; // data.vtk
  //int local_work_size = 73; // astro512.vtk

  auto k1_start = Clock::now(); 
  status = queue.enqueueNDRangeKernel(
      kernelPrefix, 
      cl::NullRange,    
      cl::NDRange(global_work_size),  
      cl::NDRange(local_work_size)); 
  queue.finish();
  auto k1_end = Clock::now(); 

  if (status == CL_SUCCESS) 
    printf("kernelPrefix success!\n");
  if (status != CL_SUCCESS) {
    printf("kernelPrefix failed: %s\n",  getErrorString(status));
    exit(0);
  }

  queue.enqueueReadBuffer(bufferOCC, CL_TRUE, 0, nCells * sizeof(int), tri_occupied);

  auto pre_start = Clock::now();
  // Create Prefix array
  int *tri_prefix = (int *) malloc(sizeof(int) * nCells);
  tri_prefix[0] = 0;
  for (int i = 1; i < nCells; ++i) {
    tri_prefix[i] = tri_prefix[i-1] + tri_occupied[i-1];
  }

  *num_tris_total = tri_occupied[nCells - 1] + tri_prefix[nCells - 1];
  auto pre_end = Clock::now();

  float *tl = (float *) malloc(sizeof(float) * (*num_tris_total) * 9);

  // Trianlge and Prefix Buffer
  cl::Buffer bufferPRE(context, CL_MEM_READ_WRITE, nCells * sizeof(int));
  cl::Buffer bufferTRI(context, CL_MEM_READ_WRITE, sizeof(float) * (*num_tris_total) * 9);

  // Write Prefix to device 
  queue.enqueueWriteBuffer(bufferPRE, CL_TRUE, 0, nCells * sizeof(int), tri_prefix);
  queue.finish();

  // Create GenerateTris Kernel 
  cl::Kernel kernelGenerateTris(program, "GenerateTris");
  kernelGenerateTris.setArg(0, bufferX);
  kernelGenerateTris.setArg(1, bufferY);
  kernelGenerateTris.setArg(2, bufferZ);
  kernelGenerateTris.setArg(3, bufferF);
  kernelGenerateTris.setArg(4, bufferD);
  kernelGenerateTris.setArg(5, bufferOCC);
  kernelGenerateTris.setArg(6, bufferPRE);
  kernelGenerateTris.setArg(7, bufferTRI);
	kernelGenerateTris.setArg(8, iso_value);

  auto k2_start = Clock::now();
  status = queue.enqueueNDRangeKernel(
      kernelGenerateTris, 
      cl::NullRange,    
      cl::NDRange(global_work_size),  
      cl::NDRange(local_work_size)); 
  queue.finish();
  auto k2_end = Clock::now();

  if (status == CL_SUCCESS) 
    printf("kernelGenerateTri success!\n");
  if (status != CL_SUCCESS) {
    printf("kernelGenerateTri failed: %s\n",  getErrorString(status));
    exit(0);
  }

  queue.enqueueReadBuffer(bufferTRI, CL_TRUE, 0,
      sizeof(float) * (*num_tris_total) * 9, tl);

  /*
  printf("Occupy Kernel: %lld microseconds\n", 
    std::chrono::duration_cast<std::chrono::microseconds>(k1_end - k1_start).count());
  printf("Prefix Sum: %lld microseconds\n", 
    std::chrono::duration_cast<std::chrono::microseconds>(pre_end - pre_start).count());
  printf("GenTri Kernel: %lld microseconds\n", 
    std::chrono::duration_cast<std::chrono::microseconds>(k2_end - k2_start).count());
  */
  printf("Occupy Kernel: %lld milliseconds\n", 
    std::chrono::duration_cast<std::chrono::milliseconds>(k1_end - k1_start).count());
  printf("Prefix Sum: %lld milliseconds\n", 
    std::chrono::duration_cast<std::chrono::milliseconds>(pre_end - pre_start).count());
  printf("GenTri Kernel: %lld milliseconds\n", 
    std::chrono::duration_cast<std::chrono::milliseconds>(k2_end - k2_start).count());
  return tl;
}
const char *getErrorString(cl_int error)
{
switch(error){
    // run-time and JIT compiler errors
    case 0: return "CL_SUCCESS";
    case -1: return "CL_DEVICE_NOT_FOUND";
    case -2: return "CL_DEVICE_NOT_AVAILABLE";
    case -3: return "CL_COMPILER_NOT_AVAILABLE";
    case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5: return "CL_OUT_OF_RESOURCES";
    case -6: return "CL_OUT_OF_HOST_MEMORY";
    case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8: return "CL_MEM_COPY_OVERLAP";
    case -9: return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";
    case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15: return "CL_COMPILE_PROGRAM_FAILURE";
    case -16: return "CL_LINKER_NOT_AVAILABLE";
    case -17: return "CL_LINK_PROGRAM_FAILURE";
    case -18: return "CL_DEVICE_PARTITION_FAILED";
    case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64: return "CL_INVALID_PROPERTY";
    case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66: return "CL_INVALID_COMPILER_OPTIONS";
    case -67: return "CL_INVALID_LINKER_OPTIONS";
    case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
    }
}
