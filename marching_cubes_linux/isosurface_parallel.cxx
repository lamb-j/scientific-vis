#include <iostream>
#include <vector>

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#include <CL/cl.hpp>
#define MAX_SOURCE_SIZE (0x100000)

int GetNumberOfCells(const int *dims);

float *GenerateIsosurfaceParallel(float *X, float *Y, float *Z, float *F,
    const int *dims, float iso_value, int *num_tris_total) {

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
  cl::Device default_device = all_devices[1];
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
  //int local_work_size = 511;
  int local_work_size = 73;

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
  if (status == CL_INVALID_WORK_GROUP_SIZE)
    printf("kernelPrefix failed: CL_INVALID_WORK_GROUP_SIZE\n");
  if (status != CL_SUCCESS) printf("kernelPrefix failed!\n");

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
  if (status == CL_INVALID_WORK_GROUP_SIZE)
    printf("kernelGenerateTri failed: CL_INVALID_WORK_GROUP_SIZE\n");
  if (status != CL_SUCCESS) printf("kernelGenerateTri failed!\n");

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
