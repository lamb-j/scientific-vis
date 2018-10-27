Build:
  cmake .
  ./isosurface file.vtk num_device

  num_device: The index of the device you want to target as the accelerator. 
    On my mac '0' points to the CPU, while '1' points to the GPU.

  NOTE: Some devices currently may fail to execute the 'GenerateTris' kernel due 
  to memory constraints. For example, this happens on my Macbook Iris graphics card
  for the larger astro512.vtk data file.

  NOTE: Currently for some devices the local_work_size (isosurface_paralle.cxx:94) 
  must be manually changed for the astro512.vtk data set, becuase the work size 
  calculated from the input dimensions is too large. (CL_INVALID_WORK_GROUP_SIZE)

  NOTE: The sequential version can be disabled by setting SEQ to 0 (main.cxx:10).

Description:
  This application implements the following:
  
    1. A serial C++ isosurface algorithm (isosurface_prefix.cxx)
    2. A parallel OpenCL isosurface algorithm (isosurface_parallel.cxx)
  
  Both the serial and parallel implementation consists of 3 steps:
  
    1. Determine how many triangles each cell contains
  
    2. Create an appropriatly sized array to store the triangles, and use a 
      prefix-sum operation to create an index array into the triangles array
      for each cell.
  
    3. Generate the triangles for each cell and store them using the prefix
      array index.

  We create OpenCL kernels for steps 1 and 3. It is also possible to parallelize the
  prefix_sum calculation, but we have not yet implemented this.

Below are some runtime results:
——- cerberus.nic.uoregon.edu (device 1)
./isosurface astro.vtk 1
dims[0]: 512, dims[1]: 512, dims[2]: 512, dims[all]: 134217728
Maching Cubes Sequential
  Triangles:	2462440
  Elapsed time:	14658 ms


Using platform: NVIDIA CUDA
Using device: Tesla K40c
kernelPrefix success!
kernelGenerateTri success!
Occupy Kernel: 38 milliseconds
Prefix Sum: 943 milliseconds
GenTri Kernel: 38 milliseconds
Maching Cubes Parallel
  Triangles:	2462440
  Elapsed time:	2966 ms

——— MacBook (device 0)
astro512.vtk
dims[0]: 512, dims[1]: 512, dims[2]: 512, dims[all]: 134217728
Maching Cubes Sequential
  Triangles:	2462440
  Elapsed time:	8027 ms

Using platform: Apple
Using device: Intel(R) Core(TM) i5-5257U CPU @ 2.70GHz
kernelPrefix success!
kernelGenerateTri success!
Occupy Kernel: 1702 milliseconds
Prefix Sum: 719 milliseconds
GenTri Kernel: 2046 milliseconds
Maching Cubes Parallel
  Triangles:	2462440
  Elapsed time:	5520 ms
