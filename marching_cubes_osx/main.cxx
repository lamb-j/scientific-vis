#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>

#include <sstream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define SEQ 1

float *GenerateIsosurfacePrefix(float *X, float *Y, float *Z, float *F, const int *dims, 
      float iso_value, int *num_tris_total);

float *GenerateIsosurfaceParallel(float *X, float *Y, float *Z, float *F, 
    const int *dims, float iso_value, int *num_tris_total, int device_num);

void RenderSurface(float *tl, int num_tris);

int main(int argc, char **argv)
{
  int  i, j;

  if (argc != 3) {
    printf("usage: astro data.vtk device_num\n");
    exit(0);
  }

  int device_num;
  std::istringstream ss1(argv[2]);
  if (!(ss1 >> device_num))
    cerr << "Invalid device number" << argv[2] << '\n';

  std::string file_name;
  std::istringstream ss2(argv[1]);
  if (!(ss2 >> file_name))
    cerr << "Invalid file name" << argv[1] << '\n';

  printf("file: %s\n", file_name.c_str());

  vtkDataSetReader *rdr = vtkDataSetReader::New();
  rdr->SetFileName(file_name.c_str());
  rdr->Update();

  int dims[3];
  vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
  rgrid->GetDimensions(dims);

  printf("dims[0]: %d, dims[1]: %d, dims[2]: %d, dims[all]: %d\n", 
      dims[0], dims[1], dims[2], dims[0]*dims[1]*dims[2]);

  float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
  float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
  float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
  float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

  float iso_value = 3.2;

#if SEQ
  auto s1 = Clock::now();
  int num_tris_seq;
  float *tl_seq = GenerateIsosurfacePrefix(X, Y, Z, F, dims, iso_value, &num_tris_seq);
  auto s2 = Clock::now();

  std::cout << "Maching Cubes Sequential\n  Triangles:\t" << num_tris_seq 
    <<  "\n  Elapsed time:\t"
    << std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1).count()
    << " ms\n\n\n";

  RenderSurface(tl_seq, num_tris_seq);
#endif

  auto p1 = Clock::now();
  int num_tris_par;
  float *tl_par = GenerateIsosurfaceParallel(X, Y, Z, F, dims, iso_value, &num_tris_par,
      device_num);
  auto p2 = Clock::now();

  std::cout << "Maching Cubes Parallel\n  Triangles:\t" << num_tris_par 
    <<  "\n  Elapsed time:\t"
    << std::chrono::duration_cast<std::chrono::milliseconds>(p2 - p1).count()
    << " ms\n";

  RenderSurface(tl_par, num_tris_par);
}
