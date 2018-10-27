#include <iostream>
#include <vector>

#include <vtkDataSetReader.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define TIME 0

int dims[3];
float *X;
float *Y;
float *Z;
float *F;

inline float InterpLin (float A, float B, float X, float fA, float fB) {
  return fA + ( (X - A) / (B - A) )*(fB - fA);
}

struct Camera
{
  double          near, far;
  double          angle;
  double          position[3];
  double          focus[3];
  double          up[3];
};


struct TransferFunction
{
  double          min;
  double          max;
  int             numBins;
  unsigned char  *colors;  // size is 3*numBins
  double         *opacities; // size is numBins

  void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
  {

    int bin = floor( InterpLin(min, max, value, 0, numBins) ); 
    if (bin < 0 || bin > 255) { 
      RGB[0] = 0; 
      RGB[1] = 0;
      RGB[2] = 0;
      opacity = 0;
      return;
    }

    RGB[0] = colors[3*bin+0];
    RGB[1] = colors[3*bin+1];
    RGB[2] = colors[3*bin+2];
    opacity = opacities[bin];
  }
};

TransferFunction SetupTransferFunction(void)
{
  int  i;

  TransferFunction rv;
  rv.min = 10;
  rv.max = 15;
  rv.numBins = 256;
  rv.colors = new unsigned char[3*256];
  rv.opacities = new double[256];
  unsigned char charOpacity[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  for (i = 0 ; i < 256 ; i++)
    rv.opacities[i] = charOpacity[i]/255.0;

  const int numControlPoints = 8;
  unsigned char controlPointColors[numControlPoints*3] = { 
    71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
    255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
  };

  double controlPointPositions[numControlPoints] = 
  { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };

  for (i = 0 ; i < numControlPoints-1 ; i++)
  {
    int start = controlPointPositions[i]*rv.numBins;
    int end   = controlPointPositions[i+1]*rv.numBins+1;
    //std::cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << std::endl;
    if (end >= rv.numBins)
      end = rv.numBins-1;
    for (int j = start ; j <= end ; j++)
    {
      double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
      if (proportion < 0 || proportion > 1.)
        continue;
      for (int k = 0 ; k < 3 ; k++)
        rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
          + controlPointColors[3*i+k];
    }
  }    

  return rv;
}

Camera SetupCamera(void)
{
  Camera rv;
  rv.focus[0] = 0;
  rv.focus[1] = 0;
  rv.focus[2] = 0;
  rv.up[0] = 0;
  rv.up[1] = -1;
  rv.up[2] = 0;
  rv.angle = 30;
  rv.near = 7.5e+7;
  rv.far = 1.4e+8;
  rv.position[0] = -8.25e+7;
  rv.position[1] = -3.45e+7;
  rv.position[2] = 3.35e+7;

  return rv;
}

void CrossProduct(double a[3], double b[3], double rv[3]) {
	rv[0] = a[1]*b[2] - a[2]*b[1];
	rv[1] = a[2]*b[0] - a[0]*b[2];
	rv[2] = a[0]*b[1] - a[1]*b[0];
}

double Magnitude(double v[3]) {
	return sqrt( (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) );
}

void FindRay(Camera C, int W, int H, int x, int y, double ray[3]) {

	double look[3];
  look[0] = C.focus[0] - C.position[0];
  look[1] = C.focus[1] - C.position[1];
  look[2] = C.focus[2] - C.position[2];

  double ru[3];
	CrossProduct(look, C.up, ru);

	double mag_ru = Magnitude(ru);
	for (int i = 0; i < 3; ++i) ru[i] /= mag_ru;

	double rv[3];
	CrossProduct(look, ru, rv);

	double mag_rv = Magnitude(rv);
	for (int i = 0; i < 3; ++i) rv[i] /= mag_rv;

	double fovx = M_PI / 6.0;
	double fovy = M_PI / 6.0;

	double rdx_tan = (2 * tan(fovx / 2.0)) / W;
	double rdy_tan = (2 * tan(fovy / 2.0)) / H;

	double rdx[3];
	for (int i = 0; i < 3; ++i) rdx[i] = rdx_tan * ru[i]; 
	double rdy[3];
	for (int i = 0; i < 3; ++i) rdy[i] = rdy_tan * rv[i]; 

	double look_norm[3];
	double look_mag = Magnitude(look);
	for (int i = 0; i < 3; ++i) look_norm[i] = look[i] / look_mag;

	for (int i = 0; i < 3; ++i) 
		ray[i] = look_norm[i] + 0.5*(2*y + 1 - W)*rdx[i] + 0.5*(2*x + 1 - H)*rdy[i];
}

// Evaluate vector field for a given sample (x,y,z)
float VectorField(double sample[3], int *last_vx, int *last_vy, int *last_vz) {

  //if (vx == dims[0] || vy == dims[1] || vz == dims[2]) return 0;
  //if (vx == 0 || vy == 0 || vz == 0) return 0;

  if (sample[0] < X[0] || sample[1] < Y[0] || sample[2] < Z[0] ) return 0;
  if (sample[0] > X[dims[0] - 1] || 
      sample[1] > Y[dims[1] - 1] || 
      sample[2] > Z[dims[2] - 1] ) return 0;

	// get X, Y, Z index
	int vx, vy, vz;
  // vx
	for (vx = *last_vx; X[vx] < sample[0] && vx < dims[0]; vx++);

  // vy
  if (sample[1] > Y[*last_vy]) 
    for (vy = *last_vy; Y[vy] < sample[1] && vy < dims[1]; vy++);
  else {
    for (vy = *last_vy - 1; Y[vy] > sample[1] && vy > 0; vy--);
    vy++;
  }

  // vz
  if (sample[2] > Z[*last_vz]) 
    for (vz = *last_vz; Z[vz] < sample[2] && vz < dims[2]; vz++);
  else {
    for (vz = *last_vz; Z[vz] > sample[2] && vz > 0; vz--);
    vz++;
  }

  //if (vx == dims[0] || vy == dims[1] || vz == dims[2]) return 0;
  //if (vx == 0 || vy == 0 || vz == 0) return 0;

  vx--; vy--; vz--;

	// get field values for cell
	int dz = (dims[0])*(dims[1]), dy = (dims[0]);
	float f[8];
	f[0] = F[vz*dz + vy*dy + vx];               // v0
	f[1] = F[vz*dz + vy*dy + (vx+1)];           // v1
	f[2] = F[vz*dz + (vy+1)*dy + vx];           // v2
	f[3] = F[vz*dz + (vy+1)*dy + (vx+1)];       // v3
	f[4] = F[(vz+1)*dz + vy*dy + vx];           // v4
	f[5] = F[(vz+1)*dz + vy*dy + (vx+1)];       // v5
	f[6] = F[(vz+1)*dz + (vy+1)*dy + vx];       // v6
	f[7] = F[(vz+1)*dz + (vy+1)*dy + (vx+1)];   // v7

	// interploate inner field value 
	float f01x = InterpLin(X[vx], X[vx+1], sample[0], f[0], f[1]);
	float f23x = InterpLin(X[vx], X[vx+1], sample[0], f[2], f[3]);
	float f45x = InterpLin(X[vx], X[vx+1], sample[0], f[4], f[5]);
	float f67x = InterpLin(X[vx], X[vx+1], sample[0], f[6], f[7]);

	float f0123y = InterpLin(Y[vy], Y[vy+1], sample[1], f01x, f23x); 
	float f4567y = InterpLin(Y[vy], Y[vy+1], sample[1], f45x, f67x); 

	float fz = InterpLin(Z[vz], Z[vz+1], sample[2], f0123y, f4567y);

  *last_vx = vx;
  *last_vy = vy;
  *last_vz = vz;
	return fz;
}

void IntersectVolume_TransferSamples(Camera C, 
                                     double ray[3], 
                                     int n_samples, 
                                     int step_size, 
                                     unsigned char final_RGB[3], 
                                     TransferFunction *tf) {

	// Scale ray to near plane
	double ray_near[3];
	for (int i = 0; i < 3; ++i) ray_near[i] = ray[i] * C.near;

	// Add ray to origin
	for (int i = 0; i < 3; ++i) ray_near[i] += C.position[i];

	double step[3];
	for (int i = 0; i < 3; ++i) { step[i] = ray[i] * step_size; }

  unsigned char sample_rgb[3];
  double sample_opacity = 0;
  double RGB[3];
  double opacity;

  // keep track of last cell
  int last_vx, last_vy, last_vz;
  last_vx = last_vy = last_vz = 0;

	// Get initial sample
	double sample[3];
	sample[0] = ray_near[0];
	sample[1] = ray_near[1];
	sample[2] = ray_near[2];
	float field_value = VectorField(sample, &last_vx, &last_vy, &last_vz);
  tf->ApplyTransferFunction(field_value, sample_rgb, sample_opacity);

  RGB[0] = sample_rgb[0];
  RGB[1] = sample_rgb[1];
  RGB[2] = sample_rgb[2];
  opacity = sample_opacity;


  // Accumulate remaining samples
	for (int i = 1; i < n_samples; ++i) {
		float field_value = VectorField(sample, &last_vx, &last_vy, &last_vz);

    if (field_value == 0) {
      // take step
      sample[0] += step[0];
      sample[1] += step[1];
      sample[2] += step[2];
      continue;
    }

    // Transfer this sample
    tf->ApplyTransferFunction(field_value, sample_rgb, sample_opacity);

    double s_rgb[3];
    s_rgb[0] = InterpLin(0, 255, sample_rgb[0], 0, 1);
    s_rgb[1] = InterpLin(0, 255, sample_rgb[1], 0, 1);
    s_rgb[2] = InterpLin(0, 255, sample_rgb[2], 0, 1);

    // Opacity correction
    sample_opacity = 1 - pow(1.0 - sample_opacity, 500.0 / (double) n_samples);

    // Composite
    RGB[0] += (1-opacity)*sample_opacity*s_rgb[0];
    RGB[1] += (1-opacity)*sample_opacity*s_rgb[1];
    RGB[2] += (1-opacity)*sample_opacity*s_rgb[2];
    opacity = opacity + (1-opacity)*sample_opacity;

    if (opacity > 1.0 - 1e-3) break;

		// take step
		sample[0] += step[0];
		sample[1] += step[1];
		sample[2] += step[2];
	}

  final_RGB[0] = floor( InterpLin(0, 1, RGB[0], 0, 255) );
  final_RGB[1] = floor( InterpLin(0, 1, RGB[1], 0, 255) );
  final_RGB[2] = floor( InterpLin(0, 1, RGB[2], 0, 255) );
	
	return;
}	

int main()
{
	// Read dataset
  vtkDataSetReader *rdr = vtkDataSetReader::New();
  rdr->SetFileName("astro512.vtk");
  rdr->Update();

  vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
  rgrid->GetDimensions(dims);

	printf("dims: %d, %d, %d\n", dims[0], dims[1], dims[2]);

  X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
  Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
  Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
  F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

  // Create pixel buffer
	int N = 1024;
	vtkSmartPointer<vtkImageData> pixelImage = vtkSmartPointer<vtkImageData>::New();
  pixelImage->SetDimensions(N, N, 1);
	pixelImage->AllocateScalars(VTK_UNSIGNED_CHAR,3);

	// Create camera
  struct Camera Cam;
	Cam = SetupCamera();

  // Calculate step size
  int n_samples = 1024;
	int step_size = ceil((Cam.far - Cam.near) / (n_samples - 1)); 

	// Create transfer function
  TransferFunction tf = SetupTransferFunction();

#if TIME
  long long rayT = 0;
  long long inter_xferT = 0;
#endif

  // For every pixel
  auto t1 = Clock::now();
	for (int x = 0; x < N; ++x) {
		for (int y = 0; y < N; ++y) {
			unsigned char* pixel = 
				static_cast<unsigned char*>(pixelImage->GetScalarPointer(x,y,0));

			// Find ray for that pixel
			double ray[3];
#if TIME
      auto ray1 = Clock::now();
#endif

			FindRay(Cam, N, N, y, x, ray);

#if TIME
      auto ray2 = Clock::now();
      rayT += 
        std::chrono::duration_cast<std::chrono::microseconds>(ray2 - ray1).count();
#endif

			// Intersect volume with ray and Calculate color from xfer function
#if TIME
      auto inter_xfer1 = Clock::now();
#endif
      unsigned char RGB[3];
 			IntersectVolume_TransferSamples(Cam, ray, n_samples, step_size, RGB, &tf);

#if TIME
      auto inter_xfer2 = Clock::now();
      inter_xferT += 
        std::chrono::duration_cast<std::chrono::microseconds>(inter_xfer2 - inter_xfer1).count();
#endif
			// Assign color to pixel
			pixel[0] = RGB[0];
			pixel[1] = RGB[1];
			pixel[2] = RGB[2];
		}
	}
  auto t2 = Clock::now();

#if TIME
  printf("Ray time: %lld ms\n", rayT);
  printf("Intersect Volume and Transfer time: %lld ms\n", inter_xferT);
#endif
  printf("Total time: %f seconds\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000.0);


  // Render pixel data
	pixelImage->Modified();

	vtkSmartPointer<vtkImageMapper> pixelMapper = vtkSmartPointer<vtkImageMapper>::New();
	pixelMapper->SetInputData(pixelImage);
  pixelMapper->SetColorWindow(255);
  pixelMapper->SetColorLevel(127.5);

	vtkSmartPointer<vtkActor2D> pixelActor = vtkSmartPointer<vtkActor2D>::New();
  pixelActor->SetMapper(pixelMapper);

  vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	ren->AddActor2D(pixelActor);

  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren); 
  renWin->SetSize(N, N);

  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);
  iren->Initialize();
  iren->Start();

	return 0;
}
