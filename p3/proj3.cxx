#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
	// 3D
	//return dims[0]*dims[1]*dims[2];
	// 2D
	return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
	// 3D
	//return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
	// 2D
	return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
	// 3D
	//return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
	// 2D
	return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
	// 3D
	//return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
	// 2D
	return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
	// 3D
	// idx[0] = pointId%dim[0];
	// idx[1] = (pointId/dims[0])%dims[1];
	// idx[2] = pointId/(dims[0]*dims[1]);

	// 2D
	idx[0] = pointId%dims[0];
	idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
	// 3D
	// idx[0] = cellId%(dims[0]-1);
	// idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
	// idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

	// 2D
	idx[0] = cellId%(dims[0]-1);
	idx[1] = cellId/(dims[0]-1);
}

void BoundingBoxForCell(const float *X, const float *Y, const int *dims,
		int cellId, float *bbox)
{
	if (cellId >= GetNumberOfCells(dims)) {
		bbox[0] = -100;
		bbox[1] = 100;
		bbox[2] = -100;
		bbox[3] = 100;

		return;
	}

	int cell_col = cellId % (dims[0] - 1);
	int cell_row = cellId / (dims[0] - 1);

	bbox[0] = X[cell_col];
	bbox[1] = X[cell_col + 1];

	bbox[2] = Y[cell_row];
	bbox[3] = Y[cell_row + 1];
}

float InterpLin (float A, float B, float X, float fA, float fB) {
	return fA + ( (X - A) / (B - A) )*(fB - fA);
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
		const float *Y, const float *F)
{
	if (pt[0] < X[0] || pt[0] > X[dims[0] - 1] ) return 0;
	if (pt[1] < Y[0] || pt[1] > Y[dims[1] - 1] ) return 0;

	// get pointId and cellId
	int idx[2];
	int i;
	for (i = 0; X[i] < pt[0]; ++i)
		; idx[0] = i - 1;
	for (i = 0; Y[i] < pt[1]; ++i)
		; idx[1] = i - 1;

	int cellId = GetCellIndex(idx, dims);
	int pointId = GetPointIndex(idx, dims);
	// get bounding box
	float bbox[4];
	BoundingBoxForCell(X, Y, dims, cellId, bbox);

	// get field values 
	float fval[4];
	fval[0] = F[pointId];
	fval[1] = F[pointId + 1];
	fval[2] = F[pointId + dims[0]];
	fval[3] = F[pointId + dims[0] + 1];

	// get top and bottom field values
	float fBot = InterpLin(bbox[0], bbox[1], pt[0], fval[0], fval[1]);
	float fTop = InterpLin(bbox[0], bbox[1], pt[0], fval[2], fval[3]);

	// get final value
	return InterpLin(bbox[2], bbox[3], pt[1], fBot, fTop);
}


void WriteImage(vtkImageData *img, const char *filename)
{
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

	vtkImageData *
NewImage(int width, int height)
{
	vtkImageData *image = vtkImageData::New();
	image->SetDimensions(width, height, 1);
	//image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
	//image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
	//image->SetNumberOfScalarComponents(3);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	//image->AllocateScalars();

	return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
	float fR = InterpLin(0, 1, F, 0, 255);
	float fG = InterpLin(0, 1, F, 0, 255);
	float fB = InterpLin(0, 1, F, 128, 255);

	RGB[0] = floor(fR);	
	RGB[1] = floor(fG);	
	RGB[2] = floor(fB);	
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
	void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
	float fR, fG, fB;
	if (F <= 0.5) {
		fR = InterpLin(0, 0.5, F, 0, 255);
		fG = InterpLin(0, 0.5, F, 0, 255);
		fB = InterpLin(0, 0.5, F, 128, 255); 
	}

	if (F > 0.5) {
		fR = InterpLin(0.5, 1, F, 255, 128);
		fG = InterpLin(0.5, 1, F, 255, 0);
		fB = InterpLin(0.5, 1, F, 255, 0);
	}

	RGB[0] = floor(fR);	
	RGB[1] = floor(fG);	
	RGB[2] = floor(fB);	
}

// ****************************************************************************
//  Function: ApplyHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
	void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
	float saturation = 1.0;
	float v = 1.0;
	float hue = InterpLin(0, 1, F, 0, 360);

	hue /= 60.f;
	// sector 0 to 5
	int i = floor (hue);
	float f = hue - i;

	// factorial part of h
	float p = 1 - saturation;
	float q = 1 - saturation * f;
	float t = 1 - saturation * ( 1 - f );

	float fR, fG, fB;

	switch (i) {
		case 0: fR = v; fG = t; fB = p; break;
		case 1: fR = q; fG = v; fB = p; break;
		case 2: fR = p; fG = v; fB = t; break;
		case 3: fR = p; fG = q; fB = v; break;
		case 4: fR = t; fG = p; fB = v; break;
		case 5: fR = v; fG = p; fB = q; break;
	}

	fR = InterpLin(0, 1, fR, 0, 255);
	fG = InterpLin(0, 1, fG, 0, 255);
	fB = InterpLin(0, 1, fB, 0, 255);

	RGB[0] = floor(fR);	
	RGB[1] = floor(fG);	
	RGB[2] = floor(fB);	
}

int main()
{
	int  i, j;

	vtkDataSetReader *rdr = vtkDataSetReader::New();
	rdr->SetFileName("proj3_data.vtk");
	rdr->Update();

	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
	rgrid->GetDimensions(dims);

	float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
	float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	int nx = 500;
	int ny = 500;

	vtkImageData *images[3];
	unsigned char *buffer[3];
	for (i = 0 ; i < 3 ; i++)
	{
		images[i] = NewImage(nx, ny);
		buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
	}

	for (i = 0 ; i < 3*nx*ny ; i++)
		for (j = 0 ; j < 3 ; j++)
			buffer[j][i] = 0;

	for (i = 0 ; i < nx ; i++)
		for (j = 0 ; j < ny ; j++)
		{
			float pt[2];
			pt[0] = InterpLin(0, 499, i, -9, 9);
			pt[1] = InterpLin(0, 499, j, -9, 9);
			float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
			float normalizedF = InterpLin(1.2, 5.02, f, 0, 1); 

			int offset = 3*(j*nx+i);
			ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
			ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
			ApplyHSVColorMap(normalizedF, buffer[2]+offset);
		}

	WriteImage(images[0], "bluehot");
	WriteImage(images[1], "difference");
	WriteImage(images[2], "hsv");
}
