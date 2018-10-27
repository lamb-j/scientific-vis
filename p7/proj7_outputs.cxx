/*=========================================================================

Program:   Visualization Toolkit
Module:    SpecularSpheres.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "TriangleList.h"
#include "tricase.cxx"

float InterpLin (float A, float B, float X, float fA, float fB) {
	return fA + ( (X - A) / (B - A) )*(fB - fA);
}

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
	return dims[0]*dims[1]*dims[2];
	// 2D
	//return dims[0]*dims[1];
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
	return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
	// 2D
	//return (dims[0]-1)*(dims[1]-1);
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
	return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
	// 2D
	//return idx[1]*dims[0]+idx[0];
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
	return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
	// 2D
	//return idx[1]*(dims[0]-1)+idx[0];
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
	idx[0] = pointId%dims[0];
	idx[1] = (pointId/dims[0])%dims[1];
	idx[2] = pointId/(dims[0]*dims[1]);

	// 2D
	//idx[0] = pointId%dims[0];
	//idx[1] = pointId/dims[0];
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
	idx[0] = cellId%(dims[0]-1);
	idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
	idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

	// 2D
	//idx[0] = cellId%(dims[0]-1);
	//idx[1] = cellId/(dims[0]-1);
}

void BoundingCubeForCell(const float *X, const float *Y, const float *Z, const int *dims,
		int cellId, float *bcube)
{
	if (cellId >= GetNumberOfCells(dims)) {
    printf("Bad cellID in BondingCubeForCell\n");
		return;
	}

  int idx[3];
  GetLogicalCellIndex(idx, cellId, dims);

	bcube[0] = X[idx[0]];
	bcube[1] = X[idx[0] + 1];

	bcube[2] = Y[idx[1]];
	bcube[3] = Y[idx[1] + 1];

	bcube[4] = Z[idx[2]];
	bcube[5] = Z[idx[2]+ 1];
}

int IdentifyCase(float *f_val, float iso_value) 
{
  int rv = 0;

  if (f_val[0] > iso_value) rv += 1;
  if (f_val[1] > iso_value) rv += 2;
  if (f_val[2] > iso_value) rv += 4;
  if (f_val[3] > iso_value) rv += 8;
  if (f_val[4] > iso_value) rv += 16;
  if (f_val[5] > iso_value) rv += 32;
  if (f_val[6] > iso_value) rv += 64;
  if (f_val[7] > iso_value) rv += 128;
	
  return rv;
};

void AddIsoSurfaces(float *X, float *Y, float *Z, float *F, const int *dims, TriangleList *tl, float iso_value) {

  static int edge_to_vertex[12][2]  = 
  { {0,1},  {1,3},  {2,3},  {0,2},
    {4,5},  {5,7},  {6,7},  {4,6},
    {0,4},  {1,5},  {2,6},  {3,7} };

	// iterate over cells
	for (int cell_id = 0; cell_id < GetNumberOfCells(dims); ++cell_id) {

		// get X, Y, Z values
		float bcube[6];   // X1, X2, Y1, Y2, Z1, Z2
		BoundingCubeForCell(X, Y, Z, dims, cell_id, bcube);

    // get Field values
    int idx[3];
    GetLogicalCellIndex(idx, cell_id, dims);

    float f_val[8];
    int vx = idx[0], vy = idx[1], vz = idx[2];
    int dz = dims[0]*dims[1], dy = dims[0];

    f_val[0] = F[vz*dz + vy*dy + vx];               // v0
    f_val[1] = F[vz*dz + vy*dy + (vx+1)];           // v1
    f_val[2] = F[vz*dz + (vy+1)*dy + vx];           // v2
    f_val[3] = F[vz*dz + (vy+1)*dy + (vx+1)];       // v3
    f_val[4] = F[(vz+1)*dz + vy*dy + vx];           // v4
    f_val[5] = F[(vz+1)*dz + vy*dy + (vx+1)];       // v5
    f_val[6] = F[(vz+1)*dz + (vy+1)*dy + vx];       // v6
    f_val[7] = F[(vz+1)*dz + (vy+1)*dy + (vx+1)];   // v7

		int icase = IdentifyCase(f_val, iso_value);

    // draw triangles 
    for (int i = 0; triCase[icase][i] != -1; i+=3) {
      printf("\n\n");
      //printf("triangle index: %d\n", i/3);
      float tri[3][3]; // (x1, y1, z1), (x2, y2, z2), (x3, y3, z3)

      for (int e = 0; e < 3; ++e) {
        int edge = triCase[icase][i + e];

        // get field values for this edge
        float f1 = f_val[edge_to_vertex[edge][0]];
        float f2 = f_val[edge_to_vertex[edge][1]];

        // determine X, Y, Z values of edge 
        float x1, x2, y1, y2, z1, z2;

        if (edge == 3 || edge == 7 || edge == 8 || edge == 10) 
        { x1 = x2 = bcube[0]; }
        if (edge == 0 || edge == 2 || edge == 4 || edge == 6) 
        { x1 = bcube[0]; x2 = bcube[1]; }
        if (edge == 1 || edge == 5 || edge == 9 || edge == 11) 
        { x1 = x2 = bcube[1]; }

        if (edge == 0 || edge == 8 || edge == 9 || edge == 4) 
        { y1 = y2 = bcube[2]; }
        if (edge == 1 || edge == 3 || edge == 5 || edge == 7) 
        { y1 = bcube[2]; y2 = bcube[3]; }
        if (edge == 2 || edge == 6 || edge == 10 || edge == 11) 
        { y1 = y2 = bcube[3]; }

        if (edge < 4)      { z1 = z2 = bcube[4]; }
        else if (edge > 7) { z1 = bcube[4]; z2 = bcube[5]; } 
        else               { z1 = z2 = bcube[5]; }

        // interpolate X, Y values of new point
        tri[e][0] = InterpLin(f1, f2, iso_value, x1, x2);
        tri[e][1] = InterpLin(f1, f2, iso_value, y1, y2);
        tri[e][2] = InterpLin(f1, f2, iso_value, z1, z2);
        //printf("edge: %d\n", edge);
      }

      tl->AddTriangle(tri[0][0], tri[0][1], tri[0][2], 
                      tri[1][0], tri[1][1], tri[1][2],
                      tri[2][0], tri[2][1], tri[2][2]); 
      //printf("Cell: %d, iso_value: %f\n", cell_id, iso_value);
      //printf("Case: %d\n", icase);
      //printf("x1: %f, x2: %f, y1: %f, y2: %f, z1: %f, z2: %f\n", 
      //    bcube[0], bcube[1], bcube[2], bcube[3], bcube[4], bcube[5]);
      //printf("f0: %f, f1: %f, f2: %f, f3: %f, f4: %f, f5: %f, f6: %f, f7: %f\n", 
      //    f_val[0], f_val[1], f_val[2], f_val[3],
      //    f_val[4], f_val[5], f_val[6], f_val[7]);
      //printf("(%f, %f, %f) -> (%f, %f, %f) -> (%f, %f, %f)\n", 
      //    tri[0][0], tri[0][1], tri[0][2], 
      //    tri[1][0], tri[1][1], tri[1][2], 
      //    tri[2][0], tri[2][1], tri[2][2] ); 
    }
	}

}

int main()
{
	int  i, j;

	vtkDataSetReader *rdr = vtkDataSetReader::New();
	rdr->SetFileName("proj7.vtk");
	rdr->Update();

	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
	rgrid->GetDimensions(dims);

	float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
	float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
	float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	TriangleList tl;

	float iso_value = 3.2;
	AddIsoSurfaces(X, Y, Z, F, dims, &tl, iso_value);

	vtkPolyData *pd = tl.MakePolyData();

	//This can be useful for debugging
	/*
		 vtkDataSetWriter *writer = vtkDataSetWriter::New();
		 writer->SetFileName("paths.vtk");
		 writer->SetInputData(pd);
		 writer->Write();
	 */

	vtkSmartPointer<vtkDataSetMapper> win1Mapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	win1Mapper->SetInputData(pd);
	win1Mapper->SetScalarRange(0, 0.15);

	vtkSmartPointer<vtkActor> win1Actor =
		vtkSmartPointer<vtkActor>::New();
	win1Actor->SetMapper(win1Mapper);

	vtkSmartPointer<vtkRenderer> ren1 =
		vtkSmartPointer<vtkRenderer>::New();

	vtkSmartPointer<vtkRenderWindow> renWin =
		vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren1);

	vtkSmartPointer<vtkRenderWindowInteractor> iren =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	ren1->AddActor(win1Actor);
	ren1->SetBackground(0.0, 0.0, 0.0);
	renWin->SetSize(800, 800);

	ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
	ren1->GetActiveCamera()->SetPosition(0,0,50);
	ren1->GetActiveCamera()->SetViewUp(0,1,0);
	ren1->GetActiveCamera()->SetClippingRange(20, 120);
	ren1->GetActiveCamera()->SetDistance(30);

	// This starts the event loop and invokes an initial render.
	//
	iren->Initialize();
	iren->Start();

	pd->Delete();
}
