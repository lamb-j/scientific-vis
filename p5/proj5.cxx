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


class SegmentList
{
	public:
		SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
		virtual      ~SegmentList() { delete [] pts; };

		void          AddSegment(float X1, float Y1, float X2, float Y2);
		vtkPolyData  *MakePolyData(void);

	protected:
		float        *pts;
		int           maxSegments;
		int           segmentIdx;
};

void SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
	pts[4*segmentIdx+0] = X1;
	pts[4*segmentIdx+1] = Y1;
	pts[4*segmentIdx+2] = X2;
	pts[4*segmentIdx+3] = Y2;
	segmentIdx++;
}

vtkPolyData* SegmentList::MakePolyData(void)
{
	int nsegments = segmentIdx;
	int numPoints = 2*(nsegments);
	vtkPoints *vtk_pts = vtkPoints::New();
	vtk_pts->SetNumberOfPoints(numPoints);
	int ptIdx = 0;
	vtkCellArray *lines = vtkCellArray::New();
	lines->EstimateSize(numPoints,2);
	for (int i = 0 ; i < nsegments ; i++)
	{
		double pt[3];
		pt[0] = pts[4*i];
		pt[1] = pts[4*i+1];
		pt[2] = 0.;
		vtk_pts->SetPoint(ptIdx, pt);
		pt[0] = pts[4*i+2];
		pt[1] = pts[4*i+3];
		pt[2] = 0.;
		vtk_pts->SetPoint(ptIdx+1, pt);
		vtkIdType ids[2] = { ptIdx, ptIdx+1 };
		lines->InsertNextCell(2, ids);
		ptIdx += 2;
	}

	vtkPolyData *pd = vtkPolyData::New();
	pd->SetPoints(vtk_pts);
	pd->SetLines(lines);
	lines->Delete();
	vtk_pts->Delete();

	return pd;
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

void LookUpTable(int lup[][4]) {

	lup[0][0] = lup[0][1] = lup[0][2] = lup[0][3] = -1;

	lup[1][0] = 0; lup[1][1] = 3; lup[1][2] = lup[1][3] = -1; 
	lup[2][0] = 0; lup[2][1] = 1; lup[2][2] = lup[2][3] = -1; 
	lup[3][0] = 1; lup[3][1] = 3; lup[3][2] = lup[3][3] = -1; 
	lup[4][0] = 2; lup[4][1] = 3; lup[4][2] = lup[4][3] = -1; 
	lup[5][0] = 0; lup[5][1] = 2; lup[5][2] = lup[5][3] = -1; 

	lup[6][0] = 0; lup[6][1] = 1; lup[6][2] = 2; lup[6][3] = 3;  

	lup[7][0] = 1; lup[7][1] = 2; lup[7][2] = lup[7][3] = -1; 
	lup[8][0] = 1; lup[8][1] = 2; lup[8][2] = lup[8][3] = -1; 

	lup[9][0] = 0; lup[9][1] = 3; lup[9][2] = 1; lup[9][3] = 2;  

	lup[10][0] = 0; lup[10][1] = 2; lup[10][2] = lup[10][3] = -1; 
	lup[11][0] = 2; lup[11][1] = 3; lup[11][2] = lup[11][3] = -1; 
	lup[12][0] = 1; lup[12][1] = 3; lup[12][2] = lup[12][3] = -1; 
	lup[13][0] = 0; lup[13][1] = 1; lup[13][2] = lup[13][3] = -1; 
	lup[14][0] = 0; lup[14][1] = 3; lup[14][2] = lup[14][3] = -1; 

	lup[15][0] = lup[15][1] = lup[15][2] = lup[15][3] = -1;
}

int IdentifyCase(float *f_val, float iso_value) 
{
  int rv = 0;

  if (f_val[0] > iso_value) rv += 1;
  if (f_val[1] > iso_value) rv += 2;
  if (f_val[3] > iso_value) rv += 4;
  if (f_val[2] > iso_value) rv += 8;
	
  return rv;
};

void AddIsoLines(float *X, float *Y, float *F, const int *dims, SegmentList *sl, float iso_value) {

	// precompute tables
	int numSegments[16] = 
	{             0, 
		1, 1, 1, 1, 1, 
		2, 1, 1, 2, 1,
		1, 1, 1, 1, 0 };

	int lup[16][4];
	LookUpTable(lup);

	// iterate over cells
	for (int cell_id = 0; cell_id < GetNumberOfCells(dims); ++cell_id) {

		// get X, Y values
		float bbox[4];
		BoundingBoxForCell(X, Y, dims, cell_id, bbox);

    // get Field values
    int idx[2];
    GetLogicalCellIndex(idx, cell_id, dims);

    float f_val[4];
    f_val[0] = F[idx[1]*dims[0] + idx[0]];
    f_val[1] = F[idx[1]*dims[0] + idx[0] + 1];
    f_val[2] = F[(idx[1] + 1)*dims[0] + idx[0] + 1];
    f_val[3] = F[(idx[1] + 1)*dims[0] + idx[0]];

		int icase = IdentifyCase(f_val, iso_value);
    int nsegments = numSegments[icase];

    // draw segments
    for (int i = 0; i < nsegments; ++i) {
      float line[2][2];

      for (int e = 0; e < 2; ++e) {
        int edge = lup[icase][2*i + e];

        // get field values for this edge
        float f1 = f_val[edge];
        float f2 = f_val[(edge + 1) % 4];

        // determine X, Y values of edge 
        float x1, x2, y1, y2;
        if (edge == 0) { x1 = bbox[0]; x2 = bbox[1]; y1 = bbox[2]; y2 = bbox[2]; }
        if (edge == 1) { x1 = bbox[1]; x2 = bbox[1]; y1 = bbox[2]; y2 = bbox[3]; }
        if (edge == 2) { x1 = bbox[1]; x2 = bbox[0]; y1 = bbox[3]; y2 = bbox[3]; }
        if (edge == 3) { x1 = bbox[0]; x2 = bbox[0]; y1 = bbox[3]; y2 = bbox[2]; }

        // interpolate X, Y values of new point
        line[e][0] = InterpLin(f1, f2, iso_value, x1, x2);
        line[e][1] = InterpLin(f1, f2, iso_value, y1, y2);
      }

      sl->AddSegment(line[0][0], line[0][1], line[1][0], line[1][1]); 
      //printf("\n\n");
      //printf("Cell: %d, iso_value: %f\n", cell_id, iso_value);
      //printf("Case: %d\n", icase);
      //printf("x1: %f, x2: %f, y1: %f, y2: %f\n", bbox[0], bbox[1], bbox[2], bbox[3]);
      //printf("f1: %f, f2: %f, f3: %f, f4: %f\n", f_val[0], f_val[1], f_val[2], f_val[3]);
      //printf("(%f, %f) -> (%f, %f)\n", line[0][0], line[0][1], line[1][0], line[1][1]);
    }
	}

}

int main()
{
	int  i, j;

	vtkDataSetReader *rdr = vtkDataSetReader::New();
	rdr->SetFileName("proj5.vtk");
	rdr->Update();

	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
	rgrid->GetDimensions(dims);

	float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
	float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	// Add 4 segments that put a frame around your isolines.  This also
	// documents how to use "AddSegment".
	SegmentList sl;
	sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
	sl.AddSegment(-10, +10, +10, +10);
	sl.AddSegment(-10, -10, -10, +10);
	sl.AddSegment(+10, -10, +10, +10);

	float iso_value = 3.2;
	AddIsoLines(X, Y, F, dims, &sl, iso_value);

	vtkPolyData *pd = sl.MakePolyData();

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
