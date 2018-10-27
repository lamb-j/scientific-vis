#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

vtkPolyData *MakePolyData(float *tl, int ntriangles) {

  int numPoints = 3*(ntriangles);

  vtkPoints *vtk_pts = vtkPoints::New();
  vtk_pts->SetNumberOfPoints(numPoints);
  int ptIdx = 0;

  vtkCellArray *tris = vtkCellArray::New();
  tris->EstimateSize(numPoints,3);

	for (int i = 0 ; i < ntriangles ; i++)
	{
		double pt[3];
		pt[0] = tl[9*i+0];
		pt[1] = tl[9*i+1];
		pt[2] = tl[9*i+2];
		vtk_pts->SetPoint(ptIdx, pt);
		pt[0] = tl[9*i+3];
		pt[1] = tl[9*i+4];
		pt[2] = tl[9*i+5];
		vtk_pts->SetPoint(ptIdx+1, pt);
		pt[0] = tl[9*i+6];
		pt[1] = tl[9*i+7];
		pt[2] = tl[9*i+8];
		vtk_pts->SetPoint(ptIdx+2, pt);
		vtkIdType ids[3] = { ptIdx, ptIdx+1, ptIdx+2 };
		tris->InsertNextCell(3, ids);
		ptIdx += 3;
	}

	vtkPolyData *pd = vtkPolyData::New();
	pd->SetPoints(vtk_pts);
	pd->SetPolys(tris);
	tris->Delete();
	vtk_pts->Delete();

	return pd;
}

void RenderSurface(float *tl, int num_tris) {
  vtkPolyData *pd = MakePolyData(tl, num_tris);

  vtkSmartPointer<vtkDataSetMapper> win1Mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  win1Mapper->SetInputData(pd);
  //win1Mapper->SetScalarRange(0, 0.15);

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
  //renWin->SetSize(800, 800);

  //ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
  //ren1->GetActiveCamera()->SetPosition(0,0,50);
  //ren1->GetActiveCamera()->SetViewUp(0,1,0);
  //ren1->GetActiveCamera()->SetClippingRange(20, 120);
  //ren1->GetActiveCamera()->SetDistance(30);

  iren->Initialize();
  iren->Start();
}
