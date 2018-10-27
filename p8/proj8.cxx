#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <vtkImageData.h>
#include <vtkContourFilter.h> 
#include <vtkHedgeHog.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkLineSource.h>
#include "vtkRuledSurfaceFilter.h"
#include "vtkRungeKutta4.h"
#include "vtkStreamTracer.h"
#include "vtkProperty.h"
#include "vtkArrowSource.h"
#include "vtkColorTransferFunction.h"
#include "vtkGlyph3D.h"

int main()
{
	int  i, j;

  double xmins[4] = {0,.5,0,.5};
  double xmaxs[4] = {0.5,1,0.5,1};
  double ymins[4] = {0,0,.5,.5};
  double ymaxs[4]= {0.5,0.5,1,1};

	vtkDataSetReader *rdr = vtkDataSetReader::New();
	rdr->SetFileName("proj8.vtk");
	rdr->Update();
  rdr->GetOutput()->GetPointData()->
    SetActiveAttribute("grad", vtkDataSetAttributes::VECTORS);
	rdr->Update();

  // -------------------------
  // Isosurface
  vtkSmartPointer<vtkColorTransferFunction> colors = 
   vtkSmartPointer<vtkColorTransferFunction>::New();
	colors->AddRGBPoint(2.5, 0,1,1);
	colors->AddRGBPoint(5.0, 1,0,1);

  vtkSmartPointer<vtkContourFilter> surface = vtkSmartPointer<vtkContourFilter>::New();
  surface->SetInputData(rdr->GetOutput());
  surface->SetValue(0, 2.5);
  surface->SetValue(1, 5.0);
  surface->Update();

	vtkSmartPointer<vtkDataSetMapper> win1Mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	win1Mapper->SetInputData(surface->GetOutput());
	win1Mapper->SetLookupTable(colors);
	win1Mapper->SetScalarRange(0, 0.15);

	vtkSmartPointer<vtkActor> win1Actor = vtkSmartPointer<vtkActor>::New();
	win1Actor->SetMapper(win1Mapper);

	vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
	ren1->AddActor(win1Actor);
  ren1->SetViewport(xmins[0], ymins[0], xmaxs[0], ymaxs[0]);

  // -------------------------
  // Slices 
  // vtkCutter

  //-- XY plane
  vtkSmartPointer<vtkPlane> XYplane =
    vtkSmartPointer<vtkPlane>::New();
  XYplane->SetOrigin(0,0,0);
  XYplane->SetNormal(0,0,1);

  vtkSmartPointer<vtkCutter> XYcutter =
    vtkSmartPointer<vtkCutter>::New();
  XYcutter->SetCutFunction(XYplane);
  XYcutter->SetInputConnection(rdr->GetOutputPort());
  XYcutter->Update();

	vtkSmartPointer<vtkDataSetMapper> XYhardyColor = 
  		vtkSmartPointer<vtkDataSetMapper>::New(); 
	XYhardyColor->SetInputConnection(XYcutter->GetOutputPort()); 
	XYhardyColor->SelectColorArray("hardyglobal"); 
	XYhardyColor->SetScalarRange(rdr->GetOutput()->GetScalarRange()); 
	XYhardyColor->Update();

	vtkNew<vtkActor> XYActor; 
	XYActor->SetMapper(XYhardyColor.Get()); 

  //-- XZ plane
  vtkSmartPointer<vtkPlane> XZplane =
    vtkSmartPointer<vtkPlane>::New();
  XZplane->SetOrigin(0,0,0);
  XZplane->SetNormal(0,1,0);

  vtkSmartPointer<vtkCutter> XZcutter =
    vtkSmartPointer<vtkCutter>::New();
  XZcutter->SetCutFunction(XZplane);
  XZcutter->SetInputConnection(rdr->GetOutputPort());
  XZcutter->Update();

	vtkSmartPointer<vtkDataSetMapper> XZhardyColor = 
  		vtkSmartPointer<vtkDataSetMapper>::New(); 
	XZhardyColor->SetInputConnection(XZcutter->GetOutputPort()); 
	XZhardyColor->SelectColorArray("hardyglobal"); 
	XZhardyColor->SetScalarRange(rdr->GetOutput()->GetScalarRange()); 
	XZhardyColor->Update();

	vtkNew<vtkActor> XZActor; 
	XZActor->SetMapper(XZhardyColor.Get()); 

  //-- YZ plane
  vtkSmartPointer<vtkPlane> YZplane =
    vtkSmartPointer<vtkPlane>::New();
  YZplane->SetOrigin(0,0,0);
  YZplane->SetNormal(1,0,0);

  vtkSmartPointer<vtkCutter> YZcutter =
    vtkSmartPointer<vtkCutter>::New();
  YZcutter->SetCutFunction(YZplane);
  YZcutter->SetInputConnection(rdr->GetOutputPort());
  YZcutter->Update();

	vtkSmartPointer<vtkDataSetMapper> YZhardyColor = 
  		vtkSmartPointer<vtkDataSetMapper>::New(); 
	YZhardyColor->SetInputConnection(YZcutter->GetOutputPort()); 
	YZhardyColor->SelectColorArray("hardyglobal"); 
	YZhardyColor->SetScalarRange(rdr->GetOutput()->GetScalarRange()); 
	YZhardyColor->Update();

	vtkNew<vtkActor> YZActor; 
	YZActor->SetMapper(YZhardyColor.Get()); 

	vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();
	ren2->AddActor(XYActor);
	ren2->AddActor(XZActor);
	ren2->AddActor(YZActor);
  ren2->SetViewport(xmins[2], ymins[2], xmaxs[2], ymaxs[2]);

  // -------------------------
  // Hedgehog 
	vtkArrowSource* arrow = vtkArrowSource::New();

  vtkSmartPointer<vtkGlyph3D> glyph3D = 
    vtkSmartPointer<vtkGlyph3D>::New();
  glyph3D->SetInputData(rdr->GetOutput());
  glyph3D->SetSourceConnection(arrow->GetOutputPort());
	glyph3D->SetScaleModeToScaleByVector();
	glyph3D->SetScaleFactor(5.0);
  glyph3D->Update();

	vtkSmartPointer<vtkDataSetMapper> win3Mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	win3Mapper->SetInputData(glyph3D->GetOutput());
	win3Mapper->SetScalarRange(0, 0.15);

	vtkSmartPointer<vtkActor> win3Actor = vtkSmartPointer<vtkActor>::New();
	win3Actor->SetMapper(win3Mapper);

	vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();
	ren3->AddActor(win3Actor);
  ren3->SetViewport(xmins[1], ymins[1], xmaxs[1], ymaxs[1]);

  // -------------------------
  // Streamlines 
  // vtkStreamTracer

  vtkLineSource *rake = vtkLineSource::New();
  rake->SetPoint1(9, 0, 0);
  rake->SetPoint2(-9, 0, 0);
  rake->SetResolution(19);

  vtkPolyDataMapper *rakeMapper = vtkPolyDataMapper::New();
  rakeMapper->SetInputConnection(rake->GetOutputPort());
  vtkActor *rakeActor = vtkActor::New();
  rakeActor->SetMapper(rakeMapper);

  vtkRungeKutta4 *integ = vtkRungeKutta4::New();
  vtkStreamTracer *sl = vtkStreamTracer::New();
  sl->SetInputConnection(rdr->GetOutputPort());
  sl->SetSourceConnection(rake->GetOutputPort());
  sl->SetIntegrator(integ);
  sl->SetMaximumPropagation(5);
  sl->SetInitialIntegrationStep(0.1);

	vtkPolyDataMapper *tracerMapper = vtkPolyDataMapper::New();
	tracerMapper->SetInputConnection(sl->GetOutputPort());
	tracerMapper->SetScalarRange(rdr->GetOutput()->GetScalarRange());

	vtkActor *tracerActor = vtkActor::New();
	tracerActor->SetMapper(tracerMapper);

	vtkSmartPointer<vtkRenderer> ren4 = vtkSmartPointer<vtkRenderer>::New();
	ren4->AddActor(tracerActor);
  ren4->SetViewport(xmins[3], ymins[3], xmaxs[3], ymaxs[3]);

	vtkSmartPointer<vtkRenderWindow> renWin =
		vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren1);  // isosurface
	renWin->AddRenderer(ren2);  // slices
	renWin->AddRenderer(ren3);  // hedgehog
	renWin->AddRenderer(ren4);  // streamline
	renWin->SetSize(800, 800);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);

	// This starts the event loop and invokes an initial render.
	//
	iren->Initialize();
	iren->Start();
}
