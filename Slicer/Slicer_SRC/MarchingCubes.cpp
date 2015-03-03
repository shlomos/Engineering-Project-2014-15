#include "MarchingCubes.h"

MarchingCubes::MarchingCubes(vtkStructuredPoints* selection) {

	this->_selection = selection;

	_surface = vtkMarchingCubes::New();
	_renderer = vtkRenderer::New();
	_renderWindow = vtkRenderWindow::New();
	_mapper = vtkPolyDataMapper::New();
	_actor = vtkActor::New();
	_mask = vtkImageMaskBits::New();
	_interactor = vtkRenderWindowInteractor::New();

	//filter segblock, only SEGMENTATION will be presented
	_mask->AddInputData(_selection);
	_mask->SetMask(static_cast<unsigned int>(FOREGROUND));
	_mask->SetOperationToAnd();

	double bounds[6];
	_selection->GetBounds(bounds);
	_surface->SetInputConnection(_mask->GetOutputPort());
	_surface->ComputeNormalsOn();
	_surface->ComputeScalarsOff();
	_surface->SetValue(0, ISO_VALUE);
	_surface->Update();
	_renderer->SetBackground(0, 0, 0.2);


	vtkSmartPointer<myVtkInteractorStyleImage3D> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage3D>::New();
	_renderWindow->SetWindowName("Mesh Viewer");

	_interactor->SetInteractorStyle(myInteractorStyle);
	_interactor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);

	_interactor->SetRenderWindow(_renderWindow);
	_interactor->Initialize();
	_renderWindow->SetInteractor(_interactor);
	
	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetNumberOfTableValues(4);
	lut->SetRange(0, 3);
	lut->SetTableValue(0, 0.9, 0.9, 0.9, 1.0);
	lut->SetTableValue(1, 0.9, 0.9, 0.9, 1.0);
	lut->SetTableValue(2, 0, 1, 0, 1.0);
	lut->SetTableValue(3, 1, 0, 0, 1.0);
	lut->Build();
	vtkPolyData* mesh = _surface->GetOutput();


	//try close holes:
	//vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
	//	vtkSmartPointer<vtkFillHolesFilter>::New();
	//fillHolesFilter->SetInputData(mesh);
	//fillHolesFilter->SetHoleSize(1000.0);
	//fillHolesFilter->Update();

	// Make the triangle windong order consistent
	//vtkSmartPointer<vtkPolyDataNormals> normals =
	//	vtkSmartPointer<vtkPolyDataNormals>::New();
	//normals->SetInputData(fillHolesFilter->GetOutput());
	//normals->ConsistencyOn();
	//normals->SplittingOff();
	//normals->Update();

	//end close holes

	_mapper->SetScalarRange(0,3);
	_mapper->SetLookupTable(lut);
	//cout << "Number of points is: " << mesh->GetNumberOfPoints() << endl;
	vtkSmartPointer<vtkIntArray> mesh_colors =
		vtkSmartPointer<vtkIntArray>::New();
	// Add the colors we created to the colors array	
	//mesh_colors->SetNumberOfValues(normals->GetOutput()->GetNumberOfPoints());
	mesh_colors->SetNumberOfValues(mesh->GetNumberOfPoints());
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++){
		mesh_colors->SetValue(i, 2);
	}
	mesh_colors->SetName("mesh_colors");
	mesh->GetPointData()->SetScalars(mesh_colors);
	mesh_colors->Modified();
	// Restore the original normals
	//normals->GetOutput()->GetPointData()->
	//	SetNormals(mesh->GetPointData()->GetNormals());
	//normals->GetOutput()->GetPointData()->
	//	SetScalars(mesh_colors);
	//_mapper->SetInputData(normals->GetOutput());
	_mapper->SetInputData(mesh);
	//_mapper->SetScalarModeToUsePointData();
	_mapper->Update();

	_actor->SetMapper(_mapper);
	_renderWindow->AddRenderer(_renderer);
	myInteractorStyle->SetDefaultRenderer(_renderer);
	myInteractorStyle->Modified();
	_renderer->AddActor(_actor);
	_renderWindow->SetCurrentCursor(VTK_CURSOR_CROSSHAIR);
	_renderer->ResetCamera();
	_renderWindow->Render();
	myInteractorStyle->Initialize("Output.obj");
	_interactor->Start();
}