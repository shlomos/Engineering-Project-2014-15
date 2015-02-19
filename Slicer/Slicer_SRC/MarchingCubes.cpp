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
	lut->SetNumberOfTableValues(3);
	lut->SetRange(0, 2);
	lut->SetTableValue((vtkIdType)NOT_ACTIVE, 0.9, 0.9, 0.9, 1.0);
	lut->SetTableValue((vtkIdType)FOREGROUND, 1, 0, 0, 1.0);
	lut->SetTableValue((vtkIdType)BACKGROUND, 0, 1, 0, 1.0);
	lut->Build();
	vtkPolyData* mesh = _surface->GetOutput();


	//try close holes:
	vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
		vtkSmartPointer<vtkFillHolesFilter>::New();
	fillHolesFilter->SetInputData(mesh);
	fillHolesFilter->SetHoleSize(1000.0);
	fillHolesFilter->Update();

	// Make the triangle windong order consistent
	vtkSmartPointer<vtkPolyDataNormals> normals =
		vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputData(fillHolesFilter->GetOutput());
	normals->ConsistencyOn();
	normals->SplittingOff();
	normals->Update();

	//end close holes

	_mapper->SetScalarRange(0,2);
	_mapper->SetLookupTable(lut);
	//cout << "Number of points is: " << mesh->GetNumberOfPoints() << endl;
	vtkSmartPointer<vtkIntArray> mesh_colors =
		vtkSmartPointer<vtkIntArray>::New();
	// Add the colors we created to the colors array	
	mesh_colors->SetNumberOfValues(normals->GetOutput()->GetNumberOfPoints());
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++){
		mesh_colors->SetValue(i, NOT_ACTIVE);
	}
	mesh_colors->SetName("mesh_colors");
	mesh->GetPointData()->SetScalars(mesh_colors);
	mesh_colors->Modified();
	// Restore the original normals
	normals->GetOutput()->GetPointData()->
		SetNormals(mesh->GetPointData()->GetNormals());
	normals->GetOutput()->GetPointData()->
		SetScalars(mesh_colors);
	_mapper->SetInputData(normals->GetOutput());
	//_mapper->SetScalarModeToUsePointData();
	_mapper->Update();

	_actor->SetMapper(_mapper);
	//_actor->GetProperty()->SetDiffuseColor(1.0, 0.3882, 0.2784);
	
	// Create a cursor (sphere)
	//vtkSmartPointer<vtkSphereSource> sphereSource =
	//	vtkSmartPointer<vtkSphereSource>::New();
	//sphereSource->SetCenter(_mapper->GetCenter());
	//sphereSource->SetRadius(5.0);

	//vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
	//	vtkSmartPointer<vtkPolyDataMapper>::New();
	//sphereMapper->SetInputData(sphereSource->GetOutput());

	//vtkSmartPointer<vtkActor> sphereActor =
	//	vtkSmartPointer<vtkActor>::New();
	//sphereActor->SetMapper(sphereMapper);
	////Set the color of the sphere
	//sphereActor->GetProperty()->SetColor(0.0, 1.0, 0.0); 
	CrosshairFactory* crossFac = CrosshairFactory::getInstance();
	vtkSmartPointer<vtkActor> crosshair= crossFac->makeCrosshair(_surface->GetOutput()->GetBounds());
	crosshair->GetMapper()->Update();
	crosshair->SetPickable(false);
	_actor->SetPickable(true);
	vtkSmartPointer<vtkRenderer> overlayRenderer = vtkSmartPointer<vtkRenderer>::New();
	//vtkInteractorStyle::SafeDownCast(_interactor->GetInteractorStyle())->AutoAdjustCameraClippingRangeOn();
	//vtkSmartPointer<vtkCamera> cameraOverlay = vtkSmartPointer<vtkCamera>::New();
	_renderWindow->SetNumberOfLayers(2);
	overlayRenderer->SetLayer(1);
	_renderer->SetLayer(0);
	_renderer->SetInteractive(1);
	overlayRenderer->SetInteractive(0);
	// sync camera
	//cameraOverlay = overlayRenderer->GetActiveCamera();
	//cameraOverlay->ShallowCopy(_renderer->GetActiveCamera());
	_renderWindow->AddRenderer(overlayRenderer);
	_renderWindow->AddRenderer(_renderer);
	overlayRenderer->AddActor(crosshair);
	_renderer->AddActor(_actor);
	_renderWindow->Render();
	_renderer->ResetCamera();
	overlayRenderer->ResetCamera();
	_renderWindow->Render();
	myInteractorStyle->Initialize("Output.obj",overlayRenderer);
	_interactor->Start();
}