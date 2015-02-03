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
	_surface->SetValue(0, ISO_VALUE);
	_renderer->SetBackground(0, 0, 0.2);


	vtkSmartPointer<myVtkInteractorStyleImage3D> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage3D>::New();
	_renderWindow->SetWindowName("Mesh Viewer");
	_renderWindow->AddRenderer(_renderer);

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
	_mapper->SetLookupTable(lut);
	_mapper->SetInputConnection(_surface->GetOutputPort());
	_mapper->SetScalarModeToUsePointData();
	_mapper->Update();

	_actor->SetMapper(_mapper);
	
	// Create a cursor (sphere)
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(_mapper->GetCenter());
	sphereSource->SetRadius(5.0);

	vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	sphereMapper->SetInputData(sphereSource->GetOutput());

	vtkSmartPointer<vtkActor> sphereActor =
		vtkSmartPointer<vtkActor>::New();
	sphereActor->SetMapper(sphereMapper);
	//Set the color of the sphere
	sphereActor->GetProperty()->SetColor(0.0, 1.0, 0.0); 

	_renderer->AddActor(sphereActor);
	_renderer->AddActor(_actor);

	_renderWindow->Render();
	myInteractorStyle->Initialize("Output.obj",sphereSource);
	_interactor->Start();
}