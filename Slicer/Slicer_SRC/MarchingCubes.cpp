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
	//_surface->ComputeScalarsOff();
	// Setup the colors array for crosshair
	//vtkSmartPointer<vtkIntArray> selection_colors =
	//	vtkSmartPointer<vtkIntArray>::New();
	//selection_colors->SetNumberOfComponents(1);//4
	//vtkPolyData* mesh = _surface->GetOutput();
	//mesh->AllocateScalars(VTK_INT, 1);
	//cout << selection->GetNumberOfCells() << endl;
	// Add the colors we created to the colors array	
	//selection_colors->SetNumberOfTuples(selection->GetNumberOfPoints());
	//selection_colors->SetName("selection_colors");
	//selection->GetPointData()->SetScalars(selection_colors);
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
	//
	//vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	//lut->SetNumberOfTableValues(3);
	//lut->SetRange(0, 2);
	//lut->SetTableValue((vtkIdType)NOT_ACTIVE, 0.9, 0.9, 0.9, 1.0);
	//lut->SetTableValue((vtkIdType)FOREGROUND, 1, 0, 0, 1.0);
	//lut->SetTableValue((vtkIdType)BACKGROUND, 0, 1, 0, 1.0);
	//lut->Build();
	//_mapper->SetLookupTable(lut);
	_mapper->SetInputConnection(_surface->GetOutputPort());
	//_mapper->SetScalarModeToUsePointData();
	_mapper->Update();
	//

	_actor->SetMapper(_mapper);

	_renderer->AddActor(_actor);

	_renderWindow->Render();
	myInteractorStyle->Initialize("Output.obj");
	_interactor->Start();
}