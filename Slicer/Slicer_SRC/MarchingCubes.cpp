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

	_renderer->SetBackground(0, 0, 0);

	vtkSmartPointer<myVtkInteractorStyleImage3D> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage3D>::New();
	myInteractorStyle->Initialize("Output.obj");
	_renderWindow->SetWindowName("Mesh Viewer");
	_renderWindow->AddRenderer(_renderer);

	_interactor->SetInteractorStyle(myInteractorStyle);
	_interactor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);

	_interactor->SetRenderWindow(_renderWindow);
	_interactor->Initialize();
	_renderWindow->SetInteractor(_interactor);

	_mapper->SetInputConnection(_surface->GetOutputPort());

	_actor->SetMapper(_mapper);

	_renderer->AddActor(_actor);

	_renderWindow->Render();
	_interactor->Start();

}