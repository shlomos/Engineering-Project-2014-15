#include "MarchingCubes.h"


MarchingCubes::MarchingCubes(vtkStructuredPoints* selection) {
	cout << "in ctor start!" << endl;
	this->_selection = selection;

	_surface = vtkMarchingCubes::New();
	_mapper = vtkPolyDataMapper::New();
	_actor = vtkActor::New();
	_mask = vtkImageMaskBits::New();

	//filter segblock, only SEGMENTATION will be presented
	_mask->AddInputData(_selection);
	_mask->SetMask(static_cast<unsigned int>(FOREGROUND));
	_mask->SetOperationToAnd();

	double bounds[6];
	_selection->GetBounds(bounds);
	_surface->SetInputConnection(_mask->GetOutputPort());
	_surface->ComputeNormalsOn();
	_surface->SetValue(0, ISO_VALUE);

	//qvtk3Ddisplayer->SetRenderWindow(renderWindow);
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(0.0, 0.0, 0.0);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(_surface->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	renderer->AddActor(actor);

	renderWindow->Render();
	interactor->Start();

	cout << "After Rendering" << endl;
	//renderWindowInteractor->Start();
	cout << "Done function!" << endl;
}