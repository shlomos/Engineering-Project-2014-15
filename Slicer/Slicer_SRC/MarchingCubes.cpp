#include "MarchingCubes.h"


MarchingCubes::MarchingCubes(vtkStructuredPoints* selection) {
	cout << "in ctor start!" << endl;
	//this->_selection = selection;

	cout << "DEBUG1" << endl;
	_surface = vtkMarchingCubes::New();
	_renderer = vtkRenderer::New();
	_renderWindow = vtkRenderWindow::New();
	_mapper = vtkPolyDataMapper::New();
	_actor = vtkActor::New();
	_mask = vtkImageMaskBits::New();

	cout << "DEBUG2" << endl;
	//filter segblock, only SEGMENTATION will be presented
	_mask->AddInputData(selection);
	_mask->SetMask(static_cast<unsigned int>(FOREGROUND));
	_mask->SetOperationToAnd();

	cout << "DEBUG3" << endl;
	double bounds[6];
	selection->GetBounds(bounds);
	_surface->SetInputConnection(_mask->GetOutputPort());
	_surface->ComputeNormalsOn();
	_surface->SetValue(0, ISO_VALUE);
	_renderer->SetBackground(0.0, 0.0, 0.0);
	_renderWindow->AddRenderer(_renderer);

	_renderer->ResetCamera(bounds);

	cout << "DEBUG4" << endl;
	//qvtk3Ddisplayer->SetRenderWindow(renderWindow);
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(.1, .2, .3);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(_surface->GetOutputPort());
	mapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	renderer->AddActor(actor);

	renderWindow->Render();
	interactor->Start();
	//_viewer->SetupInteractor();

	//vtkSmartPointer<vtkRenderWindow> renderWindow =
	// vtkSmartPointer<vtkRenderWindow>::New();
	//renderWindow->AddRenderer(renderer);*/
	//cout << "DEBUG2" << endl;

	//cout << "DEBUG3" << endl;
	////vtkSmartPointer<vtkRenderWindowInteractor> interactor =
	//// vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	// vtkSmartPointer<vtkRenderWindowInteractor>::New(;)
	//
	//_viewer->SetupInteractor(renderWindowInteractor);
	//cout << "DEBUG4" << endl;
	//renderWindowInteractor->Initialize();
	////renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	////renderWindowInteractor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);
	//cout << "DEBUG5" << endl;
	////_viewer->GetRenderer()->GetActors()->GetLastActor()->GetMapper()->ScalarVisibilityOff();
	//cout << "after long chain..." << endl;
	//
	//_viewer->Render();
	//_viewer->GetRenderer()->ResetCamera();
	//cout << "Camera reset" << endl;
	//_viewer->Render();
	//
	//qvtk3Ddisplayer->SetRenderWindow(_renderWindow);
	//mapper->SetInputConnection(_surface->GetOutputPort());
	//actor->SetMapper(mapper);
	//renderer->AddActor(actor);

	cout << "After Rendering" << endl;
	//renderWindowInteractor->Start();
	cout << "Done function!" << endl;


}

void MarchingCubes::render3D()
{
	//render
	this->_renderWindow->Render();
	//sqvtk3Ddisplayer->update();
}

void MarchingCubes::flipView(bool flip)
{
	if (flip)
		this->_renderer->GetActiveCamera()->SetRoll(180.0);
	else
		this->_renderer->GetActiveCamera()->SetRoll(0.0);
	this->_renderWindow->Render();
}



//int main(int argc, char *argv[])
//{
// vtkSmartPointer<vtkImageData> volume =
// vtkSmartPointer<vtkImageData>::New();
// double isoValue;
// if (argc < 3)
// {
// vtkSmartPointer<vtkSphereSource> sphereSource =
// vtkSmartPointer<vtkSphereSource>::New();
// sphereSource->SetPhiResolution(20);
// sphereSource->SetThetaResolution(20);
// sphereSource->Update();
//
// double bounds[6];
// sphereSource->GetOutput()->GetBounds(bounds);
// for (unsigned int i = 0; i < 6; i += 2)
// {
// double range = bounds[i + 1] - bounds[i];
// bounds[i] = bounds[i] - .1 * range;
// bounds[i + 1] = bounds[i + 1] + .1 * range;
// }
// vtkSmartPointer<vtkVoxelModeller> voxelModeller =
// vtkSmartPointer<vtkVoxelModeller>::New();
// voxelModeller->SetSampleDimensions(50, 50, 50);
// voxelModeller->SetModelBounds(bounds);
// voxelModeller->SetScalarTypeToFloat();
// voxelModeller->SetMaximumDistance(.1);
//
// voxelModeller->SetInputConnection(sphereSource->GetOutputPort());
// voxelModeller->Update();
// isoValue = 0.5;
// volume->DeepCopy(voxelModeller->GetOutput());
// }
// else
// {
// vtkSmartPointer<vtkDICOMImageReader> reader =
// vtkSmartPointer<vtkDICOMImageReader>::New();
// reader->SetDirectoryName(argv[1]);
// reader->Update();
// volume->DeepCopy(reader->GetOutput());
// isoValue = atof(argv[2]);
// }
//
// vtkSmartPointer<vtkMarchingCubes> surface =
// vtkSmartPointer<vtkMarchingCubes>::New();
//
//#if VTK_MAJOR_VERSION <= 5
// surface->SetInput(volume);
//#else
// surface->SetInputData(volume);
//#endif
// surface->ComputeNormalsOn();
// surface->SetValue(0, isoValue);
//
// vtkSmartPointer<vtkRenderer> renderer =
// vtkSmartPointer<vtkRenderer>::New();
// renderer->SetBackground(.1, .2, .3);
//
// vtkSmartPointer<vtkRenderWindow> renderWindow =
// vtkSmartPointer<vtkRenderWindow>::New();
// renderWindow->AddRenderer(renderer);
// vtkSmartPointer<vtkRenderWindowInteractor> interactor =
// vtkSmartPointer<vtkRenderWindowInteractor>::New();
// interactor->SetRenderWindow(renderWindow);
//
// vtkSmartPointer<vtkPolyDataMapper> mapper =
// vtkSmartPointer<vtkPolyDataMapper>::New();
// mapper->SetInputConnection(surface->GetOutputPort());
// mapper->ScalarVisibilityOff();
//
// vtkSmartPointer<vtkActor> actor =
// vtkSmartPointer<vtkActor>::New();
// actor->SetMapper(mapper);
//
// renderer->AddActor(actor);
//
// renderWindow->Render();
// interactor->Start();
// return EXIT_SUCCESS;
//}