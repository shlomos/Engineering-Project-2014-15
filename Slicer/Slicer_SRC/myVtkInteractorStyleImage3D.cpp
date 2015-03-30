#include "myVtkInteractorStyleImage3D.h"

myVtkInteractorStyleImage3D::myVtkInteractorStyleImage3D()
{}

void myVtkInteractorStyleImage3D::SetStatusMapper(vtkTextMapper* statusMapper) {
	_StatusMapper = statusMapper;
}

void myVtkInteractorStyleImage3D::Initialize(std::string outputName, std::string inputName, vtkStructuredPoints* selection){
	_drawSize = DEFAULT_DRAW_SIZE;
	_lal = LeapAbstractionLayer::getInstance();
	_outputName = outputName;
	_inputName = inputName;
	_hfMode = false;
	_rotLock = false;
	_currSource = -1;
	this->_selection = selection;
	//Start movement:
	this->StartRotate();
	cout << "3D Interactor iniitiated." << std::endl;
}

void myVtkInteractorStyleImage3D::LoadFromFile(){
}

void myVtkInteractorStyleImage3D::ResetAll(){
	vtkRenderer* ren = this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)ren->GetActors()->GetLastActor()->GetMapper();
	vtkPolyData* mesh = mapper->GetInput();
	vtkPointData* pd = mesh->GetPointData();
	vtkUnsignedShortArray* scalars = (vtkUnsignedShortArray*)pd->GetScalars();
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++){
		scalars->SetValue(i, NOT_ACTIVE);
	}
	scalars->Modified();
	mesh->Modified();
	pd->Modified();
	mapper->Update();
}

void myVtkInteractorStyleImage3D::WriteToFile() {
	cout << "Writing segmentation to file..." << endl;
	vtkRenderer* ren = this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)ren->GetActors()->GetLastActor()->GetMapper();
	vtkPolyData* mesh = mapper->GetInput();
	vtkSmartPointer<vtkPolyDataWriter> writer =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(mesh);
	writer->SetFileName(_outputName.c_str());
	writer->Write();
	cout << "Done." << endl;
}

typedef void(myVtkInteractorStyleImage3D::*workerFunction)();

void myVtkInteractorStyleImage3D::RemoveLeaks(){
	boost::mutex::scoped_lock scoped_lock(_canSegment_mutex);
	cout << "Started fixing mesh!" << endl;
	typedef double PixelType;
	typedef unsigned short SegPixelType;
	const unsigned char dim = 3;

	typedef itk::Image<PixelType, dim> ImageType;
	typedef itk::Image<SegPixelType, dim> SegImageType;

	MeshLeaksCorrector mlc;
	ImageType::Pointer inputImage = mlc.read3DImage<ImageType>(_inputName.c_str());
	cout << "Image reading done." << endl;
	typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;
	GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
	gradientFilter->SetInput(inputImage);
	gradientFilter->Update();
	ImageType::Pointer gradientInputImage = gradientFilter->GetOutput();

	cout << "Produced gradient image" << endl;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	typedef itk::ImageRegionIterator<SegImageType> SegIteratorType;

	IteratorType gradIt(gradientFilter->GetOutput(), gradientFilter->GetOutput()->GetLargestPossibleRegion());
	cout << "GradIt" << endl;
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
	vtkSmartPointer<vtkPolyData> vtkSegImage = vtkSmartPointer<vtkPolyData>::New();
	vtkSegImage->DeepCopy(mapper->GetInput());
	vtkSegImage->GetPointData()->SetScalars(NULL);

	cout << "Extracted mesh from actor" << endl;

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmoothPolyDataFilter::New();
	smoother->SetInputData(vtkSegImage);
	smoother->SetNumberOfIterations(MESH_SMOOTH_ITERATIONS);
	smoother->Update();
	vtkSmartPointer<vtkPolyData> mesh = smoother->GetOutput();
	cout << "Mesh smoothed" << endl;
	typedef itk::ImageToVTKImageFilter<ImageType> ConverterType;
	ConverterType::Pointer gradientConverter = ConverterType::New();
	gradientConverter->SetInput(inputImage);
	gradientConverter->Update();
	vtkSmartPointer<vtkImageData> vtkGradientImage = gradientConverter->GetOutput();

	cout << "Read CT image" << endl;

	vtkSmartPointer<vtkProbeFilter> probeFilter = vtkProbeFilter::New();
	probeFilter->SetSourceData(vtkGradientImage);
	probeFilter->SetInputData(mesh);
	probeFilter->Update();

	cout << "Probe finished" << endl;
	vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkGeometryFilter::New();
	geometryFilter->SetInputData(probeFilter->GetOutput());
	geometryFilter->Update();
	vtkSmartPointer<vtkPolyData> gradientMesh = geometryFilter->GetOutput();

	cout << "Geometric filter finished" << endl;

	vtkSmartPointer<vtkPolyData> minCurvatureMesh = mlc.polyDataToMinCurvature(mapper->GetInput());
	cout << "mlc.minCurv finished" << endl;

	//just temporary - don't forget to delete
	vtkSmartPointer<vtkPolyData> maxCurvatureMesh = mlc.polyDataToMaxCurvature(mapper->GetInput());
	//  --- up to here
	cout << "mlc.maaxCurv finished" << endl;
	vtkSmartPointer<vtkPolyData> minCutMeshLeaks = mlc.minCut(minCurvatureMesh, mapper->GetInput(), gradientMesh, MIN_CURVATURE_TAG, 1.0f);
	cout << "minCut finished" << endl;
	vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks = mlc.minCut(maxCurvatureMesh, mapper->GetInput(), gradientMesh, MAX_CURVATURE_TAG, 1.0f);
	cout << "minCut Interior finished" << endl;
	vtkSmartPointer<vtkPolyData> minCutMesh = mlc.minCutConjunction(minCutMeshLeaks, minCutMeshInteriorLeaks);
	cout << "Conjunction finished" << endl;
	mlc.attributeDilation(minCutMesh, ATTRIBUTE_DILATION_RADIUS);
	cout << "dilation finished" << endl;
	vtkSmartPointer<vtkPolyData> correctedMesh1 = mlc.laplaceInterpolation(minCutMesh);
	cout << "laplace finished" << endl;
	vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
	normals->SetInputData(correctedMesh1);
	normals->FlipNormalsOn();
	normals->Update();
	vtkSmartPointer<vtkPolyData> correctedMesh = normals->GetOutput();
	cout << "Finished fixing mesh! Wow...." << endl;
	cout << "Writing mesh to file! laplaceMesh.vtk" << endl;
	mlc.writePolyData(correctedMesh, "laplaceMesh.vtk");
	vtkSmartPointer<vtkDecimatePro> decimateFilter = vtkDecimatePro::New();
	decimateFilter->SetInputData(correctedMesh);
	decimateFilter->SetTargetReduction(DECIMATION_FACTOR);
	decimateFilter->PreserveTopologyOn();
	decimateFilter->Update();
	vtkSmartPointer<vtkPolyData> decimatedMesh = decimateFilter->GetOutput();
	cout << "Writing mesh to file! decimatedMesh.vtk" << endl;
	mlc.writePolyData(decimatedMesh, "decimatedMesh.vtk");
	//2.5D NEW code
	//Seg input image is the selection_structured_points. (upcast to vtkImageData)
	typedef itk::VTKImageToImageFilter<SegImageType> SegConverterType;
	SegConverterType::Pointer converter = SegConverterType::New();
	cout << "bedug 1" << endl;
	// selection is not NULL (checked...)
	cout << "this->_selection->GetScalarType(): " << this->_selection->GetScalarType() << endl;
	if (this->_selection->GetScalarType() == VTK_INT){
		cout << "I have INTs instead of shorts for some reason!" << endl;
	}
	converter->SetInput(this->_selection);
	cout << "deb 5" << endl;
	converter->Update(); // 
	cout << "bedug 2" << endl;
	SegImageType::Pointer segInputImage = converter->GetOutput();
	cout << "bedug 3" << endl;
	SegImageType::Pointer outputContourImage = mlc.sampleMeshOnImage<SegImageType>(/*decimatedMesh*/correctedMesh1, segInputImage);
	cout << "bedug 4" << endl;
	SegImageType::Pointer seedInputImage = mlc.sampleMeshOnImage<SegImageType>(mapper->GetInput(), segInputImage);
	cout << "bedug 5" << endl;
	SegImageType::Pointer outputImage = mlc.correctImage<SegImageType>(segInputImage, seedInputImage, outputContourImage);

	typedef itk::ImageFileWriter<SegImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	cout << "bedug 6" << endl;
	writer->SetInput(outputImage);
	cout << "Writing mesh to file! correctedImage.nii.gz" << endl;
	writer->SetFileName("correctedImage.nii.gz");
	try{
		writer->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "writing input image exception thrown" << std::endl;
		std::cerr << excp << std::endl;
		exit(1);
	}



	typedef itk::ImageToVTKImageFilter<SegImageType> SegInvConverterType;
	SegInvConverterType::Pointer correctionConverter = SegInvConverterType::New();
	cout << "bedug 7" << endl;
	correctionConverter->SetInput(outputImage);
	correctionConverter->Update();
	vtkSmartPointer<vtkImageData> vtkCorrectedImage = correctionConverter->GetOutput();
	/*vtkSmartPointer<vtkImageToStructuredPoints> convertFilter =
	vtkSmartPointer<vtkImageToStructuredPoints>::New();
	convertFilter->SetInputData(correctionConverter->GetOutput());
	convertFilter->Update();
	this->_selection->DeepCopy(convertFilter->GetOutput());*/


	cout << "bedug 8" << endl;
	vtkSmartPointer<vtkPolyData> outputMesh = mlc.createAndSmoothSurface(vtkCorrectedImage, 50);
	mlc.writePolyData(outputMesh, "final.vtk");
	cout << "bedug 9" << endl;
	// Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapperE =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapperE->SetInputData(outputMesh);

	vtkSmartPointer<vtkActor> actorE =
		vtkSmartPointer<vtkActor>::New();
	actorE->SetMapper(mapperE);

	// A renderer and render window
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// An interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	renderer->AddActor(actorE);
	renderer->SetBackground(.1, .2, .3); // Background color dark blue

	// Render
	renderWindow->SetWindowName("Corrected Mesh");
	renderWindow->Render();

	// Begin mouse interaction
	renderWindowInteractor->Start();
	cout << "We won!" << endl;
}

void myVtkInteractorStyleImage3D::OnLeftButtonDown()
{
	std::cout << "Pressed left mouse button." << std::endl;
	//MakeAnnotation(2);
	// Forward events
	vtkInteractorStyleJoystickCamera::OnLeftButtonDown();
}

void myVtkInteractorStyleImage3D::OnRightButtonDown()
{
	std::cout << "Pressed left mouse button." << std::endl;
	//MakeAnnotation(3);
	// Forward events
	vtkInteractorStyleJoystickCamera::OnRightButtonDown();
}

void myVtkInteractorStyleImage3D::OnRightButtonUp()
{
}

void myVtkInteractorStyleImage3D::OnLeftButtonUp()
{
}

void myVtkInteractorStyleImage3D::OnTimer(){

	HWND forground = GetForegroundWindow();
	if (forground) {
		char window_title[256];
		GetWindowText(forground, window_title, 256);

		if (!strcmp("Mesh Viewer", window_title)) {
			// render
			int * winSize = Interactor->GetRenderWindow()->GetSize();
			Interactor->SetEventPosition((2 * _lal->getX() + winSize[0] / 2), (_lal->getY() / LEAP_MAX_Y)*winSize[1]);
			Interactor->GetRenderWindow()->SetCursorPosition(Interactor->GetEventPosition()[0], Interactor->GetEventPosition()[1]);
		}
		vtkInteractorStyleJoystickCamera::OnTimer();
	}
}

void myVtkInteractorStyleImage3D::GetPoked2DLocation(int ijk[3]){
	int* clickPos = this->GetInteractor()->GetEventPosition();

	// Pick from this location.
	vtkSmartPointer<vtkPropPicker> picker =
		vtkSmartPointer<vtkPropPicker>::New();
	picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

	double* pos = picker->GetPickPosition();
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
	vtkPolyData* mesh = mapper->GetInput();
	vtkIdType pointId = mesh->FindPoint(pos);
	double coords[3];
	mesh->GetPoint(pointId,coords);
	double pcoords[3];
	this->_selection->ComputeStructuredCoordinates(coords, ijk, pcoords);
}

void myVtkInteractorStyleImage3D::MakeAnnotation(vtkIdType annotation){
	int* clickPos = this->GetInteractor()->GetEventPosition();

	// Pick from this location.
	vtkSmartPointer<vtkPropPicker>  picker =
		vtkSmartPointer<vtkPropPicker>::New();
	picker->Pick(clickPos[0], clickPos[1], clickPos[2], this->GetDefaultRenderer());

	double* pos = picker->GetPickPosition();
	std::cout << "Pick position (world coordinates) is: "
		<< pos[0] << " " << pos[1]
		<< " " << pos[2] << std::endl;

	std::cout << "Picked actor: " << picker->GetActor() << std::endl;
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
	vtkUnsignedShortArray* scalars = (vtkUnsignedShortArray*)mapper->GetInput()->GetPointData()->GetScalars();
	vtkPolyData* mesh = mapper->GetInput();
	vtkIdType id = mesh->FindPoint(pos);
	cout << "The Point Is: " << id << endl;
	cout << "The Pos Is: " << clickPos[0] << "," << clickPos[1] << "," << pos[2] << endl;
	cout << "Number of components is: " << scalars->GetNumberOfComponents() << endl;
	if (id > -1 && _currSource == -1){
		_currSource = id;
		cout << "Start set." << endl;
	}
	else if (id > -1 && _currSource != -1){
		//Make a line connection.
		cout << "End set." << endl;
		vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
			vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
		dijkstra->SetInputData(mesh);
		dijkstra->SetStartVertex(_currSource);
		dijkstra->SetEndVertex(id);
		dijkstra->Update();
		cout << "Found Shortest path!" << endl;
		vtkPolyData* path = dijkstra->GetOutput();
		cout << "Num. of points: " << path->GetNumberOfPoints() << endl;
		for (int i = 0; i < path->GetNumberOfPoints(); i++){
			scalars->SetValue(mesh->FindPoint(path->GetPoint(i)), annotation);
		}
		_currSource = -1;
	}
	if (id > -1){
		//cout << "component size: " << scalars->GetElementComponentSize() << endl;
		cout << scalars->GetValue(id) << endl;
		scalars->SetValue(id, annotation);
		cout << scalars->GetValue(id) << endl;
	}
	scalars->Modified();
	mesh->Modified();
	mapper->Update();
	this->Interactor->Render();
}

void myVtkInteractorStyleImage3D::OnKeyUp() {}
void myVtkInteractorStyleImage3D::OnChar() {}
void myVtkInteractorStyleImage3D::OnMouseMove(){
	vtkInteractorStyleJoystickCamera::OnMouseMove();
}
void myVtkInteractorStyleImage3D::OnKeyDown() {
	std::string key = this->GetInteractor()->GetKeySym();
	if (key.compare("Left") == 0) {
		cout << "Left arrow key was pressed." << endl;

		//set focus
		SetForegroundWindow(FindWindow(NULL, "Slicer"));

	}
	else if (key.compare("1") == 0) {
		cout << "Draw size was changed to " << _drawSize - 1 << endl;
	}
	else if (key.compare("2") == 0) {
		cout << "Draw size was changed to " << _drawSize + 1 << endl;
	}
	else if (key.compare("3") == 0) {
		Interactor->GetRenderWindow()->SetFullScreen(1);
	}
	else if (key.compare("4") == 0) {
		Interactor->GetRenderWindow()->SetFullScreen(0);
	}
	else if (key.compare("space") == 0){
		cout << "change _canSegment To FALSE" << endl;
		workerFunction f = &myVtkInteractorStyleImage3D::RemoveLeaks;
		boost::thread workerThread(f, this);
		this->Interactor->Render();
		cout << "change _canSegment To TRUE" << endl;
	}
	else if (key.compare("o") == 0) {
		cout << "Orientation key was pressed." << endl;
	}
	else if (key.compare("h") == 0) {
		_hfMode = !_hfMode;
		cout << "Hands in the air mode: " << (_hfMode ? "ON" : "OFF") << endl;
	}
	else if (key.compare("s") == 0) {
		WriteToFile();
	}
	else if (key.compare("q") == 0) {
		exit(0);
	}
	else if (key.compare("a") == 0) {
		cout << "a pressed!" << endl;
		if (_rotLock){
			_rotLock = false;
			int x = this->Interactor->GetEventPosition()[0];
			int y = this->Interactor->GetEventPosition()[1];
			this->FindPokedRenderer(x, y);
			if (this->CurrentRenderer == NULL)
			{
				cout << "Renderer is null" << endl;
				return;
			}
			this->GrabFocus(this->EventCallbackCommand);
			this->StartRotate();
		}
		else{
			this->EndRotate();
			this->Interactor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);
			if (this->Interactor)
			{
				this->ReleaseFocus();
			}
			_rotLock = true;
		}
	}
	else if (key.compare("r") == 0) {
		this->ResetAll();
	}
	else if (key.compare("l") == 0) {
		this->LoadFromFile();
	}
	else if (key.compare("f") == 0) {
		//leak
		this->MakeAnnotation(3);
	}
	else if (key.compare("b") == 0) {
		//correct
		this->MakeAnnotation(2);
	}
	else if (key.compare("g") == 0) {
		// Interior miss.
		this->MakeAnnotation(4);
	}
	else if (key.compare("v") == 0) {
		int ijk[3];
		this->GetPoked2DLocation(ijk);
		_lal->RequestUpdate(ijk[0], ijk[1], ijk[2]);
		//this->UpdateContext();
	}
	return;
}

void myVtkInteractorStyleImage3D::OnMouseWheelForward() {
	// don't forward events, otherwise the image will be zoomed
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelForward();
}
void myVtkInteractorStyleImage3D::OnMouseWheelBackward() {
	//std::cout << "Scrolled mouse wheel backward." << std::endl;
	// don't forward events, otherwise the image will be zoomed
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelBackward();
}

myVtkInteractorStyleImage3D::~myVtkInteractorStyleImage3D(){
}

vtkStandardNewMacro(myVtkInteractorStyleImage3D);
