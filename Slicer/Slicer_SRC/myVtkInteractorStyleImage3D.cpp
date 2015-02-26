#include "myVtkInteractorStyleImage3D.h"
#include <vtkRendererCollection.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkRegularPolygonSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include "constants.h"
#include <vtkDataSetMapper.h>
#include <vtkStructuredPoints.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkExtractVOI.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkCellPicker.h>
#include <vtkPointPicker.h>
#include <vtkImageMapToColors.h>
#include <vtkImageOpenClose3D.h>
#include <vtkPropPicker.h>
#include <vtkImageEuclideanDistance.h>
#include <algorithm>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkPolyDataWriter.h>
#include <vtkCamera.h>
#include <sstream>


myVtkInteractorStyleImage3D::myVtkInteractorStyleImage3D()
{}

void myVtkInteractorStyleImage3D::SetStatusMapper(vtkTextMapper* statusMapper) {
	_StatusMapper = statusMapper;
}

void myVtkInteractorStyleImage3D::Initialize(std::string outputName){
	_drawSize = DEFAULT_DRAW_SIZE;
	_lal = LeapAbstractionLayer::getInstance();
	_outputName = outputName;
	_hfMode = false;
	_rotLock = true;
	_currSource = -1;
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
	vtkIntArray* scalars = (vtkIntArray*)pd->GetScalars();
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

void myVtkInteractorStyleImage3D::doSegment() {
}

void myVtkInteractorStyleImage3D::RemoveLeaks(){
	cout << "Started fixing mesh!" << endl;
	typedef double PixelType;
	typedef unsigned short SegPixelType;
	const unsigned char dim = 3;

	typedef itk::Image<PixelType, dim> ImageType;
	typedef itk::Image<SegPixelType, dim> SegImageType;

	MeshLeaksCorrector mlc;
	ImageType::Pointer inputImage = mlc.read3DImage<ImageType>("GG.vtk");
	//SegImageType::Pointer segInputImage = read3DImage<SegImageType>(argv[SEG_INPUT_IMAGE_NAME]);
	//SegImageType::Pointer seedInputImage = read3DImage<SegImageType>(argv[SEED_INPUT_IMAGE_NAME]);
	//SegImageType::Pointer interiorSeedInputImage = read3DImage<SegImageType>(argv[INTERIOR_SEED_NAME]);

	typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;
	GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
	gradientFilter->SetInput(inputImage);
	gradientFilter->Update();
	ImageType::Pointer gradientInputImage = gradientFilter->GetOutput();

	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	typedef itk::ImageRegionIterator<SegImageType> SegIteratorType;

	IteratorType gradIt(gradientFilter->GetOutput(), gradientFilter->GetOutput()->GetLargestPossibleRegion());
	//SegIteratorType seedIt(seedInputImage, seedInputImage->GetLargestPossibleRegion());


	//typedef itk::ImageToVTKImageFilter<SegImageType> SegConverterType;
	//SegConverterType::Pointer converter = SegConverterType::New();
	//converter->SetInput(segInputImage);
	//converter->Update();
	//vtkSmartPointer<vtkImageData> vtkSegImage = converter->GetOutput();

	//SegConverterType::Pointer seedConverter = SegConverterType::New();
	//seedConverter->SetInput(seedInputImage);
	//seedConverter->Update();
	//vtkSmartPointer<vtkImageData> vtkSeedImage = seedConverter->GetOutput();
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
	vtkSmartPointer<vtkPolyData> vtkSegImage = mapper->GetInput();

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmoothPolyDataFilter::New();
	smoother->SetInputData(vtkSegImage);
	smoother->SetNumberOfIterations(MESH_SMOOTH_ITERATIONS);
	smoother->Update();
	vtkSmartPointer<vtkPolyData> mesh = smoother->GetOutput();

	vtkSmartPointer<vtkStructuredPointsReader> reader =
		vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName("GG.vtk");
	reader->Update();
	vtkSmartPointer<vtkImageData> vtkGradientImage = reader->GetOutput();
	/*typedef itk::ImageToVTKImageFilter<ImageType> ConverterType;
	ConverterType::Pointer gradientConverter = ConverterType::New();
	gradientConverter->SetInput(inputImage);
	gradientConverter->Update();
	vtkSmartPointer<vtkImageData> vtkGradientImage = gradientConverter->GetOutput();*/

	vtkSmartPointer<vtkProbeFilter> probeFilter = vtkProbeFilter::New();
	probeFilter->SetSourceData(vtkGradientImage);
	probeFilter->SetInputData(mesh);
	probeFilter->Update();
	vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkGeometryFilter::New();
	geometryFilter->SetInputData(probeFilter->GetOutput());
	geometryFilter->Update();
	vtkSmartPointer<vtkPolyData> gradientMesh = geometryFilter->GetOutput();

	vtkSmartPointer<vtkPolyData> minCurvatureMesh = mlc.polyDataToMinCurvature(mesh);

	//just temporary - don't forget to delete
	vtkSmartPointer<vtkPolyData> maxCurvatureMesh = mlc.polyDataToMaxCurvature(mesh);
	//  --- up to here

	vtkSmartPointer<vtkPolyData> minCutMeshLeaks = mlc.minCut(minCurvatureMesh, mesh, gradientMesh, MIN_CURVATURE_TAG, 1.0f);

	vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks = mlc.minCut(maxCurvatureMesh, mesh, gradientMesh, MAX_CURVATURE_TAG, 1.0f);

	vtkSmartPointer<vtkPolyData> minCutMesh = mlc.minCutConjunction(minCutMeshLeaks, minCutMeshInteriorLeaks);
	mlc.attributeDilation(minCutMesh, ATTRIBUTE_DILATION_RADIUS);
	vtkSmartPointer<vtkPolyData> correctedMesh1 = mlc.laplaceInterpolation(minCutMesh);
	vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
	normals->SetInputData(correctedMesh1);
	normals->FlipNormalsOn();
	normals->Update();
	vtkSmartPointer<vtkPolyData> correctedMesh = normals->GetOutput();
	cout << "Finished fixing mesh! Wow...." << endl;
	cout << "Writing mesh to file! laplaceMesh.vtk" << endl;
	mlc.writePolyData(correctedMesh, "laplaceMesh.vtk");
	cout << "We won!" << endl;
}

void myVtkInteractorStyleImage3D::OnLeftButtonDown()
{
	std::cout << "Pressed left mouse button." << std::endl;
	MakeAnnotation(FOREGROUND);
	// Forward events
	vtkInteractorStyleJoystickCamera::OnLeftButtonDown();
}

void myVtkInteractorStyleImage3D::OnRightButtonDown()
{
	std::cout << "Pressed left mouse button." << std::endl;
	MakeAnnotation(BACKGROUND);
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
	cout << "Got a leap event!" << endl;
	// render
	int * winSize = Interactor->GetRenderWindow()->GetSize();
	Interactor->SetEventPosition((2 * _lal->getX() + winSize[0] / 2), (_lal->getY() / LEAP_MAX_Y)*winSize[1]);
	Interactor->GetRenderWindow()->SetCursorPosition(Interactor->GetEventPosition()[0], Interactor->GetEventPosition()[1]);
	vtkInteractorStyleJoystickCamera::OnTimer();
}


void myVtkInteractorStyleImage3D::MakeAnnotation(vtkIdType annotation){
	int* clickPos = this->GetInteractor()->GetEventPosition();

	// Pick from this location.
	vtkSmartPointer<vtkPropPicker>  picker =
		vtkSmartPointer<vtkPropPicker>::New();
	picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

	double* pos = picker->GetPickPosition();
	std::cout << "Pick position (world coordinates) is: "
		<< pos[0] << " " << pos[1]
		<< " " << pos[2] << std::endl;

	std::cout << "Picked actor: " << picker->GetActor() << std::endl;
	vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
	vtkIntArray* scalars = (vtkIntArray*)mapper->GetInput()->GetPointData()->GetScalars();
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
}

void myVtkInteractorStyleImage3D::OnKeyUp() {}
void myVtkInteractorStyleImage3D::OnChar() {}
void myVtkInteractorStyleImage3D::OnMouseMove(){
	vtkInteractorStyleJoystickCamera::OnMouseMove();
}
void myVtkInteractorStyleImage3D::OnKeyDown() {
	std::string key = this->GetInteractor()->GetKeySym();
	if (key.compare("Up") == 0) {
		cout << "Up arrow key was pressed." << endl;
	}
	else if (key.compare("Down") == 0) {
		cout << "Down arrow key was pressed." << endl;
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
		this->RemoveLeaks();
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
		this->MakeAnnotation(FOREGROUND);
	}
	else if (key.compare("b") == 0) {
		this->MakeAnnotation(BACKGROUND);
	}
	else if (key.compare("v") == 0) {
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
