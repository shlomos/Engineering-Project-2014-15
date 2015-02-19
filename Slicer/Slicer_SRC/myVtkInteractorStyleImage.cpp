#include "myVtkInteractorStyleImage.h"
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
#include <vtkConnectivityFilter.h>
#include <vtkPropPicker.h>
#include <vtkImageEuclideanDistance.h>
#include <algorithm>
#include <vtkCamera.h>
#include <sstream>
#include <vtkDelaunay3D.h>


myVtkInteractorStyleImage::myVtkInteractorStyleImage()
{
	leapCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	leapCallback->SetCallback(myVtkInteractorStyleImage::ProcessLeapEvents);
	//_graph_cut = new ImageGraphCut();
	_selection_scalars = vtkSmartPointer<vtkIntArray>::New();
	this->GrabFocus(this->EventCallbackCommand);
}
void myVtkInteractorStyleImage::SetImageViewer(vtkImageViewer2* imageViewer, std::string outputName, vtkSmartPointer<vtkImageActor> selection_actor, vtkStructuredPoints* CT_image) {
	_ImageViewer = imageViewer;
	_ImageViewer->SetSliceOrientation(SLICE_ORIENTATION_XY);
	_orientation = -SLICE_ORIENTATION_XY;
	_lal = LeapAbstractionLayer::getInstance();
	_lal->setMaxSlice(imageViewer->GetSliceMax());
	_drawSize = DEFAULT_DRAW_SIZE;
	leapCallback->SetClientData(this);
	//_Slice = _MinSlice;
	_outputName = outputName;
	_hfMode = false;
	//_isSliceLocked = false;
	//_isPainting = false;
	_selection_actor = selection_actor;
	_CT_image = CT_image;
	//this->GrabFocus(this->EventCallbackCommand);
	_grow_cut = new GrowCut();
	//cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << ", Orientation: " << _orientation << std::endl;
	_selection_scalars->SetNumberOfValues(((vtkStructuredPoints*)(((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm()))->GetInput())->GetNumberOfPoints());
}
void myVtkInteractorStyleImage::SetStatusMapper(vtkTextMapper* statusMapper) {
	_StatusMapper = statusMapper;
}
//void myVtkInteractorStyleImage::setSlice(int slice){
//	if (this->Interactor->GetShiftKey() || (!_isSliceLocked && _hfMode)){
//		this->_Slice = slice;
//	}
//}
//int myVtkInteractorStyleImage::getMaxSlice(){
//	return this->_MaxSlice;
//}
void myVtkInteractorStyleImage::MoveSliceForward() {
	if (_lal->getSlice() < _lal->getMaxSlice()) {
		_lal->setSlice(_lal->getSlice() + 1);
		cout << "MoveSliceForward::Slice = " << _lal->getSlice() << std::endl;
		_ImageViewer->SetSlice(_lal->getSlice());
		int displayExtent[6];
		_ImageViewer->GetImageActor()->GetDisplayExtent(displayExtent);
		_selection_actor->SetDisplayExtent(displayExtent);
		std::string msg = StatusMessage::Format(_lal->getSlice(), _lal->getMaxSlice(), _ImageViewer->GetSliceOrientation());
		_StatusMapper->SetInput(msg.c_str());
		_ImageViewer->Render();
	}
}
void myVtkInteractorStyleImage::MoveSliceBackward() {
	if (_lal->getSlice() > 0) {
		_lal->setSlice(_lal->getSlice() - 1);
		cout << "MoveSliceBackward::Slice = " << _lal->getSlice() << std::endl;
		_ImageViewer->SetSlice(_lal->getSlice());
		int displayExtent[6];
		_ImageViewer->GetImageActor()->GetDisplayExtent(displayExtent);
		_selection_actor->SetDisplayExtent(displayExtent);
		std::string msg = StatusMessage::Format(_lal->getSlice(), _lal->getMaxSlice(), _ImageViewer->GetSliceOrientation());
		_StatusMapper->SetInput(msg.c_str());
		_ImageViewer->Render();
	}
}
void myVtkInteractorStyleImage::ToggleOrientation() {
	vtkActor* cross_actor = _ImageViewer->GetRenderer()->GetActors()->GetLastActor();
	_orientation = _orientation + 1;
	if (_orientation == 1){
		_orientation = -2;
	}
	int oooyea = std::abs(_orientation);
	_ImageViewer->SetSliceOrientation(oooyea);
	redrawCrossHair();
	int displayExtent[6];
	_ImageViewer->GetImageActor()->GetDisplayExtent(displayExtent);
	_selection_actor->SetDisplayExtent(displayExtent);
	((vtkImageActor*)cross_actor)->SetDisplayExtent(displayExtent);
	//applyCameraFixes();
	cout << "Slice GetSliceOrientation: " << _ImageViewer->GetSliceOrientation() << endl;
	_ImageViewer->SetSlice(0);
	_ImageViewer->Render();
	std::string msg = StatusMessage::Format(_ImageViewer->GetSliceMin(), _ImageViewer->GetSliceMax(), _ImageViewer->GetSliceOrientation());
	_StatusMapper->SetInput(msg.c_str());

}

void myVtkInteractorStyleImage::LoadFromFile(){
	string filename = "ref.vtk";
	cout << "Please enter the filename of the segmentation to load:... Done. Using \"ref.vtk\".." << endl;
	//cin >> filename;
	vtkSmartPointer<vtkStructuredPointsReader> sLoader =
		vtkSmartPointer<vtkStructuredPointsReader>::New();
	sLoader->SetFileName(filename.c_str());
	cout << "Reading, please wait..." << endl;
	sLoader->Update();
	vtkImageMapToColors* mapper = (vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm();
	//mapper->SetInputData(sLoader->GetOutput());
	vtkIntArray* temp = (vtkIntArray*)(sLoader->GetOutput()->GetPointData()->GetScalars());
	//for (int i = 0; i < temp->GetSize(); i++){
	//	this->_selection_scalars->SetValue(i, temp->GetValue(i));
	//}

	this->_selection_scalars->DeepCopy(temp);

	cout << "Updating view, just a little longer..." << endl;
	this->_selection_scalars->Modified();
	this->_selection_actor->Modified();
	cout << "Loaded " << filename << " succesfully!" << endl;
}

void myVtkInteractorStyleImage::ResetAll(){
	vtkStructuredPoints* selection = ((vtkStructuredPoints*)(((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm()))->GetInput());
	vtkPointData* pointData = selection->GetPointData();
	vtkIntArray* selection_scalars = (vtkIntArray*)pointData->GetScalars();
	for (int i = 0; i < selection->GetNumberOfPoints(); i++){
		selection_scalars->SetValue(i, NOT_ACTIVE);
	}
	//todo: fix this!
	//this->_graph_cut->Clean();
	selection_scalars->Modified();
}

double* myVtkInteractorStyleImage::redrawCrossHair() {
	vtkActor* cross_actor = _ImageViewer->GetRenderer()->GetActors()->GetLastActor();
	vtkSmartPointer<vtkPolyData> pd = (vtkPolyData *)((vtkPolyDataMapper*)(cross_actor->GetMapper())->GetInputAsDataSet());
	vtkSmartPointer<vtkPoints> new_pts =
		vtkSmartPointer<vtkPoints>::New();
	double* size = ((vtkStructuredPoints*)(((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm()))->GetInput())->GetBounds();
	//std::cout << "bounds: " << " [0]: " << size[0]
	//	<< " [1]: " << size[1]
	//	<< " [2]: " << size[2]
	//	<< " [3]: " << size[3]
	//	<< " [4]: " << size[4]
	//	<< " [5]: " << size[5] << std::endl;

	double cross_x = 0.0;
	double cross_y = 0.0;
	double cross_z = 0.0;
	double* temp;
	switch (_ImageViewer->GetSliceOrientation()) {
		//TODO: need to fix the "+number" factor in every set of cross_y value!!
	case SLICE_ORIENTATION_YZ:
		cross_y = std::min(size[3], std::max(size[2] + SCALE_FACTOR*((_lal->getX() - LEAP_MIN_X) / LEAP_MAX_Y)*(size[3] - size[2]) - 20, size[2]));
		cross_z = std::min(size[5], std::max(size[4] + SCALE_FACTOR*(_lal->getY() / LEAP_MAX_Y)*(size[5]-size[4]) - 20, size[4]));
		
		new_pts->InsertNextPoint(size[1], cross_y, size[4]); // vertical line
		new_pts->InsertNextPoint(size[1], cross_y, size[5]); // vertical line
		new_pts->InsertNextPoint(size[1], size[2], cross_z); // horizontal line
		new_pts->InsertNextPoint(size[1], size[3], cross_z); // horizontal line
		pd->SetPoints(new_pts);
		temp = new double[3]{ cross_z, cross_y, size[1] };
		break;
	case SLICE_ORIENTATION_XZ:
		
		cross_z = std::min(size[5], std::max(size[4] + SCALE_FACTOR*(_lal->getY() / LEAP_MAX_Y)*(abs(size[5] - size[4])) - 20, size[4]));
		cross_x = std::min(size[1], std::max(size[0] - SCALE_FACTOR*((_lal->getX() - LEAP_MIN_X) / LEAP_MAX_Y)*(abs(size[1]) - abs(size[0])) - 20, size[0]));

		new_pts->InsertNextPoint(size[0], size[2], cross_z); // vertical line
		new_pts->InsertNextPoint(size[1], size[2], cross_z); // vertical line
		new_pts->InsertNextPoint(cross_x, size[2], size[4]); // horizontal line
		new_pts->InsertNextPoint(cross_x, size[2], size[5]); // horizontal line
		pd->SetPoints(new_pts);
		temp = new double[3]{ cross_x, cross_z, size[2] };
		break;
	case SLICE_ORIENTATION_XY:
		cross_x = std::max(size[0], std::min(SCALE_FACTOR*_lal->getX(), size[1]));
		cross_y = std::min(size[3], std::max(size[2] + SCALE_FACTOR*((_lal->getY()/ LEAP_MAX_Y)*(size[3] - size[2])) - 20, size[2]));
		
		new_pts->InsertNextPoint(size[0], cross_y, size[5]);
		new_pts->InsertNextPoint(size[1], cross_y, size[5]);
		new_pts->InsertNextPoint(cross_x, size[2], size[5]);
		new_pts->InsertNextPoint(cross_x, size[3], size[5]);
		pd->SetPoints(new_pts);
		temp = new double[3]{ cross_x, cross_y, size[5] };
		break;
	}
	cross_actor->SetOrigin(_selection_actor->GetOrigin());
	cross_actor->SetScale(_selection_actor->GetScale());
	cross_actor->SetPosition(_selection_actor->GetPosition());
	int displayExtent[6];
	_ImageViewer->GetImageActor()->GetDisplayExtent(displayExtent);
	((vtkImageActor*)cross_actor)->SetDisplayExtent(displayExtent);
	_selection_actor->SetDisplayExtent(displayExtent);
	return temp;
}

void myVtkInteractorStyleImage::WriteToFile() {
	cout << "Writing segmentation to file..." << endl;
	vtkSmartPointer<vtkStructuredPointsWriter> writer =
		vtkSmartPointer<vtkStructuredPointsWriter>::New();
	writer->SetInputData((vtkStructuredPoints*)((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->GetInput());
	writer->SetFileName(_outputName.c_str());
	writer->Write();
	cout << "Done." << endl;
}

void myVtkInteractorStyleImage::translateToStructuredPoints(vtkUnstructuredGrid* component, vtkStructuredPoints* temp){
	vtkPoints* compoPoints = component->GetPoints();

	vtkPointData* pointData = temp->GetPointData();
	vtkIntArray* selection_scalars = (vtkIntArray*)pointData->GetScalars();
	vtkPointData* cPointData = component->GetPointData();
	vtkIntArray* c_selection_scalars = (vtkIntArray*)cPointData->GetScalars();

	for (int i = 0; i < compoPoints->GetNumberOfPoints(); i++){
		if (c_selection_scalars->GetValue(i) != c_selection_scalars->GetValue(max(0, i - 1))){
			cout << c_selection_scalars->GetValue(i) << endl;
		}
		selection_scalars->SetValue(i, c_selection_scalars->GetValue(i));
	}
	selection_scalars->Modified();
}

typedef void(myVtkInteractorStyleImage::*workerFunction)();

void myVtkInteractorStyleImage::doSegment() {
	boost::mutex::scoped_lock scoped_lock(_canSegment_mutex);

	vtkSmartPointer<vtkStructuredPoints> selection_structured_points = (vtkStructuredPoints*)((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->GetInput();
	
	vtkPointData* cellData1 = selection_structured_points->GetPointData();
	vtkSmartPointer<vtkIntArray> sel_scalars_copy = vtkSmartPointer<vtkIntArray>::New();
	sel_scalars_copy->DeepCopy(this->_selection_scalars);
	cellData1->SetScalars(sel_scalars_copy);

	cellData1->GetScalars()->Modified();
	cellData1->Modified();
	selection_structured_points->Modified();
	
	//vtkSmartPointer<vtkStructuredPoints> selection_structured_points_copy = vtkSmartPointer<vtkStructuredPoints>::New();
	//vtkStructuredPoints* output = vtkSmartPointer<vtkStructuredPoints>::New();
	
	// perform deep copy, in order to seperate between the segmentation given by the user, and the one computed by the machine
	//selection_structured_points_copy->DeepCopy(selection_structured_points);
	
	//_graph_cut->SetImage(selection_structured_points_copy, _CT_image);
	//output = _graph_cut->PerformSegmentation();

	_grow_cut->SetImage(_CT_image, selection_structured_points);
	_grow_cut->PerformSegmentation();
	
	//selection_structured_points = output;
	//vtkSmartPointer<vtkImageOpenClose3D> openClose =
	//	vtkSmartPointer<vtkImageOpenClose3D>::New();
	//openClose->SetInputData(selection_structured_points);
	//openClose->SetOpenValue(NOT_ACTIVE);
	//openClose->SetCloseValue(FOREGROUND);
	//openClose->SetKernelSize(XY_OPENCLOSE, XY_OPENCLOSE, Z_OPENCLOSE);
	//openClose->ReleaseDataFlagOff();
	//openClose->Update();
	
	//((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->SetInputData(selection_structured_points);
	((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->SetInputData(selection_structured_points);
}

void myVtkInteractorStyleImage::marchingCubes() {
	cout << "in marchingCubes() " << endl;
	vtkStructuredPoints* selection_structured_points = (vtkStructuredPoints*)((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->GetInput();
	_marching_cubes = new MarchingCubes(selection_structured_points);
	cout << "End of marchingCubes ctor" << endl;

}

void myVtkInteractorStyleImage::OnKeyUp() {}
void myVtkInteractorStyleImage::OnMouseMove(){
	_lal->setX(this->Interactor->GetEventPosition()[0] - 700);
	_lal->setY(this->Interactor->GetEventPosition()[1]);
}
void myVtkInteractorStyleImage::OnKeyDown() {
	std::string key = this->GetInteractor()->GetKeySym();
	if (key.compare("Up") == 0) {
		cout << "Up arrow key was pressed." << endl;
		MoveSliceForward();
	}
	else if (key.compare("Down") == 0) {
		cout << "Down arrow key was pressed." << endl;
		MoveSliceBackward();
	}
	else if (key.compare("1") == 0) {
		cout << "Draw size was changed to " << _drawSize - 1 << endl;
		this->_drawSize = std::max(MIN_DRAW_SIZE, _drawSize - 1);
	}
	else if (key.compare("2") == 0) {
		cout << "Draw size was changed to " << _drawSize + 1 << endl;
		this->_drawSize = std::min(MAX_DRAW_SIZE, _drawSize + 1);
	}
	else if (key.compare("o") == 0) {
		cout << "Orientation key was pressed." << endl;
		ToggleOrientation();
	}
	else if (key.compare("h") == 0) {
		_hfMode = !_hfMode;
		cout << "Hands free mode: " << (_hfMode ? "ON" : "OFF") << endl;
	}
	else if (key.compare("s") == 0) {
		WriteToFile();
	}
	else if (key.compare("a") == 0) {
		cout << "a was pressed!!" << endl;
	}
	else if (key.compare("r") == 0) {
		this->ResetAll();
	}
	else if (key.compare("l") == 0) {
		this->LoadFromFile();
	}
	else if (key.compare("space") == 0) {
		
		cout << "change _canSegment To FALSE" << endl;
		workerFunction f = &myVtkInteractorStyleImage::doSegment;
		boost::thread workerThread(f, this);

		// update mapper to show segmentation
		this->_selection_actor->GetMapper()->Update();
		this->_ImageViewer->Render();
		cout << "change _canSegment To TRUE" << endl;
		
	}
	else if (key.compare("m") == 0) {
		this->marchingCubes();
	}
	// forward event
}

void myVtkInteractorStyleImage::OnMouseWheelForward() {
	//std::cout << "Scrolled mouse wheel forward." << std::endl;
	MoveSliceForward();
	// don't forward events, otherwise the image will be zoomed
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelForward();
}
void myVtkInteractorStyleImage::OnMouseWheelBackward() {
	//std::cout << "Scrolled mouse wheel backward." << std::endl;
	if (_lal->getSlice() > 0) {
		MoveSliceBackward();
	}
	// don't forward events, otherwise the image will be zoomed
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelBackward();
}
void myVtkInteractorStyleImage::SetInteractor(vtkRenderWindowInteractor *i)
{
	if (i == this->Interactor)
	{
		return;
	}
	// if we already have an Interactor then stop observing it
	if (this->Interactor)
	{
		this->Interactor->RemoveObserver(this->EventCallbackCommand);
	}
	this->Interactor = i;
	// add observers for each of the events handled in ProcessEvents
	if (i)
	{
		i->AddObserver(vtkCommand::EnterEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::LeaveEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::MouseMoveEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::LeftButtonPressEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::LeftButtonReleaseEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::MiddleButtonPressEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::MiddleButtonReleaseEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::RightButtonPressEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::RightButtonReleaseEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::MouseWheelForwardEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::MouseWheelBackwardEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::ExposeEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::ConfigureEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::TimerEvent,
			this->leapCallback,
			this->Priority);
		i->AddObserver(vtkCommand::KeyPressEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::KeyReleaseEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::CharEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::DeleteEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::TDxMotionEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::TDxButtonPressEvent,
			this->EventCallbackCommand,
			this->Priority);
		i->AddObserver(vtkCommand::TDxButtonReleaseEvent,
			this->EventCallbackCommand,
			this->Priority);
	}
	this->EventForwarder->SetTarget(this->Interactor);
	if (this->Interactor)
	{
		this->AddObserver(vtkCommand::StartInteractionEvent, this->EventForwarder);
		this->AddObserver(vtkCommand::EndInteractionEvent, this->EventForwarder);
	}
	else
	{
		this->RemoveObserver(this->EventForwarder);
	}
}
void myVtkInteractorStyleImage::ProcessLeapEvents(vtkObject* object, unsigned long event, void* clientdata, void* callData){
	//cout << "Got a leap event!"<< endl;
	vtkSmartPointer<myVtkInteractorStyleImage> intStyle =
		reinterpret_cast<myVtkInteractorStyleImage*>(clientdata);
	vtkActor* cross_actor = intStyle->_ImageViewer->GetRenderer()->GetActors()->GetLastActor();
	vtkSmartPointer<vtkImageMapToColors> selection_mapper = (vtkImageMapToColors*)intStyle->_selection_actor->GetMapper()->GetInputAlgorithm();
	vtkSmartPointer<vtkStructuredPoints> selection_structured_points = (vtkStructuredPoints *)selection_mapper->GetInput();

	double* size = selection_structured_points->GetBounds();
	// changing the crosshair
	vtkSmartPointer<vtkPoints> new_pts =
		vtkSmartPointer<vtkPoints>::New();
	double* temp = intStyle->redrawCrossHair();
	double cross_1 = temp[0];
	double cross_2 = temp[1];
	double cross_3 = temp[2];
	delete temp;
	// When SHIFT key is pressed, udpate slice
	if (intStyle->Interactor->GetShiftKey() || (intStyle->_hfMode && !intStyle->_lal->getSliceLock())){
		cout << "Shift pressed" << endl;
		intStyle->_ImageViewer->SetSlice(intStyle->_lal->getSlice());
		int displayExtent[6];
		intStyle->_ImageViewer->GetImageActor()->GetDisplayExtent(displayExtent);
		intStyle->_selection_actor->SetDisplayExtent(displayExtent);
		std::string msg = StatusMessage::Format(intStyle->_lal->getSlice(), intStyle->_lal->getMaxSlice(), intStyle->_ImageViewer->GetSliceOrientation());
		intStyle->_StatusMapper->SetInput(msg.c_str());
	}

	int *selExt = selection_structured_points->GetExtent();
	

	if (intStyle->Interactor->GetControlKey() || (intStyle->_hfMode && intStyle->_lal->getPainting())) {

		vtkPointData* cellData = selection_structured_points->GetPointData();
		vtkSmartPointer<vtkIntArray> selection_scalars = vtkSmartPointer<vtkIntArray>::New();
		selection_scalars = intStyle->_selection_scalars;
		//selection_scalars->SetNumberOfValues(selection_structured_points->GetNumberOfPoints());

		double x[3];
		int ijk[3];
		double pCoord[3];
		int ijk2[3];
		int minX, maxX, minY, maxY;
		switch (intStyle->_ImageViewer->GetSliceOrientation()) {
		case SLICE_ORIENTATION_YZ:
			cout << "case 0" << endl;
			x[0] = cross_3;
			x[1] = cross_2;
			x[2] = cross_1;
			selection_structured_points->ComputeStructuredCoordinates(x, ijk, pCoord);
			ijk[0] = intStyle->_ImageViewer->GetSlice();//_lal->getSlice();
			ijk2[1] = 0;
			ijk2[2] = 0;
			ijk2[0] = ijk[0];
			minX = std::max(0, ijk[2] - intStyle->_drawSize);
			maxX = std::min(ijk[2] + intStyle->_drawSize, selExt[5]);
			minY = std::max(ijk[1] - intStyle->_drawSize, 0);
			maxY = std::min(ijk[1] + intStyle->_drawSize, selExt[3]);
			for (int i = minX; i < maxX; i++){
				ijk2[2] = i;
				for (int j = minY; j < maxY; j++){
					ijk2[1] = j;
					vtkIdType cellId = selection_structured_points->ComputePointId(ijk2);
					//cout << ijk2[0] << ":" << ijk2[1] << ":" << ijk2[2] << endl;
					selection_scalars->SetValue(cellId, FOREGROUND);
				}
			}
			break;
		case SLICE_ORIENTATION_XZ:
			cout << "case 1" << endl;
			x[0] = cross_1;
			x[1] = cross_3;
			x[2] = cross_2;
			selection_structured_points->ComputeStructuredCoordinates(x, ijk, pCoord);
			ijk[1] = intStyle->_ImageViewer->GetSlice();
			ijk2[0] = 0;
			ijk2[1] = ijk[1];
			ijk2[2] = 0;
			minX = std::max(0, ijk[0] - intStyle->_drawSize);
			maxX = std::min(ijk[0] + intStyle->_drawSize, selExt[1]);
			minY = std::max(ijk[2] - intStyle->_drawSize, 0);
			maxY = std::min(ijk[2] + intStyle->_drawSize, selExt[5]);
			for (int i = minX; i < maxX; i++){
				ijk2[0] = i;
				for (int j = minY; j < maxY; j++){
					ijk2[2] = j;
					vtkIdType cellId = selection_structured_points->ComputePointId(ijk2);
					//cout << ijk2[0] << ":" << ijk2[1] << ":" << ijk2[2] << endl;
					selection_scalars->SetValue(cellId, FOREGROUND);
				}
			}
			break;
		case SLICE_ORIENTATION_XY:
			x[0] = cross_1;
			x[1] = cross_2;
			x[2] = cross_3;
			selection_structured_points->ComputeStructuredCoordinates(x, ijk, pCoord);
			ijk[2] = intStyle->_ImageViewer->GetSlice();
			ijk2[0] = 0;
			ijk2[1] = 0;
			ijk2[2] = ijk[2];
			minX = std::max(0, ijk[0] - intStyle->_drawSize);
			maxX = std::min(ijk[0] + intStyle->_drawSize, selExt[1]);
			minY = std::max(ijk[1] - intStyle->_drawSize, 0);
			maxY = std::min(ijk[1] + intStyle->_drawSize, selExt[3]);
			for (int i = minX; i < maxX; i++){
				ijk2[0] = i;
				for (int j = minY; j < maxY; j++){
					ijk2[1] = j;
					vtkIdType cellId = selection_structured_points->ComputePointId(ijk2);
					//cout << ijk2[0] << ":" << ijk2[1] << ":" << ijk2[2] << endl;
					selection_scalars->SetValue(cellId, FOREGROUND);
				}
			}
			break;
		}

		//Update the underlying data object.
		cellData->SetScalars(selection_scalars);

		//Update the underlying data object.
		cellData->GetScalars()->Modified();
		
		////todo: remove the following lines if possible!
		cellData->Modified();
		selection_structured_points->Modified();
		selection_mapper->Modified();

		// change the crosshair's color
		vtkSmartPointer<vtkPolyData> pd = (vtkPolyData *)((vtkPolyDataMapper*)(cross_actor->GetMapper())->GetInputAsDataSet());
		vtkUnsignedCharArray* da = (vtkUnsignedCharArray*)pd->GetCellData()->GetScalars();
		da->SetValue(0, 255);
		da->SetValue(1, 0);
		da->SetValue(2, 0);
		da->SetValue(3, 255);
		da->SetValue(4, 0);
		da->SetValue(5, 0);
	}
	else if (intStyle->Interactor->GetAltKey()){

		vtkPointData* cellData = selection_structured_points->GetPointData();
		vtkSmartPointer<vtkIntArray> selection_scalars = vtkSmartPointer<vtkIntArray>::New();
		selection_scalars = intStyle->_selection_scalars; 
		//selection_scalars->SetNumberOfValues(selection_structured_points->GetNumberOfPoints());
		double x[3] = { cross_1, cross_2, cross_3 };
		int ijk[3];
		double pCoord[3];
		selection_structured_points->ComputeStructuredCoordinates(x, ijk, pCoord);
		ijk[2] = intStyle->_lal->getSlice();
		int ijk2[3] = { 0, 0, ijk[2] };
		int minX = std::max(0, ijk[0] - intStyle->_drawSize);
		int maxX = std::min(ijk[0] + intStyle->_drawSize, selExt[1]);
		int minY = std::max(ijk[1] - intStyle->_drawSize, 0);
		int maxY = std::min(ijk[1] + intStyle->_drawSize, selExt[3]);
		for (int i = minX; i < maxX; i++){
			ijk2[0] = i;
			for (int j = minY; j < maxY; j++){
				ijk2[1] = j;
				vtkIdType cellId = selection_structured_points->ComputePointId(ijk2);
				selection_scalars->SetValue(cellId, BACKGROUND);
			}
		}

		//Update the underlying data object.
		cellData->SetScalars(selection_scalars);

		//Update the underlying data object.
		cellData->GetScalars()->Modified();
		
		////todo: remove the following lines if possible!
		cellData->Modified();
		selection_structured_points->Modified();
		selection_mapper->Modified();

		// change the crosshair's color
		vtkSmartPointer<vtkPolyData> pd = (vtkPolyData *)((vtkPolyDataMapper*)(cross_actor->GetMapper())->GetInputAsDataSet());
		vtkUnsignedCharArray* da = (vtkUnsignedCharArray*)pd->GetCellData()->GetScalars();
		da->SetValue(0, 0);
		da->SetValue(1, 255);
		da->SetValue(2, 0);
		da->SetValue(3, 0);
		da->SetValue(4, 255);
		da->SetValue(5, 0);
	}
	else {
		// change the crosshair's color
		vtkSmartPointer<vtkPolyData> pd = (vtkPolyData *)((vtkPolyDataMapper*)(cross_actor->GetMapper())->GetInputAsDataSet());
		vtkUnsignedCharArray* da2 = (vtkUnsignedCharArray*)pd->GetCellData()->GetScalars();
		da2->SetValue(0, 0);
		da2->SetValue(1, 0);
		da2->SetValue(2, 255);
		da2->SetValue(3, 0);
		da2->SetValue(4, 0);
		da2->SetValue(5, 255);
	}

	// render
	cross_actor->GetMapper()->Update();
	intStyle->_selection_actor->GetMapper()->Update();
	intStyle->_ImageViewer->Render();
	
}

myVtkInteractorStyleImage::~myVtkInteractorStyleImage(){
	//delete _graph_cut;
	delete _grow_cut;
}

vtkStandardNewMacro(myVtkInteractorStyleImage);
