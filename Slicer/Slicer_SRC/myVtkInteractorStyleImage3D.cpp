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
	cout << "3D Interactor iniitiated." << std::endl;
}

void myVtkInteractorStyleImage3D::LoadFromFile(){
}

void myVtkInteractorStyleImage3D::ResetAll(){
}

void myVtkInteractorStyleImage3D::WriteToFile() {
	cout << "Writing segmentation to file..." << endl;
	/*vtkSmartPointer<vtkStructuredPointsWriter> writer =
		vtkSmartPointer<vtkStructuredPointsWriter>::New();
		writer->SetInputData((vtkStructuredPoints*)((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->GetInput());
		writer->SetFileName(_outputName.c_str());
		writer->Write();
		cout << "Done." << endl;*/
}

typedef void(myVtkInteractorStyleImage3D::*workerFunction)();

void myVtkInteractorStyleImage3D::doSegment() {
	//vtkSmartPointer<vtkImageEuclideanDistance> ed3d = vtkSmartPointer<vtkImageEuclideanDistance>::New();
	////ed3d->SetInputData((vtkImageData*)((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->GetInput());
	////ed3d->Update();
	////_graph_cut->SetImage((vtkStructuredPoints*)ed3d->GetOutput(), _CT_image);
	//vtkStructuredPoints* selection_structured_points = (vtkStructuredPoints*)((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->GetInput();
	//_graph_cut->SetImage( selection_structured_points, _CT_image );
	//_graph_cut->PerformSegmentation();
	//vtkSmartPointer<vtkImageOpenClose3D> openClose =
	//	vtkSmartPointer<vtkImageOpenClose3D>::New();
	//openClose->SetInputData(selection_structured_points);
	//openClose->SetOpenValue(NOT_ACTIVE);
	//openClose->SetCloseValue(FOREGROUND);
	//openClose->SetKernelSize(XY_OPENCLOSE, XY_OPENCLOSE, Z_OPENCLOSE);
	//openClose->ReleaseDataFlagOff();
	//openClose->Update();
	//((vtkImageMapToColors*)_selection_actor->GetMapper()->GetInputAlgorithm())->SetInputData(openClose->GetOutput());
}
void myVtkInteractorStyleImage3D::OnLeftButtonDown()
{
	std::cout << "Pressed left mouse button." << std::endl;
	// Forward events
	vtkInteractorStyleJoystickActor::OnLeftButtonDown();
}

void myVtkInteractorStyleImage3D::OnRightButtonDown()
{
	std::cout << "Pressed left mouse button." << std::endl;
	// Forward events
	vtkInteractorStyleJoystickActor::OnRightButtonDown();
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
	if (this->Interactor->GetShiftKey()){
		Interactor->SetEventPosition(_lal->getX() + winSize[0] / 2, winSize[1] / 2 - _lal->getY());
	}
	else{
		Interactor->SetEventPosition(winSize[0] / 2, winSize[1] / 2);
	}
	cout << Interactor->GetEventPosition()[0] << ":" << Interactor->GetEventPosition()[1] << endl;
	// Delegate
	vtkInteractorStyleJoystickActor::OnTimer();
	
}

void myVtkInteractorStyleImage3D::OnKeyUp() {}
void myVtkInteractorStyleImage3D::OnMouseMove(){
	vtkInteractorStyleJoystickActor::OnMouseMove();
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
	else if (key.compare("o") == 0) {
		cout << "Orientation key was pressed." << endl;
	}
	else if (key.compare("f") == 0) {
		_hfMode = !_hfMode;
		cout << "Hands free mode: " << (_hfMode ? "ON" : "OFF") << endl;
	}
	else if (key.compare("s") == 0) {
		WriteToFile();
	}
	else if (key.compare("a") == 0) {
		cout << "a pressed!" << endl;
		if (_rotLock){
			_rotLock = false;
			int x = this->Interactor->GetEventPosition()[0];
			int y = this->Interactor->GetEventPosition()[1];
			this->FindPokedRenderer(x, y);
			this->FindPickedActor(x, y);
			if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
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
	else if (key.compare("space") == 0) {

		//workerFunction f = &myVtkInteractorStyleImage::doSegment;
		//boost::thread workerThread(f, this);

		//// update mapper to show segmentation
		//this->_selection_actor->GetMapper()->Update();
		//this->_ImageViewer->Render();
	}
	else if (key.compare("m") == 0) {
		cout << "m pressed event!" << endl;
	}
	// forward event
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
