#include "myVtkInteractorStyleImage.h"


myVtkInteractorStyleImage::myVtkInteractorStyleImage()
{
	leapCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	leapCallback->SetCallback(myVtkInteractorStyleImage::ProcessLeapEvents);
}

void myVtkInteractorStyleImage::SetImageViewer(vtkImageViewer2* imageViewer) {
	_ImageViewer = imageViewer;
	_orientation=imageViewer->GetSliceOrientation();
	_MinSlice = imageViewer->GetSliceMin();
	_MaxSlice = imageViewer->GetSliceMax();
	leapCallback->SetClientData(this);
	_Slice = _MinSlice;
	_isSliceLocked = false;
	cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << ", Orientation: " << _orientation << std::endl;
}

void myVtkInteractorStyleImage::SetStatusMapper(vtkTextMapper* statusMapper) {
	_StatusMapper = statusMapper;
}

void myVtkInteractorStyleImage::setSlice(int slice){
	if(this->Interactor->GetShiftKey()){
		this->_Slice=slice;
	}
}

int myVtkInteractorStyleImage::getMaxSlice(){
	return this->_MaxSlice;
}

void myVtkInteractorStyleImage::MoveSliceForward() {
	if(_Slice < _MaxSlice) {
		_Slice += 1;
		cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
		_ImageViewer->SetSlice(_Slice);
		//std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
		//_StatusMapper->SetInput(msg.c_str());
		_ImageViewer->Render();
	}
}

void myVtkInteractorStyleImage::MoveSliceBackward() {
	if(_Slice > _MinSlice) {
		_Slice -= 1;
		cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
		_ImageViewer->SetSlice(_Slice);
		//std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
		//_StatusMapper->SetInput(msg.c_str());
		_ImageViewer->Render();
	}
}

void myVtkInteractorStyleImage::ToggleOrientation() {
	_orientation = (_orientation+1)%3;
	switch(_orientation){
	case 0:
		_ImageViewer->SetSliceOrientationToXY();
		break;
	case 1:
		_ImageViewer->SetSliceOrientationToXZ();
		break;
	case 2:
		_ImageViewer->SetSliceOrientationToYZ();
		break;
	}
	cout << "Slice orientation: " << _orientation << endl;
	_ImageViewer->SetSlice(0);
	_ImageViewer->Render();
}

void myVtkInteractorStyleImage::OnKeyDown() {
	std::string key = this->GetInteractor()->GetKeySym();
	if(key.compare("Up") == 0) {
		cout << "Up arrow key was pressed." << endl;
		MoveSliceForward();
	}
	else if(key.compare("Down") == 0) {
		cout << "Down arrow key was pressed." << endl;
		MoveSliceBackward();
	}
	else if(key.compare("o") == 0) {
		cout << "Orientation key was pressed." << endl;
		ToggleOrientation();
	}
	else if(key.compare("shift") == 0) {
		cout << "Shift key was pressed." << endl;
		lockSlice(false);
	}
	// forward event
}

void myVtkInteractorStyleImage::OnKeyUp() {
	std::string key = this->GetInteractor()->GetKeySym();
	if(key.compare("shift") == 0) {
		cout << "shift key was released." << endl;
		lockSlice(true);
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
	if(_Slice > _MinSlice) {
		MoveSliceBackward();
	}
	// don't forward events, otherwise the image will be zoomed 
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelBackward();
}

void myVtkInteractorStyleImage::lockSlice(bool state){
	_isSliceLocked = state;
}

void myVtkInteractorStyleImage::SetInteractor(vtkRenderWindowInteractor *i)
{
	if(i == this->Interactor)
	{
		return;
	}
	// if we already have an Interactor then stop observing it
	if(this->Interactor)
	{
		this->Interactor->RemoveObserver(this->EventCallbackCommand);
	}
	this->Interactor = i;
	// add observers for each of the events handled in ProcessEvents
	if(i)
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
	vtkSmartPointer<myVtkInteractorStyleImage> intStyle = 
		reinterpret_cast<myVtkInteractorStyleImage*>(clientdata);
	if(intStyle->Interactor->GetShiftKey()){
		intStyle->_ImageViewer->SetSlice(intStyle->_Slice);
		intStyle->_ImageViewer->Render();
	}
}

vtkStandardNewMacro(myVtkInteractorStyleImage);