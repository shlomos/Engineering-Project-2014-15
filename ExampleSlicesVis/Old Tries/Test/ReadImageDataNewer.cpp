#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCallbackCommand.h>
#include <vtkObjectFactory.h>
#include <vtkTextMapper.h>
#include <vtkInteractorStyleImage.h>
#include <sstream>
#include <iostream>
#include "Leap.h"


#define MY_LEAP_EVENT 864849374
// helper class to format slice status message
class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice) {
		std::stringstream tmp;
		tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
		return tmp.str();
	}
};

/******************************************************************************\
* Copyright (C) 2012-2013 Leap Motion, Inc. All rights reserved.               *
* Leap Motion proprietary and confidential. Not for distribution.              *
* Use subject to the terms of the Leap Motion SDK Agreement available at       *
* https://developer.leapmotion.com/sdk_agreement, or another agreement         *
* between Leap Motion and you, your company or other organization.             *
\******************************************************************************/

#include <iostream>
#include "Leap.h"
using namespace Leap;



// Define own interaction style
class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
	static myVtkInteractorStyleImage* New();
	vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);

protected:
	vtkImageViewer2* _ImageViewer;
	vtkTextMapper* _StatusMapper;
	int _Slice;
	int _orientation;
	int _MinSlice;


public:
	int _MaxSlice;

	void ProcessLeapEvents(vtkObject* object, unsigned long event, void* clientdata, void* callData){
		cout << "Got a leap event! yay!" << endl;
		_ImageViewer->SetSlice(121);
		_ImageViewer->Render();
	}

	void SetImageViewer(vtkImageViewer2* imageViewer) {
		_ImageViewer = imageViewer;
		_orientation=imageViewer->GetSliceOrientation();
		_MinSlice = imageViewer->GetSliceMin();
		_MaxSlice = imageViewer->GetSliceMax();
		_Slice = _MinSlice;
		cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << ", Orientation: " << _orientation << std::endl;
	}

	void SetStatusMapper(vtkTextMapper* statusMapper) {
		_StatusMapper = statusMapper;
	}

protected:
	void MoveSliceForward() {
		if(_Slice < _MaxSlice) {
			_Slice += 1;
			cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
			_ImageViewer->SetSlice(_Slice);
			//std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			//_StatusMapper->SetInput(msg.c_str());
			_ImageViewer->Render();
		}
	}

	void MoveSliceBackward() {
		if(_Slice > _MinSlice) {
			_Slice -= 1;
			cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
			_ImageViewer->SetSlice(_Slice);
			//std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			//_StatusMapper->SetInput(msg.c_str());
			_ImageViewer->Render();
		}
	}

	void ToggleOrientation() {
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

	virtual void OnEnter() {
		cout<<"We are stupid!"<<endl;
		vtkInteractorStyleImage::OnEnter();
	}


	virtual void OnKeyDown() {
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
		// forward event
	}


	virtual void OnMouseWheelForward() {
		//std::cout << "Scrolled mouse wheel forward." << std::endl;
		MoveSliceForward();
		// don't forward events, otherwise the image will be zoomed 
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		// vtkInteractorStyleImage::OnMouseWheelForward();
	}


	virtual void OnMouseWheelBackward() {
		//std::cout << "Scrolled mouse wheel backward." << std::endl;
		if(_Slice > _MinSlice) {
			MoveSliceBackward();
		}
		// don't forward events, otherwise the image will be zoomed 
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		// vtkInteractorStyleImage::OnMouseWheelBackward();
	}
};

vtkStandardNewMacro(myVtkInteractorStyleImage);

class SampleListener : public Listener,vtkObject {
public:

	vtkRenderWindowInteractor* _Interactor;
	vtkSmartPointer<vtkCallbackCommand> cb1;

	void SetImageViewer(vtkRenderWindowInteractor* interactor) {
		_Interactor = interactor;
		//cb1 = vtkSmartPointer<vtkCallbackCommand>::New();
		//cb1->SetCallback(ProcessLeapEvents);
		//this->AddObserver(vtkCommand::EnterEvent,cb1);
	}


	void onInit(const Controller& controller) {
		std::cout << "Initialized" << std::endl;
	}

	void onConnect(const Controller& controller) {
		std::cout << "Connected" << std::endl;
		controller.enableGesture(Gesture::TYPE_CIRCLE);
		controller.enableGesture(Gesture::TYPE_KEY_TAP);
		controller.enableGesture(Gesture::TYPE_SCREEN_TAP);
		controller.enableGesture(Gesture::TYPE_SWIPE);

	}

	void onDisconnect(const Controller& controller) {
		//Note: not dispatched when running in a debugger.
		std::cout << "Disconnected" << std::endl;
	}

	void onExit(const Controller& controller) {
		std::cout << "Exited" << std::endl;
	}

	void onFrame(const Controller& controller) {
		// Get the most recent frame and report some basic information
		const Frame frame = controller.frame();
		if (!frame.hands().isEmpty()) {
			// Get the first hand
			const Hand hand = frame.hands()[0];

			// Check if the hand has any fingers
			const FingerList fingers = hand.fingers();
			if (!fingers.isEmpty()) {
				// Calculate the hand's average finger tip position
				Vector avgPos;
				for (int i = 0; i < fingers.count(); ++i) {
					avgPos += fingers[i].tipPosition();
				}
				avgPos /= (float)fingers.count();
				Vector *aPoint = new Vector(0.0f, 0.0f, -150.0f);
				Vector frontmost = frame.pointables().frontmost().tipPosition();
				float realDis = (frontmost.z-aPoint->z);
				cout << "im sad!"<<endl;
				int dis = std::min(130,(int)(realDis/150.0f*130));
				cout << "im very very sad!"<<endl;
				std::cout << "try set slice to " << dis<<  std::endl;
				//Assume screen is 15 cm after leap
				std::cout << "finger is at " << frontmost.z
					<< " pointer is " << dis << "cm from screen." << std::endl;
				cout << "im desperate!"<<endl;
				cout << this->InvokeEvent(vtkCommand::UserEvent)<<endl;
				cout << this->_Interactor->InvokeEvent(vtkCommand::UserEvent) <<endl;
			}
		}

		if (!frame.hands().isEmpty()) {
			//std::cout << std::endl;
		}
	}

	void onFocusGained(const Controller& controller) {
		std::cout << "Focus Gained" << std::endl;
	}

	void onFocusLost(const Controller& controller) {
		std::cout << "Focus Lost" << std::endl;
	}
};

int main(int argc, char* argv[])
{
	// Verify input arguments
	if ( argc != 2 )
	{
		std::cout << "Usage: " << argv[0]
		<< " Filename(.vtk) as structured points." << std::endl;
		return EXIT_FAILURE;
	}

	std::string inputFilename = argv[1];

	// Read the file
	vtkSmartPointer<vtkStructuredPointsReader> reader =
		vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->Update();

	// Create a sample listener and controller
	SampleListener listener;
	Controller controller;

	// Have the sample listener receive events from the controller
	controller.addListener(listener);

	// Visualize

	//viewer
	vtkSmartPointer<vtkImageViewer2>viewer=
		vtkSmartPointer<vtkImageViewer2>::New();

	//interaction
	vtkSmartPointer<vtkRenderWindowInteractor>interactor = 
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();

	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage>::New();

	//First of all set the input for the viewer!
	viewer->SetInputConnection(reader->GetOutputPort());

	// make imageviewer2 and sliceTextMapper visible to our interactorstyle
	// to enable slice status message updates when scrolling through the slices
	myInteractorStyle->SetImageViewer(viewer);
	listener.SetImageViewer(renderWindowInteractor);
	vtkSmartPointer<vtkCallbackCommand> leapCallback = 
		vtkSmartPointer<vtkCallbackCommand>::New();
	leapCallback->SetCallback(/*Some shit...*/);
	viewer->SetupInteractor(renderWindowInteractor);
	// make the interactor use our own interactorstyle
	// cause SetupInteractor() is defining it's own default interatorstyle 
	// this must be called after SetupInteractor()
	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	renderWindowInteractor->AddObserver(vtkCommand::UserEvent,leapCallback);
	viewer->Render();
	viewer->GetRenderer()->ResetCamera();
	viewer->Render();
	renderWindowInteractor->Start();

	controller.removeListener(listener);

	return EXIT_SUCCESS;
}