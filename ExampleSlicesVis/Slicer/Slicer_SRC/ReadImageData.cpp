//Must defines before ant VTK includes if not using CMake!
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

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
#include "myVtkInteractorStyleImage.h"
#include <vtkEventForwarderCommand.h>
#include "Leap.h"
#include "SlicerLeapListener.h"
#include <sstream>
#include <iostream>


#define UPDATE_SLICE_TIMER 15
// helper class to format slice status message
class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice) {
		std::stringstream tmp;
		tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
		return tmp.str();
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
	SlicerLeapListener listener;
	Leap::Controller controller;

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
	viewer->SetupInteractor(renderWindowInteractor);
	// make the interactor use our own interactorstyle
	// cause SetupInteractor() is defining it's own default interatorstyle 
	// this must be called after SetupInteractor()
	renderWindowInteractor->Initialize();
	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	renderWindowInteractor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);
	// Have the sample listener receive events from the controller
	controller.addListener(listener);
	listener.SetImageViewer(myInteractorStyle);
	viewer->Render();
	viewer->GetRenderer()->ResetCamera();
	viewer->Render();
	renderWindowInteractor->Start();
	controller.removeListener(listener);

	return EXIT_SUCCESS;
}