//Must defines before ant VTK includes if not using CMake!
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkLine.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCellData.h>
#include <vtkObjectFactory.h>
#include <vtkCellArray.h>
#include <vtkTextMapper.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLineSource.h>
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

vtkPolyDataMapper* createCrosshair(int size){
	vtkSmartPointer<vtkPoints> pts = 
		vtkSmartPointer<vtkPoints>::New();
	pts->InsertNextPoint(-size, 0, 0);
	pts->InsertNextPoint(size, 0, 0);
	pts->InsertNextPoint(0, -size, 0);
	pts->InsertNextPoint(0, size, 0);

	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = 
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	// Add the colors we created to the colors array
	colors->InsertNextValue(255);
	colors->InsertNextValue(128);
	colors->InsertNextValue(0);

	colors->InsertNextValue(255);
	colors->InsertNextValue(128);
	colors->InsertNextValue(0);

	// Create the first line
	vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
	line0->GetPointIds()->SetId(0, 0);
	line0->GetPointIds()->SetId(1, 1);

	// Create the second line
	vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
	line1->GetPointIds()->SetId(0, 2);
	line1->GetPointIds()->SetId(1, 3);

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(line0);
	lines->InsertNextCell(line1);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	// Add the points to the dataset
	linesPolyData->SetPoints(pts);
	// Add the lines to the dataset
	linesPolyData->SetLines(lines);
	// Color the lines
	linesPolyData->GetCellData()->SetScalars(colors);
	vtkSmartPointer<vtkPolyDataMapper> crosshair = vtkSmartPointer<vtkPolyDataMapper>::New();
	crosshair->SetInputData(linesPolyData);
}

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

	// Create a leap listener and controller
	SlicerLeapListener listener;
	Leap::Controller controller;

	//crosshair
	double size = (reader->GetOutput()->GetExtent()[0]+reader->GetOutput()->GetExtent()[3]);
	vtkSmartPointer<vtkActor> crosshair = vtkSmartPointer<vtkActor>::New();
	crosshair->GetProperty()->SetLineWidth(2);
	crosshair->SetPosition(0,0,0);
	crosshair->SetMapper(createCrosshair(size));

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
	vtkSmartPointer<vtkRenderer> overlay = vtkSmartPointer<vtkRenderer>::New();
	viewer->GetRenderWindow()->SetNumberOfLayers(2);
	overlay->SetLayer(1);
	overlay->AddActor(crosshair);
	viewer->GetRenderWindow()->AddRenderer(overlay);
	viewer->GetRenderWindow()->SetSize(1280,720);
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
