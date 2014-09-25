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
#include <vtkImageBlend.h>
#include <vtkProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCellData.h>
#include <vtkObjectFactory.h>
#include <vtkRendererCollection.h>
#include <vtkExtractVOI.h>
#include <vtkImageMapper3D.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLineSource.h>
#include <vtkPointData.h>
#include "myVtkInteractorStyleImage.h"
#include <vtkEventForwarderCommand.h>
#include "Leap.h"
#include "SlicerLeapListener.h"
#include <sstream>
#include <iostream>
#include "constants.h"
#include <vtkLookupTable.h>
#include <vtkImageActor.h>
#include <vtkImageMapToColors.h>

// helper class to format slice status message
class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice) {
		std::stringstream tmp;
		tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
		return tmp.str();
	}
};

vtkSmartPointer<vtkActor> createCrosshair(double *size){
	// points for crosshair
	vtkSmartPointer<vtkPoints> pts =
		vtkSmartPointer<vtkPoints>::New();
	pts->InsertNextPoint(size[0], (size[2]+size[3])/2, size[5]);
	pts->InsertNextPoint(size[1], (size[2]+size[3])/2, size[5]);
	pts->InsertNextPoint((size[0]+size[1])/2, size[2], size[5]);
	pts->InsertNextPoint((size[0]+size[1])/2, size[3], size[5]);
	// Setup the colors array for crosshair
	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	// Add the colors we created to the colors array
	colors->InsertNextValue(0);
	colors->InsertNextValue(0);
	colors->InsertNextValue(255);

	colors->InsertNextValue(0);
	colors->InsertNextValue(0);
	colors->InsertNextValue(255);

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
	crosshair->Update();
	vtkSmartPointer<vtkActor> crosshairA = vtkSmartPointer<vtkActor>::New();
	crosshairA->GetProperty()->SetLineWidth(2);
	crosshairA->SetMapper(crosshair);
	return crosshairA;
}

int main(int argc, char* argv[])
{
	// Verify input arguments
	if (argc != 2)
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

	std::string outputName;
	cout << "Please enter a name for the output segmentation:" << endl;
	cin >> outputName;
	cout << "OK." << endl;

	// Create a leap listener and controller
	SlicerLeapListener listener;
	Leap::Controller controller;

	//selection layer
	vtkSmartPointer<vtkStructuredPoints> selection = 
		vtkSmartPointer<vtkStructuredPoints>::New();
	//reader->GetOutput()->CopyStructure(selection);
	selection->SetExtent(reader->GetOutput()->GetExtent());
	selection->SetSpacing(reader->GetOutput()->GetSpacing());
	selection->SetOrigin(reader->GetOutput()->GetOrigin());

	// Setup the colors array for crosshair
	vtkSmartPointer<vtkIntArray> selection_colors =
		vtkSmartPointer<vtkIntArray>::New();
	//selection_colors->SetNumberOfComponents(1);//4
	selection->AllocateScalars(VTK_INT, 1);
	cout << selection->GetNumberOfCells() << endl;
	// Add the colors we created to the colors array	
	selection_colors->SetNumberOfTuples(selection->GetNumberOfPoints());
	selection_colors->SetName("selection_colors");
	selection->GetPointData()->SetScalars(selection_colors);

	//crosshair
	// get window's size for crosshair creation and update
	vtkSmartPointer<vtkActor> crosshair = createCrosshair(reader->GetOutput()->GetBounds());
	crosshair->GetMapper()->Update();
	crosshair->SetPickable(false);

	// Visualize

	//viewer
	vtkSmartPointer<vtkImageViewer2>viewer =
		vtkSmartPointer<vtkImageViewer2>::New();

	//interaction
	vtkSmartPointer<vtkRenderWindowInteractor>interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();

	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage>::New();
	vtkSmartPointer<vtkImageMapToColors> mapperSel = vtkSmartPointer<vtkImageMapToColors>::New();


	//intiaite the points to be non-active
	vtkIntArray* point = (vtkIntArray*)selection->GetPointData()->GetScalars();
	for (int i = 0; i < selection->GetNumberOfPoints(); i++){
		point->SetValue(i, NOT_ACTIVE);
	}
	/*for (int i = 0; i < 20000; i++){
		point->SetValue(i, ACTIVE);
	}*/

	//Dimensions of selection can be overwritten by the CT dimensions.
	int* inputDims = selection->GetDimensions();
	std::cout << "Dims: " << " x: " << inputDims[0]
	<< " y: " << inputDims[1]
	<< " z: " << inputDims[2] << std::endl;
	std::cout << "Number of points: " << selection->GetNumberOfPoints() << std::endl;
	std::cout << "Number of cells: " << selection->GetNumberOfCells() << std::endl;

	//getting the extent, update it and set it to the extractedVOI
	int* extent = reader->GetOutput()->GetExtent();

	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetNumberOfTableValues(2);
	lut->SetRange(0,1);
	lut->SetTableValue((vtkIdType)NOT_ACTIVE, 0, 1, 0, 0.0);
	lut->SetTableValue((vtkIdType)ACTIVE, 1, 0, 0, 0.8);
	lut->Build();
	mapperSel->SetLookupTable(lut);
	mapperSel->SetInputData(selection);
	mapperSel->Update();

	//First of all set the input for the viewer!
	viewer->SetInputConnection(reader->GetOutputPort());
	// make imageviewer2 and sliceTextMapper visible to our interactorstyle
	// to enable slice status message updates when scrolling through the slices

	vtkSmartPointer<vtkImageActor> selectionA = vtkSmartPointer<vtkImageActor>::New();
	selectionA->GetMapper()->SetInputConnection(mapperSel->GetOutputPort());
	selectionA->InterpolateOff();

	myInteractorStyle->SetImageViewer(viewer, outputName, selectionA);
	viewer->SetupInteractor(renderWindowInteractor);
	//mapperSel->SetScalarModeToUseCellData();


	// make the interactor use our own interactorstyle
	// cause SetupInteractor() is defining it's own default interatorstyle 
	// this must be called after SetupInteractor()
	renderWindowInteractor->Initialize();
	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	renderWindowInteractor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);
	viewer->GetRenderer()->AddActor(selectionA);
	viewer->GetRenderer()->AddActor(crosshair);
	int displayExtent[6];
	viewer->GetImageActor()->GetDisplayExtent(displayExtent);
	selectionA->SetDisplayExtent(displayExtent);
	viewer->GetRenderWindow()->SetSize(1280, 720);
	// Have the sample listener receive events from the controller
	cout << "passed all creations " << endl;
	controller.addListener(listener);
	listener.SetImageViewer(myInteractorStyle);
	cout << "set viewer on listener" <<endl;
	viewer->Render();
	viewer->GetRenderer()->ResetCamera();
	cout << "Camera reset" <<endl;
	viewer->Render();
	renderWindowInteractor->Start();
	controller.removeListener(listener);

	return EXIT_SUCCESS;
}
