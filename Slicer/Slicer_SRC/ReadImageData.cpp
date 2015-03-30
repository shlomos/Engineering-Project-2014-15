//Must defines before ant VTK includes if not using CMake!
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLine.h>
#include <vtkActor.h>
#include <vtkImageBlend.h>
#include <vtkProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCellData.h>
#include <vtkObjectFactory.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkExtractVOI.h>
#include <vtkImageMapper3D.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkNIFTIImageReader.h>
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
#include "CrosshairFactory.h"
#include <vtkLookupTable.h>
#include <vtkImageActor.h>
#include <vtkImageMapToColors.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkUnsignedShortArray.h>


vtkSmartPointer<vtkActor> createCrosshair(double *size){
	// points for crosshair
	vtkSmartPointer<vtkPoints> pts =
		vtkSmartPointer<vtkPoints>::New();
	pts->InsertNextPoint(size[0], (size[2] + size[3]) / 2, size[5]);
	pts->InsertNextPoint(size[1], (size[2] + size[3]) / 2, size[5]);
	pts->InsertNextPoint((size[0] + size[1]) / 2, size[2], size[5]);
	pts->InsertNextPoint((size[0] + size[1]) / 2, size[3], size[5]);
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
	vtkSmartPointer<vtkNIFTIImageReader> reader1 =
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader1->SetFileName(inputFilename.c_str());
	reader1->Update();
	vtkSmartPointer<vtkImageToStructuredPoints> reader =
		vtkSmartPointer<vtkImageToStructuredPoints>::New();
	reader->SetInputConnection(reader1->GetOutputPort());
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
	vtkSmartPointer<vtkUnsignedShortArray> selection_colors =
		vtkSmartPointer<vtkUnsignedShortArray>::New();
	//selection_colors->SetNumberOfComponents(1);//4
	selection->AllocateScalars(VTK_UNSIGNED_SHORT, 1);
	cout << selection->GetNumberOfCells() << endl;
	// Add the colors we created to the colors array	
	selection_colors->SetNumberOfTuples(selection->GetNumberOfPoints());
	selection_colors->SetName("selection_colors");
	selection->GetPointData()->SetScalars(selection_colors);
	//crosshair
	CrosshairFactory* crossFac = CrosshairFactory::getInstance();
	// get window's size for crosshair creation and update
	vtkSmartPointer<vtkActor> crosshair = crossFac->makeCrosshair(reader->GetOutput()->GetBounds());//createCrosshair(reader->GetOutput()->GetBounds());
	crosshair->GetMapper()->Update();
	crosshair->SetPickable(false);

	// Visualize

	//viewer
	vtkSmartPointer<vtkImageViewer2>viewer =
		vtkSmartPointer<vtkImageViewer2>::New();

	//interaction
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();

	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage>::New();
	vtkSmartPointer<vtkImageMapToColors> mapperSel = vtkSmartPointer<vtkImageMapToColors>::New();


	//intiaite the points to be non-active
	vtkUnsignedShortArray* point = (vtkUnsignedShortArray*)selection->GetPointData()->GetScalars();
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
	lut->SetNumberOfTableValues(3);
	lut->SetRange(0, 2);
	lut->SetTableValue((vtkIdType)NOT_ACTIVE, 0, 0, 0, 0.0);
	lut->SetTableValue((vtkIdType)FOREGROUND, 1, 0, 0, 0.5);
	lut->SetTableValue((vtkIdType)BACKGROUND, 0, 0, 1, 0.5);
	lut->Build();
	mapperSel->SetLookupTable(lut);
	mapperSel->SetInputData(selection);
	mapperSel->Update();
	vtkSmartPointer<vtkImageGaussianSmooth> smoothed = vtkSmartPointer<vtkImageGaussianSmooth>::New();
	smoothed->SetInputConnection(reader->GetOutputPort());
	smoothed->SetDimensionality(3);
	smoothed->SetRadiusFactors(SMOOTHING_FACTOR_XY, SMOOTHING_FACTOR_XY, SMOOTHING_FACTOR_Z);
	//smoothed->SetStandardDeviations(0.0,0.0,0.0);
	smoothed->Update();
	//First of all set the input for the viewer!
	viewer->SetInputConnection(reader->GetOutputPort());
	// make imageviewer2 and sliceTextMapper visible to our interactorstyle
	// to enable slice status message updates when scrolling through the slices
	vtkSmartPointer<vtkImageActor> selectionA = vtkSmartPointer<vtkImageActor>::New();
	selectionA->GetMapper()->SetInputConnection(mapperSel->GetOutputPort());
	selectionA->InterpolateOff();

	// slice status message
	vtkSmartPointer<vtkTextProperty> sliceTextProp = vtkSmartPointer<vtkTextProperty>::New();
	sliceTextProp->SetFontFamilyToCourier();
	sliceTextProp->SetFontSize(20);
	sliceTextProp->SetVerticalJustificationToBottom();
	sliceTextProp->SetJustificationToLeft();

	vtkSmartPointer<vtkTextMapper> sliceTextMapper = vtkSmartPointer<vtkTextMapper>::New();
	std::string msg = StatusMessage::Format(viewer->GetSliceMin(), viewer->GetSliceMax(), viewer->GetSliceOrientation());
	sliceTextMapper->SetInput(msg.c_str());
	sliceTextMapper->SetTextProperty(sliceTextProp);

	vtkSmartPointer<vtkActor2D> sliceTextActor = vtkSmartPointer<vtkActor2D>::New();
	sliceTextActor->SetMapper(sliceTextMapper);
	sliceTextActor->SetPosition(15, 10);

	// usage hint message
	vtkSmartPointer<vtkTextProperty> usageTextProp = vtkSmartPointer<vtkTextProperty>::New();
	usageTextProp->SetFontFamilyToCourier();
	usageTextProp->SetFontSize(14);
	usageTextProp->SetVerticalJustificationToTop();
	usageTextProp->SetJustificationToLeft();

	vtkSmartPointer<vtkTextMapper> usageTextMapper = vtkSmartPointer<vtkTextMapper>::New();
	usageTextMapper->SetInput("Options:\n - Hold Shift to scroll between slices.\n - Hold Ctrl"
		" to draw segmentation.\n - Hold Alt to mark background.\n\n -- Press '1' and '2' to change"
		" brush's size\n -- Press 'o' to toggle orientation\n -- Press 's' to save segmentation\n --"
		"  Press 'l' to load segmentation\n --  Press 'r' to reset selection\n -- Press 'h' for hands-free mode.");
	usageTextMapper->SetTextProperty(usageTextProp);

	vtkSmartPointer<vtkActor2D> usageTextActor = vtkSmartPointer<vtkActor2D>::New();
	usageTextActor->SetMapper(usageTextMapper);
	usageTextActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
	usageTextActor->GetPositionCoordinate()->SetValue(0.05, 0.95);
	//Segmenter
	//Segmenter* _segmenter = new Segmenter((vtkStructuredPoints*)(((vtkImageMapToColors*)selectionA->GetMapper()->GetInputAlgorithm()))->GetInput(), reader->GetOutput());
	myInteractorStyle->SetImageViewer(viewer, outputName, inputFilename, selectionA, (vtkStructuredPoints*)smoothed->GetOutput()/*reader->GetOutput()*/);
	myInteractorStyle->SetStatusMapper(sliceTextMapper);
	viewer->SetupInteractor(renderWindowInteractor);
	//mapperSel->SetScalarModeToUseCellData();


	// make the interactor use our own interactorstyle
	// cause SetupInteractor() is defining it's own default interatorstyle 
	// this must be called after SetupInteractor()
	renderWindowInteractor->Initialize();
	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	renderWindowInteractor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);
	viewer->GetRenderer()->AddActor2D(sliceTextActor);
	viewer->GetRenderer()->AddActor2D(usageTextActor);
	viewer->GetRenderer()->AddActor(selectionA);
	viewer->GetRenderer()->AddActor(crosshair);
	int displayExtent[6];
	viewer->GetImageActor()->GetDisplayExtent(displayExtent);
	selectionA->SetDisplayExtent(displayExtent);
	viewer->GetRenderWindow()->SetPosition(0, 0);
	viewer->GetRenderWindow()->SetSize(viewer->GetRenderWindow()->GetScreenSize()[0] / 2, viewer->GetRenderWindow()->GetScreenSize()[1]-50);
	viewer->GetRenderWindow()->SetWindowName("Slicer");

	// Have the sample listener receive events from the controller
	cout << "passed all creations " << endl;
	controller.addListener(listener);
	LeapAbstractionLayer* lal = LeapAbstractionLayer::getInstance();
	listener.SetInterface(lal);
	cout << "set viewer on listener" << endl;
	viewer->Render();
	viewer->GetRenderer()->ResetCamera();
	cout << "Camera reset" << endl;
	viewer->Render();

	//set focus
	SetForegroundWindow(FindWindow(NULL, "Slicer"));

	renderWindowInteractor->Start();
	controller.removeListener(listener);

	return EXIT_SUCCESS;
}
