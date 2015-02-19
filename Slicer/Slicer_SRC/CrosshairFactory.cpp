#include <stddef.h>
#include "CrosshairFactory.h"


// Global static pointer used to ensure a single instance of the class.
CrosshairFactory* CrosshairFactory::m_pInstance = NULL;

CrosshairFactory::CrosshairFactory(){
}

CrosshairFactory* CrosshairFactory::getInstance()
{
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new CrosshairFactory;

	return m_pInstance;
}


vtkSmartPointer<vtkActor> CrosshairFactory::makeCrosshair(double* size) {
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