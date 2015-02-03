#ifndef __GROW_CUT__
#define __GROW_CUT__

// STL
#include <vector>

// vtk and stuff
#include <vtkSmartPointer.h>
#include <stdlib.h>
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
#include <vtkEventForwarderCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapToColors.h>
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
#include <vtkPropPicker.h>
#include <algorithm>
#include <vtkCamera.h>
#include <sstream>
#include <limits>
#include <vector>
#include <map>
#include "graph.h"
#include "Tumor.h"

class GrowCut {
public:
	void SetImage(vtkStructuredPoints* CT_image, vtkStructuredPoints* segmentation);
	vtkStructuredPoints* PerformSegmentation();
	void init_tumors();
	vtkIdType ComputePointId(int i, int j, int k);
	bool AddPointToTumor(Tumor::Point3D point);
	float g_function(int, int);
private:
	vtkStructuredPoints* _CT_image;
	vtkStructuredPoints* _segmentation;
	
	/** Save the tumors and their graphs*/
	typedef vector<Tumor> Tumors;
	Tumors _tumors;

	int max_color = -10000000;
	int min_color = 10000000;
};

#endif