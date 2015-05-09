#ifndef _MY_VTK_INTERCTOR_IMAGE_STYLE_
#define _MY_VTK_INTERCTOR_IMAGE_STYLE_

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
#include <vtkStructuredPointsWriter.h>
#include <vtkCellPicker.h>
#include <vtkPointPicker.h>
#include <vtkImageMapToColors.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPropPicker.h>
#include <algorithm>
#include <vtkCamera.h>
#include <sstream>
#include <vtkImageCast.h>
#include "GrowCut.h"
#include "MarchingCubes.h"
#include "LeapAbstractionLayer.h" 
#include <Windows.h>

#include <boost/thread.hpp>



// helper class to format slice status message
class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice, int orientation) {
		std::stringstream tmp;
		switch (orientation){
		case SLICE_ORIENTATION_YZ:
			tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1 << "\nOrientation " << "Sagittal";
			break;
		case SLICE_ORIENTATION_XZ:
			tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1 << "\nOrientation " << "Coronal";
			break;
		case SLICE_ORIENTATION_XY:
			tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1 << "\nOrientation " << "Axial";
			break;
		}
		return tmp.str();
	}
};


// Define own interaction style
class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
	static myVtkInteractorStyleImage* New();
	virtual ~myVtkInteractorStyleImage();
	vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);
	//float _x_position, _y_position, _z_position; // TODO: create getters and setters and move to protected.
	vtkSmartPointer<vtkUnsignedShortArray> _selection_scalars;

protected:
	vtkSmartPointer<vtkImageViewer2> _ImageViewer;
	LeapAbstractionLayer* _lal;
	vtkTextMapper* _StatusMapper;

	//ImageGraphCut* _graph_cut;
	GrowCut* _grow_cut;

	MarchingCubes* _marching_cubes;
	int _orientation;
	std::string _outputName;
	int _drawSize;
	vtkSmartPointer<vtkImageActor> _selection_actor;
	vtkStructuredPoints* _CT_image;
	std::string _inputName;
	bool _hfMode;
	//bool _isSliceLocked;
	boost::mutex _canSegment_mutex;

public:
	myVtkInteractorStyleImage();
	void SetImageViewer(vtkImageViewer2* imageViewer, std::string outputName, std::string inputName, vtkSmartPointer<vtkImageActor> selection_actor, vtkStructuredPoints* CT_image);
	void SetStatusMapper(vtkTextMapper* statusMapper);
	void doSegment();
	//void setSlice(int slice);
	//void lockSlice(bool state);
	//void SetPainting(bool state);
	int getMaxSlice();
	double* redrawCrossHair();
	void SetCrossHair(int cross_x, int cross_y);
	void applyCameraFixes();

protected:

	void MoveSliceForward();
	void MoveSliceBackward();
	void ToggleOrientation();
	void WriteToFile();
	void LoadFromFile();
	void marchingCubes();
	void ResetAll();
	virtual void OnKeyDown();
	virtual void OnKeyUp();
	virtual void OnMouseMove();
	virtual void OnTimer();
	virtual void OnMouseWheelForward();
	virtual void OnMouseWheelBackward();
};

#endif