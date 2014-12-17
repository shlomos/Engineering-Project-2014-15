#ifndef _MY_VTK_INTERCTOR_IMAGE_STYLE_3D_
#define _MY_VTK_INTERCTOR_IMAGE_STYLE_3D_

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
#include <vtkInteractorStyleJoystickActor.h>
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
#include "LeapAbstractionLayer.h"
#include "MarchingCubes.h"

#include <boost/thread.hpp>



// helper class to format slice status message
class StatusMessage3D {
public:
	static std::string Format(int slice, int maxSlice, int orientation) {
		std::stringstream tmp;
		switch (orientation){
		case SLICE_ORIENTATION_YZ:
			tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1 << "\nOrientation " << "Sagittal";
			break;
		case SLICE_ORIENTATION_XZ:
			tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1 << "\nOrientation " << "Axial";
			break;
		case SLICE_ORIENTATION_XY:
			tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1 << "\nOrientation " << "Coronal";
			break;
		}
		return tmp.str();
	}
};


// Define own interaction style
class myVtkInteractorStyleImage3D : public vtkInteractorStyleJoystickActor
{
public:
	static myVtkInteractorStyleImage3D* New();
	virtual ~myVtkInteractorStyleImage3D();
	vtkTypeMacro(myVtkInteractorStyleImage3D, vtkInteractorStyleJoystickActor);

protected:
	vtkTextMapper* _StatusMapper;
	std::string _outputName;
	LeapAbstractionLayer* _lal;
	int _drawSize;
	bool _rotLock;
	bool _hfMode;

public:
	myVtkInteractorStyleImage3D();
	void SetStatusMapper(vtkTextMapper* statusMapper);
	void Initialize(std::string outputName);

protected:
	void WriteToFile();
	void LoadFromFile();
	void doSegment();
	void ResetAll();
	virtual void OnKeyDown();
	virtual void OnLeftButtonDown();
	virtual void OnRightButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonUp();
	virtual void OnKeyUp();
	virtual void OnTimer();
	virtual void OnMouseMove();
	virtual void OnMouseWheelForward();
	virtual void OnMouseWheelBackward();
};

#endif