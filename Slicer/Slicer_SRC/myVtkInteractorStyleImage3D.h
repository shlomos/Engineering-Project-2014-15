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
#include <vtkInteractorStyleJoystickCamera.h>
#include <vtkEventForwarderCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapToColors.h>
#include <vtkRendererCollection.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkRegularPolygonSource.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <vtkImageDilateErode3D.h>
#include "constants.h"
#include <vtkDataSetMapper.h>
#include <vtkStructuredPoints.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkExtractVOI.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkCellPicker.h>
#include <vtkPointPicker.h>
#include <vtkImageMapToColors.h>
#include <vtkPropPicker.h>
#include <algorithm>
#include <vtkCamera.h>
#include <sstream>
#include "LeapAbstractionLayer.h"
#include "MarchingCubes.h"
#include <vtkImageToStructuredPoints.h>
#include <vtkNIFTIImageWriter.h>
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
class myVtkInteractorStyleImage3D : public vtkInteractorStyleJoystickCamera
{
public:
	static myVtkInteractorStyleImage3D* New();
	virtual ~myVtkInteractorStyleImage3D();
	vtkTypeMacro(myVtkInteractorStyleImage3D, vtkInteractorStyleJoystickCamera);

protected:
	vtkTextMapper* _StatusMapper;
	std::string _outputName;
	std::string _inputName;
	LeapAbstractionLayer* _lal;
	int _drawSize;
	int _numTumors;
	bool _rotLock;
	bool _hfMode;
	bool _foundLeaks;
	vtkIdType _currSource;
	vtkStructuredPoints* _selection;
	boost::mutex _canSegment_mutex;

public:
	myVtkInteractorStyleImage3D();
	void SetStatusMapper(vtkTextMapper* statusMapper);
	void SetNumTumors(int);
	void Initialize(std::string outputName, std::string inputName, vtkStructuredPoints* _selection);
	void RemoveLeaks();

protected:
	void WriteToFile();
	void LoadFromFile();
	void doSegment();
	void ResetAll();
	void MakeAnnotation(vtkIdType);
	void GetPoked2DLocation(int[3]);
	virtual void OnKeyDown();
	virtual void OnLeftButtonDown();
	virtual void OnRightButtonDown();
	virtual void OnChar();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonUp();
	virtual void OnKeyUp();
	virtual void OnTimer();
	virtual void OnMouseMove();
	virtual void OnMouseWheelForward();
	void StartInteraction3D(bool state);
	virtual void OnMouseWheelBackward();
};

#endif