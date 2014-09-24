#ifndef _MY_VTK_INTERCTOR_IMAGE_STYLE_
#define _MY_VTK_INTERCTOR_IMAGE_STYLE_

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
#include <vtkEventForwarderCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapToColors.h>

// Define own interaction style
class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
	static myVtkInteractorStyleImage* New();
	vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);
	vtkSmartPointer<vtkCallbackCommand> leapCallback;
	float _x_position, _y_position; // TODO: create getters and setters and move to protected.
	int _window_size_x;
	int _window_size_y;

protected:
	vtkSmartPointer<vtkImageViewer2> _ImageViewer;
	vtkTextMapper* _StatusMapper;
	int _orientation;
	int _MinSlice;
	int _MaxSlice;
	int _Slice;
	vtkSmartPointer<vtkImageMapToColors> _selection;
	vtkSmartPointer<vtkImageActor> _selection_actor;
	vtkSmartPointer<vtkStructuredPoints> _3D_selection;
	bool _isSliceLocked, _isPainting;
	virtual void SetInteractor(vtkRenderWindowInteractor* interactor);// method I overrode. 
	static void ProcessLeapEvents(vtkObject* object,
		unsigned long event,
		void* clientdata,
		void* calldata);
public:
	myVtkInteractorStyleImage::myVtkInteractorStyleImage();
	void SetImageViewer(vtkImageViewer2* imageViewer, int size_x, int size_y, vtkSmartPointer<vtkStructuredPoints> selection, vtkSmartPointer<vtkImageActor> selection_actor);
	void SetStatusMapper(vtkTextMapper* statusMapper);
	void setSlice(int slice);
	int getMaxSlice();

protected:
	void MoveSliceForward();
	void MoveSliceBackward();
	void ToggleOrientation();
	virtual void OnKeyDown();
	virtual void OnKeyUp();
	virtual void OnMouseWheelForward();
	virtual void OnMouseWheelBackward();
	void lockSlice(bool state);
	void startPainting(bool state);
};

#endif