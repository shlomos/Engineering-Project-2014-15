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

protected:
	vtkSmartPointer<vtkImageViewer2> _ImageViewer;
	vtkTextMapper* _StatusMapper;
	int _orientation;
	int _MinSlice;
	int _MaxSlice;
	int _Slice;
	std::string _outputName;
	int _drawSize;
	vtkSmartPointer<vtkImageActor> _selection_actor;
	bool _hfMode;
	bool _isSliceLocked;
	bool _isPainting;
	virtual void SetInteractor(vtkRenderWindowInteractor* interactor);// method I overrode. 
	static void ProcessLeapEvents(vtkObject* object,
		unsigned long event,
		void* clientdata,
		void* calldata);
public:
	myVtkInteractorStyleImage::myVtkInteractorStyleImage();
	void SetImageViewer(vtkImageViewer2* imageViewer,std::string outputName, vtkSmartPointer<vtkImageActor> selection_actor);
	void SetStatusMapper(vtkTextMapper* statusMapper);
	void setSlice(int slice);
	void lockSlice(bool state);
	void SetPainting(bool state);
	int getMaxSlice();

protected:
	void MoveSliceForward();
	void MoveSliceBackward();
	void ToggleOrientation();
	void WriteToFile();
	virtual void OnKeyDown();
	virtual void OnKeyUp();
	virtual void OnMouseWheelForward();
	virtual void OnMouseWheelBackward();
};

#endif