#ifndef __MARCHING_CUBES__
#define __MARCHING_CUBES__

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageMapToColors.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkLookupTable.h>
#include <vtkStructuredPoints.h>
#include <vtkProperty.h>
#include "CrosshairFactory.h"
#include "constants.h"
#include <vtkCursor2D.h>
#include <vtkImageMaskBits.h>
#include <vtkCamera.h>
#include <vtkFollower.h>
#include <vtkFillHolesFilter.h>
#include <vtkPolyDataNormals.h>
#include "myVtkInteractorStyleImage3D.h"
#include <vtkDijkstraGraphGeodesicPath.h>
#include "meshLeaksCorrection.h"
#include <vtkUnsignedShortArray.h>


/** This class is reponsible for creating the mesh using MarchingCubes algorithm*/
class MarchingCubes {
	
public:
	MarchingCubes(vtkStructuredPoints* selection, std::string inputName);
	~MarchingCubes();
	
private:
	vtkSmartPointer<vtkMarchingCubes>  _surface;
	vtkSmartPointer<vtkRenderer>       _renderer;
	vtkSmartPointer<vtkRenderWindow>   _renderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> _interactor;
	vtkSmartPointer<vtkPolyDataMapper> _mapper;
	vtkSmartPointer<vtkActor>          _actor;
	vtkSmartPointer<vtkImageMaskBits>  _mask;
	std::string _inputName;

	//vtkSmartPointer<vtkImageViewer2>   _viewer;
	vtkStructuredPoints* _selection;

};




#endif
