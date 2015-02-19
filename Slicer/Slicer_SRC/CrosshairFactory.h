#ifndef __CROSSHAIR_FACTORY__
#define __CROSSHAIR_FACTORY__

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
#include <vtkCellData.h>
#include <vtkLookupTable.h>
#include <vtkStructuredPoints.h>
#include <vtkProperty.h>
#include "constants.h"
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>


/** This class is reponsible for creating the mesh using MarchingCubes algorithm*/
class CrosshairFactory {
	
public:
	static CrosshairFactory* getInstance();
	CrosshairFactory();
	~CrosshairFactory();
	vtkSmartPointer<vtkActor> makeCrosshair(double* size);
	
private:
	vtkSmartPointer<vtkActor> _actor;
	CrosshairFactory(CrosshairFactory const&){};
	CrosshairFactory& operator=(CrosshairFactory const&){};
	static CrosshairFactory* m_pInstance;
};




#endif
