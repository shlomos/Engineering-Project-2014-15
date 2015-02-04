#ifndef __TUMOR__
#define __TUMOR__

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <vtkSmartPointer.h>
#include <stdlib.h>
#include <vtkpolydatamapper.h>
#include <vtkstructuredpointsreader.h>
#include <vtkimageviewer2.h>
#include <vtkcallbackcommand.h>
#include <vtkobjectfactory.h>
#include <vtktextmapper.h>
#include <vtkinteractorstyleimage.h>
#include <vtkeventforwardercommand.h>
#include <vtkdatasetmapper.h>
#include <vtkimagemaptocolors.h>
#include <vtkrenderercollection.h>
#include <vtkpolydata.h>
#include <vtkcelldata.h>
#include <vtkregularpolygonsource.h>
#include <vtkpolydata.h>
#include <vtksmartpointer.h>
#include <vtkpolydatamapper.h>
#include <vtkactor.h>
#include <vtkrenderwindow.h>
#include <vtkrenderer.h>
#include <vtkrenderwindowinteractor.h>
#include <vtkdatasetmapper.h>
#include <vtkstructuredpoints.h>
#include <vtkcell.h>
#include <vtkpointdata.h>
#include <vtkimageactor.h>
#include <vtkimagemapper3d.h>
#include <vtkextractvoi.h>
#include <vtkstructuredpointswriter.h>
#include <vtkcellpicker.h>
#include <vtkpointpicker.h>
#include <vtkimagemaptocolors.h>
#include <vtkproppicker.h>
#include <algorithm>
#include <sstream>
#include <limits>
#include "constants.h"
#include <map>

#define BB_SPACE_XY 40
#define BB_SPACE_Z  10


using namespace std;

typedef struct {
	int min_x, max_x;
	int min_y, max_y;
	int min_z, max_z;
}TumorBoundingBox;

/* stands for a single tumor and its bounding box*/
class Tumor
{
public:
	struct Point3D{
		int x, y, z;
		int nil = NULL;
		
		int point_type;
		int prev_point_type;
		float strength;
		float prev_strength;
		int color;
		
		int& operator[](int q) { switch (q){ case 0: return x; case 1: return y; case 2: return z; default: return nil; } }
		Point3D(int _x, int _y, int _z, int _point_type, float _strength, int _color) : x(_x), y(_y), z(_z), point_type(_point_type), prev_point_type(_point_type), strength(_strength), prev_strength(_strength), color(_color){

		}
	};

	map<int, Tumor::Point3D> points;
	map<int, Tumor::Point3D> active_points;

	Tumor(Tumor::Point3D limits);
	virtual ~Tumor(){
//		delete points;
	}
	

protected:
	static int numTumors;
	int t_id;
	TumorBoundingBox bBox;

public:
	int numSeeds;
	/** Add new point to this tumor, and update the bounding box accordingly*/
	void addPoint(Point3D, Point3D);

	/** Return this tumor's set of points*/
	map<int, Tumor::Point3D>* getPoints();
	
	/** Return this tumor's ID*/
	int getTId();

	/** Return this tumor's bounding box*/
	TumorBoundingBox getBoundingBox();

	/** TO BO FILLED*/
	void extendBoundingBox(int x, int y, int z);

	/** Checks wether this point is inside this particular tumor. Return TRUE if it is, and FLASE otherwise*/
	bool isInTumor(Point3D point);

	/** Returns the estimated size of the bounding box*/
	int getBBoxSize();

	/** put in neighbors an array with all the neighbors' ids*/
	void getNeighbors(Tumor::Point3D point, int* neighbors);
};

#endif