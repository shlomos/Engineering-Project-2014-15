#ifndef ImageGraphCut_H
#define ImageGraphCut_H

//
//// ITK
//#include "itkImage.h"
//#include "itkSampleToHistogramFilter.h"
//#include "itkHistogram.h"
//#include "itkListSample.h"

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


/** Perform graph cut based segmentation on an image. Image pixels can be only grayscale.
  */
class ImageGraphCut
{
public:

	~ImageGraphCut()
	{
		delete Graph;
		for (int i = 0; i < this->_map.size(); i++) {
			delete this->_map[i];
		}
	}

	/** Several initializations are done here. */
	void SetImage(vtkStructuredPoints* const image, vtkStructuredPoints* CT_image);

	/** Create and cut the graph (The main driver function). */
	vtkStructuredPoints* PerformSegmentation();

	void ImageGraphCut::Clean();

	/** Add point to the closest tumor, if there is one close enough, and create new Tumor if not*/
	bool AddPointToTumor(Tumor::Point3D point);
	vtkIdType ComputePointId(int i, int j, int k);


	/** Set the weight between the regional and boundary terms. */
	void SetLambda(const float);

	/** Set the number of bins per dimension of the foreground and background histograms. */
	void SetNumberOfHistogramBins(const int);

protected:

	/** A graph object for Kolmogorov*/
	typedef Graph<int, int, int> GraphType;

	typedef vector<Tumor> Tumors;
	typedef vector<GraphType*> Graphs;

	// map between tumors and their graphs.
	typedef map<int, GraphType*> Tumors_to_graphs_map;

	/** Save the tumors and their graphs*/
	Tumors _tumors;
	Tumors_to_graphs_map _map;

	/** The main graph object. */
	GraphType* Graph;

	/** The weighting between unary and binary terms */
	float Lambda = 0.01f;

	/** The number of bins per dimension of the foreground and background histograms */
	int NumberOfHistogramBins = 10;

	/** Create the histograms from the users selections */
	void CreateSamples();

	/** Estimate the "camera noise" */
	double ComputeNoise(Tumor tumor);

	/** Compute the distance between two neighbors pixels*/
	double ComputeDifference(vtkIdType currernt_point, vtkIdType neighbor);


	//------------------------------
	/** Create the edges between pixels and the terminals (source and sink) for the tumors*/
	void CreateTEdges_tumor(Tumor tumor);

	/** Create the edges between pixels and the terminals (source and sink) for the tumors*/
	void CreateNEdges_tumor(Tumor tumor);

	/** Create a Kolmogorov graph structure from the image and selections for every tumor */
	void CreateGraphs();

	/** Perform the s-t min cut for every tumor, and udpate the mapper accordingly*/
	void CutGraphs();

	/** Perform the s-t min cut */
	void CutGraph();

	/** Compute the probability that a point is part of a tumor using only its gray value, using prior knowledge of MU adn SIGMA*/
	int computeTumorProbability(double point_value);

	/** The selection information of the image */
	vtkStructuredPoints* _selection;
	vtkStructuredPoints* _CT_image;

	/** The gray values of the image */
	vtkStructuredPoints* _image;

};

#endif

