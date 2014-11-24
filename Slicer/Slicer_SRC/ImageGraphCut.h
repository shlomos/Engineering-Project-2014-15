/*

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ImageGraphCut_H
#define ImageGraphCut_H

//// Custom
//#include "PixelDifference.h"
//
//// Submodules
//#include "Mask/ForegroundBackgroundSegmentMask.h"
//
//// ITK
//#include "itkImage.h"
//#include "itkSampleToHistogramFilter.h"
//#include "itkHistogram.h"
//#include "itkListSample.h"

// STL
#include <vector>

//// Boost
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/one_bit_color_map.hpp>
//#include <boost/graph/stoer_wagner_min_cut.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/graph/grid_graph.hpp>

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
  //ImageGraphCut(TPixelDifferenceFunctor pixelDifferenceFunctor) :
    //PixelDifferenceFunctor(pixelDifferenceFunctor){}

  //TPixelDifferenceFunctor PixelDifferenceFunctor;

  /** This is a special type to keep track of the graph node labels. */
  //typedef itk::Image<unsigned int, 2> NodeImageType;

  /** The type of the histograms. */
  //typedef itk::Statistics::Histogram< float, itk::Statistics::DenseFrequencyContainer2 > HistogramType;

  /** The type of a list of pixels/indexes. */
  //typedef std::vector<itk::Index<2> > IndexContainer;

  /** Several initializations are done here. */
  void SetImage(vtkStructuredPoints* const image, vtkStructuredPoints* CT_image);

  /** Get the image that we are segmenting. */
  //TImage* GetImage();

  /** Create and cut the graph (The main driver function). */
  void PerformSegmentation();

  void ImageGraphCut::Clean();

  /** Add point to the closest tumor, if there is one close enough, and create new Tumor if not*/
  bool AddPointToTumor(Tumor::Point3D point);

  /** Return a list of the selected (via scribbling) pixels. */
  //IndexContainer GetSources();
 // IndexContainer GetSinks();

  /** Set the selected pixels. */
 // void SetSources(const IndexContainer& sources);
  //void SetSinks(const IndexContainer& sinks);

  /** Get the output of the segmentation. */
  //ForegroundBackgroundSegmentMask* GetSegmentMask();

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

  ///** Store the list of edges and their corresponding reverse edges. */
  //std::vector<EdgeDescriptor> ReverseEdges;

  ///** Create an edge on the graph. */
  //unsigned int AddBidirectionalEdge(unsigned int numberOfEdges, const unsigned int source,
  //                          const unsigned int target,
  //                          const float weight);

  /** The main graph object. */
  GraphType* Graph;


 // /** Maintain a list of all of the edge weights. */
 // std::vector<float> EdgeWeights;

 // /** The output segmentation */
 // //ForegroundBackgroundSegmentMask::Pointer ResultingSegments;

 // /** User specified foreground points */
 // //IndexContainer Sources;

 // /** User specified background points */
 //// IndexContainer Sinks;

  /** The weighting between unary and binary terms */
  float Lambda = 0.01f;

  /** The number of bins per dimension of the foreground and background histograms */
  int NumberOfHistogramBins = 10;

  /** An image which keeps tracks of the mapping between pixel index and graph node id */
  //NodeImageType::Pointer NodeImage;

  // Typedefs
  //typedef typename TImage::PixelType PixelType;
  //typedef itk::Statistics::ListSample<PixelType> SampleType;
  //typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType> SampleToHistogramFilterType;

  /** Create the histograms from the users selections */
  void CreateSamples();

  /** Estimate the "camera noise" */
  double ComputeNoise();

  /** Compute the distance between two neighbors pixels*/
  double ComputeDifference(vtkIdType currernt_point, vtkIdType neighbor);

  ///** Create a Kolmogorov graph structure from the image and selections */
  //void CreateGraph();

  ///** Create the edges between pixels and neighboring pixels (the grid). */
  //void CreateNEdges();

  ///** Create the edges between pixels and the terminals (source and sink). */
  //void CreateTEdges();

  //------------------------------
  /** Create the edges between pixels and the terminals (source and sink) for the tumors*/
  void CreateTEdges_tumor(Tumor tumor);

  /** Create the edges between pixels and the terminals (source and sink) for the tumors*/
  void CreateNEdges_tumor(Tumor tumor);

  /** Create a Kolmogorov graph structure from the image and selections for every tumor */
  void CreateGraphs();

  /** Perform the s-t min cut for every tumor, and udpate the mapper accordingly*/
  void CutGraphs();

  /** Attach all the graphs into*/
  //------------------------------


  /** Perform the s-t min cut */
  void CutGraph();

  /** Compute the probability that a point is part of a tumor using only its gray value, using prior knowledge of MU adn SIGMA*/
  int computeTumorProbability(double point_value);

  /** The selection information of the image */
  vtkStructuredPoints* _selection;
  vtkStructuredPoints* _CT_image;

  /** The gray values of the image */
  vtkStructuredPoints* _image;


  ///** Store the node ids of the terminals */
  //unsigned int SourceNodeId;
  //unsigned int SinkNodeId;
};

#endif

