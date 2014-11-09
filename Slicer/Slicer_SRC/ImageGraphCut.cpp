#include "ImageGraphCut.h"

// STL
#include <cmath>

void ImageGraphCut::SetImage(vtkStructuredPoints* const selection, vtkStructuredPoints* CT_image) 
{
	this->_selection = selection;
	this->_CT_image = CT_image;

}

void ImageGraphCut::CutGraph()
{
	cout << "CutGraph()..." << endl;
	// Compute mininum cut

	//boost::graph_traits<GraphType>::vertex_descriptor s = vertex(this->SourceNodeId, this->Graph);
	//boost::graph_traits<GraphType>::vertex_descriptor t = vertex(this->SinkNodeId, this->Graph);
	//std::vector<int> groups(num_vertices(this->Graph));
	//std::vector<float> residual_capacity(num_edges(this->Graph)); //this needs to be initialized to 0
	//boykov_kolmogorov_max_flow(this->Graph,
	//	boost::make_iterator_property_map(&EdgeWeights[0], get(boost::edge_index, this->Graph)),
	//	boost::make_iterator_property_map(&residual_capacity[0], get(boost::edge_index, this->Graph)),
	//	boost::make_iterator_property_map(&ReverseEdges[0], get(boost::edge_index, this->Graph)),
	//	boost::make_iterator_property_map(&groups[0], get(boost::vertex_index, this->Graph)),
	//	get(boost::vertex_index, this->Graph),
	//	s,
	//	t);

	int flow = Graph->maxflow();
	cout << "++++ FLOW IS: " << flow << endl;

	cout << "finished cutting the graph! :)" << endl;

	cout << "start update of mapper..." << endl;
	vtkIntArray* scalars = (vtkIntArray*)this->_selection->GetPointData()->GetScalars();

	for (int i = 0; i < this->_selection->GetNumberOfPoints(); i++)
	{
		
		if (Graph->what_segment(i) == GraphType::SOURCE) {
			scalars->SetValue(i, NOT_ACTIVE);
		}
		else {
			scalars->SetValue(i, FOREGROUND);
		}
	}

	cout << "finished update the mapper. Modifying..." << endl;
	scalars->Modified();
	cout << "done." << endl;

}

//
void ImageGraphCut::PerformSegmentation()
{
	cout << "perform segmentation()" << endl;
	this->CreateGraph();
	this->CutGraph();
	cout << "Finished Segmentation!" << endl;
}

double ImageGraphCut::ComputeDifference(vtkIdType current_point, vtkIdType neighbor_point)
{
	vtkPointData* pointsDataCT = _CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();
	return ((double)std::abs(scalars_CT->GetValue(neighbor_point) - scalars_CT->GetValue(current_point)));
}
// 
void ImageGraphCut::CreateNEdges()
{
	std::cout << "CreateNEdges()" << std::endl;
	// Create n-edges and set n-edge weights (links between image nodes)

	cout << "DEBUG!1" << endl;

	int ijk[3];
	int weight;
	float pixelDifference;
	vtkIdType current_point, neighbor_point;

	// Estimate the "camera noise"
	double sigma = this->ComputeNoise();

	int* extent = this->_selection->GetExtent();

	//cout << this->_selection->GetNumberOfPoints() << endl;
	//cout << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] << endl;


	for (int i = 0; i < extent[1]; i++)
	{
		for (int j = 0; j < extent[3]; j++)
		{
			for (int k = 0; k < extent[5]; k++)
			{

				ijk[0] = i;
				ijk[1] = j;
				ijk[2] = k;
				current_point = this->_selection->ComputePointId(ijk);
				//cout << "(i, j, k): (" << i << ", " << j << ", " << k << "). -->> id: " << current_point << endl;
				break;

				// r
				ijk[0] = i+1;
				ijk[1] = j;
				ijk[2] = k;
				neighbor_point = this->_selection->ComputePointId(ijk);
				//cout << "r: neighbor_point: " << neighbor_point << endl;
				pixelDifference = ComputeDifference(current_point, neighbor_point);
				//weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
				weight = std::floor(10.0*exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
				Graph->add_edge(current_point, neighbor_point, weight, weight);

				// o
				ijk[0] = i;
				ijk[1] = j;
				ijk[2] = k+1;
				neighbor_point = this->_selection->ComputePointId(ijk);
				//cout << "o: neighbor_point: " << neighbor_point << endl;
				pixelDifference = ComputeDifference(current_point, neighbor_point);
				//weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
				weight = std::floor(10.0*exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
				Graph->add_edge(current_point, neighbor_point, weight, weight);

				// u
				ijk[0] = i;
				ijk[1] = j+1;
				ijk[2] = k;
				neighbor_point = this->_selection->ComputePointId(ijk);
				//cout << "u: neighbor_point: " << neighbor_point << endl;
				pixelDifference = ComputeDifference(current_point, neighbor_point);
				//weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
				weight = std::floor(10.0*exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
				Graph->add_edge(current_point, neighbor_point, weight, weight);

				//// ru
				//ijk[0] = i+1;
				//ijk[1] = j+1;
				//ijk[2] = k;
				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());

				//// ruo
				//ijk[0] = i+1;
				//ijk[1] = j+1;
				//ijk[2] = k-1;
				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());

				//// ro
				//ijk[0] = i+1;
				//ijk[1] = j;
				//ijk[2] = k-1;
				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());


				//// uo
				//ijk[0] = i;
				//ijk[1] = j+1;
				//ijk[2] = k-1;
				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());

			}
		}
		//cout << "DEBUG!2" << endl;

	}

//	unsigned int expectedNumberOfNEdges = 2 * (
//		imageSize[0] * (imageSize[1] - 1) + // vertical edges
//		(imageSize[0] - 1) * imageSize[1]  // horizontal edges
//		);
//	std::cout << "Resizing for " << expectedNumberOfNEdges << " N-edges." << std::endl;
//	this->EdgeWeights.resize(num_edges(this->Graph) + expectedNumberOfNEdges);
//	this->ReverseEdges.resize(num_edges(this->Graph) + expectedNumberOfNEdges);
//
//	// We are only using a 4-connected structure,
//	// so the kernel (iteration neighborhood) must only be
//	// 3x3 (specified by a radius of 1)
//	itk::Size<2> radius;
//	radius.Fill(1);
//
//	typedef itk::ShapedNeighborhoodIterator<TImage> IteratorType;
//
//	// Traverse the image adding an edge between the current pixel
//	// and the pixel below it and the current pixel and the pixel to the right of it.
//	// This prevents duplicate edges (i.e. we cannot add an edge to
//	// all 4-connected neighbors of every pixel or almost every edge would be duplicated.
//	std::vector<typename IteratorType::OffsetType> neighbors;
//	typename IteratorType::OffsetType bottom = { { 0, 1 } };
//	neighbors.push_back(bottom);
//	typename IteratorType::OffsetType right = { { 1, 0 } };
//	neighbors.push_back(right);
//
//	typename IteratorType::OffsetType center = { { 0, 0 } };
//
//	IteratorType iterator(radius, this->Image, this->Image->GetLargestPossibleRegion());
//	iterator.ClearActiveList();
//	iterator.ActivateOffset(bottom);
//	iterator.ActivateOffset(right);
//	iterator.ActivateOffset(center);
//
//	// Estimate the "camera noise"
//	double sigma = this->ComputeNoise();
//
//	unsigned int currentNumberOfEdges = num_edges(this->Graph);
//	for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
//	{
//		PixelType centerPixel = iterator.GetPixel(center);
//
//		for (unsigned int i = 0; i < neighbors.size(); i++)
//		{
//			bool valid;
//			iterator.GetPixel(neighbors[i], valid);
//
//			// If the current neighbor is outside the image, skip it
//			if (!valid)
//			{
//				continue;
//			}
//
//			PixelType neighborPixel = iterator.GetPixel(neighbors[i]);
//
//			// Compute the Euclidean distance between the pixel intensities
//			float pixelDifference = PixelDifferenceFunctor.Difference(centerPixel, neighborPixel);
//
//			// Compute the edge weight
//			float weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
//			assert(weight >= 0);
//
//			// Add the edge to the graph
//			unsigned int node1 = this->NodeImage->GetPixel(iterator.GetIndex(center));
//			unsigned int node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));
//
//			currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, node1, node2, weight);
//		}
//	}
//
	std::cout << "Finished CreateNEdges()" << std::endl;
}

int ImageGraphCut::computeTumorProbability(double point_value)
{
	return std::floor(100.0*std::exp(-((pow(point_value-MU,2))/(2*pow(SIGMA,2))))/(SIGMA*sqrt(2.0*CHICKEN_PI))/(0.01638972434));
	//todo: convert to int!
}

void ImageGraphCut::createEdgesTest()
{
	cout << "testing..." << endl;
	for (int j = 0; j < this->_selection->GetNumberOfPoints(); j++) {
		this->Graph->add_node();
	}

	for (int i = 0; i < this->_selection->GetNumberOfPoints()-1; i++){
		
		// check if t_edge is to be added
		if (i == 0) {
			this->Graph->add_tweights(0, 10, 5);
		}

		if (i == this->_selection->GetNumberOfPoints() / 2) {
			this->Graph->add_edge(i, i + 1, 1, 1);
		}
		else {
			//cout << "i is: " << i << endl;
			this->Graph->add_edge(i, i + 1, 2, 2);
		}
	}
	this->Graph->add_tweights(this->_selection->GetNumberOfPoints()-1, 4, 10);
}



// Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node)
void ImageGraphCut::CreateTEdges()
{
	std::cout << "CreateTEdges()" << std::endl;

	// retreive points' scalars
	vtkPointData* pointsData = _selection->GetPointData();
	vtkIntArray* selection_scalars = (vtkIntArray*)pointsData->GetScalars();

	vtkPointData* pointsDataCT = _CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();

	//cout << "number of vertices is: " << this->_selection->GetNumberOfPoints() << endl;
	for (int i = 0; i < this->_selection->GetNumberOfPoints(); i++){
		Graph->add_node();
		if (selection_scalars->GetValue(i) == BACKGROUND)
		{
			Graph->add_tweights(i, std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
			//cout << "BACKGROUND: selection_scalars->GetValue(i) " << selection_scalars->GetValue(i) << endl;
			continue;
		}
		if (selection_scalars->GetValue(i) == FOREGROUND)
		{
			Graph->add_tweights(i, std::numeric_limits<int>::max(), std::numeric_limits<int>::min());
			//cout << "FOREGROUND:: id is" << i << endl;
			continue;
		}
		if (selection_scalars->GetValue(i) == NOT_ACTIVE){
			Graph->add_tweights(i, 1.0 - computeTumorProbability(scalars_CT->GetValue(i)), computeTumorProbability(scalars_CT->GetValue(i)));
			//Graph->add_tweights(i, std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
			continue;
		}
	}
	cout << "at the end of createTedges" << endl;
}


	

	//// Compute the histograms of the selected foreground and background pixels
	//CreateSamples();

	//itk::ImageRegionIterator<TImage> imageIterator(this->Image,	this->Image->GetLargestPossibleRegion());
	//itk::ImageRegionIterator<NodeImageType>	nodeIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
	//imageIterator.GoToBegin();
	//nodeIterator.GoToBegin();

	//while (!imageIterator.IsAtEnd())
	//{
	//	PixelType pixel = imageIterator.Get();
	//	//std::cout << "Pixels have size: " << pixel.Size() << std::endl;

	//	HistogramType::MeasurementVectorType measurementVector(pixel.Size());
	//	for (unsigned int i = 0; i < pixel.Size(); i++)
	//	{
	//		measurementVector[i] = pixel[i];
	//	}

	//	HistogramType::IndexType backgroundIndex;
	//	this->BackgroundHistogram->GetIndex(measurementVector, backgroundIndex);
	//	float sinkHistogramValue =
	//		this->BackgroundHistogram->GetFrequency(backgroundIndex);

	//	HistogramType::IndexType foregroundIndex;
	//	this->ForegroundHistogram->GetIndex(measurementVector, foregroundIndex);
	//	float sourceHistogramValue =
	//		this->ForegroundHistogram->GetFrequency(foregroundIndex);

	//	// Conver the histogram value/frequency to make it as if it came from a normalized histogram
	//	if (this->BackgroundHistogram->GetTotalFrequency() == 0 ||
	//		this->ForegroundHistogram->GetTotalFrequency() == 0)
	//	{
	//		throw std::runtime_error("The foreground or background histogram TotalFrequency is 0!");
	//	}

	//	sinkHistogramValue /= this->BackgroundHistogram->GetTotalFrequency();
	//	sourceHistogramValue /= this->ForegroundHistogram->GetTotalFrequency();

	//	if (sinkHistogramValue <= 0)
	//	{
	//		sinkHistogramValue = tinyValue;
	//	}
	//	if (sourceHistogramValue <= 0)
	//	{
	//		sourceHistogramValue = tinyValue;
	//	}

	//	// Add the edge to the graph and set its weight
	//	// log() is the natural log
	//	currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, nodeIterator.Get(), this->SinkNodeId, -this->Lambda*log(sourceHistogramValue));
	//	currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, nodeIterator.Get(), this->SourceNodeId, -this->Lambda*log(sinkHistogramValue));

	//	++imageIterator;
	//	++nodeIterator;
	//}

	// Set very high source weights for the pixels that were
	// selected as foreground by the user
//	for (unsigned int i = 0; i < this->Sources.size(); i++)
//	{
//		currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sources[i]), this->SourceNodeId,
//			this->Lambda * std::numeric_limits<float>::max());
//
//		currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sources[i]), this->SinkNodeId, 0);
//	}
//
//	// Set very high sink weights for the pixels that
//	// were selected as background by the user
//	for (unsigned int i = 0; i < this->Sinks.size(); i++)
//	{
//		currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sinks[i]), this->SourceNodeId, 0);
//
//		currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sinks[i]), this->SinkNodeId,
//			this->Lambda * std::numeric_limits<float>::max());
//	}
//
//	std::cout << "Finished CreateTEdges()" << std::endl;
//}

void ImageGraphCut::CreateGraph()
{
	cout << "CreateGraph()" << endl;

	Graph = new GraphType(this->_selection->GetNumberOfPoints(), this->_selection->GetNumberOfPoints()*16);

	CreateTEdges();
	CreateNEdges();
	//createEdgesTest();
	cout << "done creating the graph " << endl;
}

//template <typename TImage, typename TPixelDifferenceFunctor>
//void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateGraph()
//{
//	std::cout << "CreateGraph()" << std::endl;
//
//	// Add all of the nodes to the graph and store their IDs in a "node image"
//	itk::ImageRegionIterator<NodeImageType> nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
//	nodeImageIterator.GoToBegin();
//
//	unsigned int nodeId = 0;
//	while (!nodeImageIterator.IsAtEnd())
//	{
//		nodeImageIterator.Set(nodeId);
//		nodeId++;
//		++nodeImageIterator;
//	}
//
//	// Set the sink and source ids to be the two numbers immediately following the number of vertices in the grid
//	this->SinkNodeId = nodeId;
//	nodeId++;
//	this->SourceNodeId = nodeId;
//
//	CreateNEdges();
//	CreateTEdges();
//
//	{
//		itk::Size<2> imageSize = this->NodeImage->GetLargestPossibleRegion().GetSize();
//
//		std::cout << "Number of edges " << num_edges(this->Graph) << std::endl;
//		int expectedEdges = imageSize[0] * imageSize[1] * 2 * 2 + // one '2' is because there is an edge to both the source and sink from each pixel, and the other '2' is because they are double edges (bidirectional)
//			2 * (imageSize[0] - 1)*imageSize[1] + 2 * imageSize[0] * (imageSize[1] - 1); // both '2's are for the double edges (this is the number of horizontal edges + the number of vertical edges)
//		std::cout << "(Should be " << expectedEdges << " edges.)" << std::endl;
//	}
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
double ImageGraphCut::ComputeNoise()
{
//	// Compute an estimate of the "camera noise". This is used in the N-weight function.
//
//	// Since we use a 4-connected neighborhood, the kernel must be 3x3 (a rectangular radius of 1 creates a kernel side length of 3)
//	itk::Size<2> radius;
//	radius[0] = 1;
//	radius[1] = 1;
//
//	typedef itk::ShapedNeighborhoodIterator<TImage> IteratorType;
//
//	std::vector<typename IteratorType::OffsetType> neighbors;
//	typename IteratorType::OffsetType bottom = { { 0, 1 } };
//	neighbors.push_back(bottom);
//	typename IteratorType::OffsetType right = { { 1, 0 } };
//	neighbors.push_back(right);
//
//	typename IteratorType::OffsetType center = { { 0, 0 } };
//
//	IteratorType iterator(radius, this->Image, this->Image->GetLargestPossibleRegion());
//	iterator.ClearActiveList();
//	iterator.ActivateOffset(bottom);
//	iterator.ActivateOffset(right);
//	iterator.ActivateOffset(center);
//
//	double sigma = 0.0;
//	int numberOfEdges = 0;
//
//	// Traverse the image collecting the differences between neighboring pixel intensities
//	for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
//	{
//		PixelType centerPixel = iterator.GetPixel(center);
//
//		for (unsigned int i = 0; i < neighbors.size(); i++)
//		{
//			bool valid;
//			iterator.GetPixel(neighbors[i], valid);
//			if (!valid)
//			{
//				continue;
//			}
//
//			PixelType neighborPixel = iterator.GetPixel(neighbors[i]);
//
//			float colorDifference = PixelDifferenceFunctor.Difference(centerPixel, neighborPixel);
//			sigma += colorDifference;
//			numberOfEdges++;
//		}
//	}
//
//	// Normalize
//	sigma /= static_cast<double>(numberOfEdges);
//
	double sigma = 1.0;
	return sigma;
}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//typename ImageGraphCut<TImage, TPixelDifferenceFunctor>::IndexContainer ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetSources()
//{
//	return this->Sources;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetLambda(const float lambda)
//{
//	this->Lambda = lambda;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetNumberOfHistogramBins(int bins)
//{
//	this->NumberOfHistogramBins = bins;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//ForegroundBackgroundSegmentMask* ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetSegmentMask()
//{
//	return this->ResultingSegments;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//typename ImageGraphCut<TImage, TPixelDifferenceFunctor>::IndexContainer ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetSinks()
//{
//	return this->Sinks;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetSources(const IndexContainer& sources)
//{
//	this->Sources = sources;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetSinks(const IndexContainer& sinks)
//{
//	this->Sinks = sinks;
//}
//
//template <typename TImage, typename TPixelDifferenceFunctor>
//TImage* ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetImage()
//{
//	return this->Image;
//}
