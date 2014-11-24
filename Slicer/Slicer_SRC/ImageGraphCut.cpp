#include "ImageGraphCut.h"

// STL
#include <cmath>

void ImageGraphCut::SetImage(vtkStructuredPoints* const selection, vtkStructuredPoints* CT_image) 
{
	this->_selection = selection;
	this->_CT_image = CT_image;

}

void ImageGraphCut::Clean(){
	this->_tumors.clear();
	this->_map.clear();
}

void ImageGraphCut::CutGraph()
{
	cout << "CutGraph()..." << endl;

	int flow = Graph->maxflow();
	cout << "finished cutting the graph! \n" << endl;
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

void ImageGraphCut::CutGraphs() {
	cout << "CutGraphs()..." << endl;
	int flow;

	vtkIntArray* scalars = (vtkIntArray*)this->_selection->GetPointData()->GetScalars();

	for (vector<Tumor>::iterator it = this->_tumors.begin(); it != this->_tumors.end(); ++it) {
		
		//perform max flow for each graph
		flow = this->_map[(*it).getTId()]->maxflow();

		cout << "start update of mapper..." << endl;
		TumorBoundingBox bb = (*it).getBoundingBox();
		
		vtkIdType id_selection;
		int counterID = 0;
		for (int i = bb.min_x; i < bb.max_x; i++){
			for (int j = bb.min_y; j < bb.max_y; j++) {
				for (int k = bb.min_z; k < bb.max_z; k++){
					int point[3] = { i, j, k };
					id_selection = this->_selection->ComputePointId(point);

					if (this->_map[(*it).getTId()]->what_segment(counterID++) == GraphType::SOURCE) {
						scalars->SetValue(id_selection, NOT_ACTIVE);
					}
					else {
						scalars->SetValue(id_selection, FOREGROUND);
					}

				}
			}
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
	int ijk[3];
	vtkIdType current_point;
	bool createdNewTumor = true;
	int* extent = this->_selection->GetExtent();

	vtkPointData* pointsData = _selection->GetPointData();
	vtkIntArray* scalars = (vtkIntArray*)pointsData->GetScalars();
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

				if (scalars->GetValue(current_point) == FOREGROUND){
					Tumor::Point3D point = { i, j, k };
					createdNewTumor = AddPointToTumor(point);
				}
			}
		}
	}

	this->CreateGraphs();
	this->CutGraphs();


	cout << "Finished Segmentation!" << endl;
}

double ImageGraphCut::ComputeDifference(vtkIdType current_point, vtkIdType neighbor_point)
{
	vtkPointData* pointsDataCT = _CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();
	return ((double)std::abs(scalars_CT->GetValue(neighbor_point) - scalars_CT->GetValue(current_point)));
}

void ImageGraphCut::CreateNEdges_tumor(Tumor tumor){
	cout << "CreateNEdges_tumor.." << endl;

	// retreive points' scalars
	vtkPointData* pointsData = _selection->GetPointData();
	vtkIntArray* selection_scalars = (vtkIntArray*)pointsData->GetScalars();

	vtkPointData* pointsDataCT = _CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();
	
	TumorBoundingBox bb = tumor.getBoundingBox();


	int counterID = 0;
	int neighbors[3];
	int ijk[3];
	int weight;
	float pixelDifference;
	double sigma = this->ComputeNoise();
	vtkIdType current_point, neighbor_point;
	
	GraphType* g = this->_map[tumor.getTId()];

	for (int i = bb.min_x; i < bb.max_x; i++){
		for (int j = bb.min_y; j < bb.max_y; j++) {
			for (int k = bb.min_z; k < bb.max_z; k++){
				ijk[0] = i;
				ijk[1] = j;
				ijk[2] = k;
				current_point = this->_selection->ComputePointId(ijk);

				tumor.getNeighbors(counterID, neighbors);

				//R
				if (neighbors[0] != -1){
					ijk[0] = i+1;
					ijk[1] = j;
					ijk[2] = k;
					neighbor_point = this->_selection->ComputePointId(ijk);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = 5000 * std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[0], weight, weight);
				}

				//U
				if (neighbors[1] != -1){
					ijk[0] = i;
					ijk[1] = j+1;
					ijk[2] = k;
					neighbor_point = this->_selection->ComputePointId(ijk);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = 5000 * std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[1], weight, weight);
				}

				//I
				if (neighbors[2] != -1){
					ijk[0] = i;
					ijk[1] = j;
					ijk[2] = k+1;
					neighbor_point = this->_selection->ComputePointId(ijk);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = 5000 * std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[2], weight, weight);
				}
				counterID++;
			}
		}
	}

	return;
}

//void ImageGraphCut::CreateNEdges()
//{
//	std::cout << "CreateNEdges()" << std::endl;
//	// Create n-edges and set n-edge weights (links between image nodes)
//
//	//cout << "DEBUG!1" << endl;
//
//	int ijk[3];
//	int weight;
//	float pixelDifference;
//	vtkIdType current_point, neighbor_point;
//
//	// Estimate the "camera noise"
//	double sigma = this->ComputeNoise();
//
//	int* extent = this->_selection->GetExtent();
//
//	for (int i = 0; i < extent[1]; i++)
//	{
//		for (int j = 0; j < extent[3]; j++)
//		{
//			for (int k = 0; k < extent[5]; k++)
//			{
//
//				ijk[0] = i;
//				ijk[1] = j;
//				ijk[2] = k;
//				current_point = this->_selection->ComputePointId(ijk);
//				//cout << "(i, j, k): (" << i << ", " << j << ", " << k << "). -->> id: " << current_point << endl;
//				break;
//
//				// r
//				ijk[0] = i+1;
//				ijk[1] = j;
//				ijk[2] = k;
//				neighbor_point = this->_selection->ComputePointId(ijk);
//				//cout << "r: neighbor_point: " << neighbor_point << endl;
//				pixelDifference = ComputeDifference(current_point, neighbor_point);
//				//weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
//				weight = std::floor(10.0*exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
//				Graph->add_edge(current_point, neighbor_point, weight, weight);
//
//				// o
//				ijk[0] = i;
//				ijk[1] = j;
//				ijk[2] = k+1;
//				neighbor_point = this->_selection->ComputePointId(ijk);
//				//cout << "o: neighbor_point: " << neighbor_point << endl;
//				pixelDifference = ComputeDifference(current_point, neighbor_point);
//				//weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
//				weight = std::floor(10.0*exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
//				Graph->add_edge(current_point, neighbor_point, weight, weight);
//
//				// u
//				ijk[0] = i;
//				ijk[1] = j+1;
//				ijk[2] = k;
//				neighbor_point = this->_selection->ComputePointId(ijk);
//				//cout << "u: neighbor_point: " << neighbor_point << endl;
//				pixelDifference = ComputeDifference(current_point, neighbor_point);
//				//weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
//				weight = std::floor(10.0*exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
//				Graph->add_edge(current_point, neighbor_point, weight, weight);
//
//				//// ru
//				//ijk[0] = i+1;
//				//ijk[1] = j+1;
//				//ijk[2] = k;
//				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());
//
//				//// ruo
//				//ijk[0] = i+1;
//				//ijk[1] = j+1;
//				//ijk[2] = k-1;
//				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());
//
//				//// ro
//				//ijk[0] = i+1;
//				//ijk[1] = j;
//				//ijk[2] = k-1;
//				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());
//
//
//				//// uo
//				//ijk[0] = i;
//				//ijk[1] = j+1;
//				//ijk[2] = k-1;
//				//currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, current_point, this->_selection->ComputePointId(ijk), std::numeric_limits<float>::max());
//
//			}
//		}
//		//cout << "DEBUG!2" << endl;
//
//	}
//
//////	unsigned int expectedNumberOfNEdges = 2 * (
//////		imageSize[0] * (imageSize[1] - 1) + // vertical edges
//////		(imageSize[0] - 1) * imageSize[1]  // horizontal edges
//////		);
//////	std::cout << "Resizing for " << expectedNumberOfNEdges << " N-edges." << std::endl;
//////	this->EdgeWeights.resize(num_edges(this->Graph) + expectedNumberOfNEdges);
//////	this->ReverseEdges.resize(num_edges(this->Graph) + expectedNumberOfNEdges);
//////
//////	// We are only using a 4-connected structure,
//////	// so the kernel (iteration neighborhood) must only be
//////	// 3x3 (specified by a radius of 1)
//////	itk::Size<2> radius;
//////	radius.Fill(1);
//////
//////	typedef itk::ShapedNeighborhoodIterator<TImage> IteratorType;
//////
//////	// Traverse the image adding an edge between the current pixel
//////	// and the pixel below it and the current pixel and the pixel to the right of it.
//////	// This prevents duplicate edges (i.e. we cannot add an edge to
//////	// all 4-connected neighbors of every pixel or almost every edge would be duplicated.
//////	std::vector<typename IteratorType::OffsetType> neighbors;
//////	typename IteratorType::OffsetType bottom = { { 0, 1 } };
//////	neighbors.push_back(bottom);
//////	typename IteratorType::OffsetType right = { { 1, 0 } };
//////	neighbors.push_back(right);
//////
//////	typename IteratorType::OffsetType center = { { 0, 0 } };
//////
//////	IteratorType iterator(radius, this->Image, this->Image->GetLargestPossibleRegion());
//////	iterator.ClearActiveList();
//////	iterator.ActivateOffset(bottom);
//////	iterator.ActivateOffset(right);
//////	iterator.ActivateOffset(center);
//////
//////	// Estimate the "camera noise"
//////	double sigma = this->ComputeNoise();
//////
//////	unsigned int currentNumberOfEdges = num_edges(this->Graph);
//////	for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
//////	{
//////		PixelType centerPixel = iterator.GetPixel(center);
//////
//////		for (unsigned int i = 0; i < neighbors.size(); i++)
//////		{
//////			bool valid;
//////			iterator.GetPixel(neighbors[i], valid);
//////
//////			// If the current neighbor is outside the image, skip it
//////			if (!valid)
//////			{
//////				continue;
//////			}
//////
//////			PixelType neighborPixel = iterator.GetPixel(neighbors[i]);
//////
//////			// Compute the Euclidean distance between the pixel intensities
//////			float pixelDifference = PixelDifferenceFunctor.Difference(centerPixel, neighborPixel);
//////
//////			// Compute the edge weight
//////			float weight = exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
//////			assert(weight >= 0);
//////
//////			// Add the edge to the graph
//////			unsigned int node1 = this->NodeImage->GetPixel(iterator.GetIndex(center));
//////			unsigned int node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));
//////
//////			currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, node1, node2, weight);
//////		}
//////	}
//////
//	std::cout << "Finished CreateNEdges()" << std::endl;
//}


// This function returns the raw probability between 0 and 1 to be a tumor according to the gray value of the point
int ImageGraphCut::computeTumorProbability(double point_value)
{
	float max_factor = 1.0 / (MU*sqrt(2.0*CHICKEN_PI));
	return std::exp(-((pow(point_value-MU,2))/(2*pow(SIGMA,2))))/(SIGMA*sqrt(2.0*CHICKEN_PI))/(max_factor);
}


void ImageGraphCut::CreateTEdges_tumor(Tumor tumor) {

	std::cout << "CreateTEdges_tumor()" << std::endl;
	// retreive points' scalars
	vtkPointData* pointsData = _selection->GetPointData();
	vtkIntArray* selection_scalars = (vtkIntArray*)pointsData->GetScalars();

	vtkPointData* pointsDataCT = _CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();
	vtkIdType pointId;

	TumorBoundingBox bb = tumor.getBoundingBox();
	cout << "Bounding Box size:" << endl;
	cout << "bb.min_x, bb.min_x: " << bb.min_x << ", " << bb.max_x << endl;
	cout << "bb.min_y, bb.min_y: " << bb.min_y << ", " << bb.max_y << endl;
	cout << "bb.min_z, bb.min_z: " << bb.min_z << ", " << bb.max_z << endl;

	int counterID = 0;

	for (int i = bb.min_x; i < bb.max_x; i++){
		for (int j = bb.min_y; j < bb.max_y; j++) {
			for (int k = bb.min_z; k < bb.max_z; k++){
				//cout << "***tumor.getTId(): " << tumor.getTId() << endl;
				
				GraphType* g = this->_map[tumor.getTId()];
				g->add_node();

				int point[] = { i, j, k };
				pointId = _selection->ComputePointId((int*)point);
				//cout << "pointId: " << pointId << endl;
				if (selection_scalars->GetValue(pointId) == BACKGROUND)
				{
					g->add_tweights(counterID++, std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
					//cout << "BACKGROUND: selection_scalars->GetValue(pointId) is:  " << selection_scalars->GetValue(pointId) << endl;
					continue;
				}
				if (selection_scalars->GetValue(pointId) == FOREGROUND)
				{
					g->add_tweights(counterID++, std::numeric_limits<int>::max(), std::numeric_limits<int>::min());
					//cout << "FOREGROUND:: selection_scalars->GetValue(pointId) is; " << pointId << endl;
					continue;
				}
				if (selection_scalars->GetValue(pointId) == NOT_ACTIVE){					
					//todo: fix this equation!
					g->add_tweights(counterID++, std::abs(scalars_CT->GetValue(pointId) - (UPPER_BOUND + LOWER_BOUND) / 2) + MIN_POSSIBLE_VALUE, -std::abs(scalars_CT->GetValue(pointId) -(LOWER_BOUND+UPPER_BOUND) / 2 ) + MIN_POSSIBLE_VALUE + UPPER_BOUND - LOWER_BOUND);
					//g->add_tweights(counterID++, -std::abs(scalars_CT->GetValue(pointId) - (LOWER_BOUND + UPPER_BOUND) / 2) + MIN_POSSIBLE_VALUE + UPPER_BOUND - LOWER_BOUND, std::abs(scalars_CT->GetValue(pointId) - (UPPER_BOUND - 3 * LOWER_BOUND) / 2) + MIN_POSSIBLE_VALUE);
					//g->add_tweights(counterID++, std::abs(scalars_CT->GetValue(pointId) - 52), (int)((2958.0 / 42.0)*(3000 - std::abs(scalars_CT->GetValue(pointId) - 52)) - 2958.0*2958.0/42.0 + 42));
					//g->add_tweights(counterID++, std::abs(scalars_CT->GetValue(pointId) - 57) + 3000, 3094 - std::abs(scalars_CT->GetValue(pointId) - 57));
					continue;
				}
			}
		}
	}
	cout << "at the end of CreateTEdges_tumor" << endl;
}


//using ROIs
void ImageGraphCut::CreateGraphs() {

	for (int i = 0; i < this->_tumors.size(); i++){
		cout << "generating new graph for tumor number: " << i << "..." << endl;
		GraphType* g = new GraphType(this->_tumors.at(i).getBBoxSize(), this->_tumors.at(i).getBBoxSize() * 16);
		this->_map[i] = g;

		CreateTEdges_tumor(this->_tumors.at(i));
		CreateNEdges_tumor(this->_tumors.at(i));
		cout << "finished." << endl;
	}

	cout << "this->_tumors.size(): " << this->_tumors.size() << endl;
	cout << "this->_map.size()   : " << this->_map.size() << endl;
}

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

bool ImageGraphCut::AddPointToTumor(Tumor::Point3D point){
	int* extent = this->_selection->GetExtent();
	Tumor::Point3D limits = { extent[1], extent[3], extent[5] };

	for (int i = 0; i < this->_tumors.size(); i++){
		if (_tumors.at(i).isInTumor(point)){
			_tumors.at(i).addPoint(point, limits);
			return false;
		}
	}
	Tumor new_tumor(limits);
	new_tumor.addPoint(point, limits);
	_tumors.push_back(new_tumor);
	return true;
}
