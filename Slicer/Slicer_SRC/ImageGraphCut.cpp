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
					id_selection = ComputePointId(i,j,k);

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
				current_point = ComputePointId(i, j, k);

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
	// Fix overflow of unsigned type.
	double me = scalars_CT->GetValue(current_point) <= MIN_POSSIBLE_VALUE ? scalars_CT->GetValue(current_point) : (-1)*(65535 - scalars_CT->GetValue(current_point));
	double neighbor = scalars_CT->GetValue(neighbor_point) <= MIN_POSSIBLE_VALUE ? scalars_CT->GetValue(neighbor_point) : (-1)*(65535 - scalars_CT->GetValue(neighbor_point));
	//cout <<"neighboor "<< scalars_CT->GetValue(neighbor_point) << endl;
	//cout <<"me "<< scalars_CT->GetValue(current_point) << endl;
	return ((double)std::abs(neighbor - me));
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
	int weight;
	float pixelDifference;
	double sigma = this->ComputeNoise(tumor);
	vtkIdType current_point, neighbor_point;
	
	GraphType* g = this->_map[tumor.getTId()];

	for (int i = bb.min_x; i < bb.max_x; i++){
		for (int j = bb.min_y; j < bb.max_y; j++) {
			for (int k = bb.min_z; k < bb.max_z; k++){
				current_point = ComputePointId(i, j, k);

				tumor.getNeighbors(counterID, neighbors);

				//R
				if (neighbors[0] != -1){
					neighbor_point = ComputePointId(i+1, j, k);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight =  std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[0], weight, weight);
				}

				//U
				if (neighbors[1] != -1){
					neighbor_point = ComputePointId(i,j+1,k);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[1], weight, weight);
				}

				//I
				if (neighbors[2] != -1){
					neighbor_point = ComputePointId(i, j, k+1);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[2], weight, weight);
				}
				counterID++;
			}
		}
	}

	return;
}



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
				pointId = ComputePointId(i, j, k);
				//cout << "pointId: " << pointId << endl;
				if (selection_scalars->GetValue(pointId) == BACKGROUND)
				{
					g->add_tweights(counterID++, std::numeric_limits<int>::max(), 0/*std::numeric_limits<int>::min()*/);
					//cout << "BACKGROUND: selection_scalars->GetValue(pointId) is:  " << selection_scalars->GetValue(pointId) << endl;
					continue;
				}
				if (selection_scalars->GetValue(pointId) == FOREGROUND)
				{
					g->add_tweights(counterID++, /*std::numeric_limits<int>::min()*/0, std::numeric_limits<int>::max());
					//cout << "FOREGROUND:: selection_scalars->GetValue(pointId) is; " << pointId << endl;
					continue;
				}
				if (selection_scalars->GetValue(pointId) == NOT_ACTIVE){					
					//todo: fix this equation!
					int me = scalars_CT->GetValue(pointId) <= MIN_POSSIBLE_VALUE ? scalars_CT->GetValue(pointId) : (-1)*(65535 - scalars_CT->GetValue(pointId));
					g->add_tweights(counterID++, (int)(std::abs(me - (UPPER_BOUND + LOWER_BOUND) / 2) + MIN_POSSIBLE_VALUE), 
						(int)(SILENCING_FACTOR*(-std::abs(me - (LOWER_BOUND+UPPER_BOUND) / 2 ) + MIN_POSSIBLE_VALUE + UPPER_BOUND - LOWER_BOUND)));
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

vtkIdType ImageGraphCut::ComputePointId(int i, int j, int k)
{
	int ijk[3];
	int* ext = _selection->GetExtent();
	if (i>=0 && i < ext[1]){
		if (j>=0 && j < ext[3]){
			if (k>=0 && k < ext[5]){
				ijk[0] = i;
				ijk[1] = j;
				ijk[2] = k;
				return _selection->ComputePointId(ijk);
			}
		}
	}
	return -1;
}


double ImageGraphCut::ComputeNoise(Tumor tumor)
{
	vtkPointData* pointsDataCT = _CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();
	vtkIdType pointId, neighbor;

	vector<Tumor::Point3D> tumor_points = tumor.getPoints();
	double sigma = 0.0;
	int edge_counter = 0;

	for (int i = 0; i < tumor_points.size(); i++)
	{
		pointId = ComputePointId(tumor_points.at(i).x, tumor_points.at(i).y, tumor_points.at(i).z);
		// check neighbors for each point

		//r
		neighbor = ComputePointId(tumor_points.at(i).x+1, tumor_points.at(i).y, tumor_points.at(i).z);
		if (neighbor != -1){
			sigma += ComputeDifference(pointId, neighbor);
			edge_counter++;
		}

		//l
		neighbor = ComputePointId(tumor_points.at(i).x - 1, tumor_points.at(i).y, tumor_points.at(i).z);
		if (neighbor != -1){
			sigma += ComputeDifference(pointId, neighbor);
			edge_counter++;
		}

		//u
		neighbor = ComputePointId(tumor_points.at(i).x, tumor_points.at(i).y+1, tumor_points.at(i).z);
		if (neighbor != -1){
			sigma += ComputeDifference(pointId, neighbor);
			edge_counter++;
		}

		//d
		neighbor = ComputePointId(tumor_points.at(i).x, tumor_points.at(i).y-1, tumor_points.at(i).z);
		if (neighbor != -1){
			sigma += ComputeDifference(pointId, neighbor);
			edge_counter++;
		}


		//i
		neighbor = ComputePointId(tumor_points.at(i).x, tumor_points.at(i).y, tumor_points.at(i).z+1);
		if (neighbor != -1){
			sigma += ComputeDifference(pointId, neighbor);
			edge_counter++;
		}

		//o
		neighbor = ComputePointId(tumor_points.at(i).x, tumor_points.at(i).y, tumor_points.at(i).z - 1);
		if (neighbor != -1){
			sigma += ComputeDifference(pointId, neighbor);
			edge_counter++;
		}

	}
	

	// Normalize
	sigma /= static_cast<double>(edge_counter);
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
