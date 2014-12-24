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
						scalars->SetValue(id_selection, /*BACKGROUND*/  NOT_ACTIVE);
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
vtkStructuredPoints* ImageGraphCut::PerformSegmentation()
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

	return this->_selection;
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
	float unfloored_weight = 0;
	int total_weight = 0;
	int counter = 0;
	float pixelDifference;
	double sigma = this->ComputeNoise(tumor);
	cout << "+++++++++the noise sigma is: " << sigma << endl;
	int mult_constant = 1000;
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
					weight = mult_constant*std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[0], weight, weight);
					total_weight += weight;
					counter++;
					unfloored_weight += exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
				}

				//U
				if (neighbors[1] != -1){
					neighbor_point = ComputePointId(i,j+1,k);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = mult_constant*std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[1], weight, weight);
					total_weight += weight;
					counter++;
					unfloored_weight += exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
				}

				//I
				if (neighbors[2] != -1){
					neighbor_point = ComputePointId(i, j, k+1);
					pixelDifference = ComputeDifference(current_point, neighbor_point);
					weight = mult_constant*std::floor(exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma)));
					g->add_edge(counterID, neighbors[2], weight, weight);
					total_weight += weight;
					counter++;
					unfloored_weight += exp(-pow(pixelDifference, 2) / (2.0*sigma*sigma));
				}
				counterID++;
			}
		}
	}

	cout << "CreateNEdges_tumor:: average Nedges weight is: " << float(total_weight) / float(counter) << endl;
	cout << "CreateNEdges_tumor:: average Nedges weight BEFORE FLOORING is: " << float(unfloored_weight) / float(counter) << endl;

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
	
	int* bound_values;
	ComputeThresholds(tumor, bound_values);
	int tumor_lower_gray_value = bound_values[0];
	int tumor_upper_gray_value = bound_values[1];

	int counter_t = 0;
	int counter_s = 0;
	float average_s_weight = 0;
	float average_t_weight = 0;

	int counterID = 0;

	cout << "** tumor_lower_gray_value: " << tumor_lower_gray_value << endl;
	cout << "** tumor_upper_gray_value: " << tumor_upper_gray_value << endl;

	for (int i = bb.min_x; i < bb.max_x; i++){
		for (int j = bb.min_y; j < bb.max_y; j++) {
			for (int k = bb.min_z; k < bb.max_z; k++){
				
				
				GraphType* g = this->_map[tumor.getTId()];
				g->add_node();
				pointId = ComputePointId(i, j, k);

				int me = scalars_CT->GetValue(pointId) <= MIN_POSSIBLE_VALUE ? scalars_CT->GetValue(pointId) : (-1)*(65535 - scalars_CT->GetValue(pointId));

				////cout << "pointId: " << pointId << endl;
				if (selection_scalars->GetValue(pointId) == BACKGROUND)
				{
					g->add_tweights(counterID++, std::numeric_limits<int>::max(), (int)(SILENCING_FACTOR*(-std::abs(me - (tumor_lower_gray_value + tumor_upper_gray_value) / 2) + MIN_POSSIBLE_VALUE + tumor_upper_gray_value - tumor_lower_gray_value))/*std::numeric_limits<int>::min()*/);
					counter_t++;
					average_t_weight += (int)(SILENCING_FACTOR*(-std::abs(me - (tumor_lower_gray_value + tumor_upper_gray_value) / 2) + MIN_POSSIBLE_VALUE + tumor_upper_gray_value - tumor_lower_gray_value));
					continue;
				}
				if (selection_scalars->GetValue(pointId) == FOREGROUND)
				{
					g->add_tweights(counterID++, /*std::numeric_limits<int>::min()*/(int)(std::abs(me - (tumor_upper_gray_value + tumor_lower_gray_value) / 2) + MIN_POSSIBLE_VALUE), std::numeric_limits<int>::max());
					counter_s++;
					average_s_weight += (int)(std::abs(me - (tumor_upper_gray_value + tumor_lower_gray_value) / 2) + MIN_POSSIBLE_VALUE);
					continue;
				}
				if (selection_scalars->GetValue(pointId) == NOT_ACTIVE){					
					g->add_tweights(counterID++, (int)(std::abs(me - (tumor_upper_gray_value + tumor_lower_gray_value) / 2) + MIN_POSSIBLE_VALUE),
						(int)(SILENCING_FACTOR*(-std::abs(me - (tumor_lower_gray_value + tumor_upper_gray_value) / 2) + MIN_POSSIBLE_VALUE + tumor_upper_gray_value - tumor_lower_gray_value)));

					counter_s++;
					average_s_weight += (int)(std::abs(me - (tumor_upper_gray_value + tumor_lower_gray_value) / 2) + MIN_POSSIBLE_VALUE);
					counter_t++;
					average_t_weight += (int)(SILENCING_FACTOR*(-std::abs(me - (tumor_lower_gray_value + tumor_upper_gray_value) / 2) + MIN_POSSIBLE_VALUE + tumor_upper_gray_value - tumor_lower_gray_value));
					continue;
				}
			}
		}
	}
	cout << "CreateTEdges_tumor:: average t_weight is: " << float(average_t_weight) / float(counter_t) << endl;
	cout << "CreateTEdges_tumor:: average s_weight is: " << float(average_s_weight) / float(counter_s) << endl;
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

// This function computes the garpy thresholds for segmenting the tumor
void ImageGraphCut::ComputeThresholds(Tumor _tumor, int* arr){
	
	float sum = 0.0;
	int counter = 1;

	vtkPointData* pointsDataCT = this->_CT_image->GetPointData();
	vtkIntArray* scalars_CT = (vtkIntArray*)pointsDataCT->GetScalars();
	vtkIntArray* scalars_selection = (vtkIntArray*)this->_selection->GetPointData()->GetScalars();

	vector<Tumor::Point3D> tumor_points = _tumor.getPoints();
	for (int i = 0; i < tumor_points.size(); i++){
		
		vtkIdType pointId = ComputePointId(tumor_points.at(i).x, tumor_points.at(i).y, tumor_points.at(i).z);
		
		if (scalars_selection->GetValue(pointId) == FOREGROUND && scalars_CT->GetValue(pointId) < 80) {
			sum += scalars_CT->GetValue(pointId);
			counter++;
		}
	}

	float average_gray_marked = sum / counter;

	arr[0] = floor(0.01*average_gray_marked);
	arr[1] = round(2.0*average_gray_marked); //todo: should it be we relation to the compued noise per tumor?
	
	cout << "** average_gray_marked: " << average_gray_marked << endl;

}