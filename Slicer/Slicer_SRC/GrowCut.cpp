#include "GrowCut.h"


void GrowCut::SetImage(vtkStructuredPoints* CT_image, vtkStructuredPoints* segmentation) {
	this->_CT_image = CT_image;
	this->_segmentation = segmentation;
}

// returns the maximal color in the image
void GrowCut::init_tumors() {
	this->_tumors.clear();
	vtkIdType current_point;
	bool createdNewTumor = true;
	int* extent = this->_segmentation->GetExtent();
	cout << extent[1] << ":"
		<< extent[3] << ":"
		<< extent[5] << endl;
	vtkPointData* pointsData = _segmentation->GetPointData();
	vtkUnsignedShortArray* seg_scalars = (vtkUnsignedShortArray*)pointsData->GetScalars();

	vtkPointData* CT_pointsData = _CT_image->GetPointData();
	vtkIntArray* ct_scalars = (vtkIntArray*)CT_pointsData->GetScalars();
	int counter = 0;
	
	// create BB
	for (int i = 0; i < extent[1]; i++)
	{
		for (int j = 0; j < extent[3]; j++)
		{
			for (int k = 0; k < extent[5]; k++)
			{
				current_point = ComputePointId(i, j, k);
				int ct_color = ct_scalars->GetValue(current_point);
				// fix short overflow on negatives:
				ct_color = ct_color > 5000 ? (-1)*(65535 - ct_color) : ct_color;

				// add all the points in the bounding box to the tumor
				if (seg_scalars->GetValue(current_point) == FOREGROUND) {
					Tumor::Point3D point = { i, j, k, FOREGROUND, 1.0, ct_color };
					createdNewTumor = AddPointToTumor(point);
				}
			}
		}
	}


	for (int i = 0; i < extent[1]; i++)
	{
		for (int j = 0; j < extent[3]; j++)
		{
			for (int k = 0; k < extent[5]; k++)
			{
				current_point = ComputePointId(i, j, k);
				
				int ct_color = ct_scalars->GetValue(current_point);
				// fix short overflow on negatives:
				ct_color = ct_color > 5000 ? (-1)*(65535 - ct_color) : ct_color;
				max_color = max(max_color, ct_color);
				min_color = min(min_color, ct_color);
				
				// add all the points in the bounding box to the tumor
				Tumor::Point3D point = { i, j, k, seg_scalars->GetValue(current_point), 1.0f, ct_color};
				//cout << "D: " << seg_scalars->GetValue(current_point) << endl;
				if (point.point_type == NOT_ACTIVE){
					//cout << "DEBUG: 22" << endl;
					
					point.strength = 0.0f;
					point.prev_strength = 0.0f;
				}
				else{
					//cout << "DEBUG: 11" << endl;
					counter++;
				}
				createdNewTumor = AddPointToTumor(point);
			}
		}
	}
	cout << "there were num seeds: " << counter << endl;
}

vtkIdType GrowCut::ComputePointId(int i, int j, int k)
{
	int ijk[3];
	int* ext = _segmentation->GetExtent();
	if (i >= 0 && i < ext[1]){
		if (j >= 0 && j < ext[3]){
			if (k >= 0 && k < ext[5]){
				ijk[0] = i;
				ijk[1] = j;
				ijk[2] = k;
				return _segmentation->ComputePointId(ijk);
			}
		}
	}
	return -1;
}

bool GrowCut::AddPointToTumor(Tumor::Point3D point){
	int* extent = this->_segmentation->GetExtent();
	Tumor::Point3D limits = { extent[1], extent[3], extent[5], NULL, NULL, NULL};

	// check wether this point is already inside a tumor
	for (int i = 0; i < this->_tumors.size(); i++){
		if (_tumors.at(i).isInTumor(point)){
			_tumors.at(i).addPoint(point, limits);
			return false;
		}
	}

	if (point.point_type == FOREGROUND) {
		Tumor new_tumor(limits);
		new_tumor.addPoint(point, limits);
		_tumors.push_back(new_tumor);
		return true;
	}

	return false;
}

vtkStructuredPoints* GrowCut::PerformSegmentation(){

	cout << "GrowCut::PerformSegmentation()" << endl;
	
	//initialize tumors
	init_tumors();
	cout << "Tumors inited!" << endl;
	map<int, Tumor::Point3D>* points;
	int neighbors[6];
	vtkIdType point_id;

	// perform segmentation for each tumor
	cout << "this->_tumors.size(): " << this->_tumors.size() << endl;
	vtkPointData* pd = this->_segmentation->GetPointData();
	vtkUnsignedShortArray* scalars = (vtkUnsignedShortArray*)pd->GetScalars();
	int counter = 0;
	for (int i = 0; i < this->_tumors.size(); i++){
		cout << "tumor: " << i << endl;
		//perform T iterations:
		for (int t = 0; t < NUM_OF_GROW_ITERATIONS; t++){
			points = this->_tumors.at(i).getPoints();
			//cout << "iteration: " << endl;
			
			for (std::map<int, Tumor::Point3D>::iterator it = this->_tumors.at(i).active_points.begin()/*points->begin()*/; it != this->_tumors.at(i).active_points.end()/*points->end()*/; ++it) {
				Tumor::Point3D curr_point = it->second;
				this->_tumors.at(i).getNeighbors(curr_point, neighbors);
				for (int k = 0; k < 6; k++) {
					int key = neighbors[k];
					map<int, Tumor::Point3D>::iterator neighbor = points->find(key);
					if (neighbor == points->end()) {
						continue;
					}

					/*Tumor::Point3D np = points->at(key);
					if (np.prev_point_type == NOT_ACTIVE){
						continue;
					}*/

					//cout << np.color << ":" << curr_point.color << ":" << g_function(np.color, curr_point.color)*np.prev_strength << "vs" << curr_point.prev_strength << endl;
					float g_value = g_function(neighbor->second.color, curr_point.color)*curr_point.prev_strength;
					//cout << "g_value " << g_value << endl;
					//cout << "neighbor->second.prev_strength " << neighbor->second.prev_strength  << endl;
					if (g_value > neighbor->second.prev_strength){
						//cout << "DEBUG" << endl;
						neighbor->second.point_type = curr_point.prev_point_type;
						neighbor->second.strength = g_value;
						if (this->_tumors.at(i).active_points.find(key) != this->_tumors.at(i).active_points.end()) {
							this->_tumors.at(i).active_points.at(key).point_type = curr_point.prev_point_type;
							this->_tumors.at(i).active_points.at(key).strength = g_value;
						}
						else {
							this->_tumors.at(i).active_points.insert(std::pair<int, Tumor::Point3D>(key, neighbor->second));
						}

					}
					point_id = ComputePointId(it->second[0], it->second[1], it->second[2]);
					scalars->SetValue(point_id, it->second.point_type);
					scalars->Modified();
				}
			}

			for (std::map<int, Tumor::Point3D>::iterator it = this->_tumors.at(i).active_points.begin()/*points->begin()*/; it != this->_tumors.at(i).active_points.end()/*points->end()*/; ++it) {
				it->second.prev_point_type = it->second.point_type;
				it->second.prev_strength = it->second.strength;
				points->at(it->first).prev_point_type = it->second.point_type;
				points->at(it->first).prev_strength = it->second.strength;
			}
		}
		cout << "after T iterations!" << endl;

		points = this->_tumors.at(i).getPoints();
		TumorBoundingBox bBox = this->_tumors.at(i).getBoundingBox();
		//cout << "bBox.min_x :" << bBox.min_x << endl;
		//cout << "bBox.min_y :" << bBox.min_y << endl;
		//cout << "bBox.min_z :" << bBox.min_z << endl;
		//cout << "bBox.max_x :" << bBox.max_x << endl;
		//cout << "bBox.max_y :" << bBox.max_y << endl;
		//cout << "bBox.max_z :" << bBox.max_z << endl;
		
		
		//for (std::map<int, Tumor::Point3D>::iterator it = points->begin(); it != points->end(); ++it) {
		//	point_id = ComputePointId(it->second[0], it->second[1], it->second[2]);
		//	if (it->second.point_type == BACKGROUND){
		//		//counter++;
		//		it->second.point_type = NOT_ACTIVE;
		//	}
		//	counter++;
		//	scalars->SetValue(point_id, it->second.point_type);
		//}
		cout << "points->size()" << points->size() << endl;
		cout << "in tumor there are " << this->_tumors.at(i).numSeeds << " seeds!" << endl;
		cout << "this->_tumors.at(i).getBBoxSize(): " << this->_tumors.at(i).getBBoxSize() << endl;
		cout << "num of segmented is: " << counter << endl;
		cout << "" << endl; 
	}
	cout << "Done." << endl;
	return this->_segmentation;
}



float GrowCut::g_function(int Cp, int Cq) {
	//cout << "max min" << this->max_color << ":" << this->min_color << endl;
	//cout << "formula: " << 1.0-(float(abs(Cp - Cq)) / float(this->max_color - this->min_color)) << endl;
	//return 1.0 - (float(abs(Cp - Cq)) / float(this->max_color - this->min_color));
	return exp(-1.5*float(abs(Cp - Cq))*float(abs(Cp - Cq)));
}