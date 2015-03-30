#include "Tumor.h"

int Tumor::numTumors = 0;

Tumor::Tumor(Tumor::Point3D limits)
{
	this->t_id = Tumor::numTumors;
	Tumor::numTumors++;
	this->numSeeds = 0;
	this->bBox.max_x = 0;
	this->bBox.max_y = 0;
	this->bBox.max_z = 0;
	this->bBox.min_x = limits[0];
	this->bBox.min_y = limits[1];
	this->bBox.min_z = limits[2];
}

void Tumor::addPoint(Tumor::Point3D point, Tumor::Point3D limits){

	// awsome!
	int map_key = (point[0] + 2) + 1000 * (point[1] + 2) + 1000000 * (point[2] + 2);
	//cout << "map_key :" << map_key << endl;
	/*cout << "point[0]" << point[0];
	cout << " point[1] " << point[1];
	cout << " point[2] " << point[2] << endl;*/
	//map<int, Tumor::Point3D>::iterator currPoint = this->active_points.find(map_key);
	//if (currPoint != active_points.end()){
	//	currPoint->second.
	//}
	//else{
		this->points.insert(std::pair<int, Tumor::Point3D>(map_key, point));
		this->numSeeds++;
	//}

	//update active points' list
	if (point.point_type != NOT_ACTIVE) {
		this->active_points.insert(std::pair<int, Tumor::Point3D>(map_key, point));
	}

	if (point.point_type != FOREGROUND) {
		return;
	}
	//updating bounding box
	this->bBox.max_x = max(this->bBox.max_x, min(limits[0], point[0] + BB_SPACE_XY));
	//cout << "this->bBox.min_x" << this->bBox.min_x << endl;
	this->bBox.min_x = min(this->bBox.min_x, max(0, point[0] - BB_SPACE_XY));
	//cout << "this->bBox.min_x" << this->bBox.min_x << endl;

	this->bBox.max_y = max(this->bBox.max_y, min(limits[1], point[1] + BB_SPACE_XY));
	this->bBox.min_y = min(this->bBox.min_y, max(0, point[1] - BB_SPACE_XY));

	this->bBox.max_z = max(this->bBox.max_z, min(limits[2], point[2] + BB_SPACE_Z));
	this->bBox.min_z = min(this->bBox.min_z, max(0, point[2] - BB_SPACE_Z));
}

map<int, Tumor::Point3D>* Tumor::getPoints(){
	return &this->points;
}

int Tumor::getTId(){
	return this->t_id;
}

TumorBoundingBox Tumor::getBoundingBox(){
	return this->bBox;
}

void Tumor::extendBoundingBox(int x, int y, int z){
	//TBD	 
}

bool Tumor::isInTumor(Tumor::Point3D point){
	if (point[0] > bBox.max_x/* + BB_SPACE_XY*/ || point[0] < bBox.min_x/* - BB_SPACE_XY*/){
		return false;
	}

	if (point[1] > bBox.max_y/* + BB_SPACE_XY*/ || point[1] < bBox.min_y/* - BB_SPACE_XY*/){
		return false;
	}

	if (point[2] > bBox.max_z/* + BB_SPACE_Z*/ || point[2] < bBox.min_z/* - BB_SPACE_Z*/){
		return false;
	}
	return true;
}

int Tumor::getBBoxSize() {
	return (std::abs(this->bBox.max_x - this->bBox.min_x))*
		(std::abs(this->bBox.max_y - this->bBox.min_y))*
		(std::abs(this->bBox.max_z - this->bBox.min_z));
}

// should return 6/18/26 neighbors as a vector<Tumor::Point3D>
void Tumor::getNeighbors(Tumor::Point3D point, int* neighbors) {
	int x = point[0] + 2;
	int y = point[1] + 2;
	int z = point[2] + 2;

	neighbors[0] = (x - 1) + 1000 * y + 1000000 * z;
	neighbors[1] = (x + 1) + 1000 * y + 1000000 * z;
	neighbors[2] = x + 1000 * (y + 1) + 1000000 * z;
	neighbors[3] = x + 1000 * (y - 1) + 1000000 * z;
	neighbors[4] = x + 1000 * y + 1000000 * (z - 1);
	neighbors[5] = x + 1000 * y + 1000000 * (z + 1);


	/*cout << "neighbors[0]" << neighbors[0] << endl;
	cout << "neighbors[1]" << neighbors[1] << endl;
	cout << "neighbors[2]" << neighbors[2] << endl;
	cout << "neighbors[3]" << neighbors[3] << endl;
	cout << "neighbors[4]" << neighbors[4] << endl;
	cout << "neighbors[5]" << neighbors[5] << endl;*/

	//neighbors[0] = (x - 1) + 1000 * y + 1000000 * z;
	//neighbors[1] = (x - 1) + 1000 * y + 1000000 * (z + 1);
	//neighbors[2] = (x - 1) + 1000 * y + 1000000 * (z - 1);
	//neighbors[3] = (x - 1) + 1000 * (y - 1) + 1000000 * z;
	//neighbors[4] = (x - 1) + 1000 * (y - 1) + 1000000 * (z + 1);
	//neighbors[5] = (x - 1) + 1000 * (y - 1) + 1000000 * (z - 1);
	//neighbors[6] = (x - 1) + 1000 * (y + 1) + 1000000 * z;
	//neighbors[7] = (x - 1) + 1000 * (y + 1) + 1000000 * (z + 1);
	//neighbors[8] = (x - 1) + 1000 * (y + 1) + 1000000 * (z - 1);
	////
	//neighbors[9] = (x + 1) + 1000 * y + 1000000 * z;
	//neighbors[10] = (x + 1) + 1000 * y + 1000000 * (z + 1);
	//neighbors[11] = (x + 1) + 1000 * y + 1000000 * (z - 1);
	//neighbors[12] = (x + 1) + 1000 * (y - 1) + 1000000 * z;
	//neighbors[13] = (x + 1) + 1000 * (y - 1) + 1000000 * (z + 1);
	//neighbors[14] = (x + 1) + 1000 * (y - 1) + 1000000 * (z - 1);
	//neighbors[15] = (x + 1) + 1000 * (y + 1) + 1000000 * z;
	//neighbors[16] = (x + 1) + 1000 * (y + 1) + 1000000 * (z + 1);
	//neighbors[17] = (x + 1) + 1000 * (y + 1) + 1000000 * (z - 1);
	////
	//neighbors[18] = x + 1000 * y + 1000000 * (z + 1);
	//neighbors[19] = x + 1000 * y + 1000000 * (z - 1);
	//neighbors[20] = x + 1000 * (y - 1) + 1000000 * z;
	//neighbors[21] = x + 1000 * (y - 1) + 1000000 * (z + 1);
	//neighbors[22] = x + 1000 * (y - 1) + 1000000 * (z - 1);
	//neighbors[23] = x + 1000 * (y + 1) + 1000000 * z;
	//neighbors[24] = x + 1000 * (y + 1) + 1000000 * (z + 1);
	//neighbors[25] = x + 1000 * (y + 1) + 1000000 * (z - 1);








}
//int size_x = this->bBox.max_x - this->bBox.min_x;
//int size_y = this->bBox.max_y - this->bBox.min_y;
//int size_z = this->bBox.max_z - this->bBox.min_z;
//int totalSize = size_x*size_y*size_z;
//int planeSize = size_x*size_y;

//float idOnPlane = id % planeSize;
//float myRow = ceil((float)idOnPlane / (float)size_x);
//
////r
//if (idOnPlane >= max(0.0f,size_x*myRow - 1)){
//	neighbors[0] = -1;
//}
//else{
//	neighbors[0] = id + 1;
//}

////cout << "asdf movie: " << idOnPlane << "myRow: " << myRow << "size_x: " << size_x << "neighbors[0]: " << neighbors[0] << endl;
//
////u
//if (myRow >= size_y - 1){
//	neighbors[1] = -1;
//}
//else{
//	neighbors[1] = id + size_x;
//}

////o
//if (id > totalSize - size_x*size_y - 1){
//	neighbors[2] = -1;
//}
//else{
//	neighbors[2] = id + size_x*size_y;
//}