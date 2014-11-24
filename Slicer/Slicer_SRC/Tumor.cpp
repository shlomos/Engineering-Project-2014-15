#include "Tumor.h"

int Tumor::numTumors = 0;

Tumor::Tumor(Tumor::Point3D limits)
{
	this->t_id = Tumor::numTumors;
	Tumor::numTumors++;
	this->bBox.max_x = 0;
	this->bBox.max_y = 0;
	this->bBox.max_z = 0;
	this->bBox.min_x = limits[0];
	this->bBox.min_y = limits[1];
	this->bBox.min_z = limits[2];
}

void Tumor::addPoint(Tumor::Point3D point, Tumor::Point3D limits){
	//try{
	//	//assert(sizeof(point) / sizeof(*point) == 3);
	//}
	//catch (exception e){
	//	cout << "Walla Lo!" << endl;
	//	exit(1);
	//}

	this->points.push_back(point);

	//cout << "limits:" << endl;
	//cout << "limits[0]: " << *limits[0] << 
	//	"limits[1]: " << *limits[1] << 
	//	"limits[2]: " << *limits[2] << endl;
	
	//updating bounding box
	this->bBox.max_x = max(this->bBox.max_x, min(limits[0],point[0] + BB_SPACE_XY));
	this->bBox.min_x = min(this->bBox.min_x, max(0,point[0] - BB_SPACE_XY));

	this->bBox.max_y = max(this->bBox.max_y, min(limits[1], point[1] + BB_SPACE_XY));
	this->bBox.min_y = min(this->bBox.min_y, max(0, point[1] - BB_SPACE_XY));

	this->bBox.max_z = max(this->bBox.max_z, min(limits[2], point[2] + BB_SPACE_Z));
	this->bBox.min_z = min(this->bBox.min_z, max(0, point[2] - BB_SPACE_Z));
}

vector<Tumor::Point3D> Tumor::getPoints(){
	return this->points;
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
	if (point[0] > bBox.max_x + BB_SPACE_XY || point[0] < bBox.min_x - BB_SPACE_XY){
		return false;
	}

	if (point[1] > bBox.max_y + BB_SPACE_XY || point[1] < bBox.min_y - BB_SPACE_XY){
		return false;
	}

	if (point[2] > bBox.max_z + BB_SPACE_Z || point[2] < bBox.min_z - BB_SPACE_Z){
		return false;
	}
	return true;
}

int Tumor::getBBoxSize() {
	return (std::abs(this->bBox.max_x - this->bBox.min_x))*
		   (std::abs(this->bBox.max_y - this->bBox.min_y))*
		   (std::abs(this->bBox.max_z - this->bBox.min_z));
}

void Tumor::getNeighbors(vtkIdType id, int* neighbors) {
	int size_x = this->bBox.max_x - this->bBox.min_x;
	int size_y = this->bBox.max_y - this->bBox.min_y;
	int size_z = this->bBox.max_z - this->bBox.min_z;
	int totalSize = size_x*size_y*size_z;
	int planeSize = size_x*size_y;

	float idOnPlane = id % planeSize;
	float myRow = ceil((float)idOnPlane / (float)size_x);
	
	//r
	if (idOnPlane >= max(0.0f,size_x*myRow - 1)){
		neighbors[0] = -1;
	}
	else{
		neighbors[0] = id + 1;
	}

	//cout << "asdf movie: " << idOnPlane << "myRow: " << myRow << "size_x: " << size_x << "neighbors[0]: " << neighbors[0] << endl;
	
	//u
	if (myRow >= size_y - 1){
		neighbors[1] = -1;
	}
	else{
		neighbors[1] = id + size_x;
	}

	//o
	if (id > totalSize - size_x*size_y - 1){
		neighbors[2] = -1;
	}
	else{
		neighbors[2] = id + size_x*size_y;
	}

}