#include <stddef.h>
#include "LeapAbstractionLayer.h"

// Global static pointer used to ensure a single instance of the class.
LeapAbstractionLayer* LeapAbstractionLayer::m_pInstance = NULL;

LeapAbstractionLayer::LeapAbstractionLayer(){
	_isPainting = false;
	_sliceLock = false;
	_slice = 0;
	_x_position = 0;
	_y_position = 0;
	_z_position = 0;
}
/** This function is called to create an instance of the class.
Calling the constructor publicly is not allowed. The constructor
is private and is only called by this Instance function.
*/

LeapAbstractionLayer* LeapAbstractionLayer::getInstance()
{
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new LeapAbstractionLayer;

	return m_pInstance;
}

float LeapAbstractionLayer::getX(){
	return _x_position;
}

float LeapAbstractionLayer::getY(){
	return _y_position;
}

float LeapAbstractionLayer::getZ(){
	return _z_position;
}

int LeapAbstractionLayer::getSlice(){
	return _slice;
}

int LeapAbstractionLayer::getMaxSlice(){
	return _maxSlice;
}

bool LeapAbstractionLayer::getPainting(){
	return _isPainting;
}

bool LeapAbstractionLayer::getSliceLock(){
	return _sliceLock;
}

void LeapAbstractionLayer::setX(float x){
	_x_position = x;
}

void LeapAbstractionLayer::setY(float y){
	_y_position = y;
}

void LeapAbstractionLayer::setZ(float z){
	_z_position = z;
}

void LeapAbstractionLayer::setSlice(int slice){
	_slice = slice;
}

void LeapAbstractionLayer::setMaxSlice(int slice){
	_maxSlice = slice;
}

void LeapAbstractionLayer::setPainting(bool status){
	_isPainting = status;
}

void LeapAbstractionLayer::setSliceLock(bool status){
	_sliceLock = status;
}