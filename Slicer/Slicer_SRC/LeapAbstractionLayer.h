#ifndef _LEAP_ABSTRACTION_LAYER_
#define _LEAP_ABSTRACTION_LAYER_

#include <string>
#include <algorithm>

class LeapAbstractionLayer{
public:
	static LeapAbstractionLayer* getInstance();
	float getX();
	float getY();
	float getZ();
	int getSlice();
	int getMaxSlice();
	bool getPainting();
	bool getSliceLock();
	void setX(float x);
	void setY(float y);
	void setZ(float z);
	void setSlice(int slice);
	void setMaxSlice(int slice);
	void setSliceLock(bool status);
	void setPainting(bool status);

private:
	float _x_position;
	float _y_position;
	float _z_position;
	int _slice;
	int _maxSlice;
	bool _isPainting;
	bool _sliceLock;
	LeapAbstractionLayer();
	LeapAbstractionLayer(LeapAbstractionLayer const&){};
	LeapAbstractionLayer& operator=(LeapAbstractionLayer const&){};
	static LeapAbstractionLayer* m_pInstance;
};

#endif