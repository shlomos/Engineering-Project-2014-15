#ifndef _SLICER_LEAP_LISTENER_
#define _SLICER_LEAP_LISTENER_

#include "Leap.h"
#include "myVtkInteractorStyleImage.h"
#include "LeapAbstractionLayer.h"
#include <algorithm>

using namespace Leap;

class SlicerLeapListener : public Listener {
public:

	float _lastZ;
	void SetInterface(LeapAbstractionLayer* lal);
	void onInit(const Controller& controller);
	void onConnect(const Controller& controller);
	void onDisconnect(const Controller& controller);
	void onExit(const Controller& controller);
	void onFrame(const Controller& controller);
	void onFocusGained(const Controller& controller);
	void onFocusLost(const Controller& controller);

protected:
	LeapAbstractionLayer* _lal;
};

#endif