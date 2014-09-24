#ifndef _SLICER_LEAP_LISTENER_
#define _SLICER_LEAP_LISTENER_

#include "Leap.h"
#include "myVtkInteractorStyleImage.h"
#include <algorithm>

using namespace Leap;

class SlicerLeapListener : public Listener {
public:

	myVtkInteractorStyleImage* _InteractorStyle;

	void SetImageViewer(myVtkInteractorStyleImage* interactorStyle);
	void onInit(const Controller& controller);
	void onConnect(const Controller& controller);
	void onDisconnect(const Controller& controller);
	void onExit(const Controller& controller);
	void onFrame(const Controller& controller);
	void onFocusGained(const Controller& controller);
	void onFocusLost(const Controller& controller);
};

#endif