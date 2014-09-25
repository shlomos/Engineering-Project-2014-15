#include "SlicerLeapListener.h"

void SlicerLeapListener::SetImageViewer(myVtkInteractorStyleImage* interactorStyle) {
	_InteractorStyle = interactorStyle;
}


void SlicerLeapListener::onInit(const Controller& controller) {
	std::cout << "Initialized" << std::endl;
}

void SlicerLeapListener::onConnect(const Controller& controller) {
	std::cout << "Connected" << std::endl;
	controller.enableGesture(Gesture::TYPE_CIRCLE);
	controller.enableGesture(Gesture::TYPE_KEY_TAP);
	controller.enableGesture(Gesture::TYPE_SCREEN_TAP);
	controller.enableGesture(Gesture::TYPE_SWIPE);

}

void SlicerLeapListener::onDisconnect(const Controller& controller) {
	//Note: not dispatched when running in a debugger.
	std::cout << "Disconnected" << std::endl;
}

void SlicerLeapListener::onExit(const Controller& controller) {
	std::cout << "Exited" << std::endl;
}

void SlicerLeapListener::onFrame(const Controller& controller) {
	// Get the most recent frame and report some basic information
	const Frame frame = controller.frame();
	if (!frame.hands().isEmpty()) {
		// Get the first hand
		const Hand hand = frame.hands()[0];

		// Check if the hand has any fingers
		const FingerList fingers = hand.fingers();
		if (!fingers.isEmpty()) {
			// Calculate the hand's average finger tip position
			Vector avgPos;
			for (int i = 0; i < fingers.count(); ++i) {
				avgPos += fingers[i].tipPosition();
			}
			avgPos /= (float)fingers.count();
			Vector *aPoint = new Vector(0.0f, 0.0f, -150.0f);
			Vector frontmost = frame.pointables().frontmost().tipPosition();
			float realDis = (frontmost.z-aPoint->z);
			//cout << "im sad!"<<endl;
			int dis = std::max(0,std::min(_InteractorStyle->getMaxSlice(),(int)(realDis/150.0f*_InteractorStyle->getMaxSlice())));
			//cout << "im very very sad!"<<endl;
			//std::cout << "try set slice to " << dis<<  std::endl;
			//Assume screen is 15 cm after leap
			//std::cout << "finger is at " << frontmost.z
			//	<< " pointer is " << dis << "cm from screen." << std::endl;
			this->_InteractorStyle->setSlice(dis);
			this->_InteractorStyle->_x_position = frontmost.x;
			this->_InteractorStyle->_y_position = frontmost.y;
			//cout << "im desperate!"<<endl;
		}
	}

	if (!frame.hands().isEmpty()) {
		//std::cout << std::endl;
	}
}

void SlicerLeapListener::onFocusGained(const Controller& controller) {
	std::cout << "Focus Gained" << std::endl;
}

void SlicerLeapListener::onFocusLost(const Controller& controller) {
	std::cout << "Focus Lost" << std::endl;
}