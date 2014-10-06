#include "Segmenter.h"


//! 3D segmentation algorithm. Performs segmentation on a range of slices.
/*!
\param pos_x first dimension.
\param pos_y second dimension.
\param pos_z third dimension (slice number).
\param minThreshold minimum level of intensity that will pass the segmentation.
\param maxThreshold maximum level of intensity that will pass the segmentation.
\param min_Z segmentation lower slice limit.
\param max_Z segmentation upper slice limit.
*/
void Segmenter::doSegmentation3D(int pos_x, int pos_y, int pos_z, int minThreshold, int maxThreshold) 
{
	////set parameters
	//QMutexLocker locker(&mutex); //lock mutex until we go out of scope
	////check thread is not doing work
	//if (task != SLEEP)
	//{
	//	qWarning() << "Cannot do doSegmentation3D(). Thread appears to already be working!";
	//	return;
	//}

	//mSeedX = pos_x;
	//mSeedY = pos_y;
	//mSeedZ = pos_z;
	//mMinThreshold = minThreshold;
	//mMaxThreshold = maxThreshold;
	//task = SEGMENTATION_3D; //Set the thread's task
	//if (!isRunning())
	//	start();
	//else
	//	workToDo.wakeOne(); //The Thread may be sleeping, wake it up!
}

bool Segmenter::predicate3D(Node& node, int minThreshold, int maxThreshold)
{
	////Check if Node is in the current boundary
	//if (!boundaryManager->isInBoundary(node.pos_x, node.pos_y, node.pos_z))
	//	return false;
	//char* pixel_segmentation = static_cast<char*>(imagePairManager->segblock->GetScalarPointer(node.pos_x, node.pos_y, node.pos_z));
	//char pixel_visited = visited3D[node.pos_x][node.pos_y][node.pos_z];
	//short* pixel_original = static_cast<short*>(imagePairManager->original->GetScalarPointer(node.pos_x, node.pos_y, node.pos_z));
	//if (pixel_visited == 1)
	//	return false;//We've already visited this voxel
	//if ((char)pixel_segmentation[0] == imagePairManager->BLOCKING)
	//	return false;
	//if ((short)pixel_original[0] < minThreshold || (short)pixel_original[0] > maxThreshold)
	//	return false;
	////otherwise mark as visited
	//visited3D[node.pos_x][node.pos_y][node.pos_z] = 1;
	////update segmentation results
	//pixel_segmentation[0] = imagePairManager->SEGMENTATION;
	return true;
}