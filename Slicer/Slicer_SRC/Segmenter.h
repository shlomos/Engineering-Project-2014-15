#ifndef _SEGMENTER_
#define _SEGMENTER_

#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>


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

//! Node used to represent pixels (2D) and voxels(3D).
struct Node
{
	int pos_x; //!< voxel x coordinate.
	int pos_y; //!< voxel y coordinate.
	int pos_z; //!< voxel z coordinate.
	//! Struct constructor.
	Node(int x, int y, int z) :pos_x(x), pos_y(y), pos_z(z){}
};

class Segmenter //: public QThread
{
public:
	Segmenter::Segmenter(vtkSmartPointer<vtkStructuredPoints> selection_structured_points, vtkSmartPointer<vtkStructuredPoints> CT_structured_points);
	void doSegmentation3D(int pos_x, int pos_y, int pos_z, int minThreshold, int maxThreshold);
	bool predicate3D(Node& node, int minThreshold, int maxThreshold);
private:

	enum Tasks
	{
		SLEEP,
		SEGMENTATION_2D,
		SEGMENTATION_3D
	} task;

	int mSeedX, mSeedY, mSeedZ;
	int mMinThreshold, mMaxThreshold;
	unsigned int mOrientation;
};



#endif

