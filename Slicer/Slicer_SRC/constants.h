#ifndef _CONSTANTS_

#define _CONSTANTS_

#define UPDATE_SLICE_TIMER 5
#define SCALE_FACTOR 1.8
#define NUM_OF_LAYERS 1
//#define CROSS_LAYER     (NUM_OF_LAYERS>1)?NUM_OF_LAYERS-1:-1
#define SELECTION_LAYER std::max(0,NUM_OF_LAYERS-2)

#define NOT_ACTIVE 0
#define FOREGROUND 1
#define BACKGROUND 2
#define DEFAULT_DRAW_SIZE 2
#define MIN_DRAW_SIZE 1
#define MAX_DRAW_SIZE 50

// Mu and Sigma on GT. right now is mean of all cases
#define MU         58.404
#define SIGMA      4.0
#define CHICKEN_PI 3.141592

// graphCuts related values
#define LOWER_BOUND         0
#define UPPER_BOUND         104
#define MIN_POSSIBLE_VALUE  3000
#define SILENCING_FACTOR    0.99
#define NEDGES_FACTOR       1000
#define XY_OPENCLOSE        5
#define Z_OPENCLOSE         4
#define DILATION_ST_SIZE	3
#define NUM_ITER_OPENCLOSE  2
#define SMOOTHING_FACTOR_XY 25.0
#define SMOOTHING_FACTOR_Z	25.0
#define NUM_OF_GROW_ITERATIONS 50
// 
#define ISO_VALUE 0.5

#define NOWINDOW	0
#define WINDOW25D	1
#define WINDOW3D	2

#define LEAP_MAX_Y 500
#define LEAP_MIN_X -200

enum
{
	SLICE_ORIENTATION_YZ = 0,
	SLICE_ORIENTATION_XZ = 1,
	SLICE_ORIENTATION_XY = 2
};




#endif