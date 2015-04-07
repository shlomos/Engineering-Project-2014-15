#ifndef _MESH_LEAK_CORRECTOR_
#define _MESH_LEAK_CORRECTOR_


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkImage.h"
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkVTKImageImport.h>
#include <itkVTKImageToImageFilter.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkContourFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkCurvatures.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCharArray.h>
#include <vtkProbeFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkQuadricDecimation.h>
#include <vtkCellLocator.h>
#include <vtkDecimatePro.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageReader.h>
#include <itkImage.h>
#include <itkBinaryFillholeImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

#include <Eigen/Sparse>
#include <unsupported/Eigen/UmfPackSupport>
#include <unsupported/Eigen/SparseExtra>
//#include <SuiteSparse/UMFPACK/Include/umfpack.h>
#include <exception>

#include "itkImageToVTKImageFilter.h"
#include "graph.h"

#define NPRINT_LINE

#ifndef NPRINT_LINE
#define printLine() std::cout<<"line: "<<__LINE__<<"\n"
#else
#define printLine()
#endif

#define MESH_SMOOTH_ITERATIONS 200
#define GRAD_THRESH 3

#define GRADIENT_THRESHOLD 35
#define DECIMATION_FACTOR 0.5

#define ATTRIBUTE_DILATION_RADIUS 0

enum{ CORRECT_SEED = 2, LEAK_SEED = 3, INTERIOR_SEED = 4 };
enum{ INPUT_IMAGE_NAME = 1, SEG_INPUT_IMAGE_NAME, SEED_INPUT_IMAGE_NAME, INTERIOR_SEED_NAME, OUTPUT_CURVATURE_MESH, GRAPH_SIGMA };
enum{ MIN_CURVATURE_TAG = 3, MAX_CURVATURE_TAG };

class MeshLeaksCorrector{

public:

	template<class InputImageType>
	typename InputImageType::Pointer read3DImage(const char* imageName)
	{
		typedef itk::Image< double, 3 >         ImageType;
		typedef itk::ImageFileReader<InputImageType> ReaderType;

		ReaderType::Pointer reader = ReaderType::New();
		//cout << imageName << endl;
		reader->SetFileName(imageName);

		try
		{
			reader->Update();
		}
		catch (itk::ExceptionObject & excp)
		{
			std::cerr << "reading input image exception thrown" << std::endl;
			std::cerr << excp << std::endl;
			exit(1);
		}

		
		return reader->GetOutput();

	}

	template< class ItkImageType>
	vtkSmartPointer<vtkImageData> convertItkImageToVtkImage(typename ItkImageType::Pointer itkImage)
	{

		typedef itk::ImageToVTKImageFilter<ItkImageType> ConverterType;
		typename ConverterType::Pointer converter = ConverterType::New();
		converter->SetInput(itkImage);
		converter->Update();
		vtkSmartPointer<vtkImageData> vtkImage = converter->GetOutput();
		return vtkImage;

	}
	vtkSmartPointer<vtkPolyData> createAndSmoothSurface(vtkSmartPointer<vtkImageData> vtkBinaryImage, unsigned int smoothIteations);
	vtkSmartPointer<vtkPolyData> polyDataToMinCurvature(vtkSmartPointer<vtkPolyData> mesh);
	vtkSmartPointer<vtkPolyData> polyDataToMaxCurvature(vtkSmartPointer<vtkPolyData> mesh);
	void writePolyData(vtkSmartPointer< vtkPolyData> mesh, const char *fileName);
	vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);
	vtkSmartPointer<vtkIdList> GetConnectedTriangles(vtkSmartPointer<vtkPolyData> mesh, int id);
	vtkSmartPointer<vtkPolyData> minCut(vtkSmartPointer<vtkPolyData> curvatureMesh, vtkSmartPointer<vtkPolyData> seedMesh
		, vtkSmartPointer<vtkPolyData> gradientMesh, int curvatureTag, double graphSigma);
	vtkSmartPointer<vtkPolyData> minCutConjunction(vtkSmartPointer<vtkPolyData> minCutMeshLeaks,
		vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks);
	void attributeDilation(vtkSmartPointer<vtkPolyData> minCutMesh, int dilationRadius);

	//this function interpolate the mesh coordinates, given the outer surface coordinates
	vtkSmartPointer<vtkPolyData> laplaceInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh);

	//this function interpolates the normals vectors, given the outer surface known normals
	vtkSmartPointer<vtkDoubleArray> laplaceNormalsInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh);

	void computeTriangleCenter(double vertices[3][3], double center[3]);

	void computeTriangleNormal(double vertices[3][3], double normal[3]);

	//computes the weighted laplace based on the normals values
	vtkSmartPointer<vtkPolyData> laplaceWeightedInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh);

	void computeTriangleNewCoordinates(double oldVerticesCoordinates[3][3], double center[3], double oldNormal[3], double newNormal[3]
		, double newVerticesCoordinates[3][3]);

	void computeCellNormals(double vertex1Normals[3], double vertex2Normals[3], double vertex3Normals[3], double cellNormal[3]);

	void computeLaplaceOperator(vtkSmartPointer<vtkPolyData> mesh);

	//interpolates the laplacian values over unknown domain, given its values over
	//the known outer surface
	vtkSmartPointer<vtkDoubleArray> laplaceLaplaceInterpolation(vtkSmartPointer<vtkPolyData> mesh);


	vtkSmartPointer<vtkPolyData> poissonLaplaceInterpolation(vtkSmartPointer<vtkPolyData> mesh);


	//interpolate normals values, using known boundary
	vtkSmartPointer<vtkPolyData> normalsLinearInterpolation(vtkSmartPointer<vtkPolyData> mesh);

	template<class ImageType>
	vtkSmartPointer<vtkPolyData> sampleImageOnMesh(vtkSmartPointer<vtkPolyData> mesh, typename ImageType::Pointer image)
	{
		typedef typename itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
		IteratorType it(image, image->GetLargestPossibleRegion());


		//create scalar function array for gradient magnitude
		vtkSmartPointer<vtkCharArray> seedsArray =
			vtkSmartPointer<vtkCharArray>::New();
		seedsArray->SetName("seedsArray");
		seedsArray->SetNumberOfComponents(1);
		seedsArray->SetNumberOfTuples(mesh->GetNumberOfPoints());


		for (int i = 0; i < mesh->GetNumberOfPoints(); ++i)
		{
			double coordinate[3];
			mesh->GetPoint(i, coordinate);
			typename ImageType::PointType imagePoint;
			imagePoint[0] = coordinate[0];
			imagePoint[1] = coordinate[1];
			imagePoint[2] = coordinate[2];

			typename ImageType::IndexType imageIndex;
			image->TransformPhysicalPointToIndex(imagePoint, imageIndex);

			it.SetIndex(imageIndex);
			seedsArray->InsertTuple1(i, it.Value());

		}

		mesh->GetPointData()->AddArray(seedsArray);
		return mesh;

	}

	template<class ImageType>
	typename ImageType::Pointer sampleMeshOnImage(vtkSmartPointer<vtkPolyData> mesh, typename ImageType::Pointer image)
	{

		// Create the tree
		vtkSmartPointer<vtkCellLocator> cellLocator =
			vtkSmartPointer<vtkCellLocator>::New();
		cellLocator->SetDataSet(mesh);
		cellLocator->BuildLocator();


		typename ImageType::Pointer outputImage = ImageType::New();
		outputImage->SetOrigin(image->GetOrigin());
		outputImage->SetSpacing(image->GetSpacing());
		outputImage->SetRegions(image->GetLargestPossibleRegion());
		outputImage->Allocate();
		outputImage->FillBuffer(0);

		typedef typename itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
		IteratorType outIt(outputImage, outputImage->GetLargestPossibleRegion());
		IteratorType inIt(image, image->GetLargestPossibleRegion());
		for (inIt.GoToBegin(), outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt, ++inIt)
		{
			typename ImageType::PointType imagePoint;
			outputImage->TransformIndexToPhysicalPoint(outIt.GetIndex(), imagePoint);
			double testPoint[3] = { imagePoint[0], imagePoint[1], imagePoint[2] };

			if (inIt.Value() != 0)
			{
				//Find the closest points to TestPoint
				double closestPoint[3];//the coordinates of the closest point will be returned here
				double closestPointDist2; //the squared distance to the closest point will be returned here
				vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
				int subId; //this is rarely used (in triangle strips only, I believe)
				cellLocator->FindClosestPoint(testPoint, closestPoint, cellId, subId, closestPointDist2);
				if (closestPointDist2 <= 1)
				{
					outIt.Set(1);
				}
			}
		}

		///////////////

		// typename ImageType::Pointer outputImage = ImageType::New();
		// outputImage->SetOrigin(image->GetOrigin());
		// outputImage->SetSpacing(image->GetSpacing());
		// outputImage->SetRegions(image->GetLargestPossibleRegion());
		// outputImage->Allocate();
		// outputImage->FillBuffer(0);

		// typedef typename itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
		// IteratorType it(outputImage, outputImage->GetLargestPossibleRegion());


		// for( int i = 0; i < mesh->GetNumberOfPoints(); ++i)
		// {
		//     double coordinate[3];
		//     mesh->GetPoint(i , coordinate);
		//     typename ImageType::PointType imagePoint;
		//     imagePoint[0] = coordinate[0];
		//     imagePoint[1] = coordinate[1];
		//     imagePoint[2] = coordinate[2];

		//     typename ImageType::IndexType imageIndex;
		//     image->TransformPhysicalPointToIndex(imagePoint, imageIndex);

		//     it.SetIndex(imageIndex);
		//     it.Set(1);
		// }

		return outputImage;
	}
	typedef unsigned short PixelType;
	template<class ImageType>
	typename ImageType::Pointer correctImage(typename ImageType::Pointer leakyImage, typename ImageType::Pointer seedImage, typename ImageType::Pointer contourImage, int numTumors)
	{
		typedef typename itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
		IteratorType leakIt(leakyImage, leakyImage->GetLargestPossibleRegion());
		IteratorType contourIt(contourImage, contourImage->GetLargestPossibleRegion());

		for (leakIt.GoToBegin(), contourIt.GoToBegin(); !leakIt.IsAtEnd(); ++leakIt, ++contourIt)
		{
			if (contourIt.Value() == 1)
			{
				leakIt.Set(0);
			}
			if (leakIt.Value() == 2){
				leakIt.Set(0);
			}

		}

		//retain only thresholded pixels that are connected to the correction seed
		typedef typename itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectorType;

		typename ConnectorType::Pointer connector = ConnectorType::New();
		connector->SetInput(leakyImage);
		cout << "cc: bg: " << connector->GetBackgroundValue() << endl;
		connector->SetBackgroundValue(NOT_ACTIVE);
		connector->SetFullyConnected(true); // was flase
		connector->Update();

		std::cout << "Number of objects: " << connector->GetObjectCount() << std::endl;

		typedef itk::LabelShapeKeepNObjectsImageFilter< ImageType > LabelShapeKeepNObjectsImageFilterType;
		LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
		labelShapeKeepNObjectsImageFilter->SetInput(connector->GetOutput());
		labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
		labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(numTumors);
		labelShapeKeepNObjectsImageFilter->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
		labelShapeKeepNObjectsImageFilter->Update();

		return labelShapeKeepNObjectsImageFilter->GetOutput();
	}

};

#endif