#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

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
#include<vtkPolyDataNormals.h>
#include<vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <Eigen/Sparse>
#include <unsupported/Eigen/UmfPackSupport>
#include <unsupported/Eigen/SparseExtra>
//#include <SuiteSparse/UMFPACK/Include/umfpack.h>

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

enum{CORRECT_SEED = 2,LEAK_SEED=3,INTERIOR_SEED=4};
enum{INPUT_IMAGE_NAME = 1, SEG_INPUT_IMAGE_NAME,SEED_INPUT_IMAGE_NAME,INTERIOR_SEED_NAME,OUTPUT_CURVATURE_MESH,GRAPH_SIGMA};
enum{MIN_CURVATURE_TAG=3,MAX_CURVATURE_TAG};


template<class InputImageType>
typename InputImageType::Pointer read3DImage( const char* imageName ) 
{
    
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(imageName);
    
    try
    {
        reader->Update();
    }
    catch (...)
    {
        std::cerr << "reading input image exception thrown"
                  << std::endl;
        exit(1);
    }
    typename InputImageType::Pointer inputImage = reader->GetOutput();
    return inputImage;
    
}


template< class ItkImageType>
vtkSmartPointer<vtkImageData> convertItkImageToVtkImage( typename ItkImageType::Pointer itkImage )
{

    typedef itk::ImageToVTKImageFilter<ItkImageType> ConverterType;
    typename ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput(itkImage);
    converter->Update();
    vtkSmartPointer<vtkImageData> vtkImage = converter->GetOutput();
    return vtkImage;

}



vtkSmartPointer<vtkPolyData> createAndSmoothSurface( vtkSmartPointer<vtkImageData> vtkBinaryImage, unsigned int smoothIteations )
{

    vtkSmartPointer<vtkContourFilter> vtkContourFilter = vtkContourFilter::New();
    vtkContourFilter->SetInput(vtkBinaryImage);
    vtkContourFilter->SetValue(1,1);
    vtkContourFilter->Update();

    vtkSmartPointer<vtkPolyData> vtkContourMesh = vtkContourFilter->GetOutput();

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmoothPolyDataFilter::New();
    smoother->SetInput(vtkContourMesh);
    smoother->SetNumberOfIterations(smoothIteations);
    smoother->Update();
    return smoother->GetOutput();
    
}



vtkSmartPointer<vtkPolyData> polyDataToMinCurvature( vtkSmartPointer<vtkPolyData> mesh )
{

    vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
    curvaturesFilter->SetInput(mesh);
    curvaturesFilter->SetCurvatureTypeToMinimum();
    curvaturesFilter->Update();
    vtkSmartPointer<vtkPolyData> curvatureMesh = curvaturesFilter->GetOutput();
    return curvatureMesh;
    
}

vtkSmartPointer<vtkPolyData> polyDataToMaxCurvature( vtkSmartPointer<vtkPolyData> mesh )
{

    vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
    curvaturesFilter->SetInput(mesh);
    curvaturesFilter->SetCurvatureTypeToMaximum();
    curvaturesFilter->Update();
    vtkSmartPointer<vtkPolyData> curvatureMesh = curvaturesFilter->GetOutput();
    return curvatureMesh;
    
}



void writePolyData( vtkSmartPointer< vtkPolyData> mesh, const char *fileName )
{

    vtkSmartPointer<vtkPolyDataWriter> polydataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polydataWriter->SetInput(mesh);
    polydataWriter->SetFileName(fileName);
    polydataWriter->Update();

}


vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
{
    vtkSmartPointer<vtkIdList> connectedVertices =
        vtkSmartPointer<vtkIdList>::New();
    
    //get all cells that vertex 'id' is a part of
    vtkSmartPointer<vtkIdList> cellIdList =
        vtkSmartPointer<vtkIdList>::New();
    
    mesh->GetPointCells(id, cellIdList);
   
    
    for(unsigned int i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {
        //std::cout << "id " << id << " connect to " << cellIdList->GetId(i) << std::endl;
        
        vtkSmartPointer<vtkIdList> pointIdList =
            vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
        
        //std::cout << "Point id: "<<id<<std::endl;//id<<" End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;
        
        if(pointIdList->GetId(0) != id)
        {
            //std::cout << "Connected to " << pointIdList->GetId(0) << endl;
            connectedVertices->InsertNextId(pointIdList->GetId(0));
        }
        else
        {
            //std::cout << "Connected to " << pointIdList->GetId(1) << endl;
            connectedVertices->InsertNextId(pointIdList->GetId(1));
        }
    }
    //cellIdList->Delete();
    return connectedVertices;
}
vtkSmartPointer<vtkIdList> GetConnectedTriangles(vtkSmartPointer<vtkPolyData> mesh, int id)
{
    vtkSmartPointer<vtkIdList> connectedVertices =
        vtkSmartPointer<vtkIdList>::New();
    
    //get all cells that vertex 'id' is a part of
    vtkSmartPointer<vtkIdList> cellIdList =
        vtkSmartPointer<vtkIdList>::New();
    
    mesh->GetPointCells(id, cellIdList);
   
    
    for(unsigned int i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {
        
        vtkSmartPointer<vtkIdList> pointIdList =
            vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
                
        if(pointIdList->GetId(0) != id)
        {
            connectedVertices->InsertNextId(pointIdList->GetId(0));
        }
        if(pointIdList->GetId(1) != id)
        {
            connectedVertices->InsertNextId(pointIdList->GetId(1));
        }
        if(pointIdList->GetId(2) != id)
        {
            connectedVertices->InsertNextId(pointIdList->GetId(2));
        }
       
    }
    return connectedVertices;
}



// vtkSmartPointer<vtkPolyData> minCut( vtkSmartPointer<vtkPolyData> curvatureMesh ,vtkSmartPointer<vtkPolyData> seedMesh)
// {	

        
   
    
//     int num_of_nodes = curvatureMesh->GetNumberOfPoints();
//     int num_of_edges = num_of_nodes * 8;//curvatureMesh->GetNumberOfCells();
//     std::cout<<num_of_nodes<<"\n"<<num_of_edges<<"\n";

//     typedef double GraphNodeType;
    
//     GraphNodeType edge_weight = 0;
//     GraphNodeType source_weight = 0.5;
//     GraphNodeType sink_weight = 0.5;
    
//     typedef Graph<GraphNodeType,GraphNodeType,GraphNodeType> GraphType;
//     GraphType *g = new GraphType( num_of_nodes,  num_of_edges); 
//     g->add_node (num_of_nodes);

//     double currentSeedLabel = 0;
   
//     for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
//     {
//         double neighborSeedLabel = CORRECT_SEED;
        
//         seedMesh->GetPointData()->GetScalars("seedsArray")->GetTuple(meshCounter,&currentSeedLabel);
//         // special user indicated seeds
//         if( currentSeedLabel == LEAK_SEED)
//         {
//             g -> add_tweights(meshCounter, 0, 1000000 );
//         }
//         else if(currentSeedLabel == CORRECT_SEED)
//         {
      
//             g -> add_tweights(meshCounter, 100000, 0 );
          
//         }
        
//         //no user indacted seeds
//         else
//         {
//             g -> add_tweights(meshCounter, source_weight, sink_weight );
//         }
//         vtkSmartPointer<vtkIdList> neighborsList = GetConnectedVertices(curvatureMesh, meshCounter );

//         double currentMinCurvature, neighborMinCurvature;
//         curvatureMesh->GetPointData()->GetScalars()->GetTuple(meshCounter,&currentMinCurvature);
//         for(int idsCounter = 0;idsCounter < neighborsList->GetNumberOfIds(); ++ idsCounter)
//         {
//             curvatureMesh->GetPointData()->GetScalars()->GetTuple(neighborsList->GetId(idsCounter),&neighborMinCurvature);

//             double edge_weight = 1/( 1+exp(-12*(currentMinCurvature+neighborMinCurvature)));
//             g -> add_edge( meshCounter, neighborsList->GetId(idsCounter), edge_weight, edge_weight);            
//         }
        
//     }
   
//     g -> maxflow();
//     std::cout << "Done!"<< std::endl;
    
    
//     //create scalar function array for gradient magnitude
//     vtkSmartPointer<vtkCharArray> minCutArray =
//         vtkSmartPointer<vtkCharArray>::New();
//     minCutArray->SetName("minCutArray");
//     minCutArray->SetNumberOfComponents(1);
//     minCutArray->SetNumberOfTuples(curvatureMesh->GetNumberOfPoints());
    
    
    
//     for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
//     {
//         minCutArray->InsertTuple1(meshCounter,g->what_segment(meshCounter));

//     }
//     curvatureMesh->GetPointData()->AddArray(minCutArray);

//     delete g;
//     return curvatureMesh;
// }


vtkSmartPointer<vtkPolyData> minCut( vtkSmartPointer<vtkPolyData> curvatureMesh ,vtkSmartPointer<vtkPolyData> seedMesh
                                     ,vtkSmartPointer<vtkPolyData> gradientMesh, int curvatureTag,double graphSigma)
{	

        
   
    
    int num_of_nodes = curvatureMesh->GetNumberOfPoints();
    int num_of_edges = num_of_nodes * 8;//curvatureMesh->GetNumberOfCells();
    std::cout<<num_of_nodes<<"\n"<<num_of_edges<<"\n";

    typedef double GraphNodeType;
    
    GraphNodeType edge_weight = 0;
    GraphNodeType source_weight = 0.5;
    GraphNodeType sink_weight = 0.5;
    
    typedef Graph<GraphNodeType,GraphNodeType,GraphNodeType> GraphType;
    GraphType *g = new GraphType( num_of_nodes,  num_of_edges); 
    g->add_node (num_of_nodes);

    double currentSeedLabel = 0;
   
    for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
    {
        double neighborSeedLabel = CORRECT_SEED;
        
        seedMesh->GetPointData()->GetScalars("seedsArray")->GetTuple(meshCounter,&currentSeedLabel);
        // special user indicated seeds
        if( currentSeedLabel == curvatureTag)
        {
            g -> add_tweights(meshCounter, 0, 1000000 );
        }
        else if(currentSeedLabel == CORRECT_SEED)
        {
            
            g -> add_tweights(meshCounter, 100000, 0 );
            
        }
        
        //no user indacted seeds
        else
        {
            g -> add_tweights(meshCounter, source_weight, sink_weight );
        }
        vtkSmartPointer<vtkIdList> neighborsList = GetConnectedVertices(curvatureMesh, meshCounter );
        
        double currentCurvature, neighborCurvature;
        double currentGradient, neighborGradient;
        curvatureMesh->GetPointData()->GetScalars()->GetTuple(meshCounter,&currentCurvature);
        gradientMesh->GetPointData()->GetScalars()->GetTuple(meshCounter,&currentGradient);
        for(int idsCounter = 0;idsCounter < neighborsList->GetNumberOfIds(); ++ idsCounter)
        {
            curvatureMesh->GetPointData()->GetScalars()->GetTuple(neighborsList->GetId(idsCounter),&neighborCurvature);
            gradientMesh->GetPointData()->GetScalars()->GetTuple(neighborsList->GetId(idsCounter),&neighborGradient);

            double sigma = 30;
            double curvature_edge_weight;
            if(curvatureTag == MIN_CURVATURE_TAG)
            {
                curvature_edge_weight = 1/( 1+exp(-1*graphSigma*(currentCurvature+neighborCurvature)));
            }
            else
            {
                curvature_edge_weight = 1/( 1+exp(graphSigma*(currentCurvature+neighborCurvature)));

            }
            double gradient_edge_weight = 1/(   exp(  pow((currentGradient-neighborGradient),2)/sigma    )             );

            //double alpha = 0.0;
            //edge_weight = alpha * curvature_edge_weight + (1 - alpha) *
            //gradient_edge_weight;

            // if(curvature_edge_weight < gradient_edge_weight)
            //{
            edge_weight = curvature_edge_weight;
            //}
            //else
            // {
            //   edge_weight = gradient_edge_weight;
            //}
            g -> add_edge( meshCounter, neighborsList->GetId(idsCounter), edge_weight, edge_weight);            
        }
        
    }
   
    g -> maxflow();
    std::cout << "Done!"<< std::endl;
    
    
    //create scalar function array for gradient magnitude
    vtkSmartPointer<vtkCharArray> minCutArray =
        vtkSmartPointer<vtkCharArray>::New();
    if(curvatureTag == MIN_CURVATURE_TAG)
    {
        minCutArray->SetName("minCutLeaksArray");
    }
    else
    {
        minCutArray->SetName("minCutInteriorLeaksArray");
    }
    
    minCutArray->SetNumberOfComponents(1);
    minCutArray->SetNumberOfTuples(curvatureMesh->GetNumberOfPoints());
    
    
    
    for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
    {
        minCutArray->InsertTuple1(meshCounter,g->what_segment(meshCounter));

    }
    curvatureMesh->GetPointData()->AddArray(minCutArray);

    delete g;
    return curvatureMesh;
}



vtkSmartPointer<vtkPolyData> minCutConjunction(vtkSmartPointer<vtkPolyData> minCutMeshLeaks,
                                               vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks)
{
    vtkSmartPointer<vtkCharArray> minCutArray =
        vtkSmartPointer<vtkCharArray>::New();
    minCutArray->SetName("minCutArray");
   
    minCutArray->SetNumberOfComponents(1);
    minCutArray->SetNumberOfTuples(minCutMeshLeaks->GetNumberOfPoints());
    
    unsigned int num_of_nodes = minCutMeshLeaks->GetNumberOfPoints();
    for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
    {
        double label1, label2;
        minCutMeshLeaks->GetPointData()->GetScalars("minCutLeaksArray")->GetTuple(meshCounter,&label1);
        minCutMeshInteriorLeaks->GetPointData()->GetScalars("minCutInteriorLeaksArray")->GetTuple(meshCounter,&label2);

        if(label1 != 0 || label2 != 0)
        {
            minCutArray->InsertTuple1(meshCounter,1);
        }
        else
        {
            minCutArray->InsertTuple1(meshCounter,0);
        }

    }
    minCutMeshLeaks->GetPointData()->AddArray(minCutArray);
    return minCutMeshLeaks;

}

void attributeDilation( vtkSmartPointer<vtkPolyData> minCutMesh , int dilationRadius)
{
   
    vtkSmartPointer<vtkCharArray> minCutDilatedArray = vtkCharArray::SafeDownCast(minCutMesh->GetPointData()->GetArray("minCutArray"));
    
    int num_of_nodes = minCutMesh->GetNumberOfPoints();
   
    for(int i = 0;i <= dilationRadius;++i)
    {
        double  currentLabel,neighborLabel;
        for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
        {
            
            
            minCutMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&currentLabel);
            // special user indicated seeds
            if( currentLabel == 1)
            {
                vtkSmartPointer<vtkIdList> neighborsList = GetConnectedVertices(minCutMesh, meshCounter );
                for(int idsCounter = 0;idsCounter < neighborsList->GetNumberOfIds(); ++ idsCounter)
                {
                    minCutMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(neighborsList->GetId(idsCounter),&neighborLabel);
                    
                    if(neighborLabel == 0)
                    {
                        minCutDilatedArray->InsertTuple1(neighborsList->GetId(idsCounter),2);
                        
                    }
                }
                
            }
            
        }
        for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
        {
            minCutMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&currentLabel);
            
            if(currentLabel == 2)
            {
                minCutDilatedArray->InsertTuple1(meshCounter,1);
            }
        }
               
    }
}


// vtkSmartPointer<vtkPolyData> minCutConjunction(vtkSmartPointer<vtkPolyData> minCutMeshLeaks,
//                                                vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks)
// {
//     vtkSmartPointer<vtkCharArray> minCutArray =
//         vtkSmartPointer<vtkCharArray>::New();
//     minCutArray->SetName("minCutArray");
   
//     minCutArray->SetNumberOfComponents(1);
//     minCutArray->SetNumberOfTuples(minCutMeshLeaks->GetNumberOfPoints());
    
//     unsigned int num_of_nodes = minCutMeshLeaks->GetNumberOfPoints();
//     for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
//     {
//         double label1, label2;
//         minCutMeshLeaks->GetPointData()->GetScalars("minCutLeaksArray")->GetTuple(meshCounter,&label1);
//         minCutMeshInteriorLeaks->GetPointData()->GetScalars("minCutInteriorLeaksArray")->GetTuple(meshCounter,&label2);

//         if(label1 != 0 || label2 != 0)
//         {
//             minCutArray->InsertTuple1(meshCounter,1);
//         }
//         else
//         {
//             minCutArray->InsertTuple1(meshCounter,0);
//         }

//     }
//     minCutMeshLeaks->GetPointData()->AddArray(minCutArray);
//     return minCutMeshLeaks;

   


// }


//this function interpolate the mesh coordinates, given the outer surface coordinates
vtkSmartPointer<vtkPolyData> laplaceInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh)
{
    int num_of_nodes = curvatureMesh->GetNumberOfPoints();
    // Create the sparse matrix
    std::vector<int> linearSystemIDs;

    writePolyData(curvatureMesh, "minCutMesh.vtk");

    for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
    {
        double label;
        curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);

        if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
        {
            linearSystemIDs.push_back(meshCounter);
        }
    }
    unsigned int numberOfVariables = linearSystemIDs.size();
    std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
    if (numberOfVariables == 0)
    {
        std::cerr << "No masked pixels found" << std::endl;
        exit(1);
    }

    // Create the reverse mapping from pixel to variable id
    std::map<int,int> meshIDs;//(numberOfVariables);
    for(unsigned int i = 0; i < numberOfVariables; i++)
    {
        meshIDs[linearSystemIDs[i]] = i;
    }
   
    //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    
    Eigen::SparseMatrix<double> A(numberOfVariables, numberOfVariables);
    Eigen::VectorXd b_x(numberOfVariables);
    Eigen::VectorXd b_y(numberOfVariables);
    Eigen::VectorXd b_z(numberOfVariables);
   
    int currentId=linearSystemIDs[2];
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        currentId = linearSystemIDs.at(variableId);
        printLine();
        //std::cout<<"\nsize: "<<linearSystemIDs.size()<<"\ncurrent id: "<<currentId<<"\n";

       
        vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
        printLine();
        connectedVertices = GetConnectedVertices(curvatureMesh, currentId);
        printLine();
        
        // The right hand side of the equation starts equal to 0
        double bvalue_x = 0.0;
        double bvalue_y = 0.0;
        double bvalue_z = 0.0;
        
        int cellNum = connectedVertices->GetNumberOfIds();
        printLine();
        // Loop over the kernel around the current pixel
        for(unsigned int offset = 0; offset < cellNum ; ++offset)
        {
            
            
            printLine();
            double neighborLabel;
            //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
            curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
           
           
            printLine();
        
            if (neighborLabel == 1)
            {
                std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                A.insert(variableId, neighborIdP->second) = 1;
            }
            else
            {
                double coordinate[3]; 
                curvatureMesh->GetPoint(connectedVertices->GetId(offset), coordinate);             
                bvalue_x -= coordinate[0];
                bvalue_y -= coordinate[1];
                bvalue_z -= coordinate[2];
                
            }
           
            printLine();
        }
        A.insert(variableId, variableId) = -1 * cellNum;
        b_x[variableId] = bvalue_x;
        b_y[variableId] = bvalue_y;
        b_z[variableId] = bvalue_z;
      
    }
    std::cout<<"after matrix filling\n";
      // Solve the system with Eigen
    Eigen::VectorXd x(numberOfVariables);
    Eigen::VectorXd y(numberOfVariables);
    Eigen::VectorXd z(numberOfVariables);

    
   
    Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
    Eigen::VectorXd b2_x = A.adjoint() * b_x;
    Eigen::VectorXd b2_y = A.adjoint() * b_y;
    Eigen::VectorXd b2_z = A.adjoint() * b_z;
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_A(A2);
    if(!lu_of_A.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_x,&x))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_y,&y))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_z,&z))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }//     return bonesVolume;



    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        double coordinate[3];
        coordinate[0] = x[variableId];
        coordinate[1] = y[variableId];
        coordinate[2] = z[variableId];
        int currentId = linearSystemIDs.at(variableId);
        curvatureMesh->GetPoints()->SetPoint(currentId, coordinate);
    }
    return curvatureMesh;
}

//this function interpolates the normals vectors, given the outer surface known normals
vtkSmartPointer<vtkDoubleArray> laplaceNormalsInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh)
{
    std::cout<<"in laplace normals\n"<<"line: "<<__LINE__<<"\n";
    int num_of_nodes = curvatureMesh->GetNumberOfPoints();
    // Create the sparse matrix
    std::vector<int> linearSystemIDs;

    
    for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
    {
        double label;
        curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);

        if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
        {
            linearSystemIDs.push_back(meshCounter);
        }
    }
    unsigned int numberOfVariables = linearSystemIDs.size();
    std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
    if (numberOfVariables == 0)
    {
        std::cerr << "No masked pixels found" << std::endl;
        exit(1);
    }

    // Create the reverse mapping from pixel to variable id
    std::map<int,int> meshIDs;//(numberOfVariables);
    for(unsigned int i = 0; i < numberOfVariables; i++)
    {
        meshIDs[linearSystemIDs[i]] = i;
    }
   
    //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    
    Eigen::SparseMatrix<double> A(numberOfVariables, numberOfVariables);
    Eigen::VectorXd b_x(numberOfVariables);
    Eigen::VectorXd b_y(numberOfVariables);
    Eigen::VectorXd b_z(numberOfVariables);


    vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
    normals->SetInput(curvatureMesh);
    normals->ComputePointNormalsOn();
    normals->Update();

    vtkSmartPointer<vtkPolyData> normalsMesh = vtkPolyData::New();
    normalsMesh = normals->GetOutput();
    writePolyData(normalsMesh, "unChangedNormalsMesh.vtk");
    // vtkSmartPointer<vtkDoubleArray> normals = vtkDoubleArray::New();
    
    // normals = vtkDoubleArray::SafeDownCast(curvatureMesh->GetPointData()->GetNormals());
    
    int currentId=linearSystemIDs[2];
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        currentId = linearSystemIDs.at(variableId);
        //printLine();
        //std::cout<<"\nsize: "<<linearSystemIDs.size()<<"\ncurrent id: "<<currentId<<"\n";

       
        vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
        //printLine();
        connectedVertices = GetConnectedVertices(curvatureMesh, currentId);
        //printLine();
        
        // The right hand side of the equation starts equal to 0
        double bvalue_x = 0.0;
        double bvalue_y = 0.0;
        double bvalue_z = 0.0;
        
        int cellNum = connectedVertices->GetNumberOfIds();
        //printLine();
        // Loop over the kernel around the current pixel
        for(unsigned int offset = 0; offset < cellNum ; ++offset)
        {
            
            
            //printLine();
            double neighborLabel;
            //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
            curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
           

            //printLine();
        
            if (neighborLabel == 1)
            {
                std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                A.insert(variableId, neighborIdP->second) = 1;
            }
            else
            {

                double normal[3];
                //normals->GetTuple(connectedVertices->GetId(offset),normal);
                normalsMesh->GetPointData()->GetNormals()->GetTuple(connectedVertices->GetId(offset), normal);
                bvalue_x -= normal[0];
                bvalue_y -= normal[1];
                bvalue_z -= normal[2];

            }
           
            //printLine();
        }
        A.insert(variableId, variableId) = -1 * cellNum;
        b_x[variableId] = bvalue_x;
        b_y[variableId] = bvalue_y;
        b_z[variableId] = bvalue_z;
      
    }
    std::cout<<"after matrix filling\n";
      // Solve the system with Eigen
    Eigen::VectorXd x(numberOfVariables);
    Eigen::VectorXd y(numberOfVariables);
    Eigen::VectorXd z(numberOfVariables);

    
   
    Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
    Eigen::VectorXd b2_x = A.adjoint() * b_x;
    Eigen::VectorXd b2_y = A.adjoint() * b_y;
    Eigen::VectorXd b2_z = A.adjoint() * b_z;
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_A(A2);
    if(!lu_of_A.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_x,&x))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_y,&y))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_z,&z))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }//     return bonesVolume;

    // Set point normals
    vtkSmartPointer<vtkDoubleArray> pointNormalsArray = 
        vtkSmartPointer<vtkDoubleArray>::New();
    pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
    pointNormalsArray->SetNumberOfTuples(num_of_nodes);
    pointNormalsArray->SetName("interpolatedNoramlsArray");
    for(int i=0;i<num_of_nodes;++i)
    {
        double normal[3];
        normalsMesh->GetPointData()->GetNormals()->GetTuple(i, normal);
        pointNormalsArray->SetTuple(i,normal);
    }
    
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        double normal[3];
        normal[0] = x[variableId];
        normal[1] = y[variableId];
        normal[2] = z[variableId];

        double normalNorm;
        
        if (normal[0] != 0 && normal[1] != 0 && normal[2] != 0)
        {
            normalNorm = sqrt( pow(normal[0],2) + pow(normal[1], 2) + pow(normal[2],2) );
            normal[0] /= normalNorm;
            normal[1] /= normalNorm;
            normal[2] /= normalNorm;
        }

        // Add the data to the normals array
        //        pointNormalsArray->SetTuple(variableId, normal) ;
        int currentId = linearSystemIDs.at(variableId);
        pointNormalsArray->SetTuple(currentId, normal) ;
        
        //normalsMesh->GetPointData()->GetNormals()->SetTuple(currentId, normal);

    }
    curvatureMesh->GetPointData()->AddArray(pointNormalsArray);
  
    return pointNormalsArray;
}

void computeTriangleCenter(double vertices[3][3], double center[3])
{
    center[0] = 1/3 * vertices[0][0] + vertices[0][1] + vertices[0][2];
    center[1] = 1/3 * vertices[1][0] + vertices[1][1] + vertices[1][2];
    center[2] = 1/3 * vertices[2][0] + vertices[2][1] + vertices[2][2];
}

void computeTriangleNormal(double vertices[3][3],double normal[3])
{
    double vec1[3], vec2[3];
    vec1[0] = vertices[1][0] - vertices[0][0];
    vec1[1] = vertices[1][1] - vertices[0][1];
    vec1[2] = vertices[1][2] - vertices[0][2];

    vec2[0] = vertices[2][0] - vertices[0][0];
    vec2[1] = vertices[2][1] - vertices[0][1];
    vec2[2] = vertices[2][2] - vertices[0][2];

    vtkMath::Cross(vec1, vec2,normal);
    vtkMath::Normalize(normal);

    
}

//computes the weighted laplace based on the normals values
vtkSmartPointer<vtkPolyData> laplaceWeightedInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh)
{
    int num_of_nodes = curvatureMesh->GetNumberOfPoints();
    // Create the sparse matrix
    std::vector<int> linearSystemIDs;

    writePolyData(curvatureMesh, "minCutMesh.vtk");

    for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
    {
        double label;
        curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);

        if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
        {
            linearSystemIDs.push_back(meshCounter);
        }
    }
    unsigned int numberOfVariables = linearSystemIDs.size();
    std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
    if (numberOfVariables == 0)
    {
        std::cerr << "No masked pixels found" << std::endl;
        exit(1);
    }

    // Create the reverse mapping from pixel to variable id
    std::map<int,int> meshIDs;//(numberOfVariables);
    for(unsigned int i = 0; i < numberOfVariables; i++)
    {
        meshIDs[linearSystemIDs[i]] = i;
    }
   
    //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    
    Eigen::SparseMatrix<double> Ax(numberOfVariables, numberOfVariables);
    Eigen::SparseMatrix<double> Ay(numberOfVariables, numberOfVariables);
    Eigen::SparseMatrix<double> Az(numberOfVariables, numberOfVariables);
    Eigen::VectorXd b_x(numberOfVariables);
    Eigen::VectorXd b_y(numberOfVariables);
    Eigen::VectorXd b_z(numberOfVariables);
   
    int currentId=linearSystemIDs[2];
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        currentId = linearSystemIDs.at(variableId);
        printLine();
        //std::cout<<"\nsize: "<<linearSystemIDs.size()<<"\ncurrent id: "<<currentId<<"\n";

       
        vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
        printLine();
        connectedVertices = GetConnectedVertices(curvatureMesh, currentId);
        printLine();
        
        // The right hand side of the equation starts equal to 0
        double bvalue_x = 0.0;
        double bvalue_y = 0.0;
        double bvalue_z = 0.0;
        
        int cellNum = connectedVertices->GetNumberOfIds();
        printLine();

        double currentNormal[3];
        curvatureMesh->GetPointData()->GetScalars("interpolatedNoramlsArray")->GetTuple(currentId,currentNormal);

        double totalLineNormals[3];
        // Loop over the kernel around the current pixel
        for(unsigned int offset = 0; offset < cellNum ; ++offset)
        {
            
            
            printLine();
            double neighborLabel;
            double neighborNormal[3];
            //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
            curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
            curvatureMesh->GetPointData()->GetScalars("interpolatedNoramlsArray")->GetTuple
                (connectedVertices->GetId(offset),neighborNormal);

            double lineNormal[3];
            totalLineNormals[3] = {0};
            lineNormal[0] = neighborNormal[0] - currentNormal[0];
            lineNormal[1] = neighborNormal[1] - currentNormal[1];
            lineNormal[2] = neighborNormal[2] - currentNormal[2];
            vtkMath::Normalize(lineNormal);
           
            printLine();
        
            if (neighborLabel == 1)
            {
                std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                Ax.insert(variableId, neighborIdP->second) = lineNormal[0];
                Ay.insert(variableId, neighborIdP->second) = lineNormal[1];
                Az.insert(variableId, neighborIdP->second) = lineNormal[2];
            }
            else
            {
                double coordinate[3]; 
                curvatureMesh->GetPoint(connectedVertices->GetId(offset), coordinate);             
                bvalue_x -= coordinate[0];
                bvalue_y -= coordinate[1];
                bvalue_z -= coordinate[2];
                
            }

            totalLineNormals[0] += lineNormal[0];
            totalLineNormals[1] += lineNormal[1];
            totalLineNormals[2] += lineNormal[2];
            
            printLine();
        }
        Ax.insert(variableId, variableId) = -1 * totalLineNormals[0];
        Ay.insert(variableId, variableId) = -1 * totalLineNormals[1];
        Az.insert(variableId, variableId) = -1 * totalLineNormals[2];
        
        b_x[variableId] = bvalue_x;
        b_y[variableId] = bvalue_y;
        b_z[variableId] = bvalue_z;
      
    }
    std::cout<<"after matrix filling\n";
      // Solve the system with Eigen
    Eigen::VectorXd x(numberOfVariables);
    Eigen::VectorXd y(numberOfVariables);
    Eigen::VectorXd z(numberOfVariables);

    
   
    Eigen::SparseMatrix<double> A2x = Ax.adjoint() * Ax;
    Eigen::SparseMatrix<double> A2y = Ax.adjoint() * Ay;
    Eigen::SparseMatrix<double> A2z = Ax.adjoint() * Az;
    
    Eigen::VectorXd b2_x = Ax.adjoint() * b_x;
    Eigen::VectorXd b2_y = Ay.adjoint() * b_y;
    Eigen::VectorXd b2_z = Az.adjoint() * b_z;
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_Ax(A2x);
    if(!lu_of_Ax.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_Ax.solve(b2_x,&x))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_Ay(A2y);
    if(!lu_of_Ay.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_Ay.solve(b2_y,&y))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_Az(A2z);
    if(!lu_of_Az.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_Az.solve(b2_z,&z))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }//     return bonesVolume;



    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        double coordinate[3];
        coordinate[0] = x[variableId];
        coordinate[1] = y[variableId];
        coordinate[2] = z[variableId];
        int currentId = linearSystemIDs.at(variableId);
        curvatureMesh->GetPoints()->SetPoint(currentId, coordinate);
    }
    return curvatureMesh;
}

void computeTriangleNewCoordinates(double oldVerticesCoordinates[3][3], double center[3],double oldNormal[3],double newNormal[3]
                                   ,double newVerticesCoordinates[3][3])
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(oldVerticesCoordinates[0][0], oldVerticesCoordinates[0][1], oldVerticesCoordinates[0][2]);
    points->InsertNextPoint(oldVerticesCoordinates[1][0], oldVerticesCoordinates[1][1], oldVerticesCoordinates[1][2]);
    points->InsertNextPoint(oldVerticesCoordinates[2][0], oldVerticesCoordinates[2][1], oldVerticesCoordinates[2][2]);
 
    
    vtkSmartPointer<vtkPolyData> triangle = vtkPolyData::New();
    triangle->SetPoints(points);

    double axis[3];
    vtkMath::Cross(newNormal, oldNormal,axis);
    vtkMath::Normalize(axis);
    double angle = acos(vtkMath::Dot(newNormal, oldNormal)) *180 / 3.1415;
    
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->PostMultiply(); 
    transform->Translate(-center[0], -center[1], -center[2]);
    transform->RotateWXYZ(angle,axis);
    transform->Translate(center);
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInput(triangle);
    transformFilter->SetTransform(transform);
    transformFilter->Update();
    vtkSmartPointer<vtkPolyData> transformedTriangle = transformFilter->GetOutput();
    
     for( int i = 0; i < 3; ++i)
    {
        double coordinate[3];
        transformedTriangle->GetPoint(i , coordinate);
        newVerticesCoordinates[i][0] = coordinate[0];
        newVerticesCoordinates[i][1] = coordinate[1];
        newVerticesCoordinates[i][2] = coordinate[2];
    }
        
    
}

void computeCellNormals(double vertex1Normals[3],double vertex2Normals[3],double vertex3Normals[3], double cellNormal[3])
{
    cellNormal[0] = 1/3 * vertex1Normals[0] + vertex2Normals[0] + vertex3Normals[0];
    cellNormal[0] = 1/3 * vertex1Normals[1] + vertex2Normals[1] + vertex3Normals[1];
    cellNormal[0] = 1/3 * vertex1Normals[2] + vertex2Normals[2] + vertex3Normals[2];
    
    vtkMath::Normalize(cellNormal);
}

    
// vtkSmartPointer<vtkPolyData> poissonInterpolation(vtkSmartPointer<vtkPolyData> curvatureMesh, vtkSmartPointer<vtkDoubleArray> normalsArray)
// {
//     int num_of_nodes = curvatureMesh->GetNumberOfPoints();
//     // Create the sparse matrix
//     std::vector<int> linearSystemIDs;

//     writePolyData(curvatureMesh, "minCutMesh.vtk");

//     for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
//     {
//         double label;
//         curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);

//         if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
//         {
//             linearSystemIDs.push_back(meshCounter);
//         }
//     }
//     unsigned int numberOfVariables = linearSystemIDs.size();
//     std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
//     if (numberOfVariables == 0)
//     {
//         std::cerr << "No masked pixels found" << std::endl;
//         exit(1);
//     }

//     // Create the reverse mapping from pixel to variable id
//     std::map<int,int> meshIDs;//(numberOfVariables);
//     for(unsigned int i = 0; i < numberOfVariables; i++)
//     {
//         meshIDs[linearSystemIDs[i]] = i;
//     }
   
//     //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    
//     Eigen::SparseMatrix<double> A(numberOfVariables, numberOfVariables);
//     Eigen::VectorXd b_x(numberOfVariables);
//     Eigen::VectorXd b_y(numberOfVariables);
//     Eigen::VectorXd b_z(numberOfVariables);
   
//     int currentId=linearSystemIDs[2];
//     for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
//     {
//         currentId = linearSystemIDs.at(variableId);
//         printLine();
//         //std::cout<<"\nsize: "<<linearSystemIDs.size()<<"\ncurrent id: "<<currentId<<"\n";

       
//         vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
//         printLine();
//         connectedVertices = GetConnectedVertices(curvatureMesh, currentId);
//         printLine();

//         int cellNum = connectedVertices->GetNumberOfIds();
        
//         // The right hand side of the equation starts equal to 0
//         double bvalue[3];
        
//         //normalsMesh->GetPointData()->GetNormals()->GetTuple(currentId,
//         //bvalue);
//         double bvalue_x = 0, bvalue_y = 0, bvalue_z = 0;
//         // for(unsigned int offset = 0; offset < cellNum ; ++offset)
//         // {
//         //     normalsArray->GetTuple(variableId,bvalue);
            
//         //     bvalue_x = bvalue[0];
//         //     bvalue_y = bvalue[1];
//         //     bvalue_z = bvalue[2];
//         // }
       
//         printLine();
//         // Loop over the kernel around the current pixel
//         for(unsigned int offset = 0; offset < cellNum ; ++offset)
//         {
//             /// rotate triangles
//             double currentNormal[3], neighborNormal1[3], neighborNormal2[3];
//             double currentCoordinate[3], neighborCoordinate1[3], neighborCoordinate2[3];
            
//             normalsArray->GetTuple(variableId,currentNormal);
//             std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));

//             normalsArray->GetTuple(neighborIdP->second,neighborNormal1);
//             curvatureMesh->GetPoint(connectedVertices->GetId(offset), neighborCoordinate1);
//             curvatureMesh->GetPoint(currentId, currentCoordinate);
            
//             if(offset<cellNum - 1)
//             {
//                 std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset +1));
//                 normalsArray->GetTuple(neighborIdP->second,neighborNormal2);
//                 curvatureMesh->GetPoint(connectedVertices->GetId(offset+1), neighborCoordinate2);
                
//             }
//             else
//             {
//                 std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(0));
//                 normalsArray->GetTuple(neighborIdP->second,neighborNormal2);
//                 curvatureMesh->GetPoint(connectedVertices->GetId(0), neighborCoordinate2);
                
//             }

//             double triangleVertices[3][3]={
//                 {currentCoordinate[0],currentCoordinate[1],currentCoordinate[2]},
//                 {neighborCoordinate1[0],neighborCoordinate1[1],neighborCoordinate1[2]},
//                 {neighborCoordinate2[0],neighborCoordinate2[1],neighborCoordinate2[2]}
//             };

//             double triangleNormals[3][3]={
//                 {currentNormal[0],currentNormal[1],currentNormal[2]},
//                 {neighborNormal1[0],neighborNormal1[1],neighborNormal1[2]},
//                 {neighborNormal2[0],neighborNormal2[1],neighborNormal2[2]}
//             };

//             double triangleCenter[3], cellNormal[3], oldCellNormal[3];
//             double newVerticesCoordinates[3][3];
//             computeTriangleCenter(triangleVertices,triangleCenter);
//             computeCellNormals(triangleNormals,cellNormal);
//             computeTriangleNormal(triangleVertices,oldCellNormal);

//             computeTriangleNewCoordinates(triangleVertices, triangleCenter,oldCellNormal,cellNormal
//                                           ,newVerticesCoordinates);

//             bvalue_x += (newVerticesCoordinates[1][0] + newVerticesCoordinates[2][0] - newVerticesCoordinates[0][0]);
//             bvalue_y += (newVerticesCoordinates[1][1] + newVerticesCoordinates[2][1] - newVerticesCoordinates[0][1]);
//             bvalue_z += (newVerticesCoordinates[1][2] + newVerticesCoordinates[2][2] - newVerticesCoordinates[0][2]);
           
//             ///
            
//             printLine();
//             double neighborLabel;
//             //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
//             curvatureMesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
            
           
//             printLine();
        
//             if (neighborLabel == 1)
//             {
//                 std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
//                 A.insert(variableId, neighborIdP->second) = 1;
//                 //curvatureMesh->GetPoint(connectedVertices->GetId(offset), coordinate);

               
               
//             }
//             else
//             {
//                 double coordinate[3]; 
//                 curvatureMesh->GetPoint(connectedVertices->GetId(offset), coordinate);             
//                 bvalue_x -= coordinate[0];
//                 bvalue_y -= coordinate[1];
//                 bvalue_z -= coordinate[2];
                
//             }
           
//             printLine();
//         }
//         A.insert(variableId, variableId) = -1 * cellNum;
//         b_x[variableId] = bvalue_x;
//         b_y[variableId] = bvalue_y;
//         b_z[variableId] = bvalue_z;
//         std::cout<<bvalue_x<<"\t"<<bvalue_y<<"\t"<<bvalue_z<<"\n";
//     }
//     std::cout<<"after matrix filling\n";
//       // Solve the system with Eigen
//     Eigen::VectorXd x(numberOfVariables);
//     Eigen::VectorXd y(numberOfVariables);
//     Eigen::VectorXd z(numberOfVariables);

    
   
//     Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
//     Eigen::VectorXd b2_x = A.adjoint() * b_x;
//     Eigen::VectorXd b2_y = A.adjoint() * b_y;
//     Eigen::VectorXd b2_z = A.adjoint() * b_z;
    
//     Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_A(A2);
//     if(!lu_of_A.succeeded())
//     {
//         std::cerr << "decomposiiton failed!" << std::endl;
//         //exit(1);
//     }
//     if(!lu_of_A.solve(b2_x,&x))
//     {
//         std::cerr << "solving failed!" << std::endl;
//         //exit(1);
//     }
//     if(!lu_of_A.solve(b2_y,&y))
//     {
//         std::cerr << "solving failed!" << std::endl;
//         //exit(1);
//     }
//     if(!lu_of_A.solve(b2_z,&z))
//     {
//         std::cerr << "solving failed!" << std::endl;
//         //exit(1);
//     }//     return bonesVolume;



//     for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
//     {
//         double coordinate[3];
//         coordinate[0] = x[variableId];
//         coordinate[1] = y[variableId];
//         coordinate[2] = z[variableId];
//         int currentId = linearSystemIDs.at(variableId);
//         curvatureMesh->GetPoints()->SetPoint(currentId, coordinate);
//     }
//     std::cout<<"Done Poisson interpolation\n";
//     return curvatureMesh;
// }

void computeLaplaceOperator(vtkSmartPointer<vtkPolyData> mesh)
{
    
    int num_of_nodes = mesh->GetNumberOfPoints();
    
    vtkSmartPointer<vtkDoubleArray> laplaceArray = 
        vtkSmartPointer<vtkDoubleArray>::New();
    laplaceArray->SetNumberOfComponents(3);
    laplaceArray->SetNumberOfTuples(num_of_nodes);
    laplaceArray->SetName("laplaceArray");

    for(unsigned int meshCounter = 0; meshCounter < num_of_nodes; ++meshCounter)
    {
        double laplaceOperator[3] = {0};
        double currentCoordinate[3]; 

        mesh->GetPoint(meshCounter, currentCoordinate);           
        vtkSmartPointer<vtkIdList> neighborsList = GetConnectedVertices(mesh, meshCounter );

        for(int idsCounter = 0;idsCounter < neighborsList->GetNumberOfIds(); ++ idsCounter)
        {
            double neighborCoordinate[3]; 
            mesh->GetPoint(neighborsList->GetId(idsCounter), neighborCoordinate);           
            
            laplaceOperator[0] = laplaceOperator[0] - currentCoordinate[0] + neighborCoordinate[0];
            laplaceOperator[1] = laplaceOperator[1] - currentCoordinate[1] + neighborCoordinate[1];
            laplaceOperator[2] = laplaceOperator[2] - currentCoordinate[2] + neighborCoordinate[2];
        }
        laplaceArray->InsertTuple3(meshCounter,laplaceOperator[0],laplaceOperator[1],laplaceOperator[2]);
                
    }
    mesh->GetPointData()->AddArray(laplaceArray);      
}

//interpolates the laplacian values over unknown domain, given its values over
//the known outer surface
vtkSmartPointer<vtkDoubleArray> laplaceLaplaceInterpolation(vtkSmartPointer<vtkPolyData> mesh)
{
    int num_of_nodes = mesh->GetNumberOfPoints();
    // Create the sparse matrix
    std::vector<int> linearSystemIDs;
    
    for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
    {
        double label;
        mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);
        printLine();
        if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
        {
            linearSystemIDs.push_back(meshCounter);
        }
    }
    unsigned int numberOfVariables = linearSystemIDs.size();
    std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
    if (numberOfVariables == 0)
    {
        std::cerr << "No masked pixels found" << std::endl;
        exit(1);
    }

    // Create the reverse mapping from pixel to variable id
    std::map<int,int> meshIDs;//(numberOfVariables);
    for(unsigned int i = 0; i < numberOfVariables; i++)
    {
        meshIDs[linearSystemIDs[i]] = i;
    }
   
    //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    
    Eigen::SparseMatrix<double> A(numberOfVariables, numberOfVariables);
    Eigen::VectorXd b_x(numberOfVariables);
    Eigen::VectorXd b_y(numberOfVariables);
    Eigen::VectorXd b_z(numberOfVariables);


    
    int currentId=linearSystemIDs[2];
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        currentId = linearSystemIDs.at(variableId);
        //printLine();       
        vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
        printLine();
        connectedVertices = GetConnectedVertices(mesh, currentId);
        printLine();
        
        // The right hand side of the equation starts equal to 0
        double bvalue_x = 0.0;
        double bvalue_y = 0.0;
        double bvalue_z = 0.0;
        
        int cellNum = connectedVertices->GetNumberOfIds();
        printLine();
        // Loop over the kernel around the current pixel
        for(unsigned int offset = 0; offset < cellNum ; ++offset)
        {
            
            double neighborLabel;
            mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
                   
            if (neighborLabel == 1)
            {
                std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                A.insert(variableId, neighborIdP->second) = 1;
            }
            else
            {

                double laplaceOperator[3];
                //normals->GetTuple(connectedVertices->GetId(offset),normal);
                mesh->GetPointData()->GetScalars("laplaceArray")->GetTuple(connectedVertices->GetId(offset), laplaceOperator);
                bvalue_x -= laplaceOperator[0];
                bvalue_y -= laplaceOperator[1];
                bvalue_z -= laplaceOperator[2];

            }
           
            printLine();
        }
        A.insert(variableId, variableId) = -1 * cellNum;
        b_x[variableId] = bvalue_x;
        b_y[variableId] = bvalue_y;
        b_z[variableId] = bvalue_z;
      
    }
    std::cout<<"after matrix filling\n";
      // Solve the system with Eigen
    Eigen::VectorXd x(numberOfVariables);
    Eigen::VectorXd y(numberOfVariables);
    Eigen::VectorXd z(numberOfVariables);

    
   
    Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
    Eigen::VectorXd b2_x = A.adjoint() * b_x;
    Eigen::VectorXd b2_y = A.adjoint() * b_y;
    Eigen::VectorXd b2_z = A.adjoint() * b_z;
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_A(A2);
    if(!lu_of_A.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_x,&x))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_y,&y))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_z,&z))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }//     return bonesVolume;

    // Set point normals
    vtkSmartPointer<vtkDoubleArray> interpolatedLaplaceArray = 
        vtkSmartPointer<vtkDoubleArray>::New();
    interpolatedLaplaceArray->SetNumberOfComponents(3); 
    interpolatedLaplaceArray->SetName("interpolatedLaplaceArray");
    interpolatedLaplaceArray->SetNumberOfTuples(num_of_nodes);

    for(int i = 0; i<num_of_nodes;++i)
    {
        double laplaceOperator[3];
        mesh->GetPointData()->GetScalars("laplaceArray")->GetTuple(i, laplaceOperator);
        interpolatedLaplaceArray->SetTuple(i, laplaceOperator) ;
       
    }
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        double interpolatedLaplaceOperator[3];
        interpolatedLaplaceOperator[0] = x[variableId];
        interpolatedLaplaceOperator[1] = y[variableId];
        interpolatedLaplaceOperator[2] = z[variableId];

        // Add the data to the  array
        int currentId = linearSystemIDs.at(variableId);
        interpolatedLaplaceArray->SetTuple(currentId, interpolatedLaplaceOperator) ;
        

    }

   
    mesh->GetPointData()->AddArray(interpolatedLaplaceArray);

    return interpolatedLaplaceArray;
}



vtkSmartPointer<vtkPolyData> poissonLaplaceInterpolation(vtkSmartPointer<vtkPolyData> mesh)
{
    int num_of_nodes = mesh->GetNumberOfPoints();
    // Create the sparse matrix
    std::vector<int> linearSystemIDs;

    writePolyData(mesh, "minCutMesh.vtk");

    for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
    {
        double label;
        mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);

        if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
        {
            linearSystemIDs.push_back(meshCounter);
        }
    }
    unsigned int numberOfVariables = linearSystemIDs.size();
    std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
    if (numberOfVariables == 0)
    {
        std::cerr << "No masked pixels found" << std::endl;
        exit(1);
    }

    // Create the reverse mapping from pixel to variable id
    std::map<int,int> meshIDs;//(numberOfVariables);
    for(unsigned int i = 0; i < numberOfVariables; i++)
    {
        meshIDs[linearSystemIDs[i]] = i;
    }
   
    //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    
    Eigen::SparseMatrix<double> A(numberOfVariables, numberOfVariables);
    Eigen::VectorXd b_x(numberOfVariables);
    Eigen::VectorXd b_y(numberOfVariables);
    Eigen::VectorXd b_z(numberOfVariables);
   
    int currentId=linearSystemIDs[2];
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        currentId = linearSystemIDs.at(variableId);
        printLine();
        //std::cout<<"\nsize: "<<linearSystemIDs.size()<<"\ncurrent id: "<<currentId<<"\n";

       
        vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
        printLine();
        connectedVertices = GetConnectedVertices(mesh, currentId);
        printLine();

        double laplaceOperator[3];
        mesh->GetPointData()->GetScalars("interpolatedLaplaceArray")->GetTuple(currentId, laplaceOperator);

        // The right hand side of the equation starts equal to 0
        double ratio = 20;
        double bvalue_x = laplaceOperator[0]/ratio;
        double bvalue_y = laplaceOperator[1]/ratio;
        double bvalue_z = laplaceOperator[2]/ratio;
        
        int cellNum = connectedVertices->GetNumberOfIds();
        printLine();
        // Loop over the kernel around the current pixel
        for(unsigned int offset = 0; offset < cellNum ; ++offset)
        {
            
            
            printLine();
            double neighborLabel;
            //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
            mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
           
           
            printLine();
        
            if (neighborLabel == 1)
            {
                std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                A.insert(variableId, neighborIdP->second) = 1;
            }
            else
            {
                double coordinate[3]; 
                mesh->GetPoint(connectedVertices->GetId(offset), coordinate);             
                bvalue_x -= coordinate[0];
                bvalue_y -= coordinate[1];
                bvalue_z -= coordinate[2];
                
            }
           
            printLine();
        }
        A.insert(variableId, variableId) = -1 * cellNum;
        b_x[variableId] = bvalue_x;
        b_y[variableId] = bvalue_y;
        b_z[variableId] = bvalue_z;
      
    }
    std::cout<<"after matrix filling\n";
      // Solve the system with Eigen
    Eigen::VectorXd x(numberOfVariables);
    Eigen::VectorXd y(numberOfVariables);
    Eigen::VectorXd z(numberOfVariables);

    
   
    Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
    Eigen::VectorXd b2_x = A.adjoint() * b_x;
    Eigen::VectorXd b2_y = A.adjoint() * b_y;
    Eigen::VectorXd b2_z = A.adjoint() * b_z;
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_A(A2);
    if(!lu_of_A.succeeded())
    {
        std::cerr << "decomposiiton failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_x,&x))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_y,&y))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2_z,&z))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }//     return bonesVolume;



    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        double coordinate[3];
        coordinate[0] = x[variableId];
        coordinate[1] = y[variableId];
        coordinate[2] = z[variableId];
        int currentId = linearSystemIDs.at(variableId);
        mesh->GetPoints()->SetPoint(currentId, coordinate);
    }
    return mesh;
}



//interpolate normals values, using known boundary
vtkSmartPointer<vtkPolyData> normalsLinearInterpolation(vtkSmartPointer<vtkPolyData> mesh)
{
    std::cout<<"in normals linear interpolation():\n"<<"line: "<<__LINE__<<"\n";
    int num_of_nodes = mesh->GetNumberOfPoints();
    // Create the sparse matrix
    std::vector<int> linearSystemIDs;
    
    int numberOfEdges = 0;
    int numOfMarginNeighbors;

    for (unsigned int meshCounter = 0; meshCounter < num_of_nodes ;meshCounter++)
    {
        double label;
        mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(meshCounter,&label);
        
        if(label == 1 ) // The mask is non-zero representing that we want to fill this pixel
        {
            linearSystemIDs.push_back(meshCounter);
        }
    }

    // Create the reverse mapping from pixel to variable id
    std::map<int,int> meshIDs;//(numberOfVariables);
    for(unsigned int i = 0; i < linearSystemIDs.size(); i++)
    {
        meshIDs[linearSystemIDs[i]] = i;
       
    }
    
    std::vector<char>visited1(linearSystemIDs.size(),0);
    for (unsigned int counter = 0;counter<linearSystemIDs.size();++counter)
    {
        double label;
        mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(linearSystemIDs[counter],&label);
        
        vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
        connectedVertices = GetConnectedTriangles(mesh, linearSystemIDs[counter]);
        int cellNum = connectedVertices->GetNumberOfIds();
        
        for(unsigned int offset = 0; offset < cellNum ; offset++)
        {
            double neighborLabel;
            //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
            mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedVertices->GetId(offset),&neighborLabel);
            if(neighborLabel != 1)
            {
                ++numberOfEdges;
               
            }
            else
            {
                int neighborId = connectedVertices->GetId(offset);//linearSystemIDs.at(connectedTriangles->GetId(offset));            
                std::map<int,int>::iterator  neighborIdP = meshIDs.find(neighborId);
                if(visited1[neighborIdP->second] != 1)
                {
                    ++numberOfEdges;
                }
            }
        }
        visited1[counter] = 1;
        
    }
    
    
    unsigned int numberOfVariables = 3 * linearSystemIDs.size(); //three times. One per coordinate
    
    std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
    std::cout<<"numer of edges = "<<numberOfEdges<<std::endl;
    if (numberOfVariables == 0)
    {
        std::cerr << "No masked pixels found" << std::endl;
        exit(1);
    }

    
   
    //    vtkSmartPointer<vtkIdList> connectedVertices = vtkIdList::New();
    printLine();
    Eigen::SparseMatrix<double> A(numberOfEdges, numberOfVariables);
    Eigen::VectorXd b(numberOfEdges);
    printLine();

    int currentRow = 0;
    int currentId, neighbor1Id, neighbor2Id;
    numOfMarginNeighbors = 0;
    int nonMarginalEdgeNumber = 0;
    std::vector<char> vistedId(linearSystemIDs.size(),0);
    for(unsigned int variableId = 0; variableId < linearSystemIDs.size(); variableId++)
    {
        printLine();
        currentId = linearSystemIDs.at(variableId);
        std::map<int,int>::iterator  currentIdP = meshIDs.find(currentId);

       
        vtkSmartPointer<vtkIdList> connectedTriangles = vtkIdList::New();
       
        connectedTriangles = GetConnectedTriangles(mesh, currentId);

        
        
        // The right hand side of the equation starts equal to 0
        double bvalue[3] = {0};
        
        
        int cellNum = connectedTriangles->GetNumberOfIds();
       
        printLine();
        // Loop over the kernel around the current pixel
        for(unsigned int offset = 0; offset < cellNum - 1 ; offset+=2)
        {
            
            double currentNormal[3], neighbor1Normal[3], neighbor2Normal[3];
            double currentCellNormal[3];
            
            mesh->GetPointData()->GetScalars("interpolatedNoramlsArray")->GetTuple(currentId,currentNormal);
            mesh->GetPointData()->GetScalars("interpolatedNoramlsArray")->GetTuple(connectedTriangles->GetId(offset),neighbor1Normal);
            mesh->GetPointData()->GetScalars("interpolatedNoramlsArray")->GetTuple(connectedTriangles->GetId(offset+1),neighbor2Normal);
            computeCellNormals(currentNormal,neighbor1Normal,neighbor2Normal, currentCellNormal);

            printLine();
            //      std::cout<<connectedTriangles->GetId(offset)<<"\t"<<connectedTriangles->GetId(offset+1)<<"\t"<<linearSystemIDs.size()<<"\n";
            neighbor1Id = connectedTriangles->GetId(offset);//linearSystemIDs.at(connectedTriangles->GetId(offset));
            neighbor2Id = connectedTriangles->GetId(offset+1);//linearSystemIDs.at(connectedTriangles->GetId(offset+1));
            printLine();
            double neighbor1Label, neighbor2Label;
            //std::cout<<"neighbor id: "<<connectedVertices->GetId(offset)<<"\n";
            mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedTriangles->GetId(offset),&neighbor1Label);
            mesh->GetPointData()->GetScalars("minCutArray")->GetTuple(connectedTriangles->GetId(offset),&neighbor2Label);
            
           
            
            //printLine();
            std::map<int,int>::iterator  neighbor1IdP = meshIDs.find(neighbor1Id);
            if (neighbor1Label == 1)
            {
                printLine();
               

                //  std::cout<<currentRow<<"\t"<<currentIdP->second<<"\t"<<neighbor1IdP->second<<"\n";
                 
                //std::cout<<"numer of varaibles = "<<numberOfVariables<<std::endl;
                //std::cout<<"numer of edges = "<<numberOfEdges<<std::endl;
    
                //std::map<int,int>::iterator  neighborId2P =
                //meshIDs.find(neighbor2Id);
                if(vistedId[neighbor1IdP->second] != 1)
                {
                    A.insert(currentRow, 3*(currentIdP->second)) = currentCellNormal[0];
                    A.insert(currentRow, 3*(currentIdP->second)+1) = currentCellNormal[1];
                    A.insert(currentRow, 3*(currentIdP->second)+2) = currentCellNormal[2];
                    A.insert(currentRow, 3*(neighbor1IdP->second)) = -1 * currentCellNormal[0];
                    A.insert(currentRow, 3*(neighbor1IdP->second) + 1) = -1 * currentCellNormal[1];
                    A.insert(currentRow, 3*(neighbor1IdP->second) + 2) = -1 * currentCellNormal[2];
                    b[currentRow] = 0;
                    ++currentRow;
                    ++nonMarginalEdgeNumber;
                    
                }
                    printLine();
            }
            else
            {
                printLine();
                double coordinate[3];
                mesh->GetPoints()->GetPoint(neighbor1Id, coordinate);
                //std::cout<<currentRow<<"\t"<<currentId<<"\n";
                
                //std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                A.insert(currentRow, 3*(currentIdP->second))    = currentCellNormal[0];
                A.insert(currentRow, 3*( currentIdP->second)+1) = currentCellNormal[1];
                A.insert(currentRow, 3*( currentIdP->second)+2) = currentCellNormal[2];
                b[currentRow] = coordinate[0]*currentCellNormal[0]+coordinate[1]*currentCellNormal[1]+coordinate[2]*currentCellNormal[2];
                ++currentRow;
                printLine();
                ++numOfMarginNeighbors;
            }
            //printLine();
            if (neighbor2Label == 1)
            {
                printLine();
                //std::map<int,int>::iterator  neighborIdP =
                //meshIDs.find(connectedVertices->GetId(offset));
                //std::map<int,int>::iterator  neighborId1P = meshIDs.find(neighbor1Id);
                std::map<int,int>::iterator  neighborId2P = meshIDs.find(neighbor2Id);
                //std::cout<<currentRow<<"\t"<<currentIdP->second<<"\t"<<neighborId2P->second<<"\n"; 
                //std::cout<<"2 numer of varaibles = "<<numberOfVariables<<std::endl;
                //std::cout<<"2 numer of edges = "<<numberOfEdges<<std::endl;
                if(vistedId[neighborId2P->second] != 1)
                {
                    A.insert(currentRow, 3*(currentIdP->second)) = currentCellNormal[0];
                    A.insert(currentRow, 3*(currentIdP->second)+1) = currentCellNormal[1];
                    A.insert(currentRow, 3*(currentIdP->second)+2) = currentCellNormal[2];
                    A.insert(currentRow, 3*(neighborId2P->second)) = -1 * currentCellNormal[0];
                    A.insert(currentRow, 3*(neighborId2P->second) + 1) = -1 * currentCellNormal[1];
                    A.insert(currentRow, 3*(neighborId2P->second) + 2) = -1 * currentCellNormal[2];
                    b[currentRow] = 0;
                    ++currentRow;
                    ++nonMarginalEdgeNumber;
                }
                printLine();
            }
            else
            {
                printLine();
                double coordinate[3];
                mesh->GetPoints()->GetPoint(neighbor2Id, coordinate);
                
                //std::map<int,int>::iterator  neighborIdP = meshIDs.find(connectedVertices->GetId(offset));
                A.insert(currentRow, 3*(currentIdP->second)) = currentCellNormal[0];
                A.insert(currentRow, 3*(currentIdP->second)+1) = currentCellNormal[1];
                A.insert(currentRow, 3*(currentIdP->second)+2) = currentCellNormal[2];
                b[currentRow] = coordinate[0]*currentCellNormal[0]+coordinate[1]*currentCellNormal[1]+coordinate[2]*currentCellNormal[2];
                ++currentRow;
                ++numOfMarginNeighbors;
            }
                //printLine();
            }
            vistedId[currentIdP->second] = 1;
            
            
            }
    std::cout<<"number of rows: "<<currentRow-1<<"\n"<<"number of marginals = "<<numOfMarginNeighbors<<"\n"<<"number of non marginals = "<<nonMarginalEdgeNumber<<"\n";
    printLine();
    std::cout<<"after matrix filling\n";
      // Solve the system with Eigen
    Eigen::VectorXd x(numberOfVariables);
   
    
   
    Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
    Eigen::VectorXd b2 = A.adjoint() * b;
   
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::UmfPack> lu_of_A(A2);
    if(!lu_of_A.succeeded())
    {
        std::cerr << "decompositon failed!" << std::endl;
        //exit(1);
    }
    if(!lu_of_A.solve(b2,&x))
    {
        std::cerr << "solving failed!" << std::endl;
        //exit(1);
    }
   
    for(unsigned int variableId = 0; variableId < numberOfVariables; variableId+=3)
    {
        double coordinate[3];
        coordinate[0] = x[variableId];
        coordinate[1] = x[variableId+1];
        coordinate[2] = x[variableId+2];
        int currentId = linearSystemIDs.at(variableId/3);
        mesh->GetPoints()->SetPoint(currentId, coordinate);
    }
    
    return mesh;
  
}

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
    
    
    for( int i = 0; i < mesh->GetNumberOfPoints(); ++i)
    {
        double coordinate[3];
        mesh->GetPoint(i , coordinate);
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

    for(inIt.GoToBegin(),outIt.GoToBegin(); !outIt.IsAtEnd();++outIt,++inIt)
    {
        typename ImageType::PointType imagePoint;
        outputImage->TransformIndexToPhysicalPoint(outIt.GetIndex(),imagePoint);
        double testPoint[3] = {imagePoint[0], imagePoint[1], imagePoint[2]};
        
        if(inIt.Value() != 0)
        {
            //Find the closest points to TestPoint
            double closestPoint[3];//the coordinates of the closest point will be returned here
            double closestPointDist2; //the squared distance to the closest point will be returned here
            vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
            int subId; //this is rarely used (in triangle strips only, I believe)
            cellLocator->FindClosestPoint(testPoint, closestPoint, cellId, subId, closestPointDist2);
            if(closestPointDist2 <= 1)
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

template<class ImageType>
typename ImageType::Pointer correctImage(typename ImageType::Pointer leakyImage, typename ImageType::Pointer seedImage, typename ImageType::Pointer contourImage)
{
    typedef typename itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType leakIt(leakyImage,leakyImage->GetLargestPossibleRegion());
    IteratorType contourIt(contourImage,contourImage->GetLargestPossibleRegion());
    
    for (leakIt.GoToBegin(), contourIt.GoToBegin(); ! leakIt.IsAtEnd(); ++leakIt, ++contourIt)
    {
        if(contourIt.Value() == 1)
        {
            leakIt.Set(0);
        }

    }
    
    //retain only thresholded pixels that are connected to the correction seed
    typedef typename itk::ConnectedComponentImageFilter<ImageType,ImageType> ConnectorType;
    
    typename ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(leakyImage);
    connector->SetFullyConnected(false);
    connector->Update();

    IteratorType seedIt(seedImage,seedImage->GetLargestPossibleRegion());
    IteratorType connectIt(connector->GetOutput(),connector->GetOutput()->GetLargestPossibleRegion());

    int label;
    for(seedIt.GoToBegin(),connectIt.GoToBegin();!seedIt.IsAtEnd();++seedIt,++connectIt)
    {
        if(seedIt.Value() == CORRECT_SEED && connectIt.Value() != 0)
        {
            label = connectIt.Value();
        }
    }


    typedef typename itk::BinaryThresholdImageFilter<ImageType,ImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput(connector->GetOutput());
    thresholder->SetLowerThreshold(label);
    thresholder->SetUpperThreshold(label);
    thresholder->SetInsideValue(1);
    thresholder->SetOutsideValue(0);
    thresholder->Update();
    typename ImageType::Pointer correctedImage = thresholder->GetOutput();

    IteratorType correctedIt(correctedImage, correctedImage->GetLargestPossibleRegion());
    for(contourIt.GoToBegin(),correctedIt.GoToBegin();!contourIt.IsAtEnd(); ++contourIt, ++correctedIt)
    {
        if(contourIt.Value() == 1)
        {
            correctedIt.Set(1);
        }
    }
    
    return correctedImage;


}
int main( int argc, char* argv[] )
{
    typedef short PixelType;
    typedef unsigned short SegPixelType;
    const unsigned char dim = 3;

    typedef itk::Image<PixelType,dim> ImageType;
    typedef itk::Image<SegPixelType, dim> SegImageType;

    ImageType::Pointer inputImage = read3DImage<ImageType>(argv[INPUT_IMAGE_NAME]);
    SegImageType::Pointer segInputImage = read3DImage<SegImageType>(argv[SEG_INPUT_IMAGE_NAME]);
    SegImageType::Pointer seedInputImage = read3DImage<SegImageType>(argv[SEED_INPUT_IMAGE_NAME]);
    SegImageType::Pointer interiorSeedInputImage = read3DImage<SegImageType>(argv[INTERIOR_SEED_NAME]);

    typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType> GradientFilterType;
    GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
    gradientFilter->SetInput(inputImage);
    gradientFilter->Update();
    ImageType::Pointer gradientInputImage = gradientFilter->GetOutput();
    
    typedef itk::ImageRegionIterator<ImageType> IteratorType;
    typedef itk::ImageRegionIterator<SegImageType> SegIteratorType;

    IteratorType gradIt(gradientFilter->GetOutput(),gradientFilter->GetOutput()->GetLargestPossibleRegion());
    SegIteratorType seedIt(seedInputImage,seedInputImage->GetLargestPossibleRegion());


    typedef itk::ImageToVTKImageFilter<SegImageType> SegConverterType;
    SegConverterType::Pointer converter = SegConverterType::New();
    converter->SetInput(segInputImage);
    converter->Update();
    vtkSmartPointer<vtkImageData> vtkSegImage = converter->GetOutput();

    SegConverterType::Pointer seedConverter = SegConverterType::New();
    seedConverter->SetInput(seedInputImage);
    seedConverter->Update();
    vtkSmartPointer<vtkImageData> vtkSeedImage = seedConverter->GetOutput();
    
    vtkSmartPointer<vtkPolyData> mesh = createAndSmoothSurface(vtkSegImage, MESH_SMOOTH_ITERATIONS);

    typedef itk::ImageToVTKImageFilter<ImageType> ConverterType;
    ConverterType::Pointer gradientConverter = ConverterType::New();
    gradientConverter->SetInput(inputImage);
    gradientConverter->Update();
    vtkSmartPointer<vtkImageData> vtkGradientImage = gradientConverter->GetOutput();
    
    vtkSmartPointer<vtkProbeFilter> probeFilter = vtkProbeFilter::New();
    probeFilter->SetSource(vtkGradientImage);
    probeFilter->SetInput(mesh);
    probeFilter->Update();
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkGeometryFilter::New();
    geometryFilter->SetInput(probeFilter->GetOutput());
    geometryFilter->Update();
    vtkSmartPointer<vtkPolyData> gradientMesh = geometryFilter->GetOutput();
    writePolyData(gradientMesh, "gradientMesh.vtk");
    
    
    vtkSmartPointer<vtkPolyData> seedMesh = sampleImageOnMesh<SegImageType>(mesh,seedInputImage);
    
    writePolyData(seedMesh, "seedMesh.vtk");
    
    vtkSmartPointer<vtkPolyData> minCurvatureMesh = polyDataToMinCurvature(mesh);

    //just temporary - don't forget to delete
    vtkSmartPointer<vtkPolyData> maxCurvatureMesh = polyDataToMaxCurvature(mesh);
    writePolyData(maxCurvatureMesh, "maxCurvature.vtk");
    //  --- up to here

    vtkSmartPointer<vtkPolyData> minCutMeshLeaks = minCut(minCurvatureMesh,seedMesh,gradientMesh,MIN_CURVATURE_TAG,atof(argv[GRAPH_SIGMA]));

    vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks = minCut(maxCurvatureMesh,seedMesh,gradientMesh,MAX_CURVATURE_TAG,atof(argv[GRAPH_SIGMA]));

    vtkSmartPointer<vtkPolyData> minCutMesh = minCutConjunction(minCutMeshLeaks,minCutMeshInteriorLeaks);
    attributeDilation( minCutMesh , ATTRIBUTE_DILATION_RADIUS);
  
    /*computeLaplaceOperator(minCutMesh);
    vtkSmartPointer<vtkDoubleArray> normalsArray = laplaceNormalsInterpolation(minCutMesh);
    writePolyData(minCutMesh,"minCutMesh.vtk");

    vtkSmartPointer<vtkPolyData> secondOrderMesh = normalsLinearInterpolation(minCutMesh);
    writePolyData(secondOrderMesh,"secondOrderMesh.vtk");*/

        
    vtkSmartPointer<vtkPolyData> correctedMesh1 = laplaceInterpolation(minCutMesh);
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
    normals->SetInput(correctedMesh1);
    normals->FlipNormalsOn();
    normals->Update();
    vtkSmartPointer<vtkPolyData> correctedMesh = normals->GetOutput();
    writePolyData(correctedMesh,"laplaceMesh.vtk");
   
    
    vtkSmartPointer<vtkDecimatePro> decimateFilter = vtkDecimatePro::New();
    decimateFilter->SetInput(correctedMesh);
    decimateFilter->SetTargetReduction(DECIMATION_FACTOR);
    decimateFilter->PreserveTopologyOn();
    decimateFilter->Update();
    vtkSmartPointer<vtkPolyData> decimatedMesh = decimateFilter->GetOutput();
    
    
    
    writePolyData(decimatedMesh, argv[OUTPUT_CURVATURE_MESH]);
    

    SegImageType::Pointer outputContourImage = sampleMeshOnImage<SegImageType>(/*decimatedMesh*/correctedMesh1,segInputImage);

    SegImageType::Pointer outputImage = correctImage<SegImageType>(segInputImage, /*seedInputImage*/interiorSeedInputImage,outputContourImage);

    typedef itk::ImageFileWriter<SegImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(outputImage);
    writer->SetFileName("correctedImage.nii.gz");
    writer->Update();


    SegConverterType::Pointer correctionConverter = SegConverterType::New();
    correctionConverter->SetInput(outputImage);
    correctionConverter->Update();
    vtkSmartPointer<vtkImageData> vtkCorrectedImage = correctionConverter->GetOutput();
    
    vtkSmartPointer<vtkPolyData> outputMesh = createAndSmoothSurface(vtkCorrectedImage, 50);
    writePolyData(outputMesh, "final.vtk");
    
    return EXIT_SUCCESS;
    
}


