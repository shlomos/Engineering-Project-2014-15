#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
 
int main(int argc, char* argv[])
{
  // Verify input arguments
  if ( argc != 2 )
    {
    std::cout << "Usage: " << argv[0]
              << " Filename(.vtk)" << std::endl;
    return EXIT_FAILURE;
    }
 
  std::string inputFilename = argv[1];
 
  // Read the file
  vtkSmartPointer<vtkStructuredPointsReader> reader =
    vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();

 
    // Visualize
   vtkSmartPointer<vtkImageViewer2>viewer=
	vtkSmartPointer<vtkImageViewer2>::New();
	vtkSmartPointer<vtkRenderWindowInteractor>interactor = 
	vtkSmartPointer<vtkRenderWindowInteractor>::New();

interactor->SetRenderWindow( viewer->GetRenderWindow() ); 
viewer->SetupInteractor(interactor);
viewer->SetInputConnection(reader->GetOutputPort());
viewer->SetSlice(130);
viewer->GetRenderer()->ResetCamera();
viewer->Render();
interactor->Start();
 
  return EXIT_SUCCESS;
}