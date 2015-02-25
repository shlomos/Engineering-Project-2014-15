#include "MarchingCubes.h"

// Handle mouse events
class MouseInteractorStyle2 : public vtkInteractorStyleJoystickCamera
{
public:
	static MouseInteractorStyle2* New();
	vtkTypeMacro(MouseInteractorStyle2, vtkInteractorStyleJoystickCamera);
	vtkIdType _currSource = -1;
	virtual void OnRightButtonDown()
	{
		MakeAnnotation(BACKGROUND);
		vtkInteractorStyleJoystickCamera::OnRightButtonDown();
	}

	void RemoveLeaks(){
		cout << "Started fixing mesh!" << endl;
		typedef short PixelType;
		typedef unsigned short SegPixelType;
		const unsigned char dim = 3;

		typedef itk::Image<PixelType, dim> ImageType;
		typedef itk::Image<SegPixelType, dim> SegImageType;

		MeshLeaksCorrector mlc;
		ImageType::Pointer inputImage = mlc.read3DImage<ImageType>("GG.vtk");
		//SegImageType::Pointer segInputImage = read3DImage<SegImageType>(argv[SEG_INPUT_IMAGE_NAME]);
		//SegImageType::Pointer seedInputImage = read3DImage<SegImageType>(argv[SEED_INPUT_IMAGE_NAME]);
		//SegImageType::Pointer interiorSeedInputImage = read3DImage<SegImageType>(argv[INTERIOR_SEED_NAME]);

		typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;
		GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
		gradientFilter->SetInput(inputImage);
		gradientFilter->Update();
		ImageType::Pointer gradientInputImage = gradientFilter->GetOutput();

		typedef itk::ImageRegionIterator<ImageType> IteratorType;
		typedef itk::ImageRegionIterator<SegImageType> SegIteratorType;

		IteratorType gradIt(gradientFilter->GetOutput(), gradientFilter->GetOutput()->GetLargestPossibleRegion());
		//SegIteratorType seedIt(seedInputImage, seedInputImage->GetLargestPossibleRegion());


		//typedef itk::ImageToVTKImageFilter<SegImageType> SegConverterType;
		//SegConverterType::Pointer converter = SegConverterType::New();
		//converter->SetInput(segInputImage);
		//converter->Update();
		//vtkSmartPointer<vtkImageData> vtkSegImage = converter->GetOutput();

		//SegConverterType::Pointer seedConverter = SegConverterType::New();
		//seedConverter->SetInput(seedInputImage);
		//seedConverter->Update();
		//vtkSmartPointer<vtkImageData> vtkSeedImage = seedConverter->GetOutput();
		vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
		vtkSmartPointer<vtkPolyData> vtkSegImage = mapper->GetInput();

		vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmoothPolyDataFilter::New();
		smoother->SetInputData(vtkSegImage);
		smoother->SetNumberOfIterations(MESH_SMOOTH_ITERATIONS);
		smoother->Update();
		vtkSmartPointer<vtkPolyData> mesh = smoother->GetOutput();

		vtkSmartPointer<vtkStructuredPointsReader> reader =
			vtkSmartPointer<vtkStructuredPointsReader>::New();
		reader->SetFileName("GG_30_12_10_S30_liver.vtk");
		reader->Update();
		vtkSmartPointer<vtkImageData> vtkGradientImage = reader->GetOutput();
		/*typedef itk::ImageToVTKImageFilter<ImageType> ConverterType;
		ConverterType::Pointer gradientConverter = ConverterType::New();
		gradientConverter->SetInput(inputImage);
		gradientConverter->Update();
		vtkSmartPointer<vtkImageData> vtkGradientImage = gradientConverter->GetOutput();*/

		vtkSmartPointer<vtkProbeFilter> probeFilter = vtkProbeFilter::New();
		probeFilter->SetSourceData(vtkGradientImage);
		probeFilter->SetInputData(mesh);
		probeFilter->Update();
		vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkGeometryFilter::New();
		geometryFilter->SetInputData(probeFilter->GetOutput());
		geometryFilter->Update();
		vtkSmartPointer<vtkPolyData> gradientMesh = geometryFilter->GetOutput();

		vtkSmartPointer<vtkPolyData> minCurvatureMesh = mlc.polyDataToMinCurvature(mesh);

		//just temporary - don't forget to delete
		vtkSmartPointer<vtkPolyData> maxCurvatureMesh = mlc.polyDataToMaxCurvature(mesh);
		//  --- up to here

		vtkSmartPointer<vtkPolyData> minCutMeshLeaks = mlc.minCut(minCurvatureMesh, mesh, gradientMesh, MIN_CURVATURE_TAG, 1.0f);

		vtkSmartPointer<vtkPolyData> minCutMeshInteriorLeaks = mlc.minCut(maxCurvatureMesh, mesh, gradientMesh, MAX_CURVATURE_TAG, 1.0f);

		vtkSmartPointer<vtkPolyData> minCutMesh = mlc.minCutConjunction(minCutMeshLeaks, minCutMeshInteriorLeaks);
		mlc.attributeDilation(minCutMesh, ATTRIBUTE_DILATION_RADIUS);
		vtkSmartPointer<vtkPolyData> correctedMesh1 = mlc.laplaceInterpolation(minCutMesh);
		vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
		normals->SetInputData(correctedMesh1);
		normals->FlipNormalsOn();
		normals->Update();
		vtkSmartPointer<vtkPolyData> correctedMesh = normals->GetOutput();
		cout << "Finished fixing mesh! Wow...." << endl;
		cout << "Writing mesh to file! laplaceMesh.vtk" << endl;
		mlc.writePolyData(correctedMesh, "laplaceMesh.vtk");
		cout << "We won!" << endl;
	}

	virtual void OnKeyDown()
	{
		std::string key = this->GetInteractor()->GetKeySym();
		if (key.compare("g") == 0) {
			this->EndRotate();
		}else if (key.compare("h") == 0){
			this->StartRotate();
		}else if (key.compare("space") == 0){
			this->RemoveLeaks();
		}
	}

	virtual void OnLeftButtonDown()
	{
		MakeAnnotation(FOREGROUND);
		vtkInteractorStyleJoystickCamera::OnLeftButtonDown();
	}

	void MakeAnnotation(int type){
		int* clickPos = this->GetInteractor()->GetEventPosition();

		// Pick from this location.
		vtkSmartPointer<vtkPropPicker>  picker =
			vtkSmartPointer<vtkPropPicker>::New();
		picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

		double* pos = picker->GetPickPosition();
		std::cout << "Pick position (world coordinates) is: "
			<< pos[0] << " " << pos[1]
			<< " " << pos[2] << std::endl;

		std::cout << "Picked actor: " << picker->GetActor() << std::endl;
		vtkPolyDataMapper* mapper = (vtkPolyDataMapper*)(this->GetDefaultRenderer()->GetActors()->GetLastActor()->GetMapper());
		vtkIntArray* scalars = (vtkIntArray*)mapper->GetInput()->GetPointData()->GetScalars();
		vtkPolyData* mesh = mapper->GetInput();
		vtkIdType id = mesh->FindPoint(pos);
		cout << "The Point Is: " << id << endl;
		cout << "The Pos Is: " << clickPos[0] << "," << clickPos[1] << "," << pos[2] << endl;
		cout << "Number of components is: " << scalars->GetNumberOfComponents() << endl;
		if (id > -1 && _currSource == -1){
			_currSource = id;
			cout << "Start set." << endl;
		}
		else if (id > -1 && _currSource != -1){
			//Make a line connection.
			cout << "End set." << endl;
			vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
				vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
			dijkstra->SetInputData(mesh);
			dijkstra->SetStartVertex(_currSource);
			dijkstra->SetEndVertex(id);
			dijkstra->Update();
			cout << "Found Shortest path!" << endl;
			vtkPolyData* path = dijkstra->GetOutput();
			cout << "Num. of points: " << path->GetNumberOfPoints() << endl;
			for (int i = 0; i < path->GetNumberOfPoints(); i++){
				scalars->SetValue(mesh->FindPoint(path->GetPoint(i)), type);
			}
			_currSource = -1;
		}
		if (id > -1){
			//cout << "component size: " << scalars->GetElementComponentSize() << endl;
			cout << scalars->GetValue(id) << endl;
			scalars->SetValue(id, type);
			cout << scalars->GetValue(id) << endl;
		}
		scalars->Modified();
		mesh->Modified();
		mapper->Update();
	}

private:

};

vtkStandardNewMacro(MouseInteractorStyle2);

MarchingCubes::MarchingCubes(vtkStructuredPoints* selection) {

	this->_selection = selection;

	_surface = vtkMarchingCubes::New();
	_renderer = vtkRenderer::New();
	_renderWindow = vtkRenderWindow::New();
	_mapper = vtkPolyDataMapper::New();
	_actor = vtkActor::New();
	_mask = vtkImageMaskBits::New();
	_interactor = vtkRenderWindowInteractor::New();

	//filter segblock, only SEGMENTATION will be presented
	_mask->AddInputData(_selection);
	_mask->SetMask(static_cast<unsigned int>(FOREGROUND));
	_mask->SetOperationToAnd();

	double bounds[6];
	_selection->GetBounds(bounds);
	_surface->SetInputConnection(_mask->GetOutputPort());
	_surface->ComputeNormalsOn();
	_surface->ComputeScalarsOff();
	_surface->SetValue(0, ISO_VALUE);
	_surface->Update();
	_renderer->SetBackground(0, 0, 0.2);


	vtkSmartPointer<myVtkInteractorStyleImage3D> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage3D>::New();
	_renderWindow->SetWindowName("Mesh Viewer");

	//_interactor->SetInteractorStyle(myInteractorStyle);
	_interactor->CreateRepeatingTimer(UPDATE_SLICE_TIMER);

	_interactor->SetRenderWindow(_renderWindow);
	_interactor->Initialize();
	_renderWindow->SetInteractor(_interactor);
	
	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetNumberOfTableValues(3);
	lut->SetRange(0, 2);
	lut->SetTableValue((vtkIdType)NOT_ACTIVE, 0.9, 0.9, 0.9, 1.0);
	lut->SetTableValue((vtkIdType)FOREGROUND, 1, 0, 0, 1.0);
	lut->SetTableValue((vtkIdType)BACKGROUND, 0, 1, 0, 1.0);
	lut->Build();
	vtkPolyData* mesh = _surface->GetOutput();


	//try close holes:
	//vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
	//	vtkSmartPointer<vtkFillHolesFilter>::New();
	//fillHolesFilter->SetInputData(mesh);
	//fillHolesFilter->SetHoleSize(1000.0);
	//fillHolesFilter->Update();

	// Make the triangle windong order consistent
	//vtkSmartPointer<vtkPolyDataNormals> normals =
	//	vtkSmartPointer<vtkPolyDataNormals>::New();
	//normals->SetInputData(fillHolesFilter->GetOutput());
	//normals->ConsistencyOn();
	//normals->SplittingOff();
	//normals->Update();

	//end close holes

	_mapper->SetScalarRange(0,2);
	_mapper->SetLookupTable(lut);
	//cout << "Number of points is: " << mesh->GetNumberOfPoints() << endl;
	vtkSmartPointer<vtkIntArray> mesh_colors =
		vtkSmartPointer<vtkIntArray>::New();
	// Add the colors we created to the colors array	
	//mesh_colors->SetNumberOfValues(normals->GetOutput()->GetNumberOfPoints());
	mesh_colors->SetNumberOfValues(mesh->GetNumberOfPoints());
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++){
		mesh_colors->SetValue(i, NOT_ACTIVE);
	}
	mesh_colors->SetName("mesh_colors");
	mesh->GetPointData()->SetScalars(mesh_colors);
	mesh_colors->Modified();
	// Restore the original normals
	//normals->GetOutput()->GetPointData()->
	//	SetNormals(mesh->GetPointData()->GetNormals());
	//normals->GetOutput()->GetPointData()->
	//	SetScalars(mesh_colors);
	//_mapper->SetInputData(normals->GetOutput());
	_mapper->SetInputData(mesh);
	//_mapper->SetScalarModeToUsePointData();
	_mapper->Update();

	_actor->SetMapper(_mapper);

	//_actor->GetProperty()->SetDiffuseColor(1.0, 0.3882, 0.2784);
	
	// Create a cursor (sphere)
	//vtkSmartPointer<vtkSphereSource> sphereSource =
	//	vtkSmartPointer<vtkSphereSource>::New();
	//sphereSource->SetCenter(_mapper->GetCenter());
	//sphereSource->SetRadius(5.0);

	//vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
	//	vtkSmartPointer<vtkPolyDataMapper>::New();
	//sphereMapper->SetInputData(sphereSource->GetOutput());

	//vtkSmartPointer<vtkActor> sphereActor =
	//	vtkSmartPointer<vtkActor>::New();
	//sphereActor->SetMapper(sphereMapper);
	////Set the color of the sphere
	//sphereActor->GetProperty()->SetColor(0.0, 1.0, 0.0); 
	//CrosshairFactory* crossFac = CrosshairFactory::getInstance();
	//vtkSmartPointer<vtkActor> crosshair = crossFac->makeCrosshair(_actor->GetBounds());//_surface->GetOutput()->GetBounds());
	//crosshair->GetMapper()->Update();
	//crosshair->SetPickable(false);
	//_actor->SetPickable(true);
	vtkSmartPointer<vtkRenderer> overlayRenderer = vtkSmartPointer<vtkRenderer>::New();
	//vtkInteractorStyle::SafeDownCast(_interactor->GetInteractorStyle())->AutoAdjustCameraClippingRangeOn();
	//vtkSmartPointer<vtkCamera> cameraOverlay = vtkSmartPointer<vtkCamera>::New();
	_renderWindow->SetNumberOfLayers(2);
	overlayRenderer->SetLayer(1);
	_renderer->SetLayer(0);
	_renderer->SetInteractive(1);
	overlayRenderer->SetInteractive(0);
	// sync camera
	//cameraOverlay = overlayRenderer->GetActiveCamera();
	//cameraOverlay->ShallowCopy(_renderer->GetActiveCamera());
	overlayRenderer->SetInteractive(0);
	//_renderWindow->AddRenderer(overlayRenderer);
	_renderWindow->AddRenderer(_renderer);
	//overlayRenderer->AddActor(crosshair);
	//
	vtkSmartPointer<vtkCursor2D> cursor =
		vtkSmartPointer<vtkCursor2D>::New();
	cursor->SetModelBounds(-10, 10, -10, 10, 0, 0);
	cursor->AllOn();
	cursor->OutlineOff();
	cursor->Update();

	vtkSmartPointer<vtkPolyDataMapper> cursorMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	cursorMapper->SetInputData(cursor->GetOutput());
	vtkSmartPointer<vtkActor> cursorActor =
		vtkSmartPointer<vtkActor>::New();
	cursorActor->GetProperty()->SetColor(1, 0, 0);
	cursorActor->SetMapper(cursorMapper);
	overlayRenderer->AddActor(cursorActor);
	
	vtkSmartPointer<MouseInteractorStyle2> style =
		vtkSmartPointer<MouseInteractorStyle2>::New();
	style->SetDefaultRenderer(_renderer);
	style->Modified();

	_interactor->SetInteractorStyle(style);
	_renderer->AddActor(_actor);
	_renderWindow->SetCurrentCursor(VTK_CURSOR_CROSSHAIR);
	// Render and interact
	//_renderWindow->Render();
	_interactor->Initialize();
	_renderWindow->Start();
	_renderWindow->Render();
	//_renderer->ResetCamera();
	//overlayRenderer->ResetCamera();
	//_renderWindow->Render();
	//myInteractorStyle->Initialize("Output.obj",overlayRenderer);
	_interactor->Start();
}