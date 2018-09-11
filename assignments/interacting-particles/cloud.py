
if __name__ == "__main__":
    import sys

    num_proc = 0
    Np = 0
    Nt = 0
    Nint = 0
    k = 0.
    d = 0.
    dt = 0.
    s = 0
    gifname = None
    filenames = []
    for line in sys.stdin:
        import json
        import vtk
        import numpy as np

        obj = json.loads(line)
        if not num_proc:
            Np = obj['num_points']
            dt = obj['dt']
            Nt = obj['num_steps']
            Nint = obj['step_chunk']
            k = obj['k']
            d = obj['d']
            gifname = obj['gifname']

            renderer = vtk.vtkRenderer()
            renderer.SetBackground(1.,1.,1.)
            renderer.ResetCamera()
            camera = vtk.vtkCamera()
            camera.SetFocalPoint(0.,0.,0.)
            camera.SetPosition(-2.,-2.5,1.)
            camera.SetViewUp(0.,0.,1.)

            renderWindow = vtk.vtkRenderWindow()
            renderWindow.OffScreenRenderingOn()
            renderWindow.SetSize(500,370)
            renderWindow.AddRenderer(renderer)


            polyData = vtk.vtkPolyData()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(polyData)
            mapper.SetColorModeToDefault()
            mapper.SetScalarRange(0.,float(Np-1))
            mapper.SetScalarVisibility(1)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetPointSize(3)

            boxPolyData = vtk.vtkPolyData()
            box = np.array([[-0.5,-0.5,-0.5],
                            [ 0.5,-0.5,-0.5],
                            [-0.5, 0.5,-0.5],
                            [ 0.5, 0.5,-0.5],
                            [-0.5,-0.5, 0.5],
                            [ 0.5,-0.5, 0.5],
                            [-0.5, 0.5, 0.5],
                            [ 0.5, 0.5, 0.5]])
            pts = [(0,1), (2,3), (4,5), (6,7), (0,2), (1,3), (4,6), (5,7), (0,4), (1,5), (2,6), (3,7)]
            linecolor = (255, 255, 255)
            cube = vtk.vtkPolyData()
            cubepoints = vtk.vtkPoints()
            cubeedges = vtk.vtkCellArray()
            cubecolors = vtk.vtkUnsignedCharArray()
            cubecolors.SetNumberOfComponents(3)
            cubecolors.SetName("Colors")
            for i in range(8):
                cubepoints.InsertPoint(i, box[i,:])
            for i in range(12):
                edge = vtk.vtkLine()
                edge.GetPointIds().SetId(0,pts[i][0])
                edge.GetPointIds().SetId(1,pts[i][1])
                cubeedges.InsertNextCell(edge)
                cubecolors.InsertNextTuple3(0.,0.,0.)
            cube.SetPoints(cubepoints)
            cube.SetLines(cubeedges)
            cube.GetCellData().SetScalars(cubecolors)
            cubemapper = vtk.vtkPolyDataMapper()
            cubemapper.SetInputData(cube)
            cubeactor = vtk.vtkActor()
            cubeactor.SetMapper(cubemapper)
            renderer.AddActor(cubeactor)

            #renderWindowInteractor = vtk.vtkRenderWindowInteractor()
            #renderWindowInteractor.SetRenderWindow(renderWindow)

            renderer.AddActor(actor)

            cornerAnnotation = vtk.vtkCornerAnnotation()
            cornerAnnotation.SetLinearFontScaleFactor(2)
            cornerAnnotation.SetNonlinearFontScaleFactor(1)
            cornerAnnotation.SetMaximumFontSize(20)
            cornerAnnotation.GetTextProperty().SetColor(0,0,0)
            renderer.AddViewProp(cornerAnnotation)
            renderer.SetActiveCamera(camera)

            imageFilter = vtk.vtkWindowToImageFilter()
            imageFilter.SetInput(renderWindow)
            imageFilter.SetInputBufferTypeToRGB()
            imageFilter.ReadFrontBufferOff()
            imageFilter.Update()

            ogvwriter = vtk.vtkOggTheoraWriter()
            ogvwriter.SetInputConnection(imageFilter.GetOutputPort())
            ogvwriter.SetFileName(f'{gifname}.ogv')
            ogvwriter.Start()
        else:
            cornerAnnotation.SetText(0,f'k = {k}, d = {d}, dt = {dt}, T = {(num_proc - 1)* dt:9.2g}')

            try:
                X = np.array(obj['X'])
            except:
                ogvwriter.End()
                break

            X = X.T
            points = vtk.vtkPoints()
            cells = vtk.vtkCellArray()
            scalars = vtk.vtkFloatArray()
            polyData.SetPoints(points)
            polyData.SetVerts(cells)
            for i in range(X.shape[0]):
                x = X[i,:]
                for d in range(3):
                    x[d] = x[d] - round(x[d])
                pointID = points.InsertNextPoint(x)
                cells.InsertNextCell(1)
                cells.InsertCellPoint(pointID)
                scalars.InsertValue(pointID,float(pointID))
                points.Modified()
                cells.Modified()
            polyData.GetPointData().SetScalars(scalars)
            renderWindow.Render()
            imageFilter.Modified()
            ogvwriter.Write()

        num_proc = num_proc + 1



