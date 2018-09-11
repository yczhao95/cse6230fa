
if __name__ == "__main__":
    import sys
    import json
    import numpy as np
    from matplotlib import pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    from matplotlib import animation
    from mpl_toolkits.mplot3d.art3d import juggle_axes

    # from https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot#9416663
    class AnimatedScatter(object):
        """An animated scatter plot using matplotlib.animations.FuncAnimation."""
        def __init__(self, stdin=sys.stdin, numpoints=50, numsteps=1):
            self.numpoints = numpoints
            self.stream = self.data_stream()
            self.stdin = stdin
            self.numsteps = numsteps

            # Setup the figure and axes...
            self.fig = plt.figure()
            self.ax = p3.Axes3D(self.fig)
            # Then setup FuncAnimation.
            self.ani = animation.FuncAnimation(self.fig, self.update, interval=1,
                                               init_func=self.setup_plot)#, blit=True)

        def setup_plot(self):
            """Initial drawing of the scatter plot."""
            X = next(self.stream)
            self.scat = self.ax.scatter(X[0,:], X[1,:], X[2,:], c=np.array(range(X.shape[1])))#, animated=True)
            self.ax.set_xlim3d([-0.5,0.5])
            self.ax.set_ylim3d([-0.5,0.5])
            self.ax.set_zlim3d([-0.5,0.5])

            # For FuncAnimation's sake, we need to return the artist we'll be using
            # Note that it expects a sequence of artists, thus the trailing comma.
            return self.scat,

        def data_stream(self):
            while True:
                rline = self.stdin.readline()
                obj = json.loads(rline)
                try:
                    X = np.array(obj['X'])
                    for i in range(X.shape[0]):
                        for j in range(X.shape[1]):
                            X[i,j] = X[i,j] - round(X[i,j])
                    yield X
                except:
                    raise StopIteration

        def update(self, i):
            """Update the scatter plot."""
            X = next(self.stream)

            # Set x and y data...
            self.scat._offsets3d = (np.ma.ravel(X[0,:]), np.ma.ravel(X[1,:]), np.ma.ravel(X[2,:]))
            # Set colors..
            # self.scat.set_array(range(X.shape[1]))

            # We need to return the updated artist for FuncAnimation to draw..
            # Note that it expects a sequence of artists, thus the trailing comma.
            #plt.draw()
            return self.scat,

        def show(self,filename=None):
            if filename:
                self.ani.save(filename,dpi=80,writer='imagemagick')
            else:
                plt.show()

    num_proc = 0
    Np = 0
    Nt = 0
    Nint = 0
    k = 0.
    d = 0.
    dt = 0.
    s = 0
    gifname = None
    firstline = sys.stdin.readline()
    obj = json.loads(firstline)

    Np = obj['num_points']
    dt = obj['dt']
    Nt = obj['num_steps']
    Nint = obj['step_chunk']
    k = obj['k']
    d = obj['d']
    gifname = obj['gifname']

    numsteps = int(Nt) // int(Nint)
    print(Nt,Nint,numsteps)
    a = AnimatedScatter(numpoints = Np, numsteps=numsteps)

    a.show(gifname)




    #for line in sys.stdin:

    #    obj = json.loads(line)
    #    if not num_proc:
    #        Np = obj['num_points']
    #        dt = obj['dt']
    #        Nt = obj['num_steps']
    #        Nint = obj['step_chunk']
    #        k = obj['k']
    #        d = obj['d']
    #        gifname = obj['gifname']

    #        renderer = vtk.vtkRenderer()
    #        renderer.SetBackground(1.,1.,1.)
    #        renderer.ResetCamera()
    #        camera = vtk.vtkCamera()
    #        camera.SetFocalPoint(0.,0.,0.)
    #        camera.SetPosition(-2.,-2.5,1.)
    #        camera.SetViewUp(0.,0.,1.)

    #        renderWindow = vtk.vtkRenderWindow()
    #        renderWindow.OffScreenRenderingOn()
    #        renderWindow.SetSize(500,370)
    #        renderWindow.AddRenderer(renderer)


    #        polyData = vtk.vtkPolyData()

    #        mapper = vtk.vtkPolyDataMapper()
    #        mapper.SetInputData(polyData)
    #        mapper.SetColorModeToDefault()
    #        mapper.SetScalarRange(0.,float(Np-1))
    #        mapper.SetScalarVisibility(1)
    #        actor = vtk.vtkActor()
    #        actor.SetMapper(mapper)
    #        actor.GetProperty().SetPointSize(3)

    #        boxPolyData = vtk.vtkPolyData()
    #        box = np.array([[-0.5,-0.5,-0.5],
    #                        [ 0.5,-0.5,-0.5],
    #                        [-0.5, 0.5,-0.5],
    #                        [ 0.5, 0.5,-0.5],
    #                        [-0.5,-0.5, 0.5],
    #                        [ 0.5,-0.5, 0.5],
    #                        [-0.5, 0.5, 0.5],
    #                        [ 0.5, 0.5, 0.5]])
    #        pts = [(0,1), (2,3), (4,5), (6,7), (0,2), (1,3), (4,6), (5,7), (0,4), (1,5), (2,6), (3,7)]
    #        linecolor = (255, 255, 255)
    #        cube = vtk.vtkPolyData()
    #        cubepoints = vtk.vtkPoints()
    #        cubeedges = vtk.vtkCellArray()
    #        cubecolors = vtk.vtkUnsignedCharArray()
    #        cubecolors.SetNumberOfComponents(3)
    #        cubecolors.SetName("Colors")
    #        for i in range(8):
    #            cubepoints.InsertPoint(i, box[i,:])
    #        for i in range(12):
    #            edge = vtk.vtkLine()
    #            edge.GetPointIds().SetId(0,pts[i][0])
    #            edge.GetPointIds().SetId(1,pts[i][1])
    #            cubeedges.InsertNextCell(edge)
    #            cubecolors.InsertNextTuple3(0.,0.,0.)
    #        cube.SetPoints(cubepoints)
    #        cube.SetLines(cubeedges)
    #        cube.GetCellData().SetScalars(cubecolors)
    #        cubemapper = vtk.vtkPolyDataMapper()
    #        cubemapper.SetInputData(cube)
    #        cubeactor = vtk.vtkActor()
    #        cubeactor.SetMapper(cubemapper)
    #        renderer.AddActor(cubeactor)

    #        #renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    #        #renderWindowInteractor.SetRenderWindow(renderWindow)

    #        renderer.AddActor(actor)

    #        cornerAnnotation = vtk.vtkCornerAnnotation()
    #        cornerAnnotation.SetLinearFontScaleFactor(2)
    #        cornerAnnotation.SetNonlinearFontScaleFactor(1)
    #        cornerAnnotation.SetMaximumFontSize(20)
    #        cornerAnnotation.GetTextProperty().SetColor(0,0,0)
    #        renderer.AddViewProp(cornerAnnotation)
    #        renderer.SetActiveCamera(camera)

    #        imageFilter = vtk.vtkWindowToImageFilter()
    #        imageFilter.SetInput(renderWindow)
    #        imageFilter.SetInputBufferTypeToRGB()
    #        imageFilter.ReadFrontBufferOff()
    #        imageFilter.Update()

    #        ogvwriter = vtk.vtkOggTheoraWriter()
    #        ogvwriter.SetInputConnection(imageFilter.GetOutputPort())
    #        ogvwriter.SetFileName(f'{gifname}.ogv')
    #        ogvwriter.Start()
    #    else:
    #        cornerAnnotation.SetText(0,f'k = {k}, d = {d}, dt = {dt}, T = {(num_proc - 1)* dt:9.2g}')

    #        try:
    #            X = np.array(obj['X'])
    #        except:
    #            ogvwriter.End()
    #            break

    #        X = X.T
    #        points = vtk.vtkPoints()
    #        cells = vtk.vtkCellArray()
    #        scalars = vtk.vtkFloatArray()
    #        polyData.SetPoints(points)
    #        polyData.SetVerts(cells)
    #        for i in range(X.shape[0]):
    #            x = X[i,:]
    #            for d in range(3):
    #                x[d] = x[d] - round(x[d])
    #            pointID = points.InsertNextPoint(x)
    #            cells.InsertNextCell(1)
    #            cells.InsertCellPoint(pointID)
    #            scalars.InsertValue(pointID,float(pointID))
    #            points.Modified()
    #            cells.Modified()
    #        polyData.GetPointData().SetScalars(scalars)
    #        renderWindow.Render()
    #        imageFilter.Modified()
    #        ogvwriter.Write()

    #    num_proc = num_proc + 1



