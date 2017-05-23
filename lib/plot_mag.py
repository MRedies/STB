import vtk
import numpy as np
from colorsys import hsv_to_rgb


def RotX(alpha_deg):
    alpha = alpha_deg / 180 * np.pi
    
    R = np.identity(3)
    R[1,1] = np.cos(alpha)
    R[2,1] = np.sin(alpha)
    R[1,2] = -np.sin(alpha)
    R[2,2] = np.cos(alpha)

    return R

def RotZ(alpha_deg):
    alpha = alpha_deg / 180 * np.pi

    R = np.identity(3)
    R[0,0] = np.cos(alpha)
    R[1,0] = np.sin(alpha)
    R[0,1] = -np.sin(alpha)
    R[1,1] = np.cos(alpha)

    return R


def color(phi, theta):
    while(phi < 0.0):
        phi += 360

    h =  phi/360.0
    s = theta/180.0
    v = 1
    return hsv_to_rgb(h,s,v)



def make_arrow(ren, size, i_pos, phi, theta):
    pos = np.dot(RotX(-theta+90), i_pos)
    pos = np.dot(RotZ(-phi+90), pos)
    r,g,b = color(phi,theta)

    cylinder = vtk.vtkCylinderSource()
    cylinder.SetResolution(25)
    cylinder.SetHeight(0.5*size)
    cylinder.SetRadius(0.1*size)
    cylinder.SetCenter(pos[0],pos[1]-0.25*size,pos[2])


    cylinderMapper = vtk.vtkPolyDataMapper()
    cylinderMapper.SetInputConnection(cylinder.GetOutputPort())
 
    cylinderActor = vtk.vtkActor()
    cylinderActor.SetMapper(cylinderMapper)
    cylinderActor.RotateX(theta-90)
    cylinderActor.RotateZ(phi-90)
    cylinderActor.GetProperty().SetColor(r,g,b)

    ren.AddActor(cylinderActor)


    pos = np.dot(RotX(-theta+90), i_pos)
    pos = np.dot(RotZ(-phi), pos)

    cone = vtk.vtkConeSource()
    cone.SetResolution(60)
    cone.SetHeight(0.5*size)
    cone.SetCenter(pos[0]+0.25*size,pos[1],pos[2])
    cone.SetRadius(0.2*size)
    
    coneMapper = vtk.vtkPolyDataMapper()
    coneMapper.SetInputConnection(cone.GetOutputPort())


    coneActor = vtk.vtkActor()
    coneActor.SetMapper(coneMapper)
    coneActor.RotateX(theta-90)
    coneActor.RotateZ(phi)
    coneActor.GetProperty().SetColor(r,g,b)
    ren.AddActor(coneActor)

    return

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 
axes = vtk.vtkAxesActor()
ren.AddActor(axes)

phi = 0.0
theta = 0.0
pos = np.array([1.0,1.0, 1.0])
make_arrow(ren,1.0, pos, phi, theta)

ren.SetBackground(0.1, 0.2, 0.4)
renWin.SetSize(600, 600)
 
iren.Initialize()
 
ren.ResetCamera()
ren.GetActiveCamera().Zoom(1.5)
renWin.Render()
 
# Start the event loop.
iren.Start()

print("Doesn't work yet")
# print(color(0,0))
# print(color(90,0))
# print(color(0,90))
