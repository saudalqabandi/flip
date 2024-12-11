import numpy as np
import matplotlib.pyplot as plt
import vtk


with open("config.cfg") as f:
    config = f.readlines()

config = [x for x in config if x[0] != "#"]
config = [x.split("=") for x in config if "=" in x]

cfg = {}
for x in config:
    cfg[x[0].strip()] = x[1].strip()


radius = 0.5
length = float(cfg["l"])

# Create a spherocylinder
cylinder = vtk.vtkCylinderSource()
cylinder.SetRadius(radius)
cylinder.SetHeight(length)
cylinder.SetResolution(10)

sphere1 = vtk.vtkSphereSource()
sphere1.SetRadius(radius)
sphere1.SetCenter(0, length / 2, 0)
sphere1.SetThetaResolution(10)
sphere1.SetPhiResolution(10)
sphere1.SetEndPhi(180)  # Make it a hemisphere

sphere2 = vtk.vtkSphereSource()
sphere2.SetRadius(radius)
sphere2.SetCenter(0, -length / 2, 0)
sphere2.SetThetaResolution(10)
sphere2.SetPhiResolution(10)
sphere2.SetEndPhi(180)  # Make it a hemisphere

# Append the cylinder and spheres to create a spherocylinder
append_filter = vtk.vtkAppendPolyData()
append_filter.AddInputConnection(cylinder.GetOutputPort())
append_filter.AddInputConnection(sphere1.GetOutputPort())
append_filter.AddInputConnection(sphere2.GetOutputPort())
append_filter.Update()

# Write the spherocylinder to a VTK file
writer = vtk.vtkPolyDataWriter()
writer.SetFileName(f"spherocylinder_l_{cfg['l']}.vtk")
writer.SetInputData(append_filter.GetOutput())
writer.Write()
