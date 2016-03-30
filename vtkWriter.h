#ifndef vtkWriter_H
#define vtkWriter_H

#include "laplace2DKernel.h"
#include "laplace3DKernel.h"

#include <iostream>
#include <fstream>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>


using namespace std;


class VtkWriter{
private:
    const char* file_name;
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
public:
    VtkWriter(const char* file_name);
    void writePoints(
               double x_s, double x_e, int num_x,
               double y_s, double y_e, int num_y,
               double z_s, double z_e, int num_z);
    void writeScalar(
        double (*foo)(double x, double y, double z), 
        const char* dataname="scalar_field");
    
    void writeVector2(
        laplaceKern::Vector2 (*foo)(double x, double y, double z),
        const char* dataname="vector_field");

    void writeVector3(
        laplaceKern::Vector3 (*foo)(double x, double y, double z),
        const char* dataname="vector_field");

    void createFile();
};

#endif