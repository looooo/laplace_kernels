#include "vtkWriter.h"
#include <vtkPointData.h>

VtkWriter::VtkWriter(const char* file_name)
{
    this->file_name = file_name;
}

void VtkWriter::writeScalar(
    double (*foo)(double x, double y, double z),
    const char* dataname)
{
    double point[3];
    vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPointData> pointdata = this->structuredGrid->GetPointData();
    for (int i=0; i< this->points->GetNumberOfPoints(); i++)
    {
        this->points->GetPoint(i, point);
        scalars->InsertNextTuple1(foo(point[0], point[1], point[2]));
    }
    scalars->SetName(dataname);
    pointdata->AddArray(scalars);
}

void VtkWriter::writePoints(
               double x_s, double x_e, int num_x,
               double y_s, double y_e, int num_y,
               double z_s, double z_e, int num_z)
{
    double x, y, z;
    for (double zi = 0; zi <= num_z; zi++){
        if (num_z > 0)
        {
           z = z_s + double(zi) * (z_e - z_s) / double(num_z);
        }
        else
        {
            z = z_s;
        }
        for (double yi = 0; yi <= num_y; yi++){
            if (num_y > 0)
            {
                y = y_s + double(yi) * (y_e - y_s) / double(num_y);
            }
            else
            {
                y = y_s;
            }
            for (double xi = 0; xi <= num_x; xi++){
                if (num_x > 0)
                {
                    x = x_s + double(xi) * (x_e - x_s) / double(num_x);
                }
                else
                {
                    x = x_s;
                }
                this->points->InsertNextPoint(x, y, z);
            }
        }
    }
    this->structuredGrid->SetDimensions(num_x + 1, num_y + 1, num_z + 1);
    this->structuredGrid->SetPoints(this->points);
}

void VtkWriter::writeVector2(
        laplaceKern2D::Vector (*foo)(double x, double y, double z),
        const char* dataname)
{
    double point[3];
    laplaceKern2D::Vector vec;
    vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
    vectors->SetNumberOfComponents(2);
    vtkSmartPointer<vtkPointData> pointdata = this->structuredGrid->GetPointData();
    for (int i=0; i< this->points->GetNumberOfPoints(); i++)
    {
        this->points->GetPoint(i, point);
        vec = foo(point[0], point[1], point[2]);
        vectors->InsertNextTuple2(vec.x(), vec.y());
    }
    vectors->SetName(dataname);
    pointdata->AddArray(vectors);
}

void VtkWriter::writeVector3(
        laplaceKern3D::Vector (*foo)(double x, double y, double z),
        const char* dataname)
{
    double point[3];
    laplaceKern3D::Vector vec;
    vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
    vectors->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPointData> pointdata = this->structuredGrid->GetPointData();
    for (int i=0; i< this->points->GetNumberOfPoints(); i++)
    {
        this->points->GetPoint(i, point);
        vec = foo(point[0], point[1], point[2]);
        vectors->InsertNextTuple3(vec.x(), vec.y(), vec.z());
    }
    vectors->SetName(dataname);
    pointdata->AddArray(vectors);
}

void VtkWriter::createFile()
{
    this->writer->SetFileName(this->file_name);
    this->writer->SetInputData(this->structuredGrid);
    this->writer->Write();
}
