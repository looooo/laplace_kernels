#include "vtkWriter.h"
#include "laplaceKern2D.h"
#include "laplaceKern3D.h"

using namespace laplaceKern2D;
using namespace std;

VectorPtr a (new Vector(-1, 0));
VectorPtr b (new Vector( 1, 0));
Panel pan = Panel(a, b);

laplaceKern3D::Vector toVector3(Vector vec)
{
    return laplaceKern3D::Vector(vec[0], vec[1], 0);
}

double test_source(double x, double y, double z)
{
    Vector v (x, y);
    return source(v, pan.center);
}
double test_doublet(double x, double y, double z)
{
    Vector v (x, y);
    return doublet(v, pan.center, pan.t);
}
double test_vortex(double x, double y, double z)
{
    Vector v (x, y);
    return vortex(v, pan.center, pan.t);
}
double test_doublet_0(double x, double y, double z)
{
    Vector v (x, y);
    return doublet(v, pan.center, pan.t);
}
double test_source_0(double x, double y, double z)
{
    Vector v (x, y);
    return source_0(v, pan);
}
laplaceKern3D::Vector test_source_v(double x, double y, double z)
{
    Vector v (x, y);
    return toVector3(source_v(v, pan.center));
}
laplaceKern3D::Vector test_doublet_v(double x, double y, double z)
{
    Vector v (x, y);
    return toVector3(doublet_v(v, pan.center, pan.n));
}
laplaceKern3D::Vector test_vortex_v(double x, double y, double z)
{
    Vector v (x, y);
    return toVector3(vortex_v(v, pan.center));
}
laplaceKern3D::Vector test_source_0_v(double x, double y, double z)
{
    Vector v (x, y);
    return toVector3(source_0_v(v, pan));
}
laplaceKern3D::Vector test_doublet_0_v(double x, double y, double z)
{
    Vector v (x, y);
    return toVector3(doublet_0_v(v, pan));
}


int main ( int, char *[] )
{
    VtkWriter writer("elements2D.vtk");
    writer.writePoints(-3, 3, 200,
                       -3, 3, 200,
                       -0, 0, 0);
    writer.writeScalar(test_source, "source");
    writer.writeScalar(test_doublet, "doublet");
    writer.writeScalar(test_vortex, "vortex");
    writer.writeScalar(test_source_0, "source_0");
    writer.writeScalar(test_doublet_0, "doublet_0");

    writer.writeVector3(test_source_v, "source_v");
    writer.writeVector3(test_doublet_v, "doublet_v");
    writer.writeVector3(test_vortex_v, "vortex_v");
    writer.writeVector3(test_source_0_v, "source_0_v");
    writer.writeVector3(test_doublet_0_v, "doublet_0_v");
    writer.createFile();
}