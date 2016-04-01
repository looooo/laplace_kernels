#include "vtkWriter.h"
#include "laplaceKern2D.h"
#include "laplaceKern3D.h"
#include <vector>

using namespace laplaceKern3D;
using namespace std;

VectorPtr a ( new Vector(-1,-1, 0));
VectorPtr b ( new Vector(1, -1, 0));
VectorPtr c ( new Vector(1,  1, 0));
VectorPtr d ( new Vector(-1, 1, 0));

Panel pan = Panel(vector<VectorPtr> {a, b, c, d});
Edge edge = Edge(a, b);

double test_monopole(double x, double y, double z)
{
    Vector v (x, y, z);
    return monopole(v, pan.center);
}
double test_monopole_0_n0(double x, double y, double z)
{
    Vector v (x, y, z);
    return monopole_0_n0(v, pan);
}
double test_dipole(double x, double y, double z)
{
    Vector v (x, y, z);
    return dipole(v, pan.center, pan.n);
}
double test_dipole_0_n0(double x, double y, double z)
{
    Vector v (x, y, z);
    return dipole_0_n0(v, pan);
}
double test_dipole_0_sphere(double x, double y, double z)
{
    Vector v (x, y, z);
    return dipole_0_sphere(v, pan);
}
Vector test_monopole_v(double x, double y, double z)
{
    Vector v (x, y, z);
    return monopole_v(v, pan.center);
}
Vector test_monopole_0_n0_v(double x, double y, double z)
{
    Vector v (x, y, z);
    return monopole_0_n0_v(v, pan);
}
Vector test_dipole_v(double x, double y, double z)
{
    Vector v (x, y, z);
    return dipole_v(v, pan.center, pan.n);
}
Vector test_dipole_0_n0_v(double x, double y, double z)
{
    Vector v (x, y, z);
    return dipole_0_n0_v(v, pan);
}
Vector test_vortex_v(double x, double y, double z)
{
    Vector v (x, y, z);
    return vortex_v(v, edge);
}
Vector test_vortex_halfinfinity_v(double x, double y, double z)
{
    Vector v (x, y, z);
    return vortex_half_infinity_v(v, *edge.v1, edge.t);
}
Vector test_monopole_0_vsaero_v(double x, double y, double z)
{
    Vector v(x, y, z);
    return monopole_0_vsaero_v(v, pan);
}
Vector test_dipole_0_vsaero_v(double x, double y, double z)
{
    Vector v(x, y, z);
    return dipole_0_vsaero_v(v, pan);
}

int main ( int, char *[] )
{
    VtkWriter writer("elements3D.vtk");
    writer.writePoints(-3, 3, 30,
                       -3, 3, 30,
                       -3, 3, 30);
    writer.writeScalar(test_monopole, "monopole");
    writer.writeScalar(test_monopole_0_n0, "monopole_0_n0");
    writer.writeScalar(test_dipole, "dipole");
    writer.writeScalar(test_dipole_0_n0, "dipole_0_n0");
    writer.writeScalar(test_dipole_0_sphere, "dipole_0_sphere");
    
    writer.writeVector3(test_monopole_v, "monopole_v");
    writer.writeVector3(test_monopole_0_n0_v, "monopole_0_n0_v");
    writer.writeVector3(test_dipole_v, "dipole_v");
    writer.writeVector3(test_dipole_0_n0_v, "dipole_0_n0_v");
    writer.writeVector3(test_vortex_v, "vortex_v");
    writer.writeVector3(test_vortex_halfinfinity_v, "vortex_half_infinity_v");
    writer.writeVector3(test_monopole_0_vsaero_v, "monopole_0_vsaero_v");
    writer.writeVector3(test_dipole_0_vsaero_v, "dipole_0_vsaero_v");
    
    writer.createFile();
}