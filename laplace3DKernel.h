#ifndef laplace3DKernel_H
#define laplace3DKernel_H

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <sstream>

/************************-NAMING CONVENTION-************************/

/*
the influence function
    (1) element (doublet, monopol, vortex), 
    (2) dimension(2,3),
    (3) order (0,1,2), if no value -> point element
    (4) integration_sheme,
    [5] v... velocity

the arguments
    first argument is allways the target object (point in the space)
    second argument is the perturbation object (panel, line, point)

variable names
    vector in space                    point        or  pnt
    a Panel                            panel        or  pan
    the influence object               source       or  src
    the target object                  target       or  trg
    a monopol or monopol panel         monopol      or  mop
    a doublet or doublet panel         doublet      or  dip
    a velocity vector                  velocity     or  vel
    the potential                      potential    or  pot
    the jump in potential                               mue
    the jump of the pot derivative     sigma        or  sig
    a influence strength               .._influence or  .._infl

*/


using namespace std;
namespace laplaceKern
{
    typedef Eigen::Vector3d Vector3;
    double coresize = 0.0000000000000000000001;

    class Panel3{
    public:
        Panel3(vector<Vector3*>);
        vector<Vector3*> points;
        Vector3 center, n, m, l;
        double area = 0.;
        double mue;
        double sigma;
    };
    
    class Edge{;
    public:
        Edge(Vector3*, Vector3*);
        Vector3* v1;
        Vector3* v2;
        double gamma;
    };

    /************************-3D ELEMENT INFLUENCE-*********************/

    void doublet_3_0_n0(Vector3& target, Panel3* source, double& dip_infl);
    void doublet_3_0_vsaero(Vector3 & target, Panel3* source, double& dip_infl);
    double doublet_3_0_sphere(Vector3& target, Panel3* source, double& dip_infl);
    void doublet_src_3_0_n0(Vector3& target, Panel3* source, double& dip_infl, double& mop_infl);
    Vector3 doublet_src_3_0_vsaero_v(Vector3& target, Panel3* source, double mue, double sigma);
    Vector3 src_3_0_vsaero_v(Vector3& target, Panel3* source, double sigma);
    void doublet_src_3_0_vsaero(Vector3& target, Panel3* source, double& dip_infl, double& mop_infl);
    Vector3 vortex_3_0_v(Vector3& target, Vector3& source_point_0, Vector3& source_point_1);
    Vector3 vortex_3_0_v(Vector3& target, Edge& e);
    Vector3 vortex_3_0_half_infinity_v(Vector3& target, Vector3& source_direction, Vector3& source_point);
    double doublet_3(Vector3& target, Vector3& source, Vector3& normal);
    Vector3 doublet_3_v(Vector3& target, Vector3& source, Vector3& normal);

};
#endif