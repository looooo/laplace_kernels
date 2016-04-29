#ifndef laplace3DKernel_H
#define laplace3DKernel_H

#include <Eigen/Core>
#include <vector>
#include <memory>
#include <tuple>

#include <iostream>
#include <sstream>

/************************-NAMING CONVENTION-************************/

/*
the influence function
    (1) element (doublet, monopol, vortex),
    (2) order (0,1,2), if no value -> point element
    [3] integration_sheme,
    [5] v... velocity

the arguments
    first argument is allways the target object (point in the space)
    second argument is the perturbation object (panel, line, point)

variable names
    vector in space                    point        or  pnt
    a Panel                            panel        or  pan
    the influence object               source       or  src
    the target object                  target       or  trg
    a monopol or monopol panel         monopole     or  mon
    a doublet or doublet panel         dipole       or  dip
    a velocity vector                  velocity     or  vel
    the potential                      potential    or  pot
    the jump in potential                               mue
    the jump of the pot derivative     sigma        or  sig
    a influence strength               influence or or  infl
*/


namespace laplaceKern3D
{
    typedef Eigen::Vector3d Vector;
    typedef std::shared_ptr<Vector> VectorPtr;
    const double coresize = 1.0e-13;

    class Panel
    {
    public:
        Panel(std::vector<VectorPtr>);
        std::vector<VectorPtr> points;
        Vector center, normal, tangent_m, tangent_l;
        double area;
        double side_sum;
    };
    
    class Edge
    {
    public:
        Edge(VectorPtr v1, VectorPtr v2);
        VectorPtr v1;
        VectorPtr v2;
        Vector tangent;
    };

    // /************************-3D ELEMENT INFLUENCE-*********************/
    double monopole(const Vector& target, const Vector& source);
    double monopole_0_n0(const Vector& target, const Panel& source);
    double dipole(const Vector& target, const Vector& source, Vector normal);
    double dipole_0_n0(const Vector& target, const Panel& source);
    double dipole_0_sphere(const Vector& target, const Panel& source);    
    double dipole_0_vsaero(const Vector& target, const Panel& source);
    void dip_mon(const Vector& target, const Vector& source, Vector& direction, double& dip_infl, double& mon_infl);
    void dip_mon_0_n0(const Vector& target, const Panel& source, double& dip_infl, double& mon_infl);
    void dip_mon_0_vsaero(const Vector& target, const Panel& source, double& dip_infl, double& mon_infl);
    std::tuple<double, double> dip_mon(const Vector& target, const Vector& source, Vector direction);
    std::tuple<double, double> dip_mon_0_n0(const Vector& target, const Panel& source);
    std::tuple<double, double> dip_mon_0_vsaero(const Vector& target, const Panel& source);
    
    Vector monopole_v(const Vector& target, const Vector& source);
    Vector monopole_0_n0_v(const Vector& target, const Panel& source);
    Vector monopole_0_vsaero_v(const Vector& target, const Panel& source);
    Vector dipole_v(const Vector& target, const Vector& source, Vector direction);
    Vector dipole_0_n0_v(const Vector& target, const Panel& source);
    Vector dipole_0_vsaero_v(const Vector& target, const Panel& source);
    Vector vortex_v(const Vector& target, const Vector& source_point_0, const Vector& source_point_1);
    Vector vortex_v(const Vector& target, const Edge& e);
    Vector vortex_half_infinity_v(const Vector& target, const Vector& source_point, Vector source_direction);
    void dip_mon_v(const Vector& target, const Vector& source, Vector direction, Vector& dip_infl_v, Vector& mon_infl_v);
    void dip_mon_0_n0_v(const Vector& target, const Panel& source, Vector& dip_infl_v, Vector& mon_infl_v);
    void dip_mon_0_vsaero_v(const Vector& target, const Panel& source, Vector& dip_infl_v, Vector& mon_infl_v);
    std::tuple<Vector, Vector> dip_mon_v(const Vector& target, const Vector& source, Vector direction);
    std::tuple<Vector, Vector> dip_mon_0_n0_v(const Vector& target, const Panel& source);
    std::tuple<Vector, Vector> dip_mon_0_vsaero_v(const Vector& target, const Panel& source);
};
#endif