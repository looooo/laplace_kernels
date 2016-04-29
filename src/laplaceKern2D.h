#ifndef laplace2DKernel_H
#define laplace2DKernel_H

#include <Eigen/Core>
#include <vector>
#include <memory>

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


namespace laplaceKern2D
{
    typedef Eigen::Vector2d Vector;
    const double coresize = 1.0e-13;
    typedef std::shared_ptr<Vector> VectorPtr;
  
    class Panel{        
    public:
        Panel(VectorPtr p1, VectorPtr p2);
        std::vector<VectorPtr> points;
        Vector center, normal, tangent;
        double area = 0.;
    };
    



    /************************-2D ELEMENT INFLUENCE-*********************/

    double source(const Vector& target, const Vector& source);
    double doublet(const Vector& target, const Vector& source, Vector direction);
    double vortex(const Vector& target, const Vector& source, Vector direction);
    double source_0(const Vector& target, const Panel& source);
    double doublet_0(const Vector& target, const Panel& source);

    Vector source_v(const Vector& target, const Vector& source);
    Vector doublet_v(const Vector& target, const Vector& source, Vector direction);
    Vector vortex_v(const Vector& target, const Vector& source);
    Vector source_0_v(const Vector& target, const Panel& source);
    Vector doublet_0_v(const Vector& target, const Panel& source);
}
#endif