#ifndef laplace2DKernel_H
#define laplace2DKernel_H

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

namespace laplaceKern
{
    typedef Eigen::Matrix<double, 2, 1, Eigen::DontAlign> Vector2;
    double coresize = 0.0000000000000000000001;

    class Panel2{        
    public:
        Panel2(Vector2* p1, Vector2* p2)
        vector<Vector2*> points;
        Vector2 center, n, t;
        double area = 0.;
        double mue;
        double sigma;
    };


    /************************-2D ELEMENT INFLUENCE-*********************/

    double doublet_2_0(Vector2& target, Panel2& source);
    double source_2_0(Vector2& target, Panel2& source);
    Vector2 doublet_2_0_v(Vector2& target, Panel2& source);
    Vector2 source_2_0_v(Vector2& target, Panel2& source);
    double vortex_2(Vector2& target, Vector2& source, Vector2 direction);
    Vector2 vortex_2_v(Vector2& target, Vector2& source);
    double doublet_2_1(Vector2& target, Panel2& source, bool left);
    Vector2 doublet_2_1_v(Vector2& target, Panel2& source, bool left);
    double source_2_1(Vector2& target, Panel2& source, bool left);
    double source_2_1_v(Vector2& target, Panel2& source, bool left);
    double source_2(Vector2& target, Vector2& source);
    Vector2 source_2_v(Vector2& target, Vector2& source);
    double doublet_2(Vector2& target, Vector2& source, Vector2 direction);
    Vector2 doublet_2_v(Vector2& target, Vector2& source, Vector2 direction);

}
#endif