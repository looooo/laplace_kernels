# include "laplace2DKernel.h"



double doublet_2_0(Vector2& target, Panel2& source)
{
    double pn = (target - source.center).dot(source.n);
    double dx1 = (*source.points[0] - target).dot(source.t);
    double dx2 = (*source.points[1] - target).dot(source.t);
    if (pn * pn < coresize and (dx1 * dx1 + dx2 * dx2) <= (source.l * source.l))
    {
        return -0.5;
    }
    return - (atan2(pn , dx2) - atan2(pn , dx1)) / (2. * M_PI);  
}


double source_2_0(Vector2& target, Panel2& source)
{ 
    double pn, dx1, dx2, a, b, c;
    pn = (target - source.center).dot(source.n);
    dx1 = (*source.points[0] - target).dot(source.t);
    dx2 = (*source.points[1] - target).dot(source.t);
    if (pn == 0 and dx1 == 0)
    {
        a = 0;
    }
    else 
    {
        a = dx1 * log(dx1 * dx1 + pn * pn);
    }
    if (pn == 0 and dx2 == 0)
    {
        b = 0;
    }
    else 
    {
        b = - dx2 * log(dx2 * dx2 + pn * pn);
    }
    c = 2 * pn * (atan2(pn, dx2) - atan2(pn, dx1));
    return -(a + b) / (4. * M_PI);
}

Vector2 doublet_2_0_v(Vector2& target, Panel2& source)
{
    Vector2 u, w, r1, r2;
    double pn, dx1, dx2;
    pn = (target - source.center).dot(source.n);
    r1 = target - *source.points[0];
    r2 = target - *source.points[1];
    dx1 = r1.dot(source.t);
    dx2 = r2.dot(source.t);
    if (pn == 0)
    {
        if (dx1 == 0 or dx2 == 0) //the corner of the panel is singular
        {
            return r1 * 0;
        }
        u = 0 * source.t;
        w = (1 / dx1 - 1 / dx2) * source.n;
    }
    else
    {
        u = (pn / r1.dot(r1) - pn / r2.dot(r2)) * source.t;
        w = (dx1 / r1.dot(r1)- dx2 / r2.dot(r2)) * source.n;
    }
    return - (w - u) / M_PI / 2;
}

Vector2 source_2_0_v(Vector2& target, Panel2& source)
{ 
    double pn = (target - source.center).dot(source.n);
    if (pn == 0)
    {
        return -0.5 * source.n;
    }
    double dx1 = (target - *source.points[0]).dot(source.t);
    double dx2 = (target - *source.points[1]).dot(source.t);
    
    Vector2 u = log((dx1 * dx1 + pn * pn) / (dx2 * dx2 + pn * pn) ) / 2  * source.t;
    Vector2 w = (atan2(pn, dx2) - atan2(pn, dx1)) * source.n;
    return (u + w) / M_PI / 2;
}

double vortex_2(Vector2& target, Vector2& source, Vector2 direction)
{
    direction.normalize();
    double pn = (target - source).dot(normal2(direction));
    double dx = (source - target).dot(direction);
    return - atan2(pn , dx) / (2. * M_PI);  
}

Vector2 vortex_2_v(Vector2& target, Vector2& source)
{
    Vector2 r = target - source;
    if (r.norm() < coresize)
    {
        r.normalize();
        r *= coresize;
    }
    return normal2(r / r.dot(r) / 2 / M_PI);
}


double doublet_2_1(Vector2& target, Panel2& source, bool left)
{
    /* return a linear distribution doublet distribution from 1..0 if left is true
    * otherwise the distribution is from 0...1
    * for a distribution from mue2 to mue1 this function has to be called twice 
    *          mue1 * doublet_2_1(,,True) + mue2 * doublet_2_1(,,false)
    * or call a constant doublet + this linear doublet
    *          mue1 * doublet_2_0() + (mue2 - mue1) * doublet_2_1(,,false)*/
    
    double pn = (target - source.center).dot(source.n);
    Vector2 p1_t = target - *source.points[0];
    Vector2 p2_t = target - *source.points[1];
    double dx1 = p1_t.dot(source.t);
    double dx2 = p2_t.dot(source.t);
    double r1_squared = p1_t.dot(p1_t);
    double r2_squared = p2_t.dot(p2_t);
    double phi1 = atan2(pn, dx1);
    double phi2 = atan2(pn, dx2);
    if (pn * pn < coresize)
    {
        if (r1_squared < (source.l * source.l) and
        r2_squared < (source.l * source.l))
        {
            if (left){return -0.5 * p2_t.norm() / source.l;}
            else {return -0.5 * p1_t.norm() / source.l;}
        }
        return 0.;
    }
    if (left)
    {
        return  -0.5 / M_PI / source.l * (dx1 * (phi2 - phi1) + pn / 2. * log(r2_squared / r1_squared));
    }
    return 0.5 / M_PI / source.l * (dx2 * (phi2 - phi1) + pn / 2. * log(r2_squared / r1_squared)); 
}

Vector2 doublet_2_1_v(Vector2& target, Panel2& source, bool left)
{
    double pn = (target - source.center).dot(source.n);
    Vector2 p1_t = target - *source.points[0];
    Vector2 p2_t = target - *source.points[1];
    double dx1 = p1_t.dot(source.t);
    double dx2 = p2_t.dot(source.t);
    double r1_squared = p1_t.dot(p1_t);
    double r2_squared = p2_t.dot(p2_t);
    double phi1 = atan2(pn, dx1);
    double phi2 = atan2(pn, dx2);
    if (r1_squared == 0 or r2_squared == 0)
    {
        return p1_t * 0;
    }
    if (left)
    {
        Vector2 u = ((phi2 - phi1) * source.l + pn / r2_squared) * source.t;
        Vector2 w = (0.5 * source.l * log(r2_squared / r1_squared)  - dx2 / r2_squared) * source.n;
        return -0.5 / M_PI * (u + w);
    }
    else
    {
        Vector2 u = ((phi2 - phi1) * source.l + pn / r1_squared) * source.t;
        Vector2 w = (0.5 * source.l * log(r2_squared / r1_squared)  - dx1 / r1_squared) * source.n;
        return 0.5 / M_PI * (u + w);
    }
    
}


double source_2_1(Vector2& target, Panel2& source, bool left)
{
    double pn = (target - source.center).dot(source.n);
    Vector2 p1_t = target - *source.points[0];
    Vector2 p2_t = target - *source.points[1];
    double dx1 = p1_t.dot(source.t);
    double dx2 = p2_t.dot(source.t);
    double r1_squared = p1_t.dot(p1_t);
    double r2_squared = p2_t.dot(p2_t);
    double phi1 = atan2(pn, dx1);
    double phi2 = atan2(pn, dx2);    
}


double source_2_1_v(Vector2& target, Panel2& source, bool left)
{
    
}

double source_2(Vector2& target, Vector2& source)
{
    if ((target - source).norm() == 0){return 0;}
    return 1. / 2. / M_PI * log((target - source).norm());
}

Vector2 source_2_v(Vector2& target, Vector2& source)
{
    Vector2 r= target - source;
    if (r.norm() == 0){return r * 0;}
    return 1. / 2. / M_PI * r / r.dot(r);
}

double doublet_2(Vector2& target, Vector2& source, Vector2 direction)
{
    Vector2 r = target - source;
    direction.normalize();
    if (r.norm() == 0){return 0;}
    return - 1. / 2. / M_PI * direction.dot(r) / r.dot(r);
}

Vector2 doublet_2_v(Vector2& target, Vector2& source, Vector2 direction)
{
    direction.normalize();
    Vector2 r = (target - source);
    double r_dot_r = r.dot(r);
    if (r_dot_r == 0){return r;}
    r.normalize();
    Vector2 v = r.dot(direction) * r + r.dot(normal2(direction)) * normal2(r);
    return 1. / 2. / M_PI / r_dot_r * v;
}