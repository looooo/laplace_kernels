# include "laplaceKern2D.h"

using namespace std;

namespace laplaceKern2D
{
    Vector normal2(Vector t)
    {
        return Vector(-t[1], t[0]);
    }

    Panel::Panel(VectorPtr p1, VectorPtr p2)
    {
        this->points.push_back(p1);
        this->points.push_back(p2);
        this->t = *p1 - *p2;
        this->area = this->t.norm();
        this->center = (*p2 + *p1) / 2;
        this->t.normalize();
        this->n = normal2(this->t);
    }

    double source(Vector& target, Vector& source)
    {
        if ((target - source).norm() == 0){return 0;}
        return 1. / 2. / M_PI * log((target - source).norm());
    }

    double doublet(Vector& target, Vector& source, Vector& direction)
    {
        Vector r = target - source;
        direction.normalize();
        if (r.norm() == 0)
            return 0;
        return - 1. / 2. / M_PI * direction.dot(r) / r.dot(r);
    }

    double vortex(Vector& target, Vector& source, Vector& direction)
    {
        direction.normalize();
        double pn = (target - source).dot(normal2(direction));
        double dx = (source - target).dot(direction);
        return - atan2(pn , dx) / (2. * M_PI);  
    }

    double doublet_0(Vector& target, Panel& source)
    {
        double pn = (target - source.center).dot(source.n);
        double dx1 = (*source.points[0] - target).dot(source.t);
        double dx2 = (*source.points[1] - target).dot(source.t);
        if (pn * pn < coresize and (dx1 * dx1 + dx2 * dx2) <= (source.area * source.area))
        {
            return -0.5;
        }
        return - (atan2(pn , dx2) - atan2(pn , dx1)) / (2. * M_PI);  
    }

    double source_0(Vector& target, Panel& source)
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

    Vector source_v(Vector& target, Vector& source)
    {
        Vector r= target - source;
        if (r.norm() == 0){return r * 0;}
        return 1. / 2. / M_PI * r / r.dot(r);
    }

    Vector doublet_v(Vector& target, Vector& source, Vector& direction)
    {
        direction.normalize();
        Vector r = (target - source);
        double r_dot_r = r.dot(r);
        if (r_dot_r == 0){return r;}
        r.normalize();
        Vector v = r.dot(direction) * r + r.dot(normal2(direction)) * normal2(r);
        return 1. / 2. / M_PI / r_dot_r * v;
    }

    Vector vortex_v(Vector& target, Vector& source)
    {
        Vector r = target - source;
        if (r.norm() == 0)
            return Vector(0, 0);
        else if (r.norm() < coresize)
        {
            r.normalize();
            r *= coresize;
        }
        return normal2(r / r.dot(r) / 2 / M_PI);
    }

    Vector doublet_0_v(Vector& target, Panel& source)
    {
        Vector u, w, r1, r2;
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

    Vector source_0_v(Vector& target, Panel& source)
    { 
        double pn = (target - source.center).dot(source.n);

        double dx1 = (target - *source.points[0]).dot(source.t);
        double dx2 = (target - *source.points[1]).dot(source.t);
        if (pn == 0 and ((abs(dx1) + abs(dx2)) <= source.area))
        {
            return -0.5 * source.n;
        }
        Vector u = log((dx1 * dx1 + pn * pn) / (dx2 * dx2 + pn * pn) ) / 2  * source.t;
        Vector w = (atan2(pn, dx2) - atan2(pn, dx1)) * source.n;
        return (u + w) / M_PI / 2;
    }
}