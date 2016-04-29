#include "laplaceKern3D.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace std;

namespace laplaceKern3D
{
    Panel::Panel(vector<VectorPtr> points)
    {
        this->points = points;
        
        // calc geo
        Vector n_i;
        int l = this->points.size();
        this->area = 0.;
        this->side_sum = 0.;
        this->center.setZero();
        this->normal.setZero();
        for (VectorPtr point: this->points)
            this->center += *point;
        this->center /= double(l);
        for (int i=0; i < l; i++)
        {
            int j = (i == l - 1 ? 0 : i+1);     // next or first
            this->side_sum += (*this->points[j] - *this->points[i]).norm();
            n_i = (this->center - *this->points[i]).cross(this->center - *this->points[j]);
            this->normal -= n_i;
            this->area += n_i.norm();
        }
        this->area /= 2;
        this->tangent_l = (this->center - (*this->points[1] + *this->points[2]) / 2);
        this->tangent_m = this->normal.cross(this->tangent_l);
        this->normal.normalize();
        this->tangent_l.normalize();
        this->tangent_m.normalize();
    }
    
    Edge::Edge(VectorPtr v1, VectorPtr v2)
    {
        this->v1 = v1;
        this->v2 = v2;
        this->tangent = *v2 - *v1;
    }


    double monopole(const Vector& target, const Vector& source)
    {
        double diff = (target - source).norm();
        if (diff < coresize)
            return 0;
        return 1 / diff / 4 / M_PI;
    }
    
    double monopole_0_n0(const Vector& target, const Panel& source)
    {
        return source.area * monopole(target, source.center);
    }
 
    double dipole(const Vector& target, const Vector& source, Vector normal)
    {
        Vector p_v = target - source;
        if (p_v.norm() < coresize)
            return 0;
        double pn = p_v.dot(normal) / normal.norm();
        return - pn / pow(p_v.norm(), 3) / 4 / M_PI;
    }

    double dipole_0_n0(const Vector& target, const Panel& source)
    {
        return source.area * dipole(target, source.center, source.normal);
    }
    
    double dipole_0_sphere(const Vector& target, const Panel& source)
    {
        Vector u1, u2, u3;
        double dip_infl = 0;
        
        u1 = source.center - target;
        if (u1.norm() == 0.)
            return -0.5;
        else if (u1.dot(source.normal) == 0.){
            // point lies on the center of the panel
            // this will return the value of the integral, with the panel wrapped around the center point
            // in n direction.

            // to have good speeds no check if point lies on the panel is made
            return 0;
        }
        u1.normalize();
        for (int i = 0; i < source.points.size(); i++){
            int j = (i == source.points.size()-1? 0 : i+1);
            u2 = *source.points[i] - target;
            u3 = *source.points[j] - target;
            if (u2.norm())
            u2.normalize();
            u3.normalize();
            dip_infl -= atan2((u1.cross(u2)).dot(u3), (1 + u2.dot(u3) + u3.dot(u1) + u1.dot(u2))) * 0.5 / M_PI;
        }
        return dip_infl;
    }
    
    double dipole_0_vsaero(const Vector & target, const Panel& source)
    {
        Vector p_v = target - source.center;
        double pn = p_v.dot(source.normal);
        int l = source.points.size();
        double rnum;
        double dnom;
        double dip_infl = 0;

        for (int i = 0; i < l; i++){
            int j = (i == l-1? 0 : i+1);
            Vector &p1 = *source.points[i];
            Vector &p2 = *source.points[j];
            Vector s = p2 - p1;
            Vector a = target - p1;
            Vector b = target - p2;
            Vector h = a.cross(s);
            double bm = b.dot(source.tangent_m);
            double am = a.dot(source.tangent_m);
            double sl = s.dot(source.tangent_l);
            double sm = s.dot(source.tangent_m);
            double al = am * sl - a.dot(source.tangent_l) * sm;
            double pa = pow(pn, 2) * sl + al * am;
            double pb = pow(pn, 2) * sl + al * bm;
            double absa = a.norm();
            double absb = b.norm();
            double abss = s.norm();
            
            if (absa > coresize && absb > coresize && absa + absb - abss >= 0)
            {
                dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
                rnum = sm * pn * (absb * pa - absa * pb);
                if (fabs(pn) < coresize)
                {
                    double pn_sign = (pn >= 0 ? -1 : 1);
                    double al_sign = (al < 0 ? 1 : -1);
                    if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
                    else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
                }
                else
                {
                    dip_infl += atan2(rnum, dnom);
                }
            }
        }
        return dip_infl * 0.25 / M_PI;
    }
    
    void dip_mon(const Vector& target, const Vector& source, Vector direction, double& dip_infl, double& mon_infl)
    {
        dip_infl = dipole(target, source, direction);
        mon_infl = monopole(target, source);
    }


    void dip_mon_0_n0(const Vector& target, const Panel& source, double& dip_infl, double& mon_infl)
    {
        Vector p_v = target - source.center;
        double diff = p_v.norm();
        double pn = p_v.dot(source.normal);
        dip_infl = -pn * source.area / pow(diff, 3) / 4 / M_PI;
        mon_infl = source.area / diff / 4 / M_PI;
    }

    void dip_mon_0_vsaero(const Vector& target, const Panel& source, double& dip_infl, double& mon_infl)
    {
        Vector p_v = target - source.center;
        double pn = p_v.dot(source.normal);
        int l = source.points.size();
        double rnum;
        double dnom;


        for (int i = 0; i < l; i++){
            int j = (i == l-1? 0 : i+1);
            Vector &p1 = *source.points[i];
            Vector &p2 = *source.points[j];
            Vector s = p2 - p1;
            Vector a = target - p1;
            Vector b = target - p2;
            Vector h = a.cross(s);
            double bm = b.dot(source.tangent_m);
            double am = a.dot(source.tangent_m);
            double sl = s.dot(source.tangent_l);
            double sm = s.dot(source.tangent_m);
            double al = am * sl - a.dot(source.tangent_l) * sm;
            double pa = pow(pn, 2) * sl + al * am;
            double pb = pow(pn, 2) * sl + al * bm;
            double absa = a.norm();
            double absb = b.norm();
            double abss = s.norm();
            
            if (absa > coresize && absb > coresize && absa + absb - abss >= 0)
            {
            dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
            rnum = sm * pn * (absb * pa - absa * pb);

            if (fabs(pn) < coresize)
            {
                double pn_sign = (pn >= 0 ? -1 : 1);
                double al_sign = (al < 0 ? 1 : -1);
                if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
                else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
            }
            else
            {
                dip_infl += atan2(rnum, dnom); 
            }
            mon_infl -= al * log((absa + absb + abss) / (absa + absb - abss)) / abss;  
            }
        }
        mon_infl += pn * dip_infl;
        dip_infl /= 4. * M_PI;
        mon_infl /= 4. * M_PI;
    }
    
    std::tuple<double, double> dip_mon(const Vector& target, const Vector& source, Vector direction)
    {
        double dip_infl_v = dipole(target, source, direction);
        double mon_infl_v = monopole(target, source);
        return std::make_tuple(dip_infl_v, mon_infl_v);
    }
    
    std::tuple<double, double> dip_mon_0_n0(const Vector& target, const Panel& source)
    {
        double dip_infl_v;
        double mon_infl_v;
        dip_mon_0_n0(target, source, dip_infl_v, mon_infl_v);
        return std::make_tuple(dip_infl_v, mon_infl_v);
    }
    
    std::tuple<double, double> dip_mon_0_vsaero(const Vector& target, const Panel& source)
    {
        double dip_infl_v;
        double mon_infl_v;
        dip_mon_0_vsaero(target, source, dip_infl_v, mon_infl_v);
        return std::make_tuple(dip_infl_v, mon_infl_v);
    }
    
    Vector monopole_v(const Vector& target, const Vector& source)
    {
        Vector p_v = target - source;
        double diff = p_v.norm();
        if (diff < coresize)
            return Vector(0, 0, 0);
        return - 1 / pow(diff, 3) / 4 / M_PI * p_v;
    }
    
    Vector monopole_0_n0_v(const Vector& target, const Panel& source)
    {
        return source.area * monopole_v(target, source.center);
    }

    Vector dipole_v(const Vector& target, const Vector& source, Vector direction)
    {
        Vector r = target - source;
        if (r.norm() < coresize)
            return Vector(0, 0, 0);
        return (direction.cross(r).cross(r) + 2. * r * r.dot(direction)) /
               (4. * M_PI * pow(r.dot(r), 5. / 2.));
    }

    Vector dipole_0_n0_v(const Vector& target, const Panel& source)
    {
        return source.area * dipole_v(target, source.center, source.normal);
    }

    Vector vortex_v(const Vector& target, const Vector& source_point_0, const Vector& source_point_1)
    {
        Vector r0 = source_point_1 - source_point_0;
        Vector r1 = target - source_point_0;
        Vector r2 = target - source_point_1;
        Vector r_cross = r1.cross(r2);
        if (r_cross.norm() > 0.001)
        {
            return r_cross * ((r0.dot(r1) / r1.norm() - r0.dot(r2) / r2.norm()) / r_cross.dot(r_cross) / 4 / M_PI);
        }
        return Vector(0.,0.,0.);
    }

    Vector vortex_v(const Vector& target, const Edge& e)
    {
        return vortex_v(target, *e.v1, *e.v2);
    }

    Vector vortex_half_infinity_v(const Vector& target, const Vector& source_point, Vector source_direction)
    {
        Vector r1 = target - source_point;
        Vector cross = r1.cross(source_direction);
        if (cross.norm() > 0.00001)
        {
            return cross / cross.dot(cross) * (1 + r1.dot(source_direction) / r1.norm()) / 4 / M_PI;
        }
        return Vector(0.,0.,0.);
    }

    Vector monopole_0_vsaero_v(const Vector& target, const Panel& source)
    {
        Vector dip_infl_v(0, 0, 0);
        Vector mon_infl_v(0, 0, 0);
        dip_mon_0_vsaero_v(target, source, dip_infl_v, mon_infl_v);
        return mon_infl_v;
    }

    Vector dipole_0_vsaero_v(const Vector& target, const Panel& source)
    {
        Vector vel_infl(0, 0, 0);
        Vector p_v = target - source.center;
        int l = source.points.size();
        for (int i = 0; i < l; i++){
            int j = (i == l-1? 0 : i+1);
            Vector p1 = *(source.points[i]);
            Vector p2 = *(source.points[j]);
            Vector a = target - p1;
            Vector b = target - p2;
            double absa = a.norm();
            double absb = b.norm();
            if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize)
            {
                vel_infl -= (a.cross(b) * (absa + absb) / 4 / M_PI / (absa * absb * (absa * absb + a.dot(b))));
            }
        }
        return vel_infl;
    }

    void dip_mon_0_vsaero_v(const Vector& target, const Panel& source, Vector& dip_infl_v, Vector& mon_infl_v)
    {
        dip_infl_v.setZero();
        mon_infl_v.setZero();
        double dip_infl = 0;
        Vector p_v = target - source.center;
        double pn = p_v.dot(source.normal);
        int l = source.points.size();
        double rnum;
        double dnom;

        for (int i = 0; i < l; i++)
        {
            int j = (i == l-1? 0 : i+1);
            Vector &p1 = *source.points[i];
            Vector &p2 = *source.points[j];
            Vector s = p2 - p1;
            Vector a = target - p1;
            Vector b = target - p2;
            Vector h = a.cross(s);
            double bm = b.dot(source.tangent_m);
            double am = a.dot(source.tangent_m);
            double sl = s.dot(source.tangent_l);
            double sm = s.dot(source.tangent_m);
            double al = am * sl - a.dot(source.tangent_l) * sm;
            double pa = pow(pn, 2) * sl + al * am;
            double pb = pow(pn, 2) * sl + al * bm;
            double absa = a.norm();
            double absb = b.norm();
            double abss = s.norm();
            
            if (absa > coresize && absb > coresize && absa + absb - abss >= 0)
            {
                dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
                rnum = sm * pn * (absb * pa - absa * pb);

                if (fabs(pn) < coresize)
                {
                    double pn_sign = (pn >= 0 ? -1 : 1);
                    double al_sign = (al < 0 ? 1 : -1);
                    if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
                    else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
                }
                else
                {
                    dip_infl += atan2(rnum, dnom);
                }
                if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
                {
                    dip_infl_v -= a.cross(b) * (absa + absb) / (absa * absb * (absa * absb + a.dot(b)));
                    mon_infl_v -= log((absa + absb + abss) / (absa + absb - abss)) / abss * (sm * source.tangent_l - sl * source.tangent_m);
                }
            }
        }
        mon_infl_v -= source.normal * dip_infl;
        mon_infl_v /= 4. * M_PI;
        dip_infl_v /= 4. * M_PI;
    }

    void dip_mon_v(const Vector& target, const Vector& source, Vector direction, Vector& dip_infl_v, Vector& mon_infl_v)
    {
        Vector r = target - source;
        double diff = r.norm();
        if (diff > coresize)
        {
            mon_infl_v = - 1 / pow(diff, 3) / 4 / M_PI * r;
            dip_infl_v = (direction.cross(r).cross(r) + 2. * r * r.dot(direction)) /
                (4. * M_PI * pow(r.dot(r), 5. / 2.));
        }
        else
        {
            mon_infl_v.setZero();
            dip_infl_v.setZero();
        }
    }
    
    void dip_mon_0_n0_v(const Vector& target, const Panel& source, Vector& dip_infl_v, Vector& mon_infl_v)
    {
        dip_mon_v(target, source.center, source.normal, dip_infl_v, mon_infl_v);
        dip_infl_v *= source.area;
        mon_infl_v *= source.area;
    }

    std::tuple<Vector, Vector> dip_mon_v(const Vector& target, const Vector& source, Vector direction)
    {
        Vector dip_infl_v;
        Vector mon_infl_v;
        dip_mon_v(target, source, direction, dip_infl_v, mon_infl_v);
        return std::make_tuple(dip_infl_v, mon_infl_v);
    }

    std::tuple<Vector, Vector> dip_mon_0_n0_v(const Vector& target, const Panel& source)
    {
        Vector dip_infl_v;
        Vector mon_infl_v;
        dip_mon_0_n0_v(target, source, dip_infl_v, mon_infl_v);
        return std::make_tuple(dip_infl_v, mon_infl_v);
    }
    
    std::tuple<Vector, Vector> dip_mon_0_vsaero_v(const Vector& target, const Panel& source)
    {
        Vector dip_infl_v;
        Vector mon_infl_v;
        dip_mon_0_vsaero_v(target, source, dip_infl_v, mon_infl_v);
        return std::make_tuple(dip_infl_v, mon_infl_v);
    }

}