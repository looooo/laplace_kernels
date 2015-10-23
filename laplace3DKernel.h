#ifndef laplace3DKernel_H
#define laplace3DKernel_H

#include <Eigen/Core>
#include <Eigen/Geometry>
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

double coresize = 0.0000000000000000000001;

namespace laplaceKern{

    //forward decleration
    class Panel3;
    class Edge;


  
    /**
    * A 3D Vector from the Eigen library -> http://eigen.tuxfamily.org/index.php?title=Main_Page.
    */

    typedef Eigen::Vector3d Vector3;
    
    
    
    /** 
    * A Panel Vector3 is a Vector3 with some other members like a scalar potential, an vector 
    * velocity and a list of Panels which contain all the panels which use this vertex.
    */

    class PanelVector3: public Vector3{
    public:
        PanelVector3():Vector3(){}
        PanelVector3(double x, double y, double z):Vector3(x, y, z){}
        template<typename OtherDerived> PanelVector3(const Eigen::MatrixBase<OtherDerived>& other):Vector3(other){};
        template<typename OtherDerived> PanelVector3& operator=(const Eigen::MatrixBase <OtherDerived>& other){
            this->Vector3::operator=(other);
            return *this;
        }
        
        int owner = 0;
        int nr = -1;
        double mue;
        double sigma;
        bool wake_vertex = false;
        std::vector<Panel3*> panels;                //all Panels containing this vertex
        void reset_properties(){this->mue=0; this->sigma=0;}
    };
    
    class Panel3{
    private:
        int nr;
        Vector3 velocity = Vector3(0,0,0);
        double mue = 0.;
        double sigma = 0.;
    public:
        
        vector<PanelVector3*> points;
        Panel3(){};
        virtual vector<PanelVector3*> get_points()
        {
            return this->points;
        };  
        vector<Edge> get_edges();
        Panel3* get_edge_neighbour(Edge e);

        void append_point(PanelVector3* v)        
        {
            this->points.push_back(v);
            v->owner += 1;
            v->panels.push_back(this);
        };
        virtual void calc_geo()
        {
            int l;
            this->area = 0.;
            this->center = Vector3(0.,0.,0.);
            this->n = Vector3(0.,0.,0.);
            
            l = this->points.size();
            
            for (PanelVector3* point: this->points){
                this->center += *point;
            }
            this->center /= double(l);
            
            for (int i=0; i < l; i++){
                int j = (i == l - 1 ? 0 : i+1);
                this->side_sum += (*this->points[j] - *this->points[i]).norm();
                Vector3 n_i = (this->center - *this->points[i]).cross(this->center - *this->points[j]);
                this->n -= n_i;
                this->area += n_i.norm();
            }
            this->area /= 2;
            this->n.normalize();
            this->l = (this->center - (*this->points[1] + *this->points[2]) / 2);
            this->l.normalize();
            this->m = this->n.cross(this->l);
        };
        double area = 0.;
        double side_sum = 0.;

        Vector3 center, n, m, l;
        Vector3 loc_coord(Vector3);
        
        //GETTER
        virtual double get_mue(){return this->mue;};
        virtual double get_sigma(){return this->sigma;}
        virtual int get_nr(){return this->nr;}
        
        //SETTER
        void set_mue(double mue){this->mue = mue;}
        void set_sigma(double sigma){this->sigma = sigma;}
        void set_nr(int nr){this->nr = nr;}

        void reset_properties(){this->mue=0; this->sigma=0;}

        string tostring(){
            std::ostringstream out;
            out << "Panel3(";
            for (int i = 0; i < this->points.size(); i++)
            {
                out << this->points[i];
                out << (i == this->points.size() - 1 ? ")" : ", ");
            }
            return out.str();
        };
        bool is_inside(Vector3& v){
            //     check if the point v lies inside the plane
            double coresize = 0.000000000000000000000001;
            if (v == this->center)
            {
                return 1;
            }
            int l = this->points.size();
            for (int i = 0; i < l; i++)
            {
                int j = (i == l-1 ? 0 : i+1);
                if ((*this->points[j]-v).cross(v-*this->points[i]).dot(this->n) <= coresize)
                {
                    return 0;
                }
            }
            return 1;
        };
    };
    
    class Edge{
    private:
        void setup_Edge();
    public:
        Edge(){};
        Edge(PanelVector3*, PanelVector3*);
        PanelVector3* v1;
        PanelVector3* v2;
        Panel3* p1;                                         // make private + getter 
        Panel3* p2;                                         // make private + getter
        double vorticity;
    };

    /************************-3D ELEMENT INFLUENCE-*********************/

    void doublet_3_0_n0(Vector3& target, Panel3* source, double& dip_infl)
        {
            Vector3 p_v = target - source->center;
            double pn = p_v.dot(source->n);
            dip_infl = -pn * source->area / pow(p_v.norm(), 3) / 4 / M_PI;
        }

        void doublet_3_0_vsaero_v(Vector3& target, Panel3* source, Vector3& vel_infl)
        {
            Vector3 p_v = target - source->center;
            int l = source->points.size();

            for (int i = 0; i < l; i++){
                int j = (i == l-1? 0 : i+1);
                Vector3 &p1 = *source->points[i];
                Vector3 &p2 = *source->points[j];
                Vector3 a = target - p1;
                Vector3 b = target - p2;
                double absa = a.norm();
                double absb = b.norm();
                if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
                {
                    vel_infl -= (a.cross(b) * (absa + absb) / 4 / M_PI / (absa * absb * (absa * absb + a.dot(b))));
                }
            }
        }

        Vector3 doublet_3_0_vsaero_v(Vector3& target, Panel3* source)
        {
        Vector3 vel = Vector3(0, 0, 0);
        Vector3 p_v = target - source->center;
        int l = source->points.size();

        for (int i = 0; i < l; i++){
            int j = (i == l-1? 0 : i+1);
            Vector3 &p1 = *source->points[i];
            Vector3 &p2 = *source->points[j];
            Vector3 a = target - p1;
            Vector3 b = target - p2;
            double absa = a.norm();
            double absb = b.norm();
            if ((absa * absb * fabs(absa * absb + a.dot(b))) > 0.000001 )
            {
            vel -= (a.cross(b) * (absa + absb) / 4 / M_PI / (absa * absb * (absa * absb + a.dot(b))));
            }
        }
        return vel;
        }


        void doublet_3_0_vsaero(Vector3 & target, Panel3* source, double& dip_infl)
        {
        Vector3 p_v = target - source->center;
        double pn = p_v.dot(source->n);
        int l = source->points.size();
        double rnum;
        double dnom;


        for (int i = 0; i < l; i++){
            int j = (i == l-1? 0 : i+1);
            Vector3 &p1 = *source->points[i];
            Vector3 &p2 = *source->points[j];
            Vector3 s = p2 - p1;
            Vector3 a = target - p1;
            Vector3 b = target - p2;
            Vector3 h = a.cross(s);
            double bm = b.dot(source->m);
            double am = a.dot(source->m);
            double sl = s.dot(source->l);
            double sm = s.dot(source->m);
            double al = am * sl - a.dot(source->l) * sm;
            double pa = pow(pn, 2) * sl + al * am;
            double pb = pow(pn, 2) * sl + al * bm;
            double absa = a.norm();
            double absb = b.norm();
            double abss = s.norm();
            
            if (
            absa > coresize && 
            absb > coresize &&
            absa + absb - abss >= 0
            )
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
        dip_infl /= 4 * M_PI;
        }

        double doublet_3_0_sphere(Vector3& target, Panel3* source, double& dip_infl)
        {
            Vector3 u1, u2, u3;
            
            u1 = source->center - target;
            if (u1.norm() < 0.000000000000000000001){
                // point lies on the center of the panel
                // this will return the value of the integral, with the panel wrapped around the center point
                // in n direction.
                dip_infl  = -0.5;
                return dip_infl;
            }
            u1.normalize();
            for (int i = 0; i < source->points.size(); i++){
                int j = (i == source->points.size()-1? 0 : i+1);
                u2 = *source->points[i] - target;
                u3 = *source->points[j] - target;
                u2.normalize();
                u3.normalize();
                dip_infl -= atan2((u1.cross(u2)).dot(u3), (1 + u2.dot(u3) + u3.dot(u1) + u1.dot(u2))) * 0.5 / M_PI;
            }
            return dip_infl;
        }


        void doublet_src_3_0_n0(Vector3& target, Panel3* source, double& dip_infl, double& mop_infl)
        {
        Vector3 p_v = target - source->center;
        double absp = p_v.norm();
        double pn = p_v.dot(source->n);
        dip_infl = -pn * source->area / pow(absp, 3) / 4 / M_PI;
        mop_infl = source->area / absp / 4 / M_PI;
        }

        Vector3 doublet_src_3_0_vsaero_v(Vector3& target, Panel3* source, double mue, double sigma)
        {
            Vector3 vel = Vector3(0, 0, 0);
            double dip_infl = 0;
            Vector3 p_v = target - source->center;
            double pn = p_v.dot(source->n);
            int l = source->points.size();
            double rnum;
            double dnom;


            for (int i = 0; i < l; i++){
                int j = (i == l-1? 0 : i+1);
                Vector3 &p1 = *source->points[i];
                Vector3 &p2 = *source->points[j];
                Vector3 s = p2 - p1;
                Vector3 a = target - p1;
                Vector3 b = target - p2;
                Vector3 h = a.cross(s);
                double bm = b.dot(source->m);
                double am = a.dot(source->m);
                double sl = s.dot(source->l);
                double sm = s.dot(source->m);
                double al = am * sl - a.dot(source->l) * sm;
                double pa = pow(pn, 2) * sl + al * am;
                double pb = pow(pn, 2) * sl + al * bm;
                double absa = a.norm();
                double absb = b.norm();
                double abss = s.norm();
                
                if (
                absa > coresize && 
                absb > coresize &&
                absa + absb - abss >= 0
                ){
                dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
                rnum = sm * pn * (absb * pa - absa * pb);

                if (fabs(pn) < coresize){
                    double pn_sign = (pn >= 0 ? -1 : 1);
                    double al_sign = (al < 0 ? 1 : -1);
                    if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
                    else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
                }
                else{
                    dip_infl += atan2(rnum, dnom);
                }    if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
                vel -= mue * (a.cross(b) * (absa + absb) / (absa * absb * (absa * absb + a.dot(b))));
                vel -= sigma * log((absa + absb + abss) / (absa + absb - abss)) / abss * (sm * source->l - sl * source->m);
                }
            }
            vel -= sigma * source->n * dip_infl;
            return vel / 4 / M_PI;
        }

        Vector3 src_3_0_vsaero_v(Vector3& target, Panel3* source, double sigma)
        {
            Vector3 vel = Vector3(0, 0, 0);
            double dip_infl = 0;
            Vector3 p_v = target - source->center;
            double pn = p_v.dot(source->n);
            int l = source->points.size();
            double rnum;
            double dnom;


            for (int i = 0; i < l; i++){
                int j = (i == l-1? 0 : i+1);
                Vector3 &p1 = *source->points[i];
                Vector3 &p2 = *source->points[j];
                Vector3 s = p2 - p1;
                Vector3 a = target - p1;
                Vector3 b = target - p2;
                Vector3 h = a.cross(s);
                double bm = b.dot(source->m);
                double am = a.dot(source->m);
                double sl = s.dot(source->l);
                double sm = s.dot(source->m);
                double al = am * sl - a.dot(source->l) * sm;
                double pa = pow(pn, 2) * sl + al * am;
                double pb = pow(pn, 2) * sl + al * bm;
                double absa = a.norm();
                double absb = b.norm();
                double abss = s.norm();
                
                if (
                absa > coresize && 
                absb > coresize &&
                absa + absb - abss >= 0
                ){
                dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
                rnum = sm * pn * (absb * pa - absa * pb);

                if (fabs(pn) < coresize){
                    double pn_sign = (pn >= 0 ? 1 : -1);
                    double al_sign = (al < 0 ? 1 : -1);
                    if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
                    else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
                }
                else{
                    dip_infl += atan2(rnum, dnom);
                }    if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
                vel -= sigma * log((absa + absb + abss) / (absa + absb - abss)) / abss * (sm * source->l - sl * source->m);
                }
            }
            vel -= sigma * source->n * dip_infl;
            return vel / 4 / M_PI;
        }


        void doublet_src_3_0_vsaero(Vector3& target, Panel3* source, double& dip_infl, double& mop_infl)
        {
            Vector3 p_v = target - source->center;
            double pn = p_v.dot(source->n);
            int l = source->points.size();
            double rnum;
            double dnom;


            for (int i = 0; i < l; i++){
                int j = (i == l-1? 0 : i+1);
                Vector3 &p1 = *source->points[i];
                Vector3 &p2 = *source->points[j];
                Vector3 s = p2 - p1;
                Vector3 a = target - p1;
                Vector3 b = target - p2;
                Vector3 h = a.cross(s);
                double bm = b.dot(source->m);
                double am = a.dot(source->m);
                double sl = s.dot(source->l);
                double sm = s.dot(source->m);
                double al = am * sl - a.dot(source->l) * sm;
                double pa = pow(pn, 2) * sl + al * am;
                double pb = pow(pn, 2) * sl + al * bm;
                double absa = a.norm();
                double absb = b.norm();
                double abss = s.norm();
                
                if (
                absa > coresize && 
                absb > coresize &&
                absa + absb - abss >= 0
                )
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
                mop_infl -= al * log((absa + absb + abss) / (absa + absb - abss)) / abss;  
                }
            }
            mop_infl += pn * dip_infl;
            dip_infl /= 4. * M_PI;
            mop_infl /= 4. * M_PI;
        }


        Vector3 vortex_3_0_v(Vector3& target, Vector3& source_point_0, Vector3& source_point_1)
        {
        Vector3 r0 = source_point_1 - source_point_0;
        Vector3 r1 = target - source_point_0;
        Vector3 r2 = target - source_point_1;
        Vector3 r_cross = r1.cross(r2);
        if (r_cross.norm() > 0.001){
            return r_cross * ((r0.dot(r1) / r1.norm() - r0.dot(r2) / r2.norm()) / r_cross.dot(r_cross) / 4 / M_PI);
        }
        return Vector3(0.,0.,0.);
        }

        Vector3 vortex_3_0_v(Vector3& target, Edge& e)
        {
            return e.vorticity * vortex_3_0_v(target, *e.v1, *e.v2);
        }



        Vector3 vortex_3_0_half_infinity_v(Vector3& target, Vector3& source_direction, Vector3& source_point)
        {
        Vector3 r1 = target - source_point;
        Vector3 cross = r1.cross(source_direction);
        if (cross.norm() > 0.00001){
            return cross / cross.dot(cross) * (1 + r1.dot(source_direction) / r1.norm()) / 4 / M_PI;
        }
        return Vector3(0.,0.,0.);
        }

        double doublet_3(Vector3& target, Vector3& source, Vector3& normal)
        {
            Vector3 r = target - source;
            return -1. / 4 / M_PI / pow(r.dot(r), (3. / 2)) * r.dot(normal);
        }

        Vector3 doublet_3_v(Vector3& target, Vector3& source, Vector3& normal)
        {
            Vector3 r = target - source;
            return ((normal.cross(r)).cross(r) + 2. * r * (r.dot(normal))) / 4. / M_PI / pow(r.dot(r), 5. / 2.);
        }

    

}
#endif