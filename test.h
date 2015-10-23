#ifndef TEST_H
#define TEST_h

#include "vtkWriter.h"
#include "laplace3DKernel.h"
#include "laplace2DKernel.h"

class laplaceKern2DTest{
public:
    laplaceKern2DTest(){};
    ~laplaceKern2DTest(){};
    static laplaceKern::PanelVector2 a;
    static laplaceKern::PanelVector2 b;
    static laplaceKern::Panel2 pan;

    static double test_doublet_2(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::doublet_2(v, a, b - a);
    }
    static double test_doublet_2_0(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::doublet_2_0(v, pan);
    }
    static double test_doublet_2_1_left(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::doublet_2_1(v, pan, true);
    }
    static double test_doublet_2_1_right(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::doublet_2_1(v, pan, false);
    }
    static double test_source_2(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::source_2(v, a);
    }
    static double test_source_2_0(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::source_2_0(v, pan);
    }    
    static double test_source_2_1_left(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::source_2_1(v, pan, true);
    }
    static double test_source_2_1_right(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::source_2_1(v, pan, false);
    }
    static laplaceKern::Vector3 test_doublet_2_v(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::toVector3(laplaceKern::doublet_2_v(v, a, b - a));
    }
    static laplaceKern::Vector3 test_doublet_2_0_v(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::toVector3(laplaceKern::doublet_2_0_v(v, pan));
    }
    static laplaceKern::Vector3 test_doublet_2_1_left_v(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::toVector3(laplaceKern::doublet_2_1_v(v, pan, true));
    }
    static laplaceKern::Vector3 test_doublet_2_1_right_v(double x, double y, double z)
    {
        laplaceKern::Vector2 v (x, y);
        return laplaceKern::toVector3(laplaceKern::doublet_2_1_v(v, pan, false));
    }
};

#endif