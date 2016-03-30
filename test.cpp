#include "vtkWriter.h"
#include "test.h"

laplaceKern::PanelVector2 laplaceKern2DTest::a = laplaceKern::PanelVector2(1, 0);
laplaceKern::PanelVector2 laplaceKern2DTest::b = aplaceKern::PanelVector2(-1, 0);
laplaceKern::Panel2 laplaceKern2DTest::pan = laplaceKern::Panel2(&laplaceKern2DTest::a, &laplaceKern2DTest::b);

int main ( int, char *[] )
{
    VtkWriter writer("pot.vtk");
    writer.writePoints(-10, 10, 100,
                       -10, 10, 100,
                       -0, 0, 0);
    writer.writeScalar(laplaceKern2DTest::test_doublet_2, "doublet_2");
    writer.writeScalar(laplaceKern2DTest::test_doublet_2_0, "doublet_2_0");
    writer.writeScalar(laplaceKern2DTest::test_doublet_2_1_left, "doublet_2_1_l");
    writer.writeScalar(laplaceKern2DTest::test_doublet_2_1_right, "doublet_2_1_r");
    writer.writeScalar(laplaceKern2DTest::test_source_2, "source_2");
    writer.writeScalar(laplaceKern2DTest::test_source_2_0, "source_2_0");
    writer.writeScalar(laplaceKern2DTest::test_source_2_1_left, "source_2_1_l");
    writer.writeScalar(laplaceKern2DTest::test_source_2_1_right, "source_2_1_r");
    writer.writeVector3(laplaceKern2DTest::test_doublet_2_v, "doublet_2_v");
    writer.writeVector3(laplaceKern2DTest::test_doublet_2_0_v, "doublet_2_0_v");
    writer.writeVector3(laplaceKern2DTest::test_doublet_2_1_left_v, "doublet_2_1_l_v");
    writer.writeVector3(laplaceKern2DTest::test_doublet_2_1_right_v, "doublet_2_1_r_v");
    writer.createFile();
}