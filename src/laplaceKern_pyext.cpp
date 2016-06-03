#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Eigen/Core"

#include <vector>
#include <tuple>
#include <memory>

#include "laplaceKern2D.h"
#include "laplaceKern3D.h"

namespace py = pybind11;
namespace l2 = laplaceKern2D;
namespace l3 = laplaceKern3D;


// (l3::Vector (*)(const l3::Vector&, const l3::Vector&, const l3::Vector&)) &l3::vortex_v ->
// this is for overloaded functions (returnType (*)(arg1type, arg2type...) &function-name)

// PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void init_laplace_kernels(py::module &m)
{
    py::module m2 = m.def_submodule("D2", "2d elements");
    py::module m3 = m.def_submodule("D3", "3d elements");

    /************************-2D SINGULARITY ELEMENTS-************************/
    py::class_<l2::Panel>(m2, "Panel", "a 2d panel")
        .def(py::init<l2::VectorPtr, l2::VectorPtr>())
        .def_readonly("area", &l2::Panel::area)
        .def_readonly("center", &l2::Panel::center)
        .def_readonly("normal", &l2::Panel::normal)
        .def_readonly("tangent", &l2::Panel::tangent)
        .def_readonly("points", &l2::Panel::points);

    m2.def("doublet", &l2::doublet);
    m2.def("source", &l2::source);
    m2.def("vortex", &l2::vortex);
    m2.def("source_0", &l2::source_0);
    m2.def("doublet_0", &l2::doublet_0);
    m2.def("doublet_v", &l2::doublet);
    m2.def("source_v", &l2::source);
    m2.def("vortex_v", &l2::vortex);
    m2.def("source_0_v", &l2::source_0);
    m2.def("doublet_0_v", &l2::doublet_0);

    /************************-3D SINGULARITY ELEMENTS-************************/
    py::class_<l3::Panel>(m3, "Panel", "a 3d panel")
        .def(py::init<std::vector<l3::VectorPtr>>())
        .def_readonly("area", &l3::Panel::area)
        .def_readonly("normal", &l3::Panel::normal)
        .def_readonly("tangent_m", &l3::Panel::tangent_m)
        .def_readonly("tanget_l", &l3::Panel::tangent_l)
        .def_readonly("center", &l3::Panel::center)
        .def_readonly("side_sum", &l3::Panel::area);

    py::class_<l3::Edge>(m3, "Edge", "a 3d edge")
        .def(py::init<l3::VectorPtr, l3::VectorPtr>())
        .def_readonly("v1", &l3::Edge::v1)
        .def_readonly("v2", &l3::Edge::v2)
        .def_readonly("tangent", &l3::Edge::tangent);

    m3.def("monopole", &l3::monopole);
    m3.def("monopole_0_n0", &l3::monopole_0_n0);
    m3.def("dipole", &l3::dipole);
    m3.def("dipole_0_n0", &l3::dipole_0_n0);
    m3.def("dipole_0_sphere", &l3::dipole_0_sphere);
    m3.def("dipole_0_vsaero", &l3::dipole_0_vsaero);
    m3.def("dip_mon", (std::tuple<double, double> (*)
        (const l3::Vector&, const l3::Vector&, l3::Vector)) &l3::dip_mon);
    m3.def("dip_mon_0_n0", (std::tuple<double, double> (*)
        (const l3::Vector&, const l3::Panel&)) &l3::dip_mon_0_n0);
    m3.def("dip_mon_0_vsaero", (std::tuple<double, double> (*)
        (const l3::Vector&, const l3::Panel&)) &l3::dip_mon_0_vsaero);

    m3.def("monopole_v", &l3::monopole_v);
    m3.def("monopole_0_n0_v", &l3::monopole_0_n0_v);
    m3.def("monopole_0_vsaero_v", &l3::monopole_0_vsaero_v);
    m3.def("dipole_v", &l3::dipole_v);
    m3.def("dipole_0_n0_v", &l3::dipole_0_n0_v);
    m3.def("dipole_0_vsaero_v", &l3::dipole_0_vsaero_v);
    m3.def("vortex_v", (l3::Vector (*)
        (const l3::Vector&, const l3::Vector&, const l3::Vector&)) &l3::vortex_v);
    m3.def("vortex_v", (l3::Vector (*)
        (const l3::Vector&, const l3::Edge&))&l3::vortex_v);
    m3.def("vortex_half_infinity_v", &l3::vortex_half_infinity_v);
    m3.def("dip_mon_v", (std::tuple<l3::Vector, l3::Vector> (*)
        (const l3::Vector&, const l3::Vector&, l3::Vector)) &l3::dip_mon_v);
    m3.def("dip_mon_0_n0_v", (std::tuple<l3::Vector, l3::Vector> (*)
        (const l3::Vector&, const l3::Panel&)) &l3::dip_mon_0_n0_v);
    m3.def("dip_mon_0_vsaero_v", (std::tuple<l3::Vector, l3::Vector> (*)
        (const l3::Vector&, const l3::Panel&)) &l3::dip_mon_0_vsaero_v);
};


PYBIND11_PLUGIN(laplaceKern)
{
    py::module m("laplaceKern", "laplace_kernels");
    init_laplace_kernels(m);
    return m.ptr();
}

