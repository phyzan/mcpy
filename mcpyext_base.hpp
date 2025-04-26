#ifndef MCPYEXT_BASE_HPP
#define MCPYEXT_BASE_HPP

#include "ising.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using _Shape = std::vector<size_t>;


struct PySample;

struct PyBinnedSample;

//----------------functions---------------------

template<class T>
py::array_t<T> np_array(const std::vector<T>& x, const _Shape& shape={}){
    if (shape.size() == 0){
        py::array_t<T> res(x.size(), x.data());
        return res;
    }
    else{
        py::array_t<T> res(shape, x.data());
        return res;
    }
}

std::vector<double> to_vector(const py::iterable& iterable);

py::list to_pystates(const std::vector<State>& states);

py::list to_pysamples(const std::vector<Sample>& states);

Observable to_observable(py::object f);

void define_base_module(py::module& m);


struct PySample : public Sample{

    PySample(const py::iterable& array):Sample(to_vector(array)){}

    PySample(const Sample& sample) : Sample(sample){}

    PyBinnedSample pybin_it() const;
};

struct PyBinnedSample : public BinningAnalysis{

    PyBinnedSample(const py::iterable& array):BinningAnalysis(to_vector(array)){}

    PyBinnedSample(const Sample& sample) : BinningAnalysis(sample){}

    PyBinnedSample(const BinningAnalysis& obj) : BinningAnalysis(obj){}
};

void py_update_all(py::iterable obj, py::str method, const size_t& steps, const size_t& sweeps, const int& threads);

#endif