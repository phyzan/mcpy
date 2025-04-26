#include "mcpyext_base.hpp"


std::vector<double> to_vector(const py::iterable &iterable){
    std::vector<double> res;
    for (const py::handle &item : iterable)
    {
        res.push_back(py::cast<double>(item));
    }
    res.shrink_to_fit();
    return res;
}

py::list to_pystates(const std::vector<const State*>& states){
    py::list res(states.size());
    for (size_t i=0; i<states.size(); i++){
        res[i] = states[i]->safe_clone();
    }
    return res;
}

Observable to_observable(py::object f){
    return [f](const State& state){
        return f(&state).cast<double>();
    };
}

py::list to_pysamples(const std::vector<Sample>& samples){
    py::list res(samples.size());
    for (size_t i=0; i<samples.size(); i++){
        res[i] = PySample(samples[i]);
    }
    return res;
}

inline PyBinnedSample PySample::pybin_it() const{
    return this->bin_it();
}

void py_update_all(py::iterable obj, py::str method, const size_t& steps, const size_t& sweeps, const int& threads){
    std::vector<MonteCarlo*> array;
    for (const py::handle& item : obj){
        array.push_back(&item.cast<MonteCarlo&>());
    }
    update_all(array, method.cast<std::string>(), steps, sweeps, threads);
}


void define_base_module(py::module &m)
{

    py::class_<PySample>(m, "Sample", py::module_local())
        .def(py::init<py::iterable>(), py::arg("data"))
        .def_property_readonly("data", [](const PySample &self){ return np_array<double>(self.sample); })
        .def_property_readonly("N", &Sample::N)
        .def("mean", &Sample::mean)
        .def("std", &Sample::std)
        .def("popul_std", &Sample::popul_std)
        .def("error", &Sample::error)
        .def("bin_it", &PySample::pybin_it)
        .def_property_readonly("stat", &Sample::message);

    py::class_<PyBinnedSample>(m, "BinnedSample", py::module_local())
        .def(py::init<py::iterable>(), py::arg("data"))
        .def_property_readonly("max_level", &PyBinnedSample::max_level)
        .def_property_readonly("converged", &PyBinnedSample::converged)
        .def_property_readonly("tau_estimate", &PyBinnedSample::tau_estimate)
        .def_property_readonly("binned_samples", [](const PyBinnedSample& self){return to_pysamples(self.samples());});


    py::class_<State, std::unique_ptr<State>>(m, "State", py::module_local());

    py::class_<SpinState, State>(m, "SpinState", py::module_local())
        .def_property_readonly("spins", [](const SpinState& self){return np_array<int>(self.spins, {self.shape[0], self.shape[1]});})
        .def("__call__", [](const SpinState& self, long int i) {
            return self(i);
        })
        .def("__call__", [](const SpinState& self, long int i, long int j) {
            return self(i, j);
        })
        .def_property_readonly("sites", &SpinState::sites)
        .def_property_readonly("M", &SpinState::M)
        .def_property_readonly("energy", &SpinState::energy);

    py::class_<MarkovChain, std::unique_ptr<MarkovChain>>(m, "MarkovChain", py::module_local())
        .def_property_readonly("state", [](const MarkovChain& self) {return self.state().safe_clone();})
        .def("update", [](MarkovChain& self, py::str method, const size_t& steps) {return self.update(method.cast<std::string>(), steps);}, py::arg("method"), py::arg("steps")=1);


    py::class_<IsingModel2DMarkovChain, MarkovChain>(m, "IsingModel2DMarkovChain", py::module_local())
        .def(py::init<double, size_t, size_t>(), py::arg("T"), py::arg("Lx"), py::arg("Ly"))
        .def("ssf_update", &IsingModel2DMarkovChain::ssf_update)
        .def("wolff_update", &IsingModel2DMarkovChain::wolff_update);

    py::class_<MonteCarlo, std::unique_ptr<MonteCarlo>>(m, "MonteCarlo", py::module_local())
        .def(py::init<MarkovChain&>(), py::arg("markov_chain"))
        .def_property_readonly("data", [](MonteCarlo& self){return to_pystates(self.data());})
        .def_property_readonly("N", &MonteCarlo::N)
        .def("sample", [](const MonteCarlo& self, py::object obs){return PySample(self.sample(to_observable(obs)));}, py::arg("observable"))
        .def("update", &MonteCarlo::update, py::arg("method"), py::arg("steps"), py::arg("sweeps")=0)
        .def("thermalize", &MonteCarlo::thermalize, py::arg("method"), py::arg("sweeps"));
    
    py::class_<IsingModel2D, MonteCarlo>(m, "IsingModel2D", py::module_local())
        .def(py::init<double, size_t, size_t>(), py::arg("T"), py::arg("Lx"), py::arg("Ly"))
        .def_property_readonly("Temp", &IsingModel2D::T)
        .def("ssf_update", &IsingModel2D::ssf_update, py::arg("steps"), py::arg("sweeps")=0)
        .def("wolff_update", &IsingModel2D::wolff_update, py::arg("steps"), py::arg("sweeps")=0)
        .def("ssf_thermalize", &IsingModel2D::ssf_thermalize, py::arg("sweeps"))
        .def("wolff_thermalize", &IsingModel2D::wolff_thermalize, py::arg("sweeps"))
        .def("energy_sample", [](const IsingModel2D& self){return PySample(self.energy_sample());});

    m.def("update_all", py_update_all, py::arg("sims"), py::arg("method"), py::arg("steps"), py::arg("sweeps")=0, py::arg("threads")=-1);
}

