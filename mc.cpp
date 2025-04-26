#include "mc.hpp"


MarkovChain& MarkovChain::operator=(const MarkovChain& other){
    delete _state;
    _state = other._state->clone();
    _gen = other._gen;
    _uniform_dist = other._uniform_dist;
    return *this;
}

void MarkovChain::update(propagator method, const size_t& steps){
    for (size_t i=0; i<steps; i++){
        (this->*method)();
    }
}

double MarkovChain::draw_uniform(const double& a, const double& b){
    return this->_uniform_dist(_gen)*(b-a) + a;
}

MonteCarlo& MonteCarlo::operator=(const MonteCarlo& other){
    if (&other != this){
        _clear_states();
        delete this->_mc;
        this->_mc = other._mc->clone();
        this->_data = copy_states(other._data);
    }
    return *this;
}


std::vector<const State*> MonteCarlo::data() const{ //array of states used for our statistics
    std::vector<const State*> res(_data.size());
    for (size_t i=0; i<_data.size(); i++){
        res[i] = _data[i];
    }
    return res;
}

Sample MonteCarlo::sample(const Observable A) const{ //each state generates a sample element, so all monte carlo states generate an entire sample to perform statistics.
    std::vector<double> sample_array(this->N());
    for (size_t i=0; i<this->N(); i++){
        sample_array[i] = A(*this->_data[i]);
    }
    return sample_array;
}

void MonteCarlo::update(const std::string& method, const size_t& steps, const size_t& sweeps){
    propagator pr = this->_mc->method(method);
    for (size_t i=0; i<steps; i++){
        this->_mc->update(pr, sweeps+1);
        this->_data.push_back(this->_mc->state().clone());
    }
}

void MonteCarlo::_clear_states(){
    for (size_t i = 0; i < _data.size(); i++){
        delete _data[i];
    }
    _data.clear();
}

MonteCarlo::~MonteCarlo(){
    _clear_states();
    delete this->_mc;
}

void update_all(const std::vector<MonteCarlo*>& obj, const std::string& method, const size_t& steps, const size_t& sweeps, int threads){

    threads = (threads <= 0) ? omp_get_max_threads() : threads;
    #pragma omp parallel for num_threads(threads)
    for (size_t i=0; i<obj.size(); i++){
        obj[i]->update(method, steps, sweeps);
    }
}