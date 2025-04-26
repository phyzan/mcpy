#include "ising.hpp"


std::vector<int> random_spins(const size_t& Lx, const size_t& Ly){
    std::vector<int> spins(Lx*Ly);
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, 1); // generates 0 or 1

    for (int& s : spins) {
        s = dist(gen) ? 1 : -1;
    }
    return spins;
}

const int& SpinState::operator()(long int i, long int j) const{
    return this->operator()(index(i, j));
}

int& SpinState::operator()(long int i, long int j) {
    return this->operator()(index(i, j));
}

int& SpinState::operator()(const long int& i){
    return spins.at((i+spins.size()) % spins.size());
}

const int& SpinState::operator()(const long int& i) const {
    return spins.at((i+spins.size()) % spins.size());
}


size_t SpinState::sites() const{
    return this->spins.size();
}

double SpinState::M()const{
    return ::sum(this->spins);
}

double SpinState::energy() const{
    long int res = 0;
    const SpinState& s = *this;

    for (size_t i=0; i<shape[0]; i++){
        for (size_t j=0; j<shape[1]; j++){
            res -= s(i, j)*(s(i-1, j) + s(i, j-1));
        }
    }
    return res;
}

std::vector<size_t> SpinState::neighbors(const size_t& site){
    long int i = site % shape[0];
    long int j = site / shape[0];
    return {index(i-1, j), index(i+1, j), index(i, j-1), index(i, j+1)};
}

size_t SpinState::index(long int i, long int j) const{
    i = (i + shape[0]) % shape[0];
    j = (j + shape[1]) % shape[1];
    return j * shape[0] + i;
}


void IsingModel2DMarkovChain::ssf_update(){
    SpinState& S = static_cast<SpinState&>(*this->_state);
    size_t k = this->_choose_site();
    size_t i = k % S.shape[0];
    size_t j = k/S.shape[0];
    double de = 2*S(k)*(S(i-1, j)+S(i+1, j)+S(i, j-1)+S(i, j+1));
    if (this->draw_uniform(0, 1) <= exp(-de/_T)){
        S.spins[k] *= -1;
    }
}

void IsingModel2DMarkovChain::wolff_update(){
    size_t site = _choose_site();
    SpinState& S = static_cast<SpinState&>(*this->_state);
    const int s = S(site);
    const double p = 1-std::exp(-2/_T);
    std::vector<size_t> remaining = {site}; //container with sites whose neighbors we need to check

    S(site) = -s;
    size_t cluster_size = 1;
    while (remaining.size()>0){
        site = remaining.back(); remaining.pop_back();
        for (const size_t& nr : S.neighbors(site)){
            if ( (S(nr) == s ) && (this->draw_uniform(0, 1) < p)){
                S(nr) = -s;
                remaining.push_back(nr);
                cluster_size++;
            }
        }
    }
}

propagator IsingModel2DMarkovChain::method(const std::string& name) const {
    if (name == "ssf"){
        return static_cast<propagator>(&IsingModel2DMarkovChain::ssf_update);
    }
    else if (name == "wolff"){
        return static_cast<propagator>(&IsingModel2DMarkovChain::wolff_update);
    }
    else{
        throw std::runtime_error("");
    }
}
