#ifndef MC_HPP
#define MC_HPP


#include "tools.hpp"
#include <functional>
#include <random>


using Observable = std::function<double(const State&)>;

class MarkovChain;

class MonteCarlo;

using propagator = void (MarkovChain::*)();

class MarkovChain{

    //Base abstract class representing any Markov chain.
    //This class simply evolves the current state, and does not hold any data.
    //A Monte-Carlo class that implements a Markov chain is responsible for doing statistics.

public:


    inline const State& state()const{ return *_state;}

    virtual ~MarkovChain(){delete _state;}
    
    virtual propagator method(const std::string& name) const = 0;

    virtual MarkovChain* clone() const = 0;

    virtual std::unique_ptr<MarkovChain> safe_clone() const = 0;

    inline void update(const std::string& method, const size_t& steps=1){ this->update(this->method(method), steps);}

    void update(propagator method, const size_t& steps);

protected:

    MarkovChain(const State& initial_state):_state(initial_state.clone()), _gen(std::random_device()()), _uniform_dist(0, 1){}

    MarkovChain(const MarkovChain& other):_state(other._state->clone()), _gen(other._gen), _uniform_dist(other._uniform_dist){}

    MarkovChain(MarkovChain&& other):_state(std::move(other._state)), _gen(std::move(other._gen)), _uniform_dist(std::move(other._uniform_dist)){}

    MarkovChain& operator=(const MarkovChain& other);

    double draw_uniform(const double& a, const double& b);

    State* _state;
    mutable std::mt19937 _gen;

private:
    
    mutable std::uniform_real_distribution<> _uniform_dist;
};


//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


class MonteCarlo{


public:

    MonteCarlo(const MarkovChain& mc):_mc(mc.clone()){}

    MonteCarlo(const MonteCarlo& other): _mc(other._mc->clone()), _data(copy_states(other._data)){}

    MonteCarlo(MonteCarlo&& other): _mc(other._mc), _data(std::move(other._data)){}

    virtual ~MonteCarlo();

    MonteCarlo& operator=(const MonteCarlo& other);

    std::vector<const State*> data() const;

    inline const State& state(const size_t& i) const{return *_data.at(i);}

    inline size_t N() const{return _data.size();}

    Sample sample(const Observable A) const;

    void update(const std::string& method, const size_t& steps, const size_t& sweeps=0);

    inline void thermalize(const std::string& method, const size_t& sweeps){this->_mc->update(method, sweeps);}

    inline const MarkovChain& chain() const{return *this->_mc;}


protected:

    void _clear_states();

    MarkovChain* _mc = nullptr; //pointer to dynamically allocated markov chain that propagates the simulation. This is passed from derived classes to the constructor of this class
    std::vector<State*> _data = {}; //vector that holds all states obtained from the Markov chain to use for our statistics
};


void update_all(const std::vector<MonteCarlo*>& obj, const std::string& method, const size_t& steps, const size_t& sweeps, int threads);



#endif