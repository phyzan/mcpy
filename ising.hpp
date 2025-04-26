#ifndef ISING_HPP
#define ISING_HPP

#include "mc.hpp"

struct SpinState;

class IsingModel2DMarkovChain;

class IsingModel2D;

std::vector<int> random_spins(const size_t& Lx, const size_t& Ly);


struct SpinState : public State{

    /*
    Periodic 2D spin lattice
    */
    
    std::vector<int> spins;
    std::array<size_t, 2> shape;

    SpinState(const std::vector<int>& spins, const size_t& Lx, const size_t& Ly):spins(spins), shape({Lx, Ly}){
        if (Lx*Ly != spins.size()){
            throw std::runtime_error("");
        }
    }

    State* clone() const override{
        return new SpinState(*this);
    }

    std::unique_ptr<State> safe_clone() const override{
        return std::make_unique<SpinState>(*this);
    }

    const int& operator()(const long int& i) const;

    int& operator()(const long int& i);

    const int& operator()(long int i, long int j) const;

    int& operator()(long int i, long int j);

    size_t sites() const;

    double M() const;

    double energy() const;

    std::vector<size_t> neighbors(const size_t& site);

    size_t index(long int i, long int j) const;

};


class IsingModel2DMarkovChain : public MarkovChain{


public:

    const double& Temp() const{
        return _T;
    }

    IsingModel2DMarkovChain(const double& T, const size_t& Lx, const size_t& Ly) : MarkovChain(SpinState(random_spins(Lx, Ly), Lx, Ly)), _T(T), _spin_roulette(0, Lx*Ly-1){}

    inline MarkovChain* clone() const override{ return new IsingModel2DMarkovChain(*this);}

    inline std::unique_ptr<MarkovChain> safe_clone() const override{ return std::make_unique<IsingModel2DMarkovChain>(*this);}

    propagator method(const std::string& name) const override;

    void ssf_update();

    void wolff_update();

    const SpinState& ising_state() const{
        return static_cast<const SpinState&>(this->state());
    }

private:

    double _T;
    mutable std::uniform_int_distribution<size_t> _spin_roulette;

    inline size_t _choose_site() const{return _spin_roulette(this->_gen);}

};



class IsingModel2D : public MonteCarlo{


public:

    IsingModel2D(const double& T, const size_t& Lx, const size_t& Ly): MonteCarlo(IsingModel2DMarkovChain(T, Lx, Ly)){}

    void ssf_update(const size_t& steps, const size_t& sweeps = 0){
        this->update("ssf", steps, sweeps);
    }

    void wolff_update(const size_t& steps, const size_t& sweeps = 0){
        this->update("wolff", steps, sweeps);
    }

    void ssf_thermalize(const size_t& sweeps){
        this->thermalize("ssf", sweeps);
    }

    void wolff_thermalize(const size_t& sweeps){
        this->thermalize("wolff", sweeps);
    }

   inline Sample energy_sample() const{
        return this->sample([](const State& s){return static_cast<const SpinState&>(s).energy();});
   }

   inline double T() const{
        return this->chain().Temp();
   }

   inline const IsingModel2DMarkovChain& chain() const {
    return static_cast<IsingModel2DMarkovChain&>(*this->_mc);
}

private:

    inline IsingModel2DMarkovChain& _chain(){
        return static_cast<IsingModel2DMarkovChain&>(*this->_mc);
    }

};



#endif