#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include <omp.h>

std::vector<double> pow(const std::vector<double>& x, const double& p);

double mean_value(const std::vector<double>& x);

double sample_std(const std::vector<double>& x);

std::vector<double> bin_it(const std::vector<double>& data);//assumes data is divisible by 2. Returns a vector with exactly half elements, each one is the mean of 2 consecutive elements of the original vector.

template<class Scalar>
double sum(const std::vector<Scalar>& x){
    Scalar res = 0;
    for (size_t i=0; i<x.size(); i++){
        res += x[i];
    }
    return res;
}

struct Sample;

struct BinningAnalysis;

struct Sample{

    /*
    Most important values are mean() and error():
    error() is basically the error bar of our mean value.

    For large N, the mean value follows a gaussian distribution around the true population mean,
    with a 1-sigma uncertainty of .error()

    */
    std::vector<double> sample;

    Sample(const std::vector<double>& sample = {}):sample(sample){}

    inline size_t N() const{ return sample.size();} //sample size

    inline double mean() const{
        return mean_value(sample);
    } //sample mean value

    inline double std() const{
        return sample_std(sample);
    } //standard deviation <x^2> - <x>^2

    inline double popul_std() const{
        return this->std()*sqrt(this->N()/(this->N()-1.));
    } //estimate of the population std, using Bessel's correction

    inline double error() const{
        return this->std()/sqrt(this->N()-1.);
    } //estimate of the gaussian 1sigma of the sample mean's distribution
    //From the CLT, this is popul_std/sqrt(N), and we can only estimate popul_std from our sample, using the formula above.

    std::string message() const;

    BinningAnalysis bin_it() const;


};


struct BinningAnalysis{

    BinningAnalysis(std::vector<double> sample){_init(sample);}

    BinningAnalysis(const Sample& sample){_init(sample.sample);}

    size_t max_level() const{
        return _samples.size();
    }

    const bool& converged() const{
        return _converged;
    }

    double tau_estimate() const{
        return 0.5 * (pow(_samples.back().error()/_samples.front().error(), 2) - 1);
    }

    const std::vector<Sample>& samples() const{
        return _samples;
    }

protected:
    std::vector<Sample> _samples;
    bool _converged;

    void _init(std::vector<double> sample);

};


struct State{

    virtual ~State() = default;

    virtual State* clone() const = 0;

    virtual std::unique_ptr<State> safe_clone() const = 0;
};


std::vector<State*> copy_states(const std::vector<State*>& states);


#endif