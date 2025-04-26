#include "tools.hpp"

std::vector<double> pow(const std::vector<double>& x, const double& p){
    std::vector<double> res(x.size());
    for (size_t i=0; i<x.size(); i++){
        res[i] = pow(x[i], p);
    }
    return res;
}

double mean_value(const std::vector<double>& x){
    return sum(x)/x.size();
}

double sample_std(const std::vector<double>& x){
    return std::sqrt(mean_value(pow(x, 2)) - pow(mean_value(x), 2));
}


std::vector<double> bin_it(const std::vector<double>& data){
    if (data.size() % 2 != 0){
        throw std::runtime_error("Data size needs to be an even number");
    }

    std::vector<double> new_data(data.size()/2);
    for (size_t i=0; i<new_data.size(); i++){
        new_data[i] = (data.at(2*i) + data.at(2*i+1))/2.;
    }
    return new_data;
}

std::string Sample::message() const{
    return std::to_string(this->mean()) + " +/- " + std::to_string(this->error());
}

BinningAnalysis Sample::bin_it() const{
    return this->sample;
}

void BinningAnalysis::_init(std::vector<double> sample){ //passing a modifiable copy
    //We cannot keep all the samples, only a power of 2
    size_t min_samples = 64;
    size_t max_level = std::abs(log2(sample.size()/min_samples));
    size_t max_samples = min_samples * pow(2, max_level);
    //keep only the last N samples, where N is the maximum possible power of 2
    sample = std::vector<double>(sample.end() - max_samples, sample.end());
    std::vector<double> errors(max_level+1);
    std::vector<double> rel_change(max_level);
    //Now we can bin our samples
    _samples.resize(max_level+1);
    _samples[0] = sample;
    errors[0] = _samples[0].error();
    for (size_t i=1; i<max_level+1; i++){
        _samples[i] = bin_it(_samples[i-1].sample);
        errors[i] = _samples[i].error();
        rel_change[i-1] = (errors[i]-errors[i-1])/errors[i];
    }
    if (rel_change.size() > 3){
        double mean_last_change = mean_value(std::vector<double>(rel_change.end()-3, rel_change.end()));
        _converged = (mean_last_change <= 0.05);
    }
    else{
        _converged = false;
    }

}

std::vector<State*> copy_states(const std::vector<State*>& states){

    std::vector<State*> res(states.size());
    for (size_t i=0; i<states.size(); i++){
        res[i] = states[i]->clone();
    }
    return res;
}