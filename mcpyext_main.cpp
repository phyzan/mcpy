#include "mcpyext_base.hpp"

PYBIND11_MODULE(mcpy, m)
{
    define_base_module(m);
}


/*
In order to compile a python extension "mcpy", run the following command. Place the mcpy.pyi stub file next to the compile module, to assist type-hinting.
g++ -O3 -Wall -march=x86-64 -shared -std=c++20 -fopenmp -I/usr/include/python3.12 -I/usr/include/pybind11 -fPIC $(python3 -m pybind11 --includes) tools.cpp mc.cpp ising.cpp mcpyext_base.cpp mcpyext_main.cpp -o mcpy/mcpy$(python3-config --extension-suffix)
*/


//if you are compiling a pure c++ program where you run a test code in main.cpp, run this:
//g++ -O3 -Wall -march=x86-64 -std=c++20 tools.cpp mc.cpp ising.cpp main.cpp -o test