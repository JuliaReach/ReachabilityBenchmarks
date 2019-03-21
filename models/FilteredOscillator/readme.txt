### Introduction

This benchmark models a two-dimensional switched oscillator connected to
a filter that smooths the output signal. The models are parameterized 
over the order (number of 1st-order differential equations) of the filter.

The oscillator is modeled with 4 locations and 2 continuous variables,
x and y. The filter is a purely continuous system that takes as input 
the value of x, and outputs a smoothed signal z, which depends on x 
via a k-th order system of 1st-order linear differential equations.
The combined system has 4 locations and 2+k continuous variables.

### SX Model

The folder SX contains a version based on SX.

The SpaceEx files in the subfolders are extended by an additional variable
for early termination after one loop iteration in the automaton.
In the Julia version, this variable is optional.

### Options

To reproduce the settings for the [HSCC2019 Repeatability Evaluation (RE) package for the paper JuliaReach: a Toolbox for Set-Based Reachability](https://github.com/JuliaReach/HSCC2019_RE), use the option `one_loop_iteration = true`. This option amounts to restricting the number of jumps to 5.
