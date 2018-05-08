This problem is based on the ARCH 2014 Benchmark
"Networked Cooperative Platoon of Vehicles for Testing Methods and Verification Tools"
by Ibtissem Ben Makhlouf and Stefan Kowalewski
http://cps-vo.org/node/12115

The benchmarks describes the 1-dimensional positions of three autonomous vehicles with a manually driven leader. The vehicles may lose and regain communication as described by timing parameters: min. and max. time in communication, and max. time without communication. Note that the min. time without communication is zero.

The variables are 
- e_i for the gap between vehicle i and its predecessor i-1, for i=1,2,3.
- v_i is the time derivative of e_i
- a_i the acceleration of vehicle i

The parameters are
- acc_min : min. acceleration of leader vehicle
- acc_max : max. acceleration of leader vehicle
- tc : max. time before losing communication
- tb : min. time before communication is lost
- tr : max. time before regaining communication


Specification:
The vehicles should not collide, i.e., e_i >= e_{i,min}, over an unbounded time- and switching horizon.

Comments:
One can derive e_{i,min} from the bounding box of the reachable states starting from equilibrium (e_i=v_i=a_i=0). Lower values indicate higher accuracy of the reach set approximation.



