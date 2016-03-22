I. Synopsis
------------
This program simulates a synchronization algorithm of three coupled harmonic oscillators connected via a directed
strongly connected graph with an adjacency matrix given by 
            [0 1 1]
        G = [1 0 1]
            [1 0 0]
where the graphs are only able to communicate at stocastically determined isolated events.  

II. File List
------------
run.m 					%Initializes paramters, calls hybridsolver.m, and plots output
C.m						% Flow Set
f.m 					% Flow Map
D.m 					% Jump Set
g.m 					% Jump Map
hybridsolver.m 			% Hybrid equations solver

Code was built using Matlab 2015a. Other builds are untested. 