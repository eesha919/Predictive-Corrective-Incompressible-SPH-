# Predictive-Corrective-Incompressible-SPH-

The paper[1] implemented by me presents a novel, incompressible fluid simulation method based on the SPH model.
In the implemented method, incompressibility is enforced by using a prediction-correction scheme to determine the particle pressures. This approach avoids the computational expenses of solving a pressure Poisson equation, while still being able to use large time steps in the simulation. The results obtained show that predictive-corrective incompressible SPH is smoother than the common SPH model and the visual outputs are in good agreement with each other. 

=>	The project consists of the following cpps and headers:
		- Source.cpp
		- sph.cpp
		- imageloader.cpp

		- sph.h
		- particle.h
		- imageloader.h

=>	This project is built on Visual Studio 2015.
	The packages folder consists of all the libraries it is dependent on,
		namely, opengl, glew, glut, devIL
		
=>	After running the code (better in Release mode), the demo of water simulation will start.
	Following are the user controls for the demo:
	->	To vary number of particles, 
			please use number keys 1,2,3,4
	->	To swich to PCISPH demo, use character key 'p' or 'P'
	->	To vary number of particles in PCISPH mode, 
			please use number keys 8,9
	->	To switch to SPH mode or to reset the demo, use 'spacebar'
	->	To render particles without texture, use 'x' or 'X' key
	->	To render particles with textures, use 'z' or 'Z' key
	->	To quit the demo, use 'q' or 'Q' key

References 
[1] B. Solanthaler, R. Pajarola. Predictive-Corrective Incompressible SPH. ACM Trans. Graph., 28, 3 August 2009. 
