#include "particle.h"

typedef vector<int> Cell;
typedef vector<vector<int>> type1;

class sph
{
public:
	sph();
	sph(int x, int r, float step);
	vector<Particle> particleSystem;
	vector<vector<int>> neighbor_Particles;
	vector<vector<Cell>> cells;

	vec2 initParticlesCells();
	void updateSPHSystem(float step);
	void renderParticles(GLuint _textureId);
	float maxElem(float a, float b);	
	void find_Neighbor_Particles();
	void calc_dens_pres();
	void calc_forces(float step);	
	void boundaryConditions();
	void updateVelPos(float dt);
	/*PCISPH*/
	void calc_forces1(float step);
	void calc_density();
	void updatePCISPHSystem(float step);
};

