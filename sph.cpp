#include "sph.h"
#define PI 3.14159265f 
vec2 CELLS_XY;
int num_Particles;
float rest_density = 1000;
float stiffness = 10000;
int	num, radius, p_x, p_y;
float mass, H, step; 
sph::sph()
{

}
sph::sph(int x, int r, float stp)
{
	num = x;
	radius = r;
	step = stp;
	mass = (pow((1.0f/num),2)) * 1000;
	H = 2 * (1.0f/num);
	p_x = num/2.0f;	p_y = num;
	num_Particles = p_x * p_y;
	particleSystem = vector<Particle>();

	for (int i = 0; i < p_x; i++)
	{
		for (int j = 0; j < p_y; j++)
		{
			vec2 t1, t2, pos;
			t1.x = (2.0 - (2.0 / 2.5f))/2;
			t1.y = 1.5 - (2.0f*1.5 / 2.0f);
			t2.x = i*(2.0 / 2.5f) / p_x;
			t2.y = j*(2.0f*1.5 / 2.0f) / p_y;
			pos.x = t1.x + t2.x;
			pos.y = t1.y + t2.y;
			Particle p = Particle(pos);
			p.mass = mass;
			particleSystem.push_back(p);
		}
	}

	vec2 cells_XY1 = initParticlesCells();
	CELLS_XY.x = cells_XY1.x;
	CELLS_XY.y = cells_XY1.y;
	cells = vector<vector<Cell>>(CELLS_XY.x, vector<Cell>(CELLS_XY.y, Cell()));
	for (int i = 0; i < particleSystem.size(); i++)
	{
		int xx = particleSystem[i].pos.x / H;
		int yy = particleSystem[i].pos.y / H;
		cells[xx][yy].push_back(i);
	}
	cout << "Number of particles : " << num_Particles << " particles." << endl;
}

vec2 sph::initParticlesCells()
{
	int cells_X, cells_Y;
	float tempX = 2.0 / H;
	cells_X = (int)(tempX)+1;
	float tempY = 1.5 / H;
	cells_Y = (int)(tempY)+1;
	cout << "cells_X : " << cells_X << "cells_Y : " << cells_Y << endl;
	vec2 res;
	res.x = cells_X; res.y = cells_Y;
	return res;
}

void sph::renderParticles(GLuint _textureId)
{
	for (int i = 0; i < num_Particles; i++)
	{
		vec2 t;
		t.x = 0.0f; t.y = 0.0f;
		particleSystem[i].force = t;
	}
	for (int i = 0; i < num_Particles; i++)
	{
		float px = particleSystem[i].pos.x * 100;
		float py = particleSystem[i].pos.y * 100;
		float fact = 10.0f;
		if (_textureId == -100)
		{
			GLUquadric *quad;
			glPushMatrix();
			quad = gluNewQuadric();
			glColor3f(0.3f, 0.0f, 0.6f);
			glTranslatef(px, py, 0.0f);
			gluSphere(quad, radius, 100, 20);
			glPopMatrix();
			/*POINTS
			/*glEnable(GL_BLEND);
			glEnable(GL_POINT_SMOOTH);
			glPointSize(50.0f);
			glColor3f(1.0f, 1.0f, 1.0f);		
			glBegin(GL_POINTS);
			glVertex2f(px, py);
			glEnd();*/
		}
		else
		{
			glPushMatrix();
			glDisable(GL_BLEND);
			GLUquadricObj *qObj = gluNewQuadric();
			glBindTexture(GL_TEXTURE_2D, _textureId);
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			gluQuadricTexture(qObj, GL_TRUE);
			gluQuadricNormals(qObj, GLU_SMOOTH);
			glBindTexture(GL_TEXTURE_2D, _textureId);
			glTranslatef(px, py, 0.0f);
			gluSphere(qObj, radius, 100, 20);
			glDisable(GL_TEXTURE_2D);
			glPopMatrix();
		}
	}
}

void sph::updateSPHSystem(float step)
{
	//followed steps of the algorithm
	find_Neighbor_Particles();
	calc_dens_pres();
	calc_forces(step);
	updateVelPos(step);
	boundaryConditions();

	cells = vector<vector<Cell>>(CELLS_XY.x, vector<Cell>(CELLS_XY.y, Cell()));
	for (int i = 0; i < particleSystem.size(); i++)
	{
		int xx = particleSystem[i].pos.x / H;
		int yy = particleSystem[i].pos.y / H;
		cells[xx][yy].push_back(i);
	}
}

void sph::find_Neighbor_Particles()
{
	neighbor_Particles = vector<vector<int>>();
	for (Particle &p : particleSystem)
	{
		vector<int> allNeighbors_P = vector<int>();
		vector<Cell> neighbors = vector<Cell>();
		int xx = p.pos.x / H;
		int yy = p.pos.y / H;
		neighbors.push_back(cells[xx][yy]);
		if (xx > 0) neighbors.push_back(cells[xx - 1][yy]);
		if (yy > 0) neighbors.push_back(cells[xx][yy - 1]);
		if (xx > 0 && yy > 0) neighbors.push_back(cells[xx - 1][yy - 1]);
		if (xx < (CELLS_XY.x - 1)) neighbors.push_back(cells[xx + 1][yy]);
		if (yy < (CELLS_XY.y - 1)) neighbors.push_back(cells[xx][yy + 1]);
		if (xx > 0 && yy < (CELLS_XY.y - 1)) neighbors.push_back(cells[xx - 1][yy + 1]);
		if (yy > 0 && xx < (CELLS_XY.x - 1)) neighbors.push_back(cells[xx + 1][yy - 1]);
		if (xx < (CELLS_XY.x - 1) && yy < (CELLS_XY.y - 1)) neighbors.push_back(cells[xx + 1][yy + 1]);

		for (Cell &cell : neighbors)
		{
			for (int index : cell)
			{
				vec2 xx; 
				xx.x = p.pos.x - particleSystem[index].pos.x;
				xx.y = p.pos.y - particleSystem[index].pos.y;
				float dist2 = xx.x * xx.x + xx.y * xx.y;
				if (dist2 <= (H * H))
					allNeighbors_P.push_back(index);
			}
		}
		neighbor_Particles.push_back(allNeighbors_P);
	}
}

void sph::calc_dens_pres()
{
	
	for (int i = 0; i < num_Particles; i++)
	{
		vector<int> neighbors = neighbor_Particles[i];
		float Poly6_kernel, density = 0.0f;
		for (int n = 0; n < neighbors.size(); n++)
		{
			int j = neighbors[n];
			vec2 xx;
			xx.x = particleSystem[i].pos.x - particleSystem[j].pos.x;
			xx.y = particleSystem[i].pos.y - particleSystem[j].pos.y;
			float r_2 = xx.x * xx.x + xx.y * xx.y;
			float H_2 = H * H;
			if (r_2 < 0 || r_2 > H_2)
				Poly6_kernel = 0.0f;
			else
				Poly6_kernel = 315.0f / (64.0f * PI * pow(H, 9)) * pow(H_2 - r_2, 3);
			density += particleSystem[j].mass * Poly6_kernel;
		}
		particleSystem[i].den = density;
		particleSystem[i].pres = maxElem(stiffness * (particleSystem[i].den - rest_density), 0.0f);
	}
}

float sph::maxElem(float a, float b)
{
	return ((a > b) ? a : b);
}

void sph::calc_forces(float step)
{
	for (int i = 0; i < num_Particles; i++)
	{
		float r_2, visc_kernel;
		vec2 t, spiky_kernel, gravity;
		t.x = 0.0f; t.y = 0.0;
		vec2 pres_force = t; 
		vec2 grav_force = t;
		vec2 visc_force = t; 
		vector<int> neighbors = neighbor_Particles[i];
		for (int n = 0; n < neighbors.size(); n++)
		{
			int j = neighbors[n];
			vec2 xx; 
			xx.x = particleSystem[i].pos.x - particleSystem[j].pos.x;
			xx.y = particleSystem[i].pos.y - particleSystem[j].pos.y;

			//using spiky kernel for pressure (gradient)
			r_2 = sqrt(xx.x * xx.x + xx.y * xx.y);
			if (r_2 == 0.0f)
			{
				spiky_kernel.x = 0.0f;
				spiky_kernel.y = 0.0f;
			}
			else
			{
				spiky_kernel.x = (-45.0f / (PI * pow(H, 6))) * (xx.x / r_2) * pow(H - r_2, 2);
				spiky_kernel.y = (-45.0f / (PI * pow(H, 6))) * (xx.y / r_2) * pow(H - r_2, 2); 
			}
			
			pres_force.x += particleSystem[j].mass * (particleSystem[i].pres + particleSystem[j].pres) / (2.0f * particleSystem[j].den) * spiky_kernel.x;
			pres_force.y += particleSystem[j].mass * (particleSystem[i].pres + particleSystem[j].pres) / (2.0f * particleSystem[j].den) * spiky_kernel.y;
		
			//using viscosity kernel (Laplacian)
			r_2 = sqrt(xx.x * xx.x + xx.y * xx.y);
			visc_kernel = 45.0f / (PI * pow(H, 6)) * (H - r_2);
			visc_force.x += particleSystem[j].mass * (particleSystem[j].vel.x - particleSystem[i].vel.x) / particleSystem[j].den * visc_kernel;
			visc_force.y += particleSystem[j].mass * (particleSystem[j].vel.y - particleSystem[i].vel.y) / particleSystem[j].den * visc_kernel;
		}

		//Gravity
		gravity.x = 0.0f; gravity.y = -13000;
		grav_force.x = particleSystem[i].den * gravity.x;
		grav_force.y = particleSystem[i].den * gravity.y;

		float visc = 18000;
		particleSystem[i].force.x += (-1 * pres_force.x) + (visc_force.x*visc) + grav_force.x;
		particleSystem[i].force.y += (-1 * pres_force.y) + (visc_force.y*visc) + grav_force.y;

	}
}

void sph::updateVelPos(float step)
{
	for (int i = 0; i < num_Particles; i++)
	{
		particleSystem[i].vel.x += step * particleSystem[i].force.x / particleSystem[i].den;
		particleSystem[i].vel.y += step * particleSystem[i].force.y / particleSystem[i].den;
		particleSystem[i].pos.x += step * particleSystem[i].vel.x;
		particleSystem[i].pos.y += step * particleSystem[i].vel.y;
	}
}

void sph::boundaryConditions()
{
	float px, py, vx, vy;
	for (int i = 0; i < num_Particles; i++)
	{
		px = particleSystem[i].pos.x;
		vx = particleSystem[i].vel.x;
		py = particleSystem[i].pos.y;
		vy = particleSystem[i].vel.y;
		if (px < 0.0f){
			particleSystem[i].pos.x = 0.0f;
			particleSystem[i].vel.x = -1*(vx/2);
		}
		if (py < 0.0f){
			particleSystem[i].pos.y = 0.0f;
			particleSystem[i].vel.y = -1 * (vy / 2);
		}
		if (px > 2.0){
			particleSystem[i].pos.x = 2.0f;
			particleSystem[i].vel.x = -1 * (vx/2);
		}
		if (py > 1.5){
			particleSystem[i].pos.y = 1.5f;
			particleSystem[i].vel.y = -1 * (vy/2);
		}
	}
}

/*PCISPH components*/
void sph::updatePCISPHSystem(float step)
{
	//followed steps of the algorithm
	find_Neighbor_Particles();
	calc_density();
	calc_forces1(step);

	for (int i = 0; i < num_Particles; i++)
	{
		particleSystem[i].pres = 0.0f;
	}

	float errD = 0.0f;
	int iter = 0, minIteration = 5;
	while (iter <= minIteration)
	{
		iter += 1;
		for (int i = 0; i < num_Particles; i++)
		{
			particleSystem[i].vel_p.x += step * particleSystem[i].force.x / particleSystem[i].den;
			particleSystem[i].vel_p.y += step * particleSystem[i].force.y / particleSystem[i].den;
			particleSystem[i].pos_p.x += step * particleSystem[i].vel_p.x;
			particleSystem[i].pos_p.y += step * particleSystem[i].vel_p.y;
		}

		for (int i = 0; i < num_Particles; i++)
		{
			vector<int> neighbors = neighbor_Particles[i];
			float Poly6_kernel, density_p = 0.0f;
			for (int n = 0; n < neighbors.size(); n++)
			{
				int j = neighbors[n];
				vec2 xx;
				xx.x = particleSystem[i].pos.x - particleSystem[j].pos.x;
				xx.y = particleSystem[i].pos.y - particleSystem[j].pos.y;
				float r_2 = xx.x * xx.x + xx.y * xx.y;
				float H_2 = H * H;

				if (r_2 < 0 || r_2 > H_2)
					Poly6_kernel = 0.0f;
				else
					Poly6_kernel = 315.0f / (64.0f * PI * pow(H, 9)) * pow(H_2 - r_2, 3);
				density_p += particleSystem[j].mass * Poly6_kernel;
			}
			particleSystem[i].den_p = density_p;
		}
	}

	for (int i = 0; i < num_Particles; i++)
	{
		particleSystem[i].pres = maxElem(stiffness * (particleSystem[i].den_p - rest_density), 0.0f);
	}
	calc_forces(step);
	updateVelPos(step);
	boundaryConditions();

	cells = vector<vector<Cell>>(CELLS_XY.x, vector<Cell>(CELLS_XY.y, Cell()));
	for (int i = 0; i < particleSystem.size(); i++)
	{
		int xx = particleSystem[i].pos.x / H;
		int yy = particleSystem[i].pos.y / H;
		cells[xx][yy].push_back(i);
	}
}

void sph::calc_density()
{
	for (int i = 0; i < num_Particles; i++)
	{
		vector<int> neighbors = neighbor_Particles[i];
		float Poly6_kernel, density = 0.0f;
		for (int n = 0; n < neighbors.size(); n++)
		{
			int j = neighbors[n];
			vec2 xx;
			xx.x = particleSystem[i].pos.x - particleSystem[j].pos.x;
			xx.y = particleSystem[i].pos.y - particleSystem[j].pos.y;
			float r_2 = xx.x * xx.x + xx.y * xx.y;
			float H_2 = H * H;
			if (r_2 < 0 || r_2 > H_2)
				Poly6_kernel = 0.0f;
			else
				Poly6_kernel = 315.0f / (64.0f * PI * pow(H, 9)) * pow(H_2 - r_2, 3);
			density += particleSystem[j].mass * Poly6_kernel;
		}
		particleSystem[i].den = density;
	}
}

void sph::calc_forces1(float step)
{
	for (int i = 0; i < num_Particles; i++)
	{
		float r_2, visc_kernel;
		vec2 t, spiky_kernel, gravity;
		t.x = 0.0f; t.y = 0.0;
		vec2 pres_force = t;
		vec2 grav_force = t;
		vec2 visc_force = t;
		vector<int> neighbors = neighbor_Particles[i];
		for (int n = 0; n < neighbors.size(); n++)
		{
			int j = neighbors[n];
			vec2 xx;
			xx.x = particleSystem[i].pos.x - particleSystem[j].pos.x;
			xx.y = particleSystem[i].pos.y - particleSystem[j].pos.y;

			//using spiky kernel for pressure (gradient)
			r_2 = sqrt(xx.x * xx.x + xx.y * xx.y);
			if (r_2 == 0.0f)
			{
				spiky_kernel.x = 0.0f;
				spiky_kernel.y = 0.0f;
			}
			else
			{
				spiky_kernel.x = (-45.0f / (PI * pow(H, 6))) * (xx.x / r_2) * pow(H - r_2, 2);
				spiky_kernel.y = (-45.0f / (PI * pow(H, 6))) * (xx.y / r_2) * pow(H - r_2, 2);
			}

			pres_force.x += particleSystem[j].mass * (particleSystem[i].pres + particleSystem[j].pres) / (2.0f * particleSystem[j].den) * spiky_kernel.x;
			pres_force.y += particleSystem[j].mass * (particleSystem[i].pres + particleSystem[j].pres) / (2.0f * particleSystem[j].den) * spiky_kernel.y;


			//using viscosity kernel (Laplacian)
			r_2 = sqrt(xx.x * xx.x + xx.y * xx.y);
			visc_kernel = 45.0f / (PI * pow(H, 6)) * (H - r_2);
			visc_force.x += particleSystem[j].mass * (particleSystem[j].vel.x - particleSystem[i].vel.x) / particleSystem[j].den * visc_kernel;
			visc_force.y += particleSystem[j].mass * (particleSystem[j].vel.y - particleSystem[i].vel.y) / particleSystem[j].den * visc_kernel;
		}

		//Gravity
		gravity.x = 0.0f; gravity.y = -12000;
		grav_force.x = particleSystem[i].den * gravity.x;
		grav_force.y = particleSystem[i].den * gravity.y;

		float visc = 12000;
		particleSystem[i].force.x += (visc_force.x*visc) + grav_force.x;
		particleSystem[i].force.y += (visc_force.y*visc) + grav_force.y;

	}
}