#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <cstdio>

/////////Forward declarations
void PlaneCollisionCalculus(float dt, int index);

namespace LilSpheres {
	extern const int maxParticles;
	extern int firstParticleIdx;
	float lifeExpectancy;
	extern int particleCount;
	extern void updateParticles(int startIdx, int count, float* array_data);
	extern glm::vec3 directorVector;
	extern glm::vec3 velocity;
}
namespace Utils 
{
	float standardDirectorVector[3] = { 1, 1, 1 };
	float standardVelocity[3] = { 4, -1, 1 };
	glm::vec3 floatToVec(float* values)
	{
		glm::vec3 temp = glm::vec3(values[0], values[1], values[2]);
		return temp;
	}
	//3 vertex of the cube defined in render.cpp file class Box()
	glm::vec3 points[3] = { glm::vec3(-5,0,-5),glm::vec3(5,0,-5),glm::vec3(5,0,5) };
}

namespace {
	static struct PhysParams {
		glm::vec3 acceleration;
		float min = 0.f;
		float max = 10.f;
	} p_pars;

	static struct ParticleSystem {
		glm::vec3 *position;
		float *timeLeft;
		int numParticles;
		glm::vec3 *directorVector;
		glm::vec3 *velocity;
	} s_PS;
}

void Exemple_GUI() {
	ImGui::SliderFloat("Min Position Range", &p_pars.min, 0.f, 4.f);
	ImGui::SliderFloat("Max Position Ramge", &p_pars.max, 6.f, 10.f);
	ImGui::SliderFloat("Life Expectancy in seconds", &LilSpheres::lifeExpectancy, 0.1f, 10.f);
	ImGui::SliderFloat3("Director Vector", Utils::standardDirectorVector, 0.f, 1.f);
	ImGui::SliderFloat3("Starting Velocity", Utils::standardVelocity, -10.f, 10.f);
}

void Exemple_PhysicsInit() 
{
	p_pars.acceleration = glm::vec3(0, -9.81, 0);
	s_PS.numParticles = 100;
	s_PS.directorVector = new glm::vec3[s_PS.numParticles];
	s_PS.velocity = new glm::vec3[s_PS.numParticles];
	s_PS.position = new glm::vec3[s_PS.numParticles];
	s_PS.timeLeft = new float[s_PS.numParticles];
	extern bool renderParticles; renderParticles = true;
	LilSpheres::firstParticleIdx = 0;
	LilSpheres::lifeExpectancy = 1.f;
	LilSpheres::particleCount = s_PS.numParticles;
	for (int i = 0; i < s_PS.numParticles; i++)
	{
		s_PS.directorVector[i] = Utils::floatToVec(Utils::standardDirectorVector);
		s_PS.velocity[i] = Utils::floatToVec(Utils::standardVelocity);
		s_PS.timeLeft[i] = LilSpheres::lifeExpectancy;
		float x = -5 + p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
		float y = p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
		float z = -5 + p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
		printf("Creating particle with position %.2f, %.2f, %.2f\n", x, y, z);
		s_PS.position[i] = glm::vec3(x, y, z);
	}
	LilSpheres::updateParticles(0, s_PS.numParticles, &(s_PS.position[0].x));
}

void Exemple_PhysicsUpdate(float dt) {
	for (int i = 0; i < s_PS.numParticles; i++) {		
		s_PS.timeLeft[i] -= dt;
		PlaneCollisionCalculus(dt, i);
		if (s_PS.timeLeft[i] < 0)
		{
			s_PS.directorVector[i] = Utils::floatToVec(Utils::standardDirectorVector);
			s_PS.velocity[i] = Utils::floatToVec(Utils::standardVelocity);
			s_PS.timeLeft[i] = LilSpheres::lifeExpectancy;
			float x = -5 + p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
			float y = p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
			float z = -5 + p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
			printf("Particle Died, creating new particle with position %.2f, %.2f, %.2f\n", x, y, z);
			s_PS.position[i] = glm::vec3(x, y, z);
		}
	}
	LilSpheres::updateParticles(0, s_PS.numParticles, &(s_PS.position[0].x));
}

void PlaneCollisionCalculus(float dt, int index)
{
	glm::vec3 tempPos = s_PS.position[index] + (dt * (s_PS.velocity[index] * s_PS.directorVector[index]));
	glm::vec3 tempVel = s_PS.velocity[index] + (dt * p_pars.acceleration);
	glm::vec3 PQ = glm::vec3(Utils::points[1] - Utils::points[0]);
	glm::vec3 PR = glm::vec3(Utils::points[2] - Utils::points[0]);
	glm::vec3 normal = glm::normalize(glm::cross(PR, PQ));
	if (glm::dot(normal, s_PS.position[index]) < 0 || glm::dot(normal, tempPos) < 0)
	{
		//Si colision
		s_PS.position[index] = tempPos - 2 * (glm::dot(normal, tempPos)) *  normal;
		s_PS.velocity[index] = tempVel - 2 * (glm::dot(normal, tempVel)) * normal;
	}
	else
	{
		//No colision
		s_PS.position[index] = s_PS.position[index] + (dt * (s_PS.velocity[index] * s_PS.directorVector[index]));
		s_PS.velocity[index] = s_PS.velocity[index] + (dt * p_pars.acceleration);
	}
}

void Exemple_PhysicsCleanup() {
	delete[] s_PS.position;
}

