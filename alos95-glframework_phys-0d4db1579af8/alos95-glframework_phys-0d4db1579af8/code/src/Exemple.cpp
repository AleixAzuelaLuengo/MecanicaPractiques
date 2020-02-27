#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <cstdio>

/////////Forward declarations

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
	class Plane
	{
		public:
		float x, y, z;
		float A, B, C, D;
		glm::vec3 normal;
		Plane() {
			x = y = z = A = B = C = D = 0;
			normal = glm::vec3(0);
		};
		Plane(glm::vec3 point1, glm::vec3 point2, glm::vec3 point3)
		{
			glm::vec3 PQ = glm::vec3(point2 - point1);
			glm::vec3 PR = glm::vec3(point3 - point1);
			normal = glm::normalize(glm::cross(PR, PQ));
			A = normal.x;
			B = normal.y;
			C = normal.z;
			D = (-A * point1.x - B * point1.y - C * point1.z);

			printf("%.2f, %.2f, %.2f, %.2f", A, B, C, D);
		};
		void constructorPlane(glm::vec3 point1, glm::vec3 rectPoint, glm::vec3 center)
		{
			normal = glm::normalize((rectPoint - center));
			A = normal.x;
			B = normal.y;
			C = normal.z;
			D = (-A * point1.x - B * point1.y - C * point1.z);
		};
		void constructPlanefromTwoPositions(glm::vec3 prevPos, glm::vec3 tempPos)
		{
			
		}
		bool hasCollisioned(glm::vec3 position, glm::vec3 tempPos)
		{
			float uno = (glm::dot(glm::normalize(normal), position) + D);
			float dos = (glm::dot(glm::normalize(normal), tempPos) + D);
			return (uno * dos <= 0);
		};
	};

	class Sphere
	{
	public:
		glm::vec3 center;
		float radius;
		Sphere();
		Sphere(glm::vec3 pos, float rad)
		{
			center = pos;
			radius = rad;
		}
	};

	Plane cubePlaneCollision[6];
	float standardDirectorVector[3] = { 0, 0, 0 };
	float standardVelocity[3] = { 0, -10, 0 };
	glm::vec3 floatToVec(float* values)
	{
		glm::vec3 temp = glm::vec3(values[0], values[1], values[2]);
		return temp;
	}
	//3 vertex of the cube defined in render.cpp file class Box()
	float cubeVerts[] = {
		// -5,0,-5 -- 5, 10, 5
		-5.f,  0.f, -5.f,
		 5.f,  0.f, -5.f,
		 5.f,  0.f,  5.f,
		-5.f,  0.f,  5.f,
		-5.f, 10.f, -5.f,
		 5.f, 10.f, -5.f,
		 5.f, 10.f,  5.f,
		-5.f, 10.f,  5.f,
	};
	glm::vec3 pointsPlane1[3] = { glm::vec3(-5,0,-5), glm::vec3(5,0,-5),glm::vec3(5,0,5) };  // down
	glm::vec3 pointsPlane2[3] = { glm::vec3(-5,0,5),  glm::vec3(-5,10,5), glm::vec3(-5,0,-5)};  // front
	glm::vec3 pointsPlane3[3] = { glm::vec3(-5,0,-5), glm::vec3(-5,10,-5),glm::vec3(5,0,-5) }; //right
	glm::vec3 pointsPlane4[3] = { glm::vec3(-5,0,5), glm::vec3(5,0,5),glm::vec3(5.f, 10.f,5.f) };//left
	glm::vec3 pointsPlane5[3] = { glm::vec3(5.f, 10.f,  5.f), glm::vec3(5,0,-5),glm::vec3(5,0,5) }; //back
	glm::vec3 pointsPlane6[3] = { glm::vec3(-5.f, 10.f, -5.f),glm::vec3(5.f, 10.f, -5.f),glm::vec3(5.f, 10.f,  5.f) };//up

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
	p_pars.acceleration = glm::vec3(0, -0, 0);
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
		//printf("Creating particle with position %.2f, %.2f, %.2f\n", x, y, z);
		s_PS.position[i] = glm::vec3(x, y, z);
	}
	Utils::cubePlaneCollision[0] = Utils::Plane(Utils::pointsPlane1[0], Utils::pointsPlane1[1], Utils::pointsPlane1[2]);
	Utils::cubePlaneCollision[1] = Utils::Plane(Utils::pointsPlane2[0], Utils::pointsPlane2[1], Utils::pointsPlane2[2]);
	Utils::cubePlaneCollision[2] = Utils::Plane(Utils::pointsPlane3[0], Utils::pointsPlane3[1], Utils::pointsPlane3[2]);
	Utils::cubePlaneCollision[3] = Utils::Plane(Utils::pointsPlane4[0], Utils::pointsPlane4[1], Utils::pointsPlane4[2]);
	Utils::cubePlaneCollision[4] = Utils::Plane(Utils::pointsPlane5[0], Utils::pointsPlane5[1], Utils::pointsPlane5[2]);
	Utils::cubePlaneCollision[5] = Utils::Plane(Utils::pointsPlane6[0], Utils::pointsPlane6[1], Utils::pointsPlane6[2]);
	LilSpheres::updateParticles(0, s_PS.numParticles, &(s_PS.position[0].x));
}

void PlaneCollisionCalculus(float dt, int index, Utils::Plane plane)
{
	glm::vec3 tempPos = s_PS.position[index] + (dt * (s_PS.velocity[index]));
	glm::vec3 tempVel = s_PS.velocity[index] + (dt * p_pars.acceleration);

	if (plane.hasCollisioned(s_PS.position[index], tempPos))
	{
		//Si colision
		s_PS.position[index] = tempPos - 2 * (glm::dot(plane.normal, tempPos) + plane.D) *  glm::normalize(plane.normal);
		s_PS.velocity[index] = tempVel - 2 * (glm::dot(plane.normal, tempVel)) * glm::normalize(plane.normal);
	}
	else
	{
		//No colision
		s_PS.position[index] = s_PS.position[index] + (dt * (s_PS.velocity[index]));
		s_PS.velocity[index] = s_PS.velocity[index] + (dt * p_pars.acceleration);
	}
}

void SphereCollisionCalculus(float dt, int index, Utils::Sphere sphere)
{
	Utils::Plane spherePlane;
	glm::vec3 tempPos = s_PS.position[index] + (dt * (s_PS.velocity[index] * s_PS.directorVector[index]));
	glm::vec3 tempVel = s_PS.velocity[index] + (dt * p_pars.acceleration);
	glm::vec3 rectVec = tempPos - s_PS.position[index];
	//Second Grade Equation Solved
	float landaSquared = (glm::pow(rectVec.x, 2)) + (glm::pow(rectVec.y, 2)) + (glm::pow(rectVec.z, 2));
	float landa = (2 * (s_PS.position[index].x * rectVec.x)) + (2 * (s_PS.position[index].y * rectVec.y)) + (2 * (s_PS.position[index].z * rectVec.z));
	float number = (glm::pow(s_PS.position[index].x, 2)) + (glm::pow(s_PS.position[index].y, 2)) + (glm::pow(s_PS.position[index].z, 2));
	float solutionLandaPlus = (-landa + glm::sqrt(glm::pow(landa, 2) - 4 * (landaSquared)* number)) / (2 * landaSquared);
	float solutionLandaMinus = (-landa - glm::sqrt(glm::pow(landa, 2) - 4 * (landaSquared)* number)) / (2 * landaSquared);
	//Calculate the 2 point to see which one is the closest one
	glm::vec3 recPointPlus = s_PS.position[index] + (solutionLandaPlus * (rectVec));
	glm::vec3 recPointMinus = s_PS.position[index] + (solutionLandaMinus * (rectVec));
	if (glm::distance(s_PS.position[index], recPointPlus) > glm::distance(s_PS.position[index], recPointMinus))
	{
		spherePlane.constructorPlane(recPointMinus, recPointMinus, sphere.center);
	}
	else
	{
		spherePlane.constructorPlane(recPointPlus, recPointMinus, sphere.center);
	}
	PlaneCollisionCalculus(dt, index, spherePlane);
}

void Exemple_PhysicsUpdate(float dt) {
	for (int i = 0; i < s_PS.numParticles; i++) {		
		s_PS.timeLeft[i] -= dt;

		//for(int j=0; j<5; j++)
			PlaneCollisionCalculus(dt, i, Utils::cubePlaneCollision[6]);
			if(i ==0)printf("Particle has position %.2f, %.2f, %.2f\n", s_PS.position[0].x, s_PS.position[0].y, s_PS.position[0].z);
		if (s_PS.timeLeft[i] < 0)
		{
			s_PS.directorVector[i] = Utils::floatToVec(Utils::standardDirectorVector);
			s_PS.velocity[i] = Utils::floatToVec(Utils::standardVelocity);
			s_PS.timeLeft[i] = LilSpheres::lifeExpectancy;
			float x = -5 + p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
			float y = p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
			float z = -5 + p_pars.min + (float)rand() / (RAND_MAX / (p_pars.max - p_pars.min));
			//printf("Particle Died, creating new particle with position %.2f, %.2f, %.2f\n", x, y, z);
			s_PS.position[i] = glm::vec3(x, y, z);
		}
	}
	LilSpheres::updateParticles(0, s_PS.numParticles, &(s_PS.position[0].x));
}

void Exemple_PhysicsCleanup() {
	delete[] s_PS.position;
}

