#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <cstdio>
#include <stdlib.h>
#include <random>
#include <vector>

/////////Forward declarations
extern bool renderSphere;
extern bool renderCapsule;
extern bool renderParticles;
extern bool renderCloth;
extern bool renderCube;
glm::vec3 *particlePosition;
namespace {
	static struct PhysParams {
		glm::vec3 acceleration;
		float pPos[3] = { 2.5f, 1.f, 2.5f };
		float min = 0.f;
		float max = 10.f;
	} p_pars;

	static struct ParticleSystem {
		glm::vec3 spawnPosDefault = glm::vec3(0, 5, 0);
		glm::vec3 startingSpawnPoint = glm::vec3(0, 3, 0);
		glm::vec3 spawnVelDefalut = glm::vec3(1, 2, 0);
		float displacement = -5 + startingSpawnPoint.x;
		int spawnVelYMin = 2;
		int spawnVelYMax = 5;
		int spawnVelZMin = -2;
		int spawnVelZMax = 2;
		char* items[2] = { "Fountain","Waterfall" };
		char* current_item = items[0];
		float emissionRate = 100;
		std::vector<glm::vec3> position;
		std::vector<float> timeLeft;
		int numParticles;
		int particlesPerFrame;
		std::vector<glm::vec3> velocity;
		void updateParticlesPerFrame(int newN) { particlesPerFrame = newN; };
	} s_PS;
}

namespace LilSpheres {
	extern const int maxParticles;
	extern int firstParticleIdx;
	float lifeExpectancy;
	extern int particleCount;
	extern void updateParticles(int startIdx, int count, float* array_data);
	extern glm::vec3 directorVector;
	extern glm::vec3 velocity;
}

namespace Sphere {
	extern void updateSphere(glm::vec3 pos, float radius = 1.f);
}
namespace Capsule {
	extern void updateCapsule(glm::vec3 posA, glm::vec3 posB, float radius = 1.f);
}
namespace Utils
{
	class Plane
	{
	public:
		float A, B, C, D;
		glm::vec3 normal;
		Plane() {
			A = B = C = D = 0;
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
			return (((glm::dot((normal), position) + D) * (glm::dot((normal), tempPos) + D)) <= 0);
		};

		void CalculusPostPlaneCollision(glm::vec3 &tempPos, glm::vec3 &tempVel)
		{
			//Get the next position and velocity and calculate where the particle should be when it bounces
			tempPos -= 2 * (glm::dot(normal, tempPos) + D) * glm::normalize(normal);
			tempVel -= 2 * (glm::dot(normal, tempVel)) * glm::normalize(normal);
		}
	};

	class Sphere
	{
	public:
		glm::vec3 center;
		float radius;
		Sphere()
		{
			center = glm::vec3(0);
			radius = 0;
		}
		Sphere(glm::vec3 pos, float rad)
		{
			center = pos;
			radius = rad;
		}

		bool hasCollisioned(glm::vec3 tempPos)
		{
			return glm::distance(tempPos, center) < radius;
		}


		void SphereCollisionCalculus(glm::vec3 &tempPos, glm::vec3 &tempVel, int index, ParticleSystem particle)
		{
			glm::vec3 normal;
			float D;
			float d = glm::distance(tempPos, center);
			//landa squared
			float landaSquared = (glm::pow(particle.velocity[index].x, 2)) + (glm::pow(particle.velocity[index].y, 2)) + (glm::pow(particle.velocity[index].z, 2));
			//landa
			float landaX = (2 * (particle.position[index].x * -particle.velocity[index].x)) + (2 * (particle.velocity[index].x) * -center.x);
			float landaY = (2 * (particle.position[index].y * -particle.velocity[index].y)) + (2 * (particle.velocity[index].y) * -center.y);
			float landaZ = (2 * (particle.position[index].z * -particle.velocity[index].z)) + (2 * (particle.velocity[index].z) * -center.z);
			float landa = landaX + landaY + landaZ;
			//independent number
			float numberX = glm::pow(particle.position[index].x, 2) + glm::pow(-center.x, 2) + (2 * (particle.position[index].x * -center.x));
			float numberY = glm::pow(particle.position[index].y, 2) + glm::pow(-center.y, 2) + (2 * (particle.position[index].y * -center.y));
			float numberZ = glm::pow(particle.position[index].z, 2) + glm::pow(-center.z, 2) + (2 * (particle.position[index].z * -center.z));
			float number = numberX + numberY + numberZ - glm::pow(radius, 2);
			//2nd grade equations solved
			float solutionLandaPlus = (-landa + glm::sqrt(glm::pow(landa, 2) + (-4 * (landaSquared)* number))) / (2 * landaSquared);
			float solutionLandaMinus = (-landa - glm::sqrt(glm::pow(landa, 2) + (-4 * (landaSquared)* number))) / (2 * landaSquared);
			//Calculate the 2 point to see which one is the closest one
			glm::vec3 recPointPlus = particle.position[index] + (solutionLandaPlus * particle.velocity[index]);
			glm::vec3 recPointMinus = particle.position[index] + (solutionLandaMinus * particle.velocity[index]);
			if (glm::distance(particle.position[index], recPointPlus) > glm::distance(particle.position[index], recPointMinus))
			{
				normal = recPointMinus - center;
				normal = glm::normalize(normal);
				D = -((normal.x * recPointMinus.x) + (normal.y * recPointMinus.y) + (normal.z * recPointMinus.z));
			}
			else
			{
				normal = recPointPlus - center;
				normal = glm::normalize(normal);
				D = -((normal.x * recPointPlus.x) + (normal.y * recPointPlus.y) + (normal.z * recPointPlus.z));
			}
			tempPos -= 2 * (glm::dot(normal, tempPos) + D) * (normal);
			tempVel -= 2 * (glm::dot(normal, tempVel)) * (normal);
		}
	};

	class Capsule
	{
		enum capsuleCollisions { TOPSPHERE, BOTTOMSPHERE, CYLLINDER, NONE };
	public:
		Sphere topSemiSphere;
		Sphere bottomSemiSphere;

		Capsule()
		{
			topSemiSphere = Sphere();
			bottomSemiSphere = Sphere();
		}

		Capsule(Sphere top, Sphere bottom)
		{
			topSemiSphere = top;
			bottomSemiSphere = bottom;
		}

		bool hasCollisioned(glm::vec3 tempPos)
		{
			float distToTopSphere;
			float distToBottomSphere;
			bool collisionTop = false;
			bool collisionBottom = false;
			bool collisionCyllinder = false;
			//OPTION A
			////Calculate the distances to see if we colldie
			glm::vec3 AB = topSemiSphere.center - bottomSemiSphere.center;
			glm::vec3 cross = glm::cross((tempPos - bottomSemiSphere.center), AB);
			float ABmodule = (glm::sqrt(glm::pow(AB.x, 2) + glm::pow(AB.y, 2) + glm::pow(AB.z, 2)));
			float crossModule = glm::sqrt(glm::pow(cross.x, 2) + glm::pow(cross.y, 2) + glm::pow(cross.z, 2));
			float distToRect = crossModule / ABmodule;

			/*//OPTION B
			glm::vec3 closestPoint = getClosestPoint(topSemiSphere.center, bottomSemiSphere.center, tempPos);
			float distToRect = glm::distance(closestPoint,tempPos);*/
			if (topSemiSphere.hasCollisioned(tempPos))
			{
				collisionTop = true;
				distToTopSphere = glm::distance(tempPos, topSemiSphere.center);
			}
			if (bottomSemiSphere.hasCollisioned(tempPos))
			{
				collisionBottom = true;
				distToBottomSphere = glm::distance(tempPos, bottomSemiSphere.center);
			}
			if (distToRect < topSemiSphere.radius || distToRect < bottomSemiSphere.radius)
			{
				collisionCyllinder = true;
			}
			//Calculate which distance is the lowest
			if (collisionTop && collisionCyllinder)
			{
				if (distToTopSphere < distToRect)
					return true;
				else if (distToRect < distToTopSphere)
					return true;
				else
					return false;
			}
			else if (collisionBottom && collisionCyllinder)
			{
				if (distToBottomSphere < distToRect)
					return true;
				else if (distToRect < distToBottomSphere)
					return true;
				else
					return false;
			}
			else if (collisionTop)
				return true;
			else if (collisionBottom)
				return true;
			else if (collisionCyllinder)
				return true;
			else return false;
		}

		void CapsuleCollisionCalculus(glm::vec3 &tempPos, glm::vec3 &tempVel, int index, ParticleSystem particle)
		{
			if (hasCollisioned(tempPos))
			{
				glm::vec3 closestPoint = getClosestPoint(topSemiSphere.center, bottomSemiSphere.center, tempPos);
				glm::vec3 normal;
				float D;
				float d = glm::distance(tempPos, closestPoint);
				//landa squared
				float landaSquared = (glm::pow(particle.velocity[index].x, 2)) + (glm::pow(particle.velocity[index].y, 2)) + (glm::pow(particle.velocity[index].z, 2));
				//landa
				float landaX = (2 * (particle.position[index].x * -particle.velocity[index].x)) + (2 * (particle.velocity[index].x) * -closestPoint.x);
				float landaY = (2 * (particle.position[index].y * -particle.velocity[index].y)) + (2 * (particle.velocity[index].y) * -closestPoint.y);
				float landaZ = (2 * (particle.position[index].z * -particle.velocity[index].z)) + (2 * (particle.velocity[index].z) * -closestPoint.z);
				float landa = landaX + landaY + landaZ;
				//independent number
				float numberX = glm::pow(particle.position[index].x, 2) + glm::pow(-closestPoint.x, 2) + (2 * (particle.position[index].x * -closestPoint.x));
				float numberY = glm::pow(particle.position[index].y, 2) + glm::pow(-closestPoint.y, 2) + (2 * (particle.position[index].y * -closestPoint.y));
				float numberZ = glm::pow(particle.position[index].z, 2) + glm::pow(-closestPoint.z, 2) + (2 * (particle.position[index].z * -closestPoint.z));
				float number = numberX + numberY + numberZ - glm::pow(topSemiSphere.radius, 2);
				//2nd grade equations solved
				float solutionLandaPlus = (-landa + glm::sqrt(glm::pow(landa, 2) + (-4 * (landaSquared)* number))) / (2 * landaSquared);
				float solutionLandaMinus = (-landa - glm::sqrt(glm::pow(landa, 2) + (-4 * (landaSquared)* number))) / (2 * landaSquared);
				//Calculate the 2 point to see which one is the closest one
				glm::vec3 recPointPlus = particle.position[index] + (solutionLandaPlus * particle.velocity[index]);
				glm::vec3 recPointMinus = particle.position[index] + (solutionLandaMinus * particle.velocity[index]);
				if (glm::distance(particle.position[index], recPointPlus) > glm::distance(particle.position[index], recPointMinus))
				{
					normal = recPointMinus - closestPoint;
					normal = glm::normalize(normal);
					D = -((normal.x * recPointMinus.x) + (normal.y * recPointMinus.y) + (normal.z * recPointMinus.z));
				}
				else
				{
					normal = recPointPlus - closestPoint;
					normal = glm::normalize(normal);
					D = -((normal.x * recPointPlus.x) + (normal.y * recPointPlus.y) + (normal.z * recPointPlus.z));
				}
				tempPos -= 2 * (glm::dot(normal, tempPos) + D) * (normal);
				tempVel -= 2 * (glm::dot(normal, tempVel)) * (normal);
			}
		}

		glm::vec3 getClosestPoint(glm::vec3 A, glm::vec3 B, glm::vec3 particlePosition)
		{
			glm::vec3 AP = particlePosition - A;       //Vector from A to P   
			glm::vec3 AB = B - A;       //Vector from A to B  

			float magnitudeAB = pow(AB.x, 2) + pow(AB.y, 2) + pow(AB.z, 2);     //Magnitude of AB vector (it's length squared)     
			float ABAPproduct = glm::dot(AP, AB);		//The DOT product of AP and AB     
			float distance = ABAPproduct / magnitudeAB; //The normalized "distance" from a to your closest point  

			if (distance < 0)     //Check if P projection is over vectorAB     
				return A;
			else if (distance > 1)
				return B;
			else
				return A + AB * distance;
		}
	};

	Plane cubePlaneCollision[6];
	float standardDirectorVector[3] = { 0, 0, 0 };
	float standardVelocity[3] = { 0, -5, 0 };
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
	//Capsula
	Capsule capsule = Utils::Capsule(Utils::Sphere(glm::vec3(-3.f, 2.f, 0.f), 1.f), Utils::Sphere(glm::vec3(3.f, 2.f, 0.f), 1.f));
	//Esfera
	Sphere sphere = Utils::Sphere(glm::vec3(0.f, 1.f, 0.f), 1.f);
	//Cub
	glm::vec3 pointsPlane1[3] = { glm::vec3(-5,0,-5), glm::vec3(5,0,-5),glm::vec3(5,0,5) };  // down
	glm::vec3 pointsPlane2[3] = { glm::vec3(-5,0,5),  glm::vec3(-5,10,5), glm::vec3(-5,0,-5) };  // front
	glm::vec3 pointsPlane3[3] = { glm::vec3(-5,0,-5), glm::vec3(-5,10,-5),glm::vec3(5,0,-5) }; //right
	glm::vec3 pointsPlane4[3] = { glm::vec3(-5,0,5), glm::vec3(5,0,5),glm::vec3(5.f, 10.f,5.f) };//left
	glm::vec3 pointsPlane5[3] = { glm::vec3(5.f, 10.f,  5.f), glm::vec3(5,0,-5),glm::vec3(5,0,5) }; //back
	glm::vec3 pointsPlane6[3] = { glm::vec3(-5.f, 10.f, -5.f),glm::vec3(5.f, 10.f, -5.f),glm::vec3(5.f, 10.f,  5.f) };//up

}




void Exemple_GUI()
{
	if (ImGui::CollapsingHeader("PARTICLE VARIABLES"))
	{
		ImGui::SliderFloat("Emission Rate (particles / second)", &s_PS.emissionRate, 100.f, 400.f);
		ImGui::SliderFloat("Min Position Range", &p_pars.min, 0.f, 4.f);
		ImGui::SliderFloat("Max Position Range", &p_pars.max, 6.f, 10.f);
		ImGui::SliderFloat("Life Expectancy in seconds", &LilSpheres::lifeExpectancy, 1.f, 10.f);
	}
	if (ImGui::CollapsingHeader("EMISSION TYPE VARIABLES"))
	{
		if (ImGui::BeginCombo("Spawning Type", s_PS.current_item))
		{
			for (int n = 0; n < IM_ARRAYSIZE(s_PS.items); n++)
			{
				bool is_selected = (s_PS.current_item == s_PS.items[n]);
				if (ImGui::Selectable(s_PS.items[n], is_selected))
					s_PS.current_item = s_PS.items[n];
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}

			ImGui::EndCombo();
		}
		if (s_PS.current_item == s_PS.items[0]) // fountain
		{
			ImGui::Text("FOUNTAIN VARIABLES");
			ImGui::SliderFloat3("Starting position", &s_PS.spawnPosDefault.x, -5.f, 5.f);
			ImGui::SliderInt("Min Direction", &s_PS.spawnVelYMin, -5.f, 5.f);
			ImGui::SliderInt("Max Direction", &s_PS.spawnVelYMax, 1, 5.f);
		}
		else if (s_PS.current_item == s_PS.items[1]) //waterfall
		{
			ImGui::Text("WATERFALL VARIABLES");
			ImGui::SliderFloat3("Starting direction", &s_PS.spawnVelDefalut.x, -5.f, 5.f);
			ImGui::SliderFloat3("Starting Point", &s_PS.startingSpawnPoint.x, -5.f, 5.f);
			ImGui::SliderFloat("Max Displament on X axis", &s_PS.displacement, -0, 5 - s_PS.startingSpawnPoint.x);
		}
	}
	if (ImGui::CollapsingHeader("SPHERE VARIABLES"))
	{
		ImGui::Checkbox("Render Sphere", &renderSphere);
		if (ImGui::SliderFloat3("Sphere center", &Utils::sphere.center.x, -5.f, 5.f))
		{
			Sphere::updateSphere(Utils::sphere.center, Utils::sphere.radius);
		}
		if (ImGui::SliderFloat("Sphere radius", &Utils::sphere.radius, 1.f, 3.f))
		{
			Sphere::updateSphere(Utils::sphere.center, Utils::sphere.radius);
		}
	}
	if (ImGui::CollapsingHeader("CAPSULE VARIABLES"))
	{
		ImGui::Checkbox("Render Capsule", &renderCapsule);
		if (ImGui::SliderFloat3("Capsule Top Sphere center", &Utils::capsule.topSemiSphere.center.x, -5.f, 5.f))
		{
			Capsule::updateCapsule(Utils::capsule.topSemiSphere.center, Utils::capsule.bottomSemiSphere.center, Utils::capsule.topSemiSphere.radius);
		}
		if (ImGui::SliderFloat3("Capsule Bottom Sphere center", &Utils::capsule.bottomSemiSphere.center.x, -5.f, 5.f))
		{
			Capsule::updateCapsule(Utils::capsule.topSemiSphere.center, Utils::capsule.bottomSemiSphere.center, Utils::capsule.topSemiSphere.radius);
		}
		if (ImGui::SliderFloat("Capsule Radius", &Utils::capsule.topSemiSphere.radius, 1.f, 5.f))
		{
			Capsule::updateCapsule(Utils::capsule.topSemiSphere.center, Utils::capsule.bottomSemiSphere.center, Utils::capsule.topSemiSphere.radius);
		}

	}
}

void Exemple_PhysicsInit()
{
	std::mt19937 rng(0.154687f);
	Capsule::updateCapsule(Utils::capsule.topSemiSphere.center, Utils::capsule.bottomSemiSphere.center, Utils::capsule.topSemiSphere.radius);
	p_pars.acceleration = glm::vec3(0, -9.81, 0);
	s_PS.particlesPerFrame = glm::round(s_PS.emissionRate / 1000 * 33.3f);
	s_PS.numParticles = 0;
	renderParticles = true;
	LilSpheres::firstParticleIdx = 0;
	LilSpheres::lifeExpectancy = 1.f;
	Utils::cubePlaneCollision[0] = Utils::Plane(Utils::pointsPlane1[0], Utils::pointsPlane1[1], Utils::pointsPlane1[2]);
	Utils::cubePlaneCollision[1] = Utils::Plane(Utils::pointsPlane2[0], Utils::pointsPlane2[1], Utils::pointsPlane2[2]);
	Utils::cubePlaneCollision[2] = Utils::Plane(Utils::pointsPlane3[0], Utils::pointsPlane3[1], Utils::pointsPlane3[2]);
	Utils::cubePlaneCollision[3] = Utils::Plane(Utils::pointsPlane4[0], Utils::pointsPlane4[1], Utils::pointsPlane4[2]);
	Utils::cubePlaneCollision[4] = Utils::Plane(Utils::pointsPlane5[0], Utils::pointsPlane5[1], Utils::pointsPlane5[2]);
	Utils::cubePlaneCollision[5] = Utils::Plane(Utils::pointsPlane6[0], Utils::pointsPlane6[1], Utils::pointsPlane6[2]);
}

void Exemple_PhysicsUpdate(float dt) {
	s_PS.particlesPerFrame = glm::round(s_PS.emissionRate / 1000 * 33.3f);
	int temp = s_PS.numParticles;
	for (int i = s_PS.numParticles; i < s_PS.particlesPerFrame + temp; i++)
	{
		float x;
		float y;
		float z;
		s_PS.numParticles++;
		s_PS.timeLeft.push_back(LilSpheres::lifeExpectancy);
		if (s_PS.current_item == s_PS.items[0])
		{
			glm::vec3 randVelocity = glm::vec3((rand()) % (s_PS.spawnVelZMax - s_PS.spawnVelZMin + 1) + s_PS.spawnVelZMin, (rand()) % (s_PS.spawnVelYMax - s_PS.spawnVelYMin + 1) + s_PS.spawnVelYMin, 5);
			s_PS.velocity.push_back((randVelocity));
			x = s_PS.spawnPosDefault.x;
			y = s_PS.spawnPosDefault.y;
			z = s_PS.spawnPosDefault.z;
		}
		else if ((s_PS.current_item == s_PS.items[1]))
		{
			glm::vec3 randPos = glm::vec3((rand()) % (int)((s_PS.startingSpawnPoint.x + s_PS.displacement) - s_PS.startingSpawnPoint.x + 1) + (s_PS.startingSpawnPoint.x), s_PS.startingSpawnPoint.y, s_PS.startingSpawnPoint.z);
			s_PS.velocity.push_back((s_PS.spawnVelDefalut));
			x = randPos.x;
			y = randPos.y;
			z = randPos.z;
		}
		s_PS.position.push_back(glm::vec3(x, y, z));
	}
	for (int i = 0; i < s_PS.numParticles; i++) {
		s_PS.timeLeft[i] -= dt;
		//Calculate next position and next velocity
		glm::vec3 tempPos = s_PS.position[i] + (dt * (s_PS.velocity[i]));
		glm::vec3 tempVel = s_PS.velocity[i] + (dt * p_pars.acceleration);
		//Check collisions with all 6 planes in the cube 
		for (int j = 0; j < 6; j++) {
			if (Utils::cubePlaneCollision[j].hasCollisioned(s_PS.position[i], tempPos))
			{
				//Compute the collisions and save them in tempPos and tempVel
				Utils::cubePlaneCollision[j].CalculusPostPlaneCollision(tempPos, tempVel);
			}
		}
		//Collision with Sphere
		if (renderSphere)
		{
			if (Utils::sphere.hasCollisioned(tempPos))
			{
				Utils::sphere.SphereCollisionCalculus(tempPos, tempVel, i, s_PS);
			}
		}
		//Collision with Capsule
		if (renderCapsule)
		{
			Utils::capsule.CapsuleCollisionCalculus(tempPos, tempVel, i, s_PS);
		}

		//Change the position.
		//When we compute the collisions, all these changes are saved in tempPos and tempVel,
		//so these are the final positions after all collisions are computed.
		s_PS.position[i] = tempPos;
		s_PS.velocity[i] = tempVel;


		//Check life of the particle
		if (s_PS.timeLeft[i] < 0)
		{
			s_PS.numParticles--;
			s_PS.position.erase(s_PS.position.begin() + i);
			s_PS.velocity.erase(s_PS.velocity.begin() + i);
			s_PS.timeLeft.erase(s_PS.timeLeft.begin() + i);
			i--;
		}
	}
	particlePosition = new glm::vec3[s_PS.numParticles];
	for (int i = 0; i < s_PS.numParticles; i++)
	{
		particlePosition[i] = s_PS.position[i];
	}
	LilSpheres::particleCount = s_PS.numParticles;
	LilSpheres::updateParticles(0, s_PS.numParticles, &(particlePosition[0].x));
	delete[] particlePosition;
}

void Exemple_PhysicsCleanup() {
}

