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
long nanCrashes = 0;
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
		//Default spawning positions and velocity
		glm::vec3 spawnPosDefault = glm::vec3(0, 5, 0);
		glm::vec3 startingSpawnPoint = glm::vec3(0, 3, 0);
		glm::vec3 spawnVelDefalut = glm::vec3(1, 2, 0);

		//Displacement X for waterfall
		float displacement = -5 + startingSpawnPoint.x;

		//Min and Max spawn velocities to use random velocities when spawned
		int spawnVelYMin = 2;
		int spawnVelYMax = 5;
		int spawnVelZMin = -2;
		int spawnVelZMax = 2;

		//Items for drop field to choose between these two types
		char* items[2] = { "Fountain","Waterfall" };
		char* current_item = items[0];

		//Default emission rate
		float emissionRate = 100;

		//Particle system variables
		std::vector<glm::vec3> position;
		std::vector<glm::vec3> prevPosition;
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
			//Create null plane at 0
			A = B = C = D = 0;
			normal = glm::vec3(0);
		};
		Plane(glm::vec3 point1, glm::vec3 point2, glm::vec3 point3)
		{
			//Create a plane with 3 points.

			//We use 1 point and 2 vectors, and then find the normal to find the 3rd vector
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
		bool hasCollisioned(glm::vec3 position, glm::vec3 tempPos)
		{
			//Check if a particle has collisioned with the plane.
			//We get the particle previous position and the tempPos (the one that might have gone through the plane)
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
			//Define a sphere by its center and the radius
			center = pos;
			radius = rad;
		}

		bool hasCollisioned(glm::vec3 tempPos)
		{
			//Check if the distance between the position of the particle 
			//and the center of the sphere is smaller than the radius
			//to see if the particle is inside.
			return glm::distance(tempPos, center) < radius;
		}


		void SphereCollisionCalculus(glm::vec3 &tempPos, glm::vec3 &velocity, glm::vec3 &position)
		{
			//Check collision of the particle with the sphere


			glm::vec3 normal;
			float D = 0.f;

			//Get the distance between the particle position and the center
			//landa squared
			float landaSquared = (glm::pow(velocity.x, 2)) + (glm::pow(velocity.y, 2)) + (glm::pow(velocity.z, 2));
			//landa
			float landaX = -2 * (velocity.x * center.x) + 2 * (position.x * velocity.x);
			float landaY = -2 * (velocity.y * center.y) + 2 * (position.y * velocity.y);
			float landaZ = -2 * (velocity.z * center.z) + 2 * (position.z * velocity.z);
			float landa = landaX + landaY + landaZ;
			//independent number
			float numberX = glm::pow(position.x, 2) + glm::pow(center.x, 2) - 2 * (position.x * center.x);
			float numberY = glm::pow(position.y, 2) + glm::pow(center.y, 2) - 2 * (position.y * center.y);
			float numberZ = glm::pow(position.z, 2) + glm::pow(center.z, 2) - 2 * (position.z * center.z);
			float number = numberX + numberY + numberZ - glm::pow(radius, 2);
			//2nd grade equations solved
			float det = (glm::pow(landa, 2) - (4 * (landaSquared)* number));
			if (det < 0)
			{
				printf("Arrel Negativa \n");
				system("pause");
			}
			float solutionLandaPlus = (-landa + glm::sqrt((det)))/ (2 * landaSquared);
			float solutionLandaMinus = (-landa - glm::sqrt((det))) / (2 * landaSquared);

			//Calculate the 2 point to see which one is the closest one
			glm::vec3 recPointPlus = position + (solutionLandaPlus * velocity);
			glm::vec3 recPointMinus = position + (solutionLandaMinus * velocity);
			if (glm::distance(position, recPointPlus) > glm::distance(position, recPointMinus))
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

			//Update the position to bounce from the surface
			position -= 2 * (glm::dot(normal, position) + D) * (normal);
			tempPos -= 2 * (glm::dot(normal, tempPos) + D) * (normal);
			velocity -= 2 * (glm::dot(normal, velocity)) * (normal);
		}
	};

	class Capsule
	{
		//Define the parts of the capsule
		enum capsuleCollisions { TOPSPHERE, BOTTOMSPHERE, CYLLINDER, NONE };

	public:
		//Define the 2 spheres (top and bottom) of the capsule
		Sphere topSemiSphere;
		Sphere bottomSemiSphere;

		Capsule()
		{
			topSemiSphere = Sphere();
			bottomSemiSphere = Sphere();
		}

		Capsule(Sphere top, Sphere bottom)
		{
			//Create a capsule with two spheres
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
			//Calculate the distances to see if we collide
			glm::vec3 AB = topSemiSphere.center - bottomSemiSphere.center;
			glm::vec3 cross = glm::cross((tempPos - bottomSemiSphere.center), AB);
			float ABmodule = (glm::sqrt(glm::pow(AB.x, 2) + glm::pow(AB.y, 2) + glm::pow(AB.z, 2)));
			float crossModule = glm::sqrt(glm::pow(cross.x, 2) + glm::pow(cross.y, 2) + glm::pow(cross.z, 2));
			float distToRect = crossModule / ABmodule;

			//OPTION B
			if (topSemiSphere.hasCollisioned(tempPos))			//Check collision with top sphere
			{
				collisionTop = true;
				distToTopSphere = glm::distance(tempPos, topSemiSphere.center);
			}
			if (bottomSemiSphere.hasCollisioned(tempPos))		//Check collision with bottom sphere
			{
				collisionBottom = true;
				distToBottomSphere = glm::distance(tempPos, bottomSemiSphere.center);
			}
			//If any distance is smaller than the radius, it's in collision
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

	#pragma region DEFINITION OF GEOMETRICAL FIGURES
	//Define the cube with 6 planes
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
	//Capsule
	Capsule capsule = Utils::Capsule(Utils::Sphere(glm::vec3(-3.f, 2.f, 0.f), 1.f), Utils::Sphere(glm::vec3(3.f, 2.f, 0.f), 1.f));
	//Sphere
	Sphere sphere = Utils::Sphere(glm::vec3(0.f, 1.f, 0.f), 1.f);
	//Cube
	glm::vec3 pointsPlane1[3] = { glm::vec3(-5,0,-5), glm::vec3(5,0,-5),glm::vec3(5,0,5) };  // down
	glm::vec3 pointsPlane2[3] = { glm::vec3(-5,0,5),  glm::vec3(-5,10,5), glm::vec3(-5,0,-5) };  // front
	glm::vec3 pointsPlane3[3] = { glm::vec3(-5,0,-5), glm::vec3(-5,10,-5),glm::vec3(5,0,-5) }; //right
	glm::vec3 pointsPlane4[3] = { glm::vec3(-5,0,5), glm::vec3(5,0,5),glm::vec3(5.f, 10.f,5.f) };//left
	glm::vec3 pointsPlane5[3] = { glm::vec3(5.f, 10.f,  5.f), glm::vec3(5,0,-5),glm::vec3(5,0,5) }; //back
	glm::vec3 pointsPlane6[3] = { glm::vec3(-5.f, 10.f, -5.f),glm::vec3(5.f, 10.f, -5.f),glm::vec3(5.f, 10.f,  5.f) };//up
#pragma endregion

	float moduleOfAVector(glm::vec3 vec)
	{
		return glm::sqrt(glm::pow(vec.x, 2) + glm::pow(vec.y, 2) + glm::pow(vec.z, 2));
	}
}

namespace ClothMesh {
	float timeSimulation = 0;
	extern const int numCols;
	extern const int numRows;
	const int numWaves = 2;
	float amplitude[2] = { 0.3f, 0.5f};
	float omega[2] = { 1 , 0.5f };
	glm::vec3 waveVector[2] = { glm::vec3(1, 0, 0) , glm::vec3(1, 0, 1) };
	extern void updateClothMesh(float* array_data);
	glm::vec3 clothPositions[18][14];
	glm::vec3 fToBeApplied[18][14];
	glm::vec3 **initialPos;
	glm::vec3 **velocity;
	float mass = 1;
	float originalSpringRowsLenght;
	float diagonalOriginalLenght;
	float originalSpringColsLenght; 
	float kElasticity = 1000.f;
	float kDamping = 50.f;

	void WaveUpdate(float dt)
	{
		
	}
	void resetMesh()
	{
		ClothMesh::diagonalOriginalLenght = glm::sqrt(glm::pow(ClothMesh::originalSpringRowsLenght, 2) + glm::pow(ClothMesh::originalSpringColsLenght, 2));
		float currentRowPos = -5;
		float currentColPos = -5;
		for (int i = 0; i < ClothMesh::numRows; i++)
		{
			for (int j = 0; j < ClothMesh::numCols; j++)
			{
				ClothMesh::clothPositions[i][j] = glm::vec3((float)currentRowPos, 5, (float)currentColPos);
				ClothMesh::initialPos[i][j] = glm::vec3((float)currentRowPos, 5, (float)currentColPos);
				ClothMesh::velocity[i][j] = glm::vec3(0.f, 0.f, 0.f);
				currentColPos += ClothMesh::originalSpringColsLenght;
			}
			currentRowPos += ClothMesh::originalSpringRowsLenght;
			currentColPos = -5;
		}
	}
	glm::vec3 calculusF(glm::vec3 clothPos, glm::vec3 nextPos, glm::vec3 velocity, glm::vec3 nextVelocity, float originalSpringLenght)
	{
		glm::vec3 retorno = glm::vec3(0.f, 0.f, 0.f);
		glm::vec3 vector = clothPos - nextPos;
		float mod = Utils::moduleOfAVector(vector);
		glm::vec3 normal = ((clothPos - nextPos) / mod);
		float Elasticity = (ClothMesh::kElasticity * (mod - originalSpringLenght));
		float dot = glm::dot<float>((velocity - nextVelocity), normal);
		float Damping = ClothMesh::kDamping * dot;
		if (!std::isnan(Damping) || !std::isnan(Elasticity))
		{
			retorno = glm::vec3(-1.f, -1.f, -1.f) * ((Elasticity + Damping) * normal);
		}
		else
		{
			nanCrashes++;
			retorno = glm::vec3(0.f, 0.f, 0.f);
		}
		return retorno;
	}
}

void Exemple_GUI()
{
	
}

void Exemple_PhysicsInit()
{
	renderSphere = true;
	renderCloth = true;
	Utils::cubePlaneCollision[0] = Utils::Plane(Utils::pointsPlane1[0], Utils::pointsPlane1[1], Utils::pointsPlane1[2]);
	Utils::cubePlaneCollision[1] = Utils::Plane(Utils::pointsPlane2[0], Utils::pointsPlane2[1], Utils::pointsPlane2[2]);
	Utils::cubePlaneCollision[2] = Utils::Plane(Utils::pointsPlane3[0], Utils::pointsPlane3[1], Utils::pointsPlane3[2]);
	Utils::cubePlaneCollision[3] = Utils::Plane(Utils::pointsPlane4[0], Utils::pointsPlane4[1], Utils::pointsPlane4[2]);
	Utils::cubePlaneCollision[4] = Utils::Plane(Utils::pointsPlane5[0], Utils::pointsPlane5[1], Utils::pointsPlane5[2]);
	Utils::cubePlaneCollision[5] = Utils::Plane(Utils::pointsPlane6[0], Utils::pointsPlane6[1], Utils::pointsPlane6[2]);
	Utils::sphere.center = glm::vec3((-4) + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (-2 - (-4)))), (0) + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (8 - (0)))), (-2) + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 - (-2)))));
	Utils::sphere.radius = (1.0f) + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 - (1.0f))));
	Sphere::updateSphere(Utils::sphere.center, Utils::sphere.radius);
	ClothMesh::initialPos = new glm::vec3*[ClothMesh::numRows];
	ClothMesh::velocity = new glm::vec3*[ClothMesh::numRows];
	for (int i = 0; i < ClothMesh::numRows; i++)
	{
		ClothMesh::initialPos[i] = new glm::vec3[ClothMesh::numCols];
		ClothMesh::velocity[i] = new glm::vec3[ClothMesh::numCols];
	}
	float incrementRows = 10.f / ClothMesh::numRows;
	float incrementCols = 10.f / ClothMesh::numCols;
	ClothMesh::originalSpringColsLenght = incrementCols;
	ClothMesh::originalSpringRowsLenght = incrementRows;
	ClothMesh::resetMesh();
}

void Exemple_PhysicsUpdate(float dt) 
{
	//Verlet
	glm::vec3 tempPos;
	glm::vec3 tempVel;
	glm::vec3 FTotal = glm::vec3(0.f, 0.f, 0.f);
	if (renderCloth)
	{
		ClothMesh::timeSimulation += dt;
		for (int i = 0; i < ClothMesh::numRows; i++)
		{
			for (int j = 0; j < ClothMesh::numCols; j++)
			{
				glm::vec3 tempPos = ClothMesh::initialPos[i][j];
				glm::vec3 aux(0,0,0);
				for (int k = 0; k < ClothMesh::numWaves; k++)
				{
					aux = (ClothMesh::waveVector[k] / ((2 * glm::pi<float>()) / glm::normalize(ClothMesh::waveVector[k]))) * ClothMesh::amplitude[k] * glm::sin(glm::dot(ClothMesh::waveVector[k], ClothMesh::clothPositions[i][j]) - ClothMesh::omega[k] * ClothMesh::timeSimulation);
					tempPos.x += aux.x;
					tempPos.z += aux.z;
					tempPos.y += ClothMesh::amplitude[k] * glm::cos(glm::dot(ClothMesh::waveVector[k], ClothMesh::initialPos[i][j]) - ClothMesh::omega[k] * ClothMesh::timeSimulation);
				}
				ClothMesh::clothPositions[i][j] = tempPos;
			}
		}
		ClothMesh::updateClothMesh(&(ClothMesh::clothPositions[0][0].x));
	}
}

void Exemple_PhysicsCleanup() {
}
