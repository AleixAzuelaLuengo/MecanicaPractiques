#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <cstdio>
#include <stdlib.h>
#include <random>
#include <vector>
#include <ctime>
#include <iostream>

/////////Forward declarations
extern bool renderSphere;
extern bool renderCapsule;
extern bool renderParticles;
extern bool renderCloth;
extern bool renderCube;
int SIMULATION_TIME = 15;
std::clock_t start;
double duration;

namespace Sphere {
	extern void updateSphere(glm::vec3 pos, float radius = 1.f);
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
	float originalSpringColsLenght; 
	float fluidDesnity = 1.f;
	int positionX;
	int positionY;
	void resetMesh()
	{
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
}

namespace Utils
{
	
	bool timer()
	{
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		return (duration >= SIMULATION_TIME);
	}

	class Sphere
	{
	public:
		glm::vec3 velocity = glm::vec3(0, 0, 0);
		glm::vec3 center;

		float radius;
		float mass = 1;
		float friction = 1;
		float acum = 0;
		float acum2 = 0;
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
		void reset_sphere()
		{
			center = glm::vec3(std::rand() % 8 - 4, 9, std::rand() % 8 - 4);
			radius = std::rand() % 1 + 1;
			mass = std::rand() % 3 + 1;
			velocity = glm::vec3(0, 0, 0);
			friction = 1;
			acum = 0;
			acum2 = 0;

			duration = 0;
			start = clock();
		}
		
	};

	#pragma region DEFINITION OF GEOMETRICAL FIGURES
	//Define the cube with 6 planes
	Sphere sphere = Utils::Sphere(glm::vec3(0.f, 1.f, 0.f), 1.f);

}


void Exemple_GUI()
{
	ImGui::InputFloat("Sphere Mass", &Utils::sphere.mass);
	if (ImGui::Button("Reset sphere", ImVec2(100, 30)))
	{
		Utils::sphere.reset_sphere();
	}

}

void Exemple_PhysicsInit()
{
	std::srand(time(NULL));
	renderSphere = true;
	renderCloth = true;
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
	Utils::sphere.reset_sphere();
	Sphere::updateSphere(Utils::sphere.center, Utils::sphere.radius);

	glm::vec3 tempPos = Utils::sphere.center;
	tempPos.y = 5;
	std::vector<std::pair<glm::vec3, glm::vec2>>clothP;
	for (int i = 0; i < ClothMesh::numRows; i++)
	{
		for (int j = 0; j < ClothMesh::numCols; j++)
		{
			if ((glm::abs(ClothMesh::clothPositions[i][j].x - tempPos.x) < Utils::sphere.radius) && (glm::abs(ClothMesh::clothPositions[i][j].z - tempPos.z) < Utils::sphere.radius))
			{
				clothP.push_back(std::pair<glm::vec3, glm::vec2>(glm::vec3(ClothMesh::clothPositions[i][j]), glm::vec2(i, j)));
			}
		}
	}
	glm::vec2 minorPos = clothP[0].second;
	float minorVec = glm::distance(clothP[0].first, tempPos);
	for (int i = 1; i < clothP.size(); i++)
	{
		float temp = glm::distance(clothP[i].first, tempPos);
		if (minorVec > temp)
		{
			minorVec = temp;
			minorPos = clothP[i].second;
		}
	}
	ClothMesh::positionX = minorPos.x;
	ClothMesh::positionY = minorPos.y;
}

void Exemple_PhysicsUpdate(float dt) 
{
	//Verlet
	glm::vec3 tempPos;
	glm::vec3 spherePos;
	glm::vec3 tempVel;
	glm::vec3 buoyancy;
	glm::vec3 FTotal = glm::vec3(0.f, 0.f, 0.f);
	if (renderCloth)
	{
		if (Utils::timer())
		{
			Utils::sphere.reset_sphere();
			Sphere::updateSphere(Utils::sphere.center, Utils::sphere.radius);
			glm::vec3 tempPos = Utils::sphere.center;
			tempPos.y = 5;
			std::vector<std::pair<glm::vec3, glm::vec2>>clothP;
			for (int i = 0; i < ClothMesh::numRows; i++)
			{
				for (int j = 0; j < ClothMesh::numCols; j++)
				{
					if ((glm::abs(ClothMesh::clothPositions[i][j].x - tempPos.x) < Utils::sphere.radius) && (glm::abs(ClothMesh::clothPositions[i][j].z - tempPos.z) < Utils::sphere.radius))
					{
						clothP.push_back(std::pair<glm::vec3, glm::vec2>(glm::vec3(ClothMesh::clothPositions[i][j]), glm::vec2(i, j)));
					}
				}
			}
			glm::vec2 minorPos = clothP[0].second;
			glm::vec3 minorVec = clothP[0].first - Utils::sphere.center;
			for (int i = 1; i < clothP.size(); i++)
			{
				glm::vec3 temp = clothP[i].first - Utils::sphere.center;
				if (minorVec.length > temp.length)
				{
					minorVec = temp;
					minorPos = clothP[i].second;
				}
			}
			ClothMesh::positionX = minorPos.x;
			ClothMesh::positionY = minorPos.y;
			start = std::clock();
			duration = 0;
		}

		#pragma region WAVE UPDATE
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
					tempPos.x -= aux.x;
					tempPos.z -= aux.z;
					tempPos.y += ClothMesh::amplitude[k] * glm::cos(glm::dot(ClothMesh::waveVector[k], ClothMesh::initialPos[i][j]) - ClothMesh::omega[k] * ClothMesh::timeSimulation);
				}
				ClothMesh::clothPositions[i][j] = tempPos;
			}
		}
		#pragma endregion

		#pragma region BUOYANCY
		//BUOYANCY
		float height;
		float intensityBuoyancy;
		//in the air
		if (Utils::sphere.center.y - ClothMesh::clothPositions[ClothMesh::positionX][ClothMesh::positionY].y >= Utils::sphere.radius)
		{
			intensityBuoyancy = 0;
			Utils::sphere.friction = 1;
			printf("---\n");



		}
		else if(ClothMesh::clothPositions[ClothMesh::positionX][ClothMesh::positionY].y - Utils::sphere.center.y  >= Utils::sphere.radius)
		{
			intensityBuoyancy = 1;


			if (Utils::sphere.velocity.y > 0)
			{
				//going up in the water
				Utils::sphere.friction = 0.3;
			}
			else
			{
				//going down in the water
				Utils::sphere.friction = 0.6;
			}
		}
		else
		{
			if (Utils::sphere.velocity.y > 0)
			{
				//going up in the water
				Utils::sphere.friction = 0.3 - Utils::sphere.acum;
				if (0.3 - Utils::sphere.acum>0)
					Utils::sphere.acum += 0.001;

			}
			else
			{
				//going down in the water
				Utils::sphere.friction = 0.6 - Utils::sphere.acum2;
				if (0.3 - Utils::sphere.acum2 > 0)
					Utils::sphere.acum2 += 0.001;

			}

			float heightDiference = 0;
			for (int k = 0; k < ClothMesh::numWaves; k++)
			{
				heightDiference += ClothMesh::amplitude[k] * glm::cos(glm::dot(ClothMesh::waveVector[k], Utils::sphere.center));
			}
			height = ClothMesh::clothPositions[ClothMesh::positionX][ClothMesh::positionY].y + heightDiference;
			intensityBuoyancy = (2 * Utils::sphere.radius - (Utils::sphere.center.y + Utils::sphere.radius - height)) / (2 * Utils::sphere.radius);
		}

		buoyancy = ClothMesh::fluidDesnity * (glm::vec3(0, 9.81f, 0)) * ((ClothMesh::initialPos[ClothMesh::positionX][ClothMesh::positionY] - tempPos.y) * intensityBuoyancy);
		tempVel += dt * ((glm::vec3(0, -9.81f, 0)) + buoyancy / Utils::sphere.mass);
		Utils::sphere.velocity += tempVel* Utils::sphere.friction;
		Utils::sphere.center += Utils::sphere.velocity * dt;
		#pragma endregion

		Sphere::updateSphere(Utils::sphere.center, Utils::sphere.radius);
		ClothMesh::updateClothMesh(&(ClothMesh::clothPositions[0][0].x));
	}
}

void Exemple_PhysicsCleanup() {
}
