#include "SphereSystemSimulator.h"
#include <time.h>       /* time */

std::function<float(float)> SphereSystemSimulator::m_Kernels[5] = {
	[](float x) {return 1.0f; },              // Constant, m_iKernel = 0
	[](float x) {return 1.0f - x; },          // Linear, m_iKernel = 1, as given in the exercise Sheet, x = d/2r
	[](float x) {return (1.0f - x)*(1.0f - x); }, // Quadratic, m_iKernel = 2
	[](float x) {return 1.0f / (x)-1.0f; },     // Weak Electric Charge, m_iKernel = 3
	[](float x) {return 1.0f / (x*x) - 1.0f; },   // Electric Charge, m_iKernel = 4
};

// SphereSystemSimulator member functions
SphereSystemSimulator::SphereSystemSimulator()
{
	m_iKernel = 1;
	m_iTestCase = 0;
	m_iNumSpheres = 100;
	m_fRadius = 0.05;
	m_fMass = 1;
	m_fDamping = 5;
	m_fStiffness = 1000;
}

tuple<int, int, int> SphereSystemSimulator::intToTuple(int cell)
{
	int m = get<1>(gridSize);
	int k = get<2>(gridSize);
	return make_tuple(cell / m / k, (cell % (m * k)) % k, cell % k);
}

void SphereSystemSimulator::addRandomSystem() {
	default_random_engine generator;
	generator.seed(time(NULL));
	uniform_int_distribution<int> dist(0, 1);
	uniform_real_distribution<float> pos(-0.5 + m_fRadius, 0.5 - m_fRadius);
	uniform_real_distribution<float> length(3 * m_fRadius, 0.5);

	int fig = dist(generator);
	switch (fig) {
	case 0: {
		Vec3 posVec = Vec3(pos(generator), pos(generator), pos(generator));
		Sphere s1 = Sphere(posVec, Vec3(0, 0, 0), true);
		posVec = Vec3(pos(generator), pos(generator), pos(generator));
		Sphere s2 = Sphere(posVec, Vec3(0, 0, 0), true);
		posVec = Vec3(pos(generator), pos(generator), pos(generator));
		Sphere s3 = Sphere(posVec, Vec3(0, 0, 0), true);
		posVec = Vec3(pos(generator), pos(generator), pos(generator));
		Sphere s4 = Sphere(posVec, Vec3(0, 0, 0), true);
		posVec = Vec3(pos(generator), pos(generator), pos(generator));
		Sphere s5 = Sphere(posVec, Vec3(0, 0, 0), true);
		posVec = Vec3(pos(generator), pos(generator), pos(generator));
		Sphere s6 = Sphere(posVec, Vec3(0, 0, 0), true);
		int nsprings = SpringSystem.size();
		int nspheres = SphereSystem.size();
		SphereSystem.push_back(s1);
		SphereSystem.push_back(s2);
		SphereSystem.push_back(s3);
		SphereSystem.push_back(s4);
		SphereSystem.push_back(s5);
		SphereSystem.push_back(s6);

		SpringSystem.push_back(&SphereSystem.at(nspheres));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 1));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 2));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 3));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 4));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 5));

		Spring sp1 = Spring(nsprings, nsprings + 1, sqrt(s1.position.squaredDistanceTo(s2.position)));
		Spring sp2 = Spring(nsprings, nsprings + 2, sqrt(s1.position.squaredDistanceTo(s3.position)));
		Spring sp3 = Spring(nsprings, nsprings + 3, sqrt(s1.position.squaredDistanceTo(s4.position)));
		Spring sp4 = Spring(nsprings, nsprings + 4, sqrt(s1.position.squaredDistanceTo(s5.position)));
		Spring sp5 = Spring(nsprings + 1, nsprings + 2, sqrt(s2.position.squaredDistanceTo(s3.position)));
		Spring sp6 = Spring(nsprings + 2, nsprings + 3, sqrt(s3.position.squaredDistanceTo(s4.position)));

		Spring sp13 = Spring(nsprings + 1, nsprings + 3, sqrt(s2.position.squaredDistanceTo(s4.position)));
		Spring sp14 = Spring(nsprings + 2, nsprings + 4, sqrt(s3.position.squaredDistanceTo(s5.position)));

		Spring sp7 = Spring(nsprings + 3, nsprings + 4, sqrt(s4.position.squaredDistanceTo(s5.position)));
		Spring sp8 = Spring(nsprings + 1, nsprings + 4, sqrt(s2.position.squaredDistanceTo(s5.position)));
		Spring sp9 = Spring(nsprings + 5, nsprings + 1, sqrt(s6.position.squaredDistanceTo(s2.position)));
		Spring sp10 = Spring(nsprings + 5, nsprings + 2, sqrt(s6.position.squaredDistanceTo(s3.position)));
		Spring sp11 = Spring(nsprings + 5, nsprings + 3, sqrt(s6.position.squaredDistanceTo(s4.position)));
		Spring sp12 = Spring(nsprings + 5, nsprings + 4, sqrt(s6.position.squaredDistanceTo(s5.position)));

		springs.push_back(sp1);
		springs.push_back(sp2);
		springs.push_back(sp3);
		springs.push_back(sp4);
		springs.push_back(sp5);
		springs.push_back(sp6);
		springs.push_back(sp7);
		springs.push_back(sp8);
		springs.push_back(sp9);
		springs.push_back(sp10);
		springs.push_back(sp11);
		springs.push_back(sp12);
		springs.push_back(sp13);
		springs.push_back(sp14);
		break;
	}
	case 1: {
		float l = length(generator) / 2;
		uniform_real_distribution<float> posCube(-0.5 + m_fRadius + l, 0.5 - m_fRadius - l);
		Vec3 posVec = Vec3(posCube(generator), posCube(generator), posCube(generator));
		Sphere s1 = Sphere(Vec3(-l, -l, -l) + posVec, Vec3(0, 0, 0), true);
		Sphere s2 = Sphere(Vec3(l,  -l, -l) + posVec, Vec3(0, 0, 0), true);
		Sphere s3 = Sphere(Vec3(l, -l, l) + posVec, Vec3(0, 0, 0), true);
		Sphere s4 = Sphere(Vec3(-l, -l, l) + posVec, Vec3(0, 0, 0), true);
		Sphere s5 = Sphere(Vec3(-l, l, -l) + posVec, Vec3(0, 0, 0), true);
		Sphere s6 = Sphere(Vec3(l, l, -l) + posVec, Vec3(0, 0, 0), true);
		Sphere s7 = Sphere(Vec3(l, l, l) + posVec, Vec3(0, 0, 0), true);
		Sphere s8 = Sphere(Vec3(-l, l, l) + posVec, Vec3(0, 0, 0), true);

		int nsprings = SpringSystem.size();
		int nspheres = SphereSystem.size();
		SphereSystem.push_back(s1);
		SphereSystem.push_back(s2);
		SphereSystem.push_back(s3);
		SphereSystem.push_back(s4);
		SphereSystem.push_back(s5);
		SphereSystem.push_back(s6);
		SphereSystem.push_back(s7);
		SphereSystem.push_back(s8);

		SpringSystem.push_back(&SphereSystem.at(nspheres));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 1));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 2));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 3));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 4));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 5));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 6));
		SpringSystem.push_back(&SphereSystem.at(nspheres + 7));

		Spring sp1 = Spring(nsprings, nsprings + 1, sqrt(s1.position.squaredDistanceTo(s2.position)));
		Spring sp2 = Spring(nsprings + 1, nsprings + 2, sqrt(s2.position.squaredDistanceTo(s3.position)));
		Spring sp3 = Spring(nsprings + 2, nsprings + 3, sqrt(s3.position.squaredDistanceTo(s4.position)));
		Spring sp4 = Spring(nsprings + 3, nsprings, sqrt(s1.position.squaredDistanceTo(s4.position)));

		Spring sp5 = Spring(nsprings + 4, nsprings + 5, sqrt(s5.position.squaredDistanceTo(s6.position)));
		Spring sp6 = Spring(nsprings + 5, nsprings + 6, sqrt(s6.position.squaredDistanceTo(s7.position)));
		Spring sp7 = Spring(nsprings + 6, nsprings + 7, sqrt(s7.position.squaredDistanceTo(s8.position)));
		Spring sp8 = Spring(nsprings + 4, nsprings + 7, sqrt(s5.position.squaredDistanceTo(s8.position)));

		Spring sp9 = Spring(nsprings, nsprings + 4, sqrt(s1.position.squaredDistanceTo(s5.position)));
		Spring sp10 = Spring(nsprings + 1, nsprings + 5, sqrt(s2.position.squaredDistanceTo(s6.position)));
		Spring sp11 = Spring(nsprings + 2, nsprings + 6, sqrt(s3.position.squaredDistanceTo(s7.position)));
		Spring sp12 = Spring(nsprings + 3, nsprings + 7, sqrt(s4.position.squaredDistanceTo(s8.position)));

		Spring s1p1 = Spring(nsprings, nsprings + 2, sqrt(s1.position.squaredDistanceTo(s3.position)));
		Spring s1p2 = Spring(nsprings + 1, nsprings + 4, sqrt(s2.position.squaredDistanceTo(s4.position)));
		Spring s1p3 = Spring(nsprings + 4, nsprings + 6, sqrt(s5.position.squaredDistanceTo(s7.position)));
		Spring s1p4 = Spring(nsprings + 5, nsprings + 7, sqrt(s6.position.squaredDistanceTo(s8.position)));

		Spring s1p5 = Spring(nsprings, nsprings + 5, sqrt(s1.position.squaredDistanceTo(s6.position)));
		Spring s1p6 = Spring(nsprings + 1, nsprings + 4, sqrt(s2.position.squaredDistanceTo(s5.position)));
		Spring s1p7 = Spring(nsprings + 2, nsprings + 7, sqrt(s3.position.squaredDistanceTo(s8.position)));
		Spring s1p8 = Spring(nsprings + 3, nsprings + 6, sqrt(s4.position.squaredDistanceTo(s7.position)));

		Spring s1p9 = Spring(nsprings, nsprings + 7, sqrt(s1.position.squaredDistanceTo(s8.position)));
		Spring s1p10 = Spring(nsprings + 3, nsprings + 4, sqrt(s4.position.squaredDistanceTo(s5.position)));
		Spring s1p11 = Spring(nsprings + 1, nsprings + 6, sqrt(s2.position.squaredDistanceTo(s7.position)));
		Spring s1p12 = Spring(nsprings + 2, nsprings + 5, sqrt(s3.position.squaredDistanceTo(s6.position)));

		Spring s2p1 = Spring(nsprings, nsprings + 6, sqrt(s1.position.squaredDistanceTo(s7.position)));
		Spring s2p2 = Spring(nsprings + 1, nsprings + 7, sqrt(s2.position.squaredDistanceTo(s8.position)));
		Spring s2p3 = Spring(nsprings + 2, nsprings + 4, sqrt(s3.position.squaredDistanceTo(s5.position)));
		Spring s2p4 = Spring(nsprings + 3, nsprings + 5, sqrt(s4.position.squaredDistanceTo(s6.position)));

		springs.push_back(sp1);
		springs.push_back(sp2);
		springs.push_back(sp3);
		springs.push_back(sp4);
		springs.push_back(sp5);
		springs.push_back(sp6);
		springs.push_back(sp7);
		springs.push_back(sp8);
		springs.push_back(sp9);
		springs.push_back(sp10);
		springs.push_back(sp11);
		springs.push_back(sp12);

		springs.push_back(s1p1);
		springs.push_back(s1p2);
		springs.push_back(s1p3);
		springs.push_back(s1p4);
		springs.push_back(s1p5);
		springs.push_back(s1p6);
		springs.push_back(s1p7);
		springs.push_back(s1p8);
		springs.push_back(s1p9);
		springs.push_back(s1p10);
		springs.push_back(s1p11);
		springs.push_back(s1p12);

		springs.push_back(s2p1);
		springs.push_back(s2p2);
		springs.push_back(s2p3);
		springs.push_back(s2p4);
		break;
	}
	}
}


int SphereSystemSimulator::tupleToInt(tuple<int, int, int> cell) 
{
	int m = get<1>(gridSize);
	int k = get<2>(gridSize);
	return get<0>(cell) * m * k + get<1>(cell) * k + get<2>(cell);
}

const char * SphereSystemSimulator::getTestCasesStr()
{
	return "Project";
}

void SphereSystemSimulator::initUI(DrawingUtilitiesClass * DUC)
{
	this->DUC = DUC;
}

void SphereSystemSimulator::FillSpheres() 
{
	int sN = 0;
	for (int i = 0; i < sqrt(m_iNumSpheres * 1. / 2); i++) {
		for (int j = 0; j < sqrt(m_iNumSpheres * 1. / 2); j++) {
			Sphere sphere = Sphere(Vec3(-0.5 + m_fRadius + j * 2 * m_fRadius, 0.5 - m_fRadius, -0.5 + m_fRadius + i * 2 * m_fRadius), Vec3(0, 0, 0));
			Sphere sphere1 = Sphere(Vec3(-0.5 + m_fRadius + j * 2 * m_fRadius, 0.5 - m_fRadius * 3, -0.5 + m_fRadius + i * 2 * m_fRadius), Vec3(0, 0, 0));

			int x = (int)round((sphere.position.x + 0.5 - m_fRadius) / 2 / m_fRadius);
			int y = (int)round((sphere.position.y + 0.5 - m_fRadius) / 2 / m_fRadius);
			int z = (int)round((sphere.position.z + 0.5 - m_fRadius) / 2 / m_fRadius);
			sphere.oldCell = make_tuple(x, y, z);
			sphere.newCell = make_tuple(x, y, z);
			SphereSystem.push_back(sphere);
			SphereSystem1.push_back(sphere);
			sN++;
			if (sN == m_iNumSpheres)
				break;

			x = (int)round((sphere1.position.x + 0.5 - m_fRadius) / 2 / m_fRadius);
			y = (int)round((sphere1.position.y + 0.5 - m_fRadius) / 2 / m_fRadius);
			z = (int)round((sphere1.position.z + 0.5 - m_fRadius) / 2 / m_fRadius);
			sphere1.oldCell = make_tuple(x, y, z);
			sphere1.newCell = make_tuple(x, y, z);
			SphereSystem.push_back(sphere1);
			SphereSystem1.push_back(sphere1);
			sN++;
			if (sN == m_iNumSpheres)
				break;
		}
		if (sN == m_iNumSpheres)
			break;
	}
}

void SphereSystemSimulator::reset()
{
	m_mouse.x = m_mouse.y = 0;
	m_trackmouse.x = m_trackmouse.y = 0;
	m_oldtrackmouse.x = m_oldtrackmouse.y = 0;

	SphereSystem.clear();
	SphereSystem1.clear();
	SpringSystem.clear();
	springs.clear();
	grid.clear();

	gridSize = make_tuple((int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius));
	for (int i = 0; i < 1. / 2 / m_fRadius; i++)
		for (int j = 0; j < 1. / 2 / m_fRadius; j++)
			for (int k = 0; k < 1. / 2 / m_fRadius; k++) {
				set<Sphere*> list;
				grid.push_back(list);
			}

	FillSpheres();

	auto neededSystem = &SphereSystem;
	for (int i = 0; i < neededSystem->size(); i++) {
		Sphere* sphere = &neededSystem->at(i);
		int cell = tupleToInt(sphere->oldCell);
		grid.at(cell).insert(sphere);
	}

}

void SphereSystemSimulator::drawFrame(ID3D11DeviceContext * pd3dImmediateContext)
{
	for each (Spring spring in springs)
	{
		DUC->beginLine();
		DUC->drawLine(SpringSystem.at(spring.masspoint1)->position, Vec3(1, 1, 1), SpringSystem.at(spring.masspoint2)->position, Vec3(1, 1, 1));
		DUC->endLine();
	}

	for each (Sphere sphere in SphereSystem)
	{
		if (!sphere.isSpring)
			DUC->setUpLighting(Vec3(), 0.4*Vec3(1, 1, 1), 100, 0.6*Vec3(0.97, 0.86, 1));
		else
			DUC->setUpLighting(Vec3(), 0.4*Vec3(1, 1, 1), 100, 0.6*Vec3(0.97, 0, 0));
		DUC->drawSphere(sphere.position, m_fRadius);
	}
	
}

void SphereSystemSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	reset();
}

void SphereSystemSimulator::externalForcesCalculations(float timeElapsed)
{
	Point2D mouseDiff;
	mouseDiff.x = m_trackmouse.x - m_oldtrackmouse.x;
	mouseDiff.y = m_trackmouse.y - m_oldtrackmouse.y;
	if (mouseDiff.x != 0 || mouseDiff.y != 0)
	{
		Mat4 worldViewInv = Mat4(DUC->g_camera.GetWorldMatrix() * DUC->g_camera.GetViewMatrix());
		worldViewInv = worldViewInv.inverse();
		Vec3 inputView = Vec3((float)mouseDiff.x, (float)-mouseDiff.y, 0);

		Vec3 inputWorld = worldViewInv.transformVectorNormal(inputView);
		// find a proper scale!
		float inputScale = 3;
		inputWorld = inputWorld * inputScale;
		extForce = inputWorld;
	}
	else {
		if (extForce.X != 0 || extForce.Y != 0 || extForce.Z != 0) {
			externalForce = externalForce1 = extForce;
		}
		extForce = Vec3(0, 0, 0);
	}
}

void SphereSystemSimulator::simulateTimestep(float timeStep)
{
	if (m_iTestCase == 0)
		midpointGrid(timeStep);
}

void SphereSystemSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void SphereSystemSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

int SphereSystemSimulator::changeNumber(int n)
{
	int layer = floor(1. / 2 / m_fRadius + 0.001) * floor(1. / 2 / m_fRadius + 0.001);
	if (n > 2 * layer)
		m_iNumSpheres = 2 * layer;
	else
		m_iNumSpheres = n;
	reset();
	return m_iNumSpheres;
}

void SphereSystemSimulator::changeMass(float m)
{
	m_fMass = m;
}

int SphereSystemSimulator::changeRadius(float r)
{
	int layer = floor(1. / 2 / r + 0.001) * floor(1. / 2 / r + 0.001);
	if (m_iNumSpheres > 2 * layer)
		m_iNumSpheres = 2 * layer;
	m_fRadius = r;
	reset();
	return m_iNumSpheres;
}

Vec3 SphereSystemSimulator::calculateForces(int i) 
{
	Vec3 force = Vec3(0, 0, 0);
	function<float(float)> kernel;
	Sphere sphere1 = SphereSystem.at(i);
	for (int j = 0; j < SphereSystem.size(); j++)
		if (i != j) {
			Sphere sphere2 = SphereSystem.at(j);
			if (sphere2.isSpring || sphere1.isSpring)
				kernel = m_Kernels[m_iKernelSpring];
			else
				kernel = m_Kernels[m_iKernel];
			float distance = sqrt(sphere1.oldPosition.squaredDistanceTo(sphere2.oldPosition));
			if (distance < 2 * m_fRadius && distance != 0) {
				Vec3 direction = (sphere1.oldPosition - sphere2.oldPosition) / distance;
				force += direction * kernel(distance / 2 / m_fRadius) * 1500;
			}
		}

	return force;
}

Vec3 SphereSystemSimulator::calculateForcesGrid(Sphere sphere1)
{
	Vec3 force = Vec3(0, 0, 0);
	auto cell = sphere1.oldCell;
	int n = get<0>(gridSize);
	int m = get<1>(gridSize);
	int o = get<2>(gridSize);
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <=1; j++)
			for (int k = -1; k <= 1; k++) {
				int x = get<0>(cell) + i;
				int y = get<1>(cell) + j;
				int z = get<2>(cell) + k;
				auto nearCell = x * m * o + y * o + z;
				if (x >= 0 && x < n && y >= 0 && y < m && z >= 0 && z < o) {
					for each (Sphere* sphere2 in grid.at(nearCell)) {
						float distance = sqrt(sphere1.oldPosition.squaredDistanceTo(sphere2->oldPosition));
						if (distance < 2 * m_fRadius && distance != 0) {
							Vec3 direction = (sphere1.oldPosition - sphere2->oldPosition) / distance;
							force += direction * m_Kernels[m_iKernel](distance / 2 / m_fRadius) * 2000;
						}
					}
				}
			}
	return force;
}

Vec3 SphereSystemSimulator::calculateSpringForces(int index) {
	// this thing calculates forces for the point x_i, use it every time you need it
	Vec3 force = Vec3(0, 0, 0);
	for each (Spring spring in springs)
	{
		if (spring.masspoint1 == index) {
			float length = SpringSystem.at(spring.masspoint1)->oldPosition.squaredDistanceTo(SpringSystem.at(spring.masspoint2)->oldPosition);
			length = sqrt(length);
			force -= m_fStiffness * (length - spring.initialLength) / length * (SpringSystem.at(spring.masspoint1)->oldPosition - SpringSystem.at(spring.masspoint2)->oldPosition);
		}
		else if (spring.masspoint2 == index) {
			float length = SpringSystem.at(spring.masspoint1)->oldPosition.squaredDistanceTo(SpringSystem.at(spring.masspoint2)->oldPosition);
			length = sqrt(length);
			force += m_fStiffness * (length - spring.initialLength) / length * (SpringSystem.at(spring.masspoint1)->oldPosition - SpringSystem.at(spring.masspoint2)->oldPosition);
		}
	}
	return force;
}

void SphereSystemSimulator::midpointGrid(float timestep)
{
	auto neededSystem = &SphereSystem;
	default_random_engine generator;
	uniform_real_distribution<float> distribution(1.f, 1.5f);
	vector<Sphere*> changedSpheres;
	for (int i = 0; i < SpringSystem.size(); i++) {
		Sphere* sphere = SpringSystem.at(i);
		sphere->force += calculateSpringForces(i);
	}

	for (int i = 0; i < (*neededSystem).size(); i++) {
		Sphere* sphere = &(*neededSystem).at(i);
		Vec3 forces = calculateForcesGrid(*sphere) + externalForce1 * distribution(generator) + sphere->force;
		sphere->force = Vec3(0, 0, 0);
		Vec3 gravity = Vec3(0, 50, 0);
		Vec3 y_x2 = sphere->oldVelocity + (forces - gravity - sphere->oldVelocity * m_fDamping) / m_fMass * timestep / 2;
		Vec3 y1_x2 = sphere->oldVelocity + timestep * (forces - gravity - y_x2 * m_fDamping) / m_fMass;
		Vec3 y1_x1 = sphere->oldPosition + timestep * y1_x2;

		sphere->position = y1_x1;
		sphere->velocity = y1_x2;

		if (sphere->position.y < m_fRadius - 0.5) {
			sphere->position.y = m_fRadius - 0.5;
			sphere->velocity.y = 0;
		}
		if (sphere->position.y > 0.5 - m_fRadius) {
			sphere->position.y = 0.5 - m_fRadius;
			sphere->velocity.y = 0;
		}

		if (sphere->position.x < m_fRadius - 0.5) {
			sphere->position.x = m_fRadius - 0.5;
			sphere->velocity.x = 0;
		}
		if (sphere->position.x > 0.5 - m_fRadius) {
			sphere->position.x = 0.5 - m_fRadius;
			sphere->velocity.x = 0;
		}
		
		if (sphere->position.z < m_fRadius - 0.5) {
			sphere->position.z = m_fRadius - 0.5;
			sphere->velocity.z = 0;
		}
		if (sphere->position.z > 0.5 - m_fRadius) {
			sphere->position.z = 0.5 - m_fRadius;
			sphere->velocity.z = 0;
		}

		int x = (int)round((sphere->position.x + 0.5 - m_fRadius) / 2 / m_fRadius);
		int y = (int)round((sphere->position.y + 0.5 - m_fRadius) / 2 / m_fRadius);
		int z = (int)round((sphere->position.z + 0.5 - m_fRadius) / 2 / m_fRadius);

		if (make_tuple(x, y, z) != sphere->oldCell) {
			sphere->newCell = make_tuple(x, y, z);
			changedSpheres.push_back(sphere);
		}
	}
	for (int i = 0; i < changedSpheres.size(); i++) {
		Sphere* sphere = changedSpheres.at(i);
		grid.at(tupleToInt(sphere->oldCell)).erase(sphere);
		grid.at(tupleToInt(sphere->newCell)).insert(sphere);
		sphere->oldCell = sphere->newCell;
	}

	for (int i = 0; i < neededSystem->size(); i++) {
		Sphere* sphere = &neededSystem->at(i);
		sphere->oldPosition = sphere->position;
		sphere->oldVelocity = sphere->velocity;
	}
	externalForce1 = Vec3(0, 0, 0);
}
