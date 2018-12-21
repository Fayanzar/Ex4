#include "SphereSystemSimulator.h"

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
	m_fDamping = 25;
}

tuple<int, int, int> SphereSystemSimulator::intToTuple(int cell)
{
	int m = get<1>(gridSize);
	int k = get<2>(gridSize);
	return make_tuple(cell / m / k, (cell % (m * k)) % k, cell % k);
}

int SphereSystemSimulator::tupleToInt(tuple<int, int, int> cell) 
{
	int m = get<1>(gridSize);
	int k = get<2>(gridSize);
	return get<0>(cell) * m * k + get<1>(cell) * k + get<2>(cell);
}

const char * SphereSystemSimulator::getTestCasesStr()
{
	return "Demo1,Demo2,Demo3,DemoPerf";
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
	grid.clear();

	gridSize = make_tuple((int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius));
	for (int i = 0; i < 1. / 2 / m_fRadius; i++)
		for (int j = 0; j < 1. / 2 / m_fRadius; j++)
			for (int k = 0; k < 1. / 2 / m_fRadius; k++) {
				set<Sphere*> list;
				grid.push_back(list);
			}

	FillSpheres();

	auto neededSystem = m_iTestCase == 2 ? &SphereSystem1 : &SphereSystem;
	for (int i = 0; i < neededSystem->size(); i++) {
		Sphere* sphere = &neededSystem->at(i);
		int cell = tupleToInt(sphere->oldCell);
		grid.at(cell).insert(sphere);
	}

	if (m_iTestCase == 3) {
		MuTime myTimer;
		for (int k = 0; k < 3; k++) {
			switch (k) {
			case 1: 
				m_fRadius = 0.02;
				break;
			case 2:
				m_fRadius = 0.005;
				break;
			}
			SphereSystem.clear();
			SphereSystem1.clear();
			grid.clear();

			gridSize = make_tuple((int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius));
			for (int i = 0; i < 1. / 2 / m_fRadius; i++)
				for (int j = 0; j < 1. / 2 / m_fRadius; j++)
					for (int k = 0; k < 1. / 2 / m_fRadius; k++) {
						set<Sphere*> list;
						grid.push_back(list);
					}
			FillSpheres();
			cout << "n = " << m_iNumSpheres << endl;
			cout << "r = " << m_fRadius << endl;
			myTimer.get();
			for (int i = 0; i < 10; i++)
				midpoint(0.001);
			cout << "naïve: " << myTimer.update().time << endl;

			SphereSystem.clear();
			SphereSystem1.clear();
			grid.clear();

			gridSize = make_tuple((int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius), (int)ceil(1. / 2 / m_fRadius));
			for (int i = 0; i < 1. / 2 / m_fRadius; i++)
				for (int j = 0; j < 1. / 2 / m_fRadius; j++)
					for (int k = 0; k < 1. / 2 / m_fRadius; k++) {
						set<Sphere*> list;
						grid.push_back(list);
					}

			FillSpheres();
			myTimer.get();
			for (int i = 0; i < 10; i++)
				midpointGrid(0.001);
			cout << "accelerated: " << myTimer.update().time << endl;
			m_iNumSpheres *= 10;
		}
	}
}

void SphereSystemSimulator::drawFrame(ID3D11DeviceContext * pd3dImmediateContext)
{
	if (m_iTestCase != 2 || showWhite) {
		DUC->setUpLighting(Vec3(), 0.4*Vec3(1, 1, 1), 100, 0.6*Vec3(0.97, 0.86, 1));
		for each (Sphere sphere in SphereSystem)
		{
			DUC->drawSphere(sphere.position, m_fRadius);
		}
	}
	if (m_iTestCase == 2 && showRed) {
		DUC->setUpLighting(Vec3(), 0.4*Vec3(1, 1, 1), 100, 0.6*Vec3(1, 0, 0));
		for (Sphere sphere : SphereSystem1) {
			DUC->drawSphere(sphere.position, m_fRadius);
		}
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
	if (m_iTestCase == 1)
		midpointGrid(timeStep);
	if (m_iTestCase == 0)
		midpoint(timeStep);
	if (m_iTestCase == 2) {
		midpoint(timeStep);
		midpointGrid(timeStep);
	}
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
	r = r <= 0.25 ? r : 0.25;
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
	Sphere sphere1 = SphereSystem.at(i);
	for (int j = 0; j < SphereSystem.size(); j++)
		if (i != j) {
			Sphere sphere2 = SphereSystem.at(j);
			float distance = sqrt(sphere1.oldPosition.squaredDistanceTo(sphere2.oldPosition));
			if (distance < 2 * m_fRadius && distance != 0) {
				Vec3 direction = (sphere1.oldPosition - sphere2.oldPosition) / distance;
				force += direction * m_Kernels[m_iKernel](distance / 2 / m_fRadius) * 2000;
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

void SphereSystemSimulator::midpointGrid(float timestep)
{
	auto neededSystem = m_iTestCase == 2 ? &SphereSystem1 : &SphereSystem;
	default_random_engine generator;
	uniform_real_distribution<float> distribution(1.f, 1.5f);
	vector<Sphere*> changedSpheres;
	for (int i = 0; i < (*neededSystem).size(); i++) {
		Sphere* sphere = &(*neededSystem).at(i);
		Vec3 forces = calculateForcesGrid(*sphere) + externalForce1 * distribution(generator);
		Vec3 gravity = Vec3(0, 100, 0);
		Vec3 y_x2 = sphere->oldVelocity + (forces - gravity - sphere->oldVelocity * m_fDamping) / m_fMass * timestep / 2;
		Vec3 y1_x2 = sphere->oldVelocity + timestep * (forces - gravity - y_x2 * m_fDamping) / m_fMass;
		Vec3 y1_x1 = sphere->oldPosition + timestep * y1_x2;

		sphere->position = y1_x1;
		if (sphere->position.y < m_fRadius - 0.5)
			sphere->position.y = m_fRadius - 0.5;
		if (sphere->position.y > 0.5 - m_fRadius)
			sphere->position.y = 0.5 - m_fRadius;

		if (sphere->position.x < m_fRadius - 0.5)
			sphere->position.x = m_fRadius - 0.5;
		if (sphere->position.x > 0.5 - m_fRadius)
			sphere->position.x = 0.5 - m_fRadius;

		if (sphere->position.z < m_fRadius - 0.5)
			sphere->position.z = m_fRadius - 0.5;
		if (sphere->position.z > 0.5 - m_fRadius)
			sphere->position.z = 0.5 - m_fRadius;

		int x = (int)round((sphere->position.x + 0.5 - m_fRadius) / 2 / m_fRadius);
		int y = (int)round((sphere->position.y + 0.5 - m_fRadius) / 2 / m_fRadius);
		int z = (int)round((sphere->position.z + 0.5 - m_fRadius) / 2 / m_fRadius);

		if (make_tuple(x, y, z) != sphere->oldCell) {
			sphere->newCell = make_tuple(x, y, z);
			changedSpheres.push_back(sphere);
		}
		sphere->velocity = y1_x2;
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

void SphereSystemSimulator::midpoint(float timestep) {
	default_random_engine generator;
	uniform_real_distribution<float> distribution(1.f, 1.5f);
	for (int i = 0; i < SphereSystem.size(); i++) {
		Vec3 forces = calculateForces(i) + externalForce * distribution(generator);
		Sphere* sphere = &SphereSystem.at(i);
		Vec3 gravity = Vec3(0, 100, 0);
		Vec3 y_x2 = sphere->oldVelocity + (forces - gravity - sphere->oldVelocity * m_fDamping) / m_fMass * timestep / 2;
		Vec3 y1_x2 = sphere->oldVelocity + timestep * (forces - gravity - y_x2 * m_fDamping) / m_fMass;
		Vec3 y1_x1 = sphere->oldPosition + timestep * y1_x2;

		sphere->position = y1_x1;
		if (sphere->position.y < m_fRadius - 0.5)
			sphere->position.y = m_fRadius - 0.5;
		if (sphere->position.y > 0.5 - m_fRadius)
			sphere->position.y = 0.5 - m_fRadius;

		if (sphere->position.x < m_fRadius - 0.5)
			sphere->position.x = m_fRadius - 0.5;
		if (sphere->position.x > 0.5 - m_fRadius)
			sphere->position.x = 0.5 - m_fRadius;

		if (sphere->position.z < m_fRadius - 0.5)
			sphere->position.z = m_fRadius - 0.5;
		if (sphere->position.z > 0.5 - m_fRadius)
			sphere->position.z = 0.5 - m_fRadius;
		sphere->velocity = y1_x2;
	}
	for (int i = 0; i < SphereSystem.size(); i++) {
		Sphere* sphere = &SphereSystem.at(i);
		sphere->oldPosition = sphere->position;
		sphere->oldVelocity = sphere->velocity;
	}
	externalForce = Vec3(0, 0, 0);

}
