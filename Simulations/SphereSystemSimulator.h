#ifndef SPHSYSTEMSIMULATOR_h
#define SPHSYSTEMSIMULATOR_h
#include "Simulator.h"
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <random>
#endif
//#include "spheresystem.h", add your sphere system header file

#define NAIVEACC 0
#define GRIDACC 1

using namespace std;

struct Sphere {
	Sphere(Vec3 p, Vec3 v) {
		position = p;
		oldPosition = p;
		velocity = v;
		oldVelocity = v;
		oldCell = tuple<int, int, int>(0, 0, 0);
		newCell = tuple<int, int, int>(0, 0, 0);
		force = Vec3(0, 0, 0);
		isSpring = false;
	}
	Sphere(Vec3 p, Vec3 v, bool iS) {
		position = p;
		oldPosition = p;
		velocity = v;
		oldVelocity = v;
		oldCell = tuple<int, int, int>(0, 0, 0);
		newCell = tuple<int, int, int>(0, 0, 0);
		force = Vec3(0, 0, 0);
		isSpring = iS;
	}
	Vec3 position;
	Vec3 oldPosition;
	Vec3 velocity;
	Vec3 oldVelocity;
	tuple<int, int, int> oldCell;
	tuple<int, int, int> newCell;
	Vec3 force;
	bool isSpring;
};

struct Spring {
	Spring(int mp1, int mp2, float iL) {
		masspoint1 = mp1;
		masspoint2 = mp2;
		initialLength = iL;
	}
	int masspoint1;
	int masspoint2;
	float initialLength;
};

class SphereSystemSimulator:public Simulator{
public:
	// Construtors
	SphereSystemSimulator();
	// Functions
	tuple<int, int, int> intToTuple(int cell);
	int tupleToInt(tuple<int, int, int> cell);
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void externalForcesCalculations(float timeElapsed);
	void simulateTimestep(float timeStep);
	void onClick(int x, int y);
	void onMouse(int x, int y);
	int changeNumber(int n);
	void changeMass(float m);
	int changeRadius(float r);
	void FillSpheres();
	Vec3 calculateForces(int i);
	void addRandomSystem();
	Vec3 calculateForcesGrid(Sphere sphere1);

	Vec3 calculateSpringForces(int index);

	void midpointGrid(float timestep);
	
protected:
	// Attributes
	Vec3 externalForce;
	Vec3 externalForce1;
	Vec3 extForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
	float m_fMass;
	float m_fRadius;
	float m_fStiffness;
	float m_fForceScaling;
	float m_fDamping;
	int m_iNumSpheres;
	
	int m_iKernel; // index of the m_Kernels[5], more detials in SphereSystemSimulator.cpp
	int m_iKernelSpring;
	static std::function<float(float)> m_Kernels[5];
	
	int m_iAccelerator; // switch between NAIVEACC and GRIDACC, (optionally, KDTREEACC, 2)
	
	vector<Sphere> SphereSystem;
	vector<Sphere> SphereSystem1;
	tuple<int, int, int> gridSize;
	vector<set<Sphere*>> grid;

	vector<Sphere*> SpringSystem;
	vector<Spring> springs;
	// for Demo 3 only:
	// you will need multiple SphereSystem objects to do comparisons in Demo 3
	// m_iAccelerator should be ignored.
	// SphereSystem * m_pSphereSystemGrid; 

};