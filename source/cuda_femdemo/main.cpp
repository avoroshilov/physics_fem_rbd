#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include <stdio.h>		// Simple I/O
#include <math.h>		// Simple math Routines
#include <windows.h>	// Header for Windows
#include <vector>
#include <algorithm>
#include <functional>

#include "helpers/use_cuda.h"

#include "all.h"
#include "drawer/gldrawer.h"
#include "drawer/camera.h"
#include "math/m_Quaternion.h"
#include "coldet/gjk.h"
#include "coldet/simplealgs.h"
#include "helpers/globalallocator.h"
#include "dynamics/joints.h"
#include "dynamics/solver_base.h"
#include "dynamics/solver_pgs.h"
#include "dynamics/solver_cg.h"
#include "dynamics/solver_lcpcg.h"
#ifdef USE_CUDA
#	include "cuda_solver/solver_lcpcg_cuda.h"
#	include "cuda_solver/solver_jacobi_cuda.h"
#	pragma comment (lib, "cuda_solver.lib")
#endif
#include "FEM/FEM.h"
#include "FEM/joints.h"
#include "helpers/timer.h"


bool UpdateFEM = false;
bool gbFreeze = false;
bool gbStep = false;

/*
Defines description

// Defines demo-scene
//////////////////////////////////////////////////////////////////////////
//		1: FEM
//		2: Simple ball-sockets
//		3: Large-scale ball-sockets
//		4: Chandeliare
//////////////////////////////////////////////////////////////////////////
#define DEMO_FEMCUBE				1

// What object should we simulate
#define DEMO_FEMCUBE_MODEL			0

// Is FEM-Object in test scene fixed by the left side
#define DEMOCUBE_FIXED				0
// Is Ball-Sockets chain attached to FEM-Object, could be 0, 1 or 2
#define DEMO_FEMCUBE_NUMCHAINS		2
// Should stress testing longbox be generated
#define DEMO_FEMCUBE_LONGBOX		0
// Is Plane constraint acting on FEM-Object
#define DEMO_FEMCUBE_PLANE			1

*/

#define TEST_CUSTOM					0
#define TEST_FEMBOX_CHAIN			1
#define TEST_FEMBOX_2CHAINS			2
#define TEST_FEMLONGBOX_FIXED		3
#define TEST_FEMBUNNY_HEAVY			4
#define TEST_FEMBUNNY_HEAVY_CHAIN	5
#define TEST_PRESET					TEST_FEMBOX_2CHAINS

#if (TEST_PRESET == TEST_FEMBOX_CHAIN)

	// FEM box attached to a single chain, falling on a plane with friction
#	define DEMO_FEMCUBE				1
#	define DEMO_FEMCUBE_MODEL		0
#	define DEMOCUBE_FIXED			0
#	define DEMO_FEMCUBE_NUMCHAINS	1
#	define DEMO_FEMCUBE_LONGBOX		0
#	define DEMO_FEMCUBE_PLANE		1

#elif (TEST_PRESET == TEST_FEMBOX_2CHAINS)

	// FEM box attached to 2 chains
#	define DEMO_FEMCUBE				1
#	define DEMO_FEMCUBE_MODEL		0
#	define DEMOCUBE_FIXED			0
#	define DEMO_FEMCUBE_NUMCHAINS	2
#	define DEMO_FEMCUBE_LONGBOX		0
#	define DEMO_FEMCUBE_PLANE		1

#elif (TEST_PRESET == TEST_FEMLONGBOX_FIXED)

// FEM box attached to 2 chains
#	define DEMO_FEMCUBE				1
#	define DEMO_FEMCUBE_MODEL		0
#	define DEMOCUBE_FIXED			1
#	define DEMO_FEMCUBE_NUMCHAINS	0
#	define DEMO_FEMCUBE_LONGBOX		1
#	define DEMO_FEMCUBE_PLANE		0

#elif (TEST_PRESET == TEST_FEMBUNNY_HEAVY)

	// Heavy FEM bunny
#	define DEMO_FEMCUBE				1
#	define DEMO_FEMCUBE_MODEL		1
#	define DEMOCUBE_FIXED			0
#	define DEMO_FEMCUBE_NUMCHAINS	0
#	define DEMO_FEMCUBE_LONGBOX		0
#	define DEMO_FEMCUBE_PLANE		1

#elif (TEST_PRESET == TEST_FEMBUNNY_HEAVY_CHAIN)

	// Heavy FEM bunny, attached to a single chain
#	define DEMO_FEMCUBE				1
#	define DEMO_FEMCUBE_MODEL		1
#	define DEMOCUBE_FIXED			0
#	define DEMO_FEMCUBE_NUMCHAINS	1
#	define DEMO_FEMCUBE_LONGBOX		0
#	define DEMO_FEMCUBE_PLANE		1

#else

	// Custom
#	define DEMO_FEMCUBE				3
#	define DEMO_FEMCUBE_MODEL		0
#	define DEMOCUBE_FIXED			0
#	define DEMO_FEMCUBE_NUMCHAINS	0
#	define DEMO_FEMCUBE_LONGBOX		0
#	define DEMO_FEMCUBE_PLANE		0

#endif

// Colorizes FE-strain if enables, otherwise - lambda vector magnitude (CPU only)
//	WARNING: requires tweaking based on the material parameters to get proper visualization (look for FEstrain)
//			 currently tweaked for longbox test
#define COLORING_STRAIN		1

#define HW_SOLVER_NONE		0
#define HW_SOLVER_LCPCG		1
#define HW_SOLVER_JACOBI	2
#define HARDWARE_SOLVER		HW_SOLVER_NONE

#if ((!defined(USE_CUDA)) && (HARDWARE_SOLVER != HW_SOLVER_NONE))
#error Using hardware solver without CUDA (check helpers/use_cuda.h)
#endif

#define SW_SOLVER_LCPCG		0
#define SW_SOLVER_PGS		1
#define SW_SOLVER_CG		2
#define SOFTWARE_SOLVER		SW_SOLVER_LCPCG

#define LIGHT_RENDER 0

#define SOLVER_GOLD 0

#if (DEMO_FEMCUBE == 4)
namespace ChandelierDemo
{
	BallSocket SimpleSocket;
	Vec3 Anchor;
	float JointERP = 0.6f, JointCFM = 0.005f;
	float invMass = 1.0f / 1.0f;
	CMatrix3 invInertiaMatrix;
};

template <class SolverType>
void ChandeliereAddLevel(int depth, int baseNodePosIdx, float step, SolverType & FEMTestSolver, std::vector<BallSocket, GlobalAllocator<BallSocket>> & ballsockets)
{
	if (depth > 4)
		return;

	step /= 1.5f;
	
	Vec3 baseNodePos = Vec3(FEMTestSolver.m_NodePosRot_x[baseNodePosIdx],
							FEMTestSolver.m_NodePosRot_y[baseNodePosIdx],
							FEMTestSolver.m_NodePosRot_z[baseNodePosIdx]);
	int baseNodeRotIdx = baseNodePosIdx + 1;

	float stepY = -step;//  2.0f;
	Vec3 nx(-step, stepY, 0.0f), px( step, stepY, 0.0f),
		 nz(0.0f, stepY, -step), pz(0.0f, stepY,  step);

#define ADD_JOINT(nodeIndex) \
	{ \
		int NodeIdx1 = baseNodePosIdx; \
		int NodeIdx2 = nodeIndex; \
		\
		ChandelierDemo::Anchor.x = FEMTestSolver.m_NodePosRot_x[NodeIdx2]; \
		ChandelierDemo::Anchor.y = FEMTestSolver.m_NodePosRot_y[NodeIdx1]; \
		ChandelierDemo::Anchor.z = FEMTestSolver.m_NodePosRot_z[NodeIdx2]; \
		\
		ChandelierDemo::SimpleSocket.Init(FEMTestSolver, ChandelierDemo::JointERP, ChandelierDemo::JointCFM, ChandelierDemo::Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1); \
		ballsockets.push_back(ChandelierDemo::SimpleSocket); \
	}

	int Node0Idx = FEMTestSolver.m_NumNodes;
	FEMTestSolver.AddTranslationalNode(	ChandelierDemo::invMass,
										baseNodePos + nx,
										Vec3(0.0f, 0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ChandelierDemo::invInertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
	ADD_JOINT(Node0Idx);

	int Node1Idx = FEMTestSolver.m_NumNodes;
	FEMTestSolver.AddTranslationalNode(	ChandelierDemo::invMass,
										baseNodePos + px,
										Vec3(0.0f, 0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ChandelierDemo::invInertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
	ADD_JOINT(Node1Idx);

	int Node2Idx = FEMTestSolver.m_NumNodes;
	FEMTestSolver.AddTranslationalNode(	ChandelierDemo::invMass,
										baseNodePos + nz,
										Vec3(0.0f, 0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ChandelierDemo::invInertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
	ADD_JOINT(Node2Idx);

	int Node3Idx = FEMTestSolver.m_NumNodes;
	FEMTestSolver.AddTranslationalNode(	ChandelierDemo::invMass,
										baseNodePos + pz,
										Vec3(0.0f, 0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ChandelierDemo::invInertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
	ADD_JOINT(Node3Idx);

	ChandeliereAddLevel(depth + 1, Node0Idx, step, FEMTestSolver, ballsockets);
	ChandeliereAddLevel(depth + 1, Node1Idx, step, FEMTestSolver, ballsockets);
	ChandeliereAddLevel(depth + 1, Node2Idx, step, FEMTestSolver, ballsockets);
	ChandeliereAddLevel(depth + 1, Node3Idx, step, FEMTestSolver, ballsockets);
}


template <class SolverType>
void InitChandeliere(SolverType & FEMTestSolver, std::vector<BallSocket, GlobalAllocator<BallSocket>> & ballsockets)
{
	FEMTestSolver.SetPrecision(0.0001f);
	FEMTestSolver.SetMaxIterations(1000);

	CMatrix3 ZeroMatrix;
	
	ChandelierDemo::invInertiaMatrix.MakeIdentity();

	// r = 0.5; r*r = 0.25;
	// mr*r * (2/5) = 0.4 * mass * 0.25
	float mass = 1.0f / ChandelierDemo::invMass;
	ChandelierDemo::invInertiaMatrix.mMatrix[0][0] = 1.0f / (0.1f * mass);
	ChandelierDemo::invInertiaMatrix.mMatrix[1][1] = 1.0f / (0.1f * mass);
	ChandelierDemo::invInertiaMatrix.mMatrix[2][2] = 1.0f / (0.1f * mass);

	ZeroMatrix.MakeZero();

	// Node 1 - Fixed
	FEMTestSolver.AddTranslationalNode(	0.0f,
										Vec3(0.0f, 0.0f, 0.0f),
										Vec3(0.0f, 0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ZeroMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);

	FEMTestSolver.AddTranslationalNode(	ChandelierDemo::invMass,
										Vec3(0.0f, -1.0f, 0.0f),
										Vec3(0.0f,  0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ChandelierDemo::invInertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);

	int NodeIdx1 = 1;
	int NodeIdx2 = 3;

	ChandelierDemo::Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
	ChandelierDemo::Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
	ChandelierDemo::Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

	ChandelierDemo::SimpleSocket.Init(FEMTestSolver, ChandelierDemo::JointERP, ChandelierDemo::JointCFM, ChandelierDemo::Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

	ballsockets.push_back(ChandelierDemo::SimpleSocket);

	ChandeliereAddLevel(0, 3, 8.0f, FEMTestSolver, ballsockets);
}
#endif

class demoscene
{
public:

	std::vector<FEMJoint, GlobalAllocator<FEMJoint>> tetrajoints;
	std::vector<BallSocket, GlobalAllocator<BallSocket>> ballsockets;
	std::vector<BallSocket_FEM, GlobalAllocator<BallSocket_FEM>> ballsockets_fem;

	std::vector<PlaneConstraint, GlobalAllocator<PlaneConstraint>> planejoints;

#if (HARDWARE_SOLVER == HW_SOLVER_LCPCG)

	SolverLCPCG_CUDA FEMTestSolver;

#if (SOLVER_GOLD == 1)
 	SolverLCPCG FEMGoldSolver;
#endif

#elif (HARDWARE_SOLVER == HW_SOLVER_JACOBI)

	SolverJacobi_CUDA FEMTestSolver;

#if (SOLVER_GOLD == 1)
 	SolverLCPCG FEMGoldSolver;
#endif

#else

#	if (SOFTWARE_SOLVER == SW_SOLVER_CG)
#	pragma message(__FILE__" : warning: Non-LCP solver is selected! Constraint limits won't work!")
	SolverCG FEMTestSolver;
#	elif (SOFTWARE_SOLVER == SW_SOLVER_PGS)
	SolverPGS FEMTestSolver;
#	else
	SolverLCPCG FEMTestSolver;
#endif

#if (SOLVER_GOLD == 1)
	SolverPGS FEMGoldSolver;
#endif

#endif

	// Key to determine which lambda should we render
	int LambdaRender;

	demoscene():
#if (SOLVER_GOLD == 1)
				FEMGoldSolver(15, 0.001f),
#endif
#if (DEMO_FEMCUBE_LONGBOX == 1)
				CubeSegWidth(96), CubeSegHeight(2), CubeSegDepth(1),
#else
				CubeSegWidth(8), CubeSegHeight(2), CubeSegDepth(1),
#endif
//				CubeSegWidth(8), CubeSegHeight(8), CubeSegDepth(8),
				Time(0.0f), FEMTestSolver(15, 0.001f),
				m_FE_TetraNum(0), m_FE_NodeNum(0)
	{
	}

	void Init(unsigned int deviceID)
	{
 		FEMTestSolver.SetInitialGuessCutFlag(false);
 		FEMTestSolver.SetInitialGuessCutThreshold(0.0f);
// 		FEMGoldSolver.SetInitialGuessCutFlag(true);
// 		FEMGoldSolver.SetInitialGuessCutThreshold(0.0f);

		int i;

#if (HARDWARE_SOLVER != HW_SOLVER_NONE)

#if (HARDWARE_SOLVER == HW_SOLVER_LCPCG)
		HRSRC hRSrc_PTX = FindResource(NULL, MAKEINTRESOURCE(IDR_CUDA_PTX1), "CUDA_PTX");
#elif (HARDWARE_SOLVER == HW_SOLVER_JACOBI)
		HRSRC hRSrc_PTX = FindResource(NULL, MAKEINTRESOURCE(IDR_CUDA_PTX2), "CUDA_PTX");
#endif

		int RSrc_PTX_Size = SizeofResource(NULL, hRSrc_PTX);
		char *ptr_ptx = new char[RSrc_PTX_Size + 1];//(char *)LockResource(LoadResource(NULL, hRSrc_PTX));

		memcpy( ptr_ptx, LockResource(LoadResource(NULL, hRSrc_PTX)), RSrc_PTX_Size * sizeof(char) );
		ptr_ptx[RSrc_PTX_Size] = '\0';

		FEMTestSolver.m_MainPTX_Image = ptr_ptx;
//		FEMTestSolver.m_MainPTX_Image = (char *)LockResource(LoadResource(NULL, hRSrc_PTX));
		FEMTestSolver.InitCUDA(deviceID);

		delete [] ptr_ptx;

#endif

#if (DEMO_FEMCUBE == 1)


		LambdaRender = 0;

		FEMTestSolver.SetPrecision(0.0001f);
//		FEMTestSolver.SetPrecision(0.0f);
		FEMTestSolver.SetMaxIterations(1000);

#if (SOLVER_GOLD == 1)
		FEMGoldSolver.SetPrecision(0.0001f);
//		FEMGoldSolver.SetPrecision(0.0f);
		FEMGoldSolver.SetMaxIterations(500);
#endif

#if (DEMO_FEMCUBE_MODEL == 1)

		CMatrix3 CubeRotation1;
		CubeRotation1.MakeIdentity();

		FEMCube.FEMNodes.push_back(Vec3(0.0f, 0.0f, 0.0f));
//		BuildCube(FEMCube, Vec3(0.0f, 10.5f, 0.0f), Vec3(4.0f, 1.0f, 0.5f), CubeRotation1, CubeSegWidth, CubeSegHeight, CubeSegDepth);
//		BuildCube(FEMCube, Vec3(0.0f, 10.5f, 0.0f), Vec3(4.0f, 4.0f, 4.0f), CubeRotation1, CubeSegWidth, CubeSegHeight, CubeSegDepth);
		BuildFromFile("resources\\bunny.tet", FEMCube, Vec3(0.0f, 7.5f, 0.0f), Vec3(25.0f, 25.0f, 25.0f), CubeRotation1);
//		BuildFromFile("resources\\dragon.tet", FEMCube, Vec3(0.0f, 1.5f, 0.0f), Vec3(1.0f, 1.0f, 1.0f), CubeRotation1);

		CheckSelfIntersections(FEMCube);

//		BuildFromFile_Tetgen("resources\\sph.1.node", "resources\\sph.1.ele", FEMCube, Vec3(0.0f, 1.5f, 0.0f), Vec3(0.025f, 0.025f, 0.025f), CubeRotation1);

		// Set cow position is needed for beautiful chain attaching
//		BuildFromFile("resources\\cow_out.tet", FEMCube, Vec3(-0.8f, 9.7f, -0.5f), Vec3(1.0f, 1.0f, 1.0f), CubeRotation1);
		
		gLog.Print("Number of FE: %d;\n", FEMCube.FEMTetras.size());

		int size = FEMCube.FEMNodes.size();
		m_FE_NodeNum = size;

		// based on density of rubber == 1100 kg / m^3
		float TotalMass = 2200.0f;

	// First vertex should be zero-vertex
	//Object.FEMNodes.push_back(Vec3(0.0f, 0.0f, 0.0f));

		//TODO: remember the offset for the cube in the solver
		for (int node = 1; node < size; ++node)
		{
/*
			float mass = TotalMass / (float)size;

#if (DEMOCUBE_FIXED == 1)
// 			float invmass = node < 2 ? 0.0f : (1.0f / mass);
// 			if (node == 5 || node == 13) invmass = 0.0f;

 			float invmass = (node == 13) ? 0.0f : (1.0f / mass);
#else
			float invmass = 1.0f / mass;
#endif
*/
			float invmass = 0.0f;


			unsigned idx = FEMTestSolver.AddTranslationalNode(invmass, FEMCube.FEMNodes[node], Vec3(0.0f, 0.0f, 0.0f));
		}

		tetrajoints.resize(FEMCube.FEMTetras.size());
		
		std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>>::iterator itTetraEnd = FEMCube.FEMTetras.end(), itTetra;
		std::vector<FEMJoint, GlobalAllocator<FEMJoint>>::iterator itJointEnd = tetrajoints.end(), itJoint;

		// Elasticity
//		float Young_Modulus = 1000.0f * 1000.0f * 1000.0f * 1000.0f * 1.0f;
		float Young_Modulus = 1000.0f * 1000.0f * 1.0f * 1.0f;
		float Poisson_Ratio = 0.33f;

		// Plasticity
		// c_creep [0..1/dt]
		// c_yeld - threshold for plastic strain
		// c_max - maximum plastic strain for FE
//		float c_yeld = 0.05f, c_creep = 0.01f, c_max = 10.0f;
		float c_yeld = 0.05f, c_creep = 0.01f, c_max = 10.0f;
		float beta = 0.05f;										// Beta (C.O.)

//		float Regularization = 0.5f;
		float Regularization = 0.005f;

		unsigned FE_cnt = 0;
		for (itTetra = FEMCube.FEMTetras.begin(), itJoint = tetrajoints.begin(); itTetra != itTetraEnd; ++itTetra, ++itJoint)
		{
			++FE_cnt;

			// Regularization = 0.0001
#if (SOLVER_GOLD == 1)
			itJoint->Init(FEMGoldSolver, FEMCube.FEMNodes, itTetra->ind1, itTetra->ind2, itTetra->ind3, itTetra->ind4,
							Young_Modulus, Poisson_Ratio, c_yeld, c_creep, c_max, beta, 0.005f);
#endif

			float mass = TotalMass / (float)size;


			FEMTestSolver.m_NodeInvMass0_00[itTetra->ind1] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_11[itTetra->ind1] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_22[itTetra->ind1] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_00[itTetra->ind2] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_11[itTetra->ind2] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_22[itTetra->ind2] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_00[itTetra->ind3] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_11[itTetra->ind3] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_22[itTetra->ind3] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_00[itTetra->ind4] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_11[itTetra->ind4] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass0_22[itTetra->ind4] = 1.0f / mass;

			FEMTestSolver.m_NodeInvMass_00[itTetra->ind1] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_11[itTetra->ind1] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_22[itTetra->ind1] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_00[itTetra->ind2] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_11[itTetra->ind2] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_22[itTetra->ind2] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_00[itTetra->ind3] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_11[itTetra->ind3] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_22[itTetra->ind3] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_00[itTetra->ind4] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_11[itTetra->ind4] = 1.0f / mass;
			FEMTestSolver.m_NodeInvMass_22[itTetra->ind4] = 1.0f / mass;


			itJoint->Init(FEMTestSolver, FEMCube.FEMNodes, itTetra->ind1, itTetra->ind2, itTetra->ind3, itTetra->ind4,
							Young_Modulus, Poisson_Ratio, c_yeld, c_creep, c_max, beta, Regularization);
		}

		m_FE_TetraNum = FE_cnt;

#elif (DEMO_FEMCUBE_MODEL == 0)

		CMatrix3 CubeRotation1;
		CubeRotation1.MakeIdentity();

		FEMCube.FEMNodes.push_back(Vec3(0.0f, 0.0f, 0.0f));
#if (DEMO_FEMCUBE_LONGBOX == 1)
		BuildCube(FEMCube, Vec3(0.0f, 10.5f, 0.0f), Vec3(48.0f, 1.0f, 0.5f), CubeRotation1, CubeSegWidth, CubeSegHeight, CubeSegDepth);
#else
		BuildCube(FEMCube, Vec3(0.0f, 10.5f, 0.0f), Vec3(4.0f, 1.0f, 0.5f), CubeRotation1, CubeSegWidth, CubeSegHeight, CubeSegDepth);
#endif
		
		gLog.Print("Number of FE: %d;\n", FEMCube.FEMTetras.size());

		int size = FEMCube.FEMNodes.size();

		m_FE_NodeNum = size;

		// based on density of rubber == 1100 kg / m^3
#if (DEMO_FEMCUBE_LONGBOX == 1)
		float TotalMass = 2200.0f;
#else
		float TotalMass = 20.0f;
#endif

		//TODO: remember the offset for the cube in the solver
		for (int node = 1; node < size; ++node)
		{
			float mass = TotalMass / (float)size;

#if (DEMOCUBE_FIXED == 1)
			float invmass = node < 2 ? 0.0f : (1.0f / mass);
			if (node == 5 || node == 13) invmass = 0.0f;
#else
			float invmass = 1.0f / mass;
#endif

			unsigned idx = FEMTestSolver.AddTranslationalNode(invmass, FEMCube.FEMNodes[node], Vec3(0.0f, 0.0f, 0.0f));

#if (SOLVER_GOLD == 1)
			FEMGoldSolver.AddTranslationalNode(invmass, FEMCube.FEMNodes[node], Vec3(0.0f, 0.0f, 0.0f));
#endif
		}

		tetrajoints.resize(FEMCube.FEMTetras.size());
		
		std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>>::iterator itTetraEnd = FEMCube.FEMTetras.end(), itTetra;
		std::vector<FEMJoint, GlobalAllocator<FEMJoint>>::iterator itJointEnd = tetrajoints.end(), itJoint;

		// Elasticity
//		float Young_Modulus = 1000.0f * 1000.0f * 1000.0f * 1000.0f * 1.0f;
#if (DEMO_FEMCUBE_LONGBOX == 1)
		float Young_Modulus = 1000.0f * 1000.0f * 1000.0f * 1.0f * 1.0f;
#else
		float Young_Modulus = 1000.0f * 5.0f * 1.0f * 1.0f * 1.0f;
#endif
		float Poisson_Ratio = 0.35f;

		// Plasticity
		// c_creep [0..1/dt]
		// c_yeld - threshold for plastic strain
		// c_max - maximum plastic strain for FE
		float c_yeld = 0.05f, c_creep = 0.01f, c_max = 0.3f;
		float beta = 0.03f;										// Beta (C.O.)

		unsigned FE_cnt = 0;
		for (itTetra = FEMCube.FEMTetras.begin(), itJoint = tetrajoints.begin(); itTetra != itTetraEnd; ++itTetra, ++itJoint)
		{
			++FE_cnt;

			// Regularization = 0.0001
#if (SOLVER_GOLD == 1)
			itJoint->Init(FEMGoldSolver, FEMCube.FEMNodes, itTetra->ind1, itTetra->ind2, itTetra->ind3, itTetra->ind4,
							Young_Modulus, Poisson_Ratio, c_yeld, c_creep, c_max, beta, 0.005f);
#endif

			itJoint->Init(FEMTestSolver, FEMCube.FEMNodes, itTetra->ind1, itTetra->ind2, itTetra->ind3, itTetra->ind4,
							Young_Modulus, Poisson_Ratio, c_yeld, c_creep, c_max, beta, 0.000005f);
		}

		m_FE_TetraNum = FE_cnt;

#endif

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////


#if (HARDWARE_SOLVER != HW_SOLVER_NONE)

		FEMTestSolver.SendFEMInitParams_GPU(tetrajoints);

#endif

#if (DEMO_FEMCUBE_PLANE == 1)

		planejoints.resize(size);

		std::vector<PlaneConstraint, GlobalAllocator<PlaneConstraint>>::iterator itPlaneEnd = planejoints.end(), itPlane;

		Vec3 PlaneNrm(0, 1, 0), PlanePnt(0, -3, 0);
		unsigned Node_cnt = 0;

		for (itPlane = planejoints.begin(); itPlane != itPlaneEnd; ++itPlane)
		{
//			itPlane->Init(FEMTestSolver, 0.05f, 0.0f, PlanePnt, PlaneNrm, Node_cnt++, 0.05f);
			itPlane->Init(FEMTestSolver, 0.08f, 0.0f, PlanePnt, PlaneNrm, Node_cnt++, 0.01f);

#if (SOLVER_GOLD == 1)
			itPlane->Init(FEMGoldSolver, 0.08f, 0.0f, PlanePnt, PlaneNrm, Node_cnt, 0.01f);
#endif
		}

#endif

#if (DEMO_FEMCUBE_NUMCHAINS >= 1)

		// CHAIN
		int StartNumNodes = FEMTestSolver.m_NumNodes;

		float mass = 2.0f;

		CMatrix3 InertiaMatrix, ZeroMatrix;
		InertiaMatrix.MakeIdentity();
		InertiaMatrix.mMatrix[0][0] = 1.0f / (0.1f * mass);
		InertiaMatrix.mMatrix[1][1] = 1.0f / (0.1f * mass);
		InertiaMatrix.mMatrix[2][2] = 1.0f / (0.1f * mass);

		ZeroMatrix.MakeZero();

		int NumLinks = 10;

		// Add rigid bodies
		// 1st
		for (i = 0; i < NumLinks + 1; ++i)
		{
#if (SOLVER_GOLD == 1)
			FEMGoldSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(25 - i*2.0f, 10.3f, -0.3f), NullVec3);
			FEMGoldSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
#endif

			FEMTestSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(25 - i*2.0f, 10.3f, -0.3f), NullVec3);
			FEMTestSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
		}

		BallSocket SimpleSocket;
		Vec3 Anchor;

		// 1st
		for (i = 0; i < NumLinks; ++i)
		{
			Anchor.x =	 (FEMTestSolver.m_NodePosRot_x[StartNumNodes + (i) * 2]
						+ FEMTestSolver.m_NodePosRot_x[StartNumNodes + (i + 1) * 2]) * 0.5f;
			Anchor.y =	 (FEMTestSolver.m_NodePosRot_y[StartNumNodes + (i) * 2]
						+ FEMTestSolver.m_NodePosRot_y[StartNumNodes + (i + 1) * 2]) * 0.5f;
			Anchor.z =	 (FEMTestSolver.m_NodePosRot_z[StartNumNodes + (i) * 2]
						+ FEMTestSolver.m_NodePosRot_z[StartNumNodes + (i + 1) * 2]) * 0.5f;

#if (SOLVER_GOLD == 1)
			SimpleSocket.Init(FEMGoldSolver, 0.5f, 0.01f, Anchor, StartNumNodes + (i)*2, StartNumNodes + (i)*2+1, StartNumNodes + (i + 1)*2, StartNumNodes + (i + 1)*2+1);
#endif

			SimpleSocket.Init(FEMTestSolver, 0.5f, 0.01f, Anchor, StartNumNodes + (i)*2, StartNumNodes + (i)*2+1, StartNumNodes + (i + 1)*2, StartNumNodes + (i + 1)*2+1);

			ballsockets.push_back(SimpleSocket);
		}

		int LastLinearTriple = FEMTestSolver.m_NumNodes - 2;

		Anchor.x = FEMTestSolver.m_NodePosRot_x[LastLinearTriple] - 0.5f;
		Anchor.y = FEMTestSolver.m_NodePosRot_y[LastLinearTriple];
		Anchor.z = FEMTestSolver.m_NodePosRot_z[LastLinearTriple];

		int trIdx[3];

 		trIdx[0] = 49;
 		trIdx[1] = 50;
 		trIdx[2] = 52;

/*
		// Short FEMCube
		trIdx[0] = 49;
		trIdx[1] = 51;
		trIdx[2] = 52;

		// Short Cow
		trIdx[0] = 48;
		trIdx[1] = 50;
		trIdx[2] = 44;
*/

		BallSocket_FEM SimpleFemSocket;


#if (SOLVER_GOLD == 1)
		SimpleFemSocket.Init(FEMGoldSolver, 0.5f, 0.01f, Anchor, LastLinearTriple, LastLinearTriple + 1, trIdx[0], trIdx[1], trIdx[2]);
#endif

		SimpleFemSocket.Init(FEMTestSolver, 0.5f, 0.01f, Anchor, LastLinearTriple, LastLinearTriple + 1, trIdx[0], trIdx[1], trIdx[2]);
		ballsockets_fem.push_back(SimpleFemSocket);


#if (DEMO_FEMCUBE_NUMCHAINS > 1)
		// CHAIN
		StartNumNodes = FEMTestSolver.m_NumNodes;

		NumLinks = 10;

		// Add rigid bodies
		// 1st
		for (i = 0; i < NumLinks + 1; ++i)
		{
#if (SOLVER_GOLD == 1)
			FEMGoldSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(-25 + i*2.0f, 10.3f, -0.3f), NullVec3);
			FEMGoldSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
#endif

			FEMTestSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(-25 + i*2.0f, 10.3f, -0.3f), NullVec3);
			FEMTestSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
		}

		// 1st
		for (i = 0; i < NumLinks; ++i)
		{
			Anchor.x =	 (FEMTestSolver.m_NodePosRot_x[StartNumNodes + (i) * 2]
						+ FEMTestSolver.m_NodePosRot_x[StartNumNodes + (i + 1) * 2]) * 0.5f;
			Anchor.y =	 (FEMTestSolver.m_NodePosRot_y[StartNumNodes + (i) * 2]
						+ FEMTestSolver.m_NodePosRot_y[StartNumNodes + (i + 1) * 2]) * 0.5f;
			Anchor.z =	 (FEMTestSolver.m_NodePosRot_z[StartNumNodes + (i) * 2]
						+ FEMTestSolver.m_NodePosRot_z[StartNumNodes + (i + 1) * 2]) * 0.5f;

#if (SOLVER_GOLD == 1)
			SimpleSocket.Init(FEMGoldSolver, 0.5f, 0.01f, Anchor, StartNumNodes + (i)*2, StartNumNodes + (i)*2+1, StartNumNodes + (i + 1)*2, StartNumNodes + (i + 1)*2+1);
#endif

			SimpleSocket.Init(FEMTestSolver, 0.5f, 0.01f, Anchor, StartNumNodes + (i)*2, StartNumNodes + (i)*2+1, StartNumNodes + (i + 1)*2, StartNumNodes + (i + 1)*2+1);

			ballsockets.push_back(SimpleSocket);
		}

		LastLinearTriple = FEMTestSolver.m_NumNodes - 2;

		Anchor.x = FEMTestSolver.m_NodePosRot_x[LastLinearTriple] - 0.5f;
		Anchor.y = FEMTestSolver.m_NodePosRot_y[LastLinearTriple];
		Anchor.z = FEMTestSolver.m_NodePosRot_z[LastLinearTriple];

 		trIdx[0] = 1;
 		trIdx[1] = 2;
 		trIdx[2] = 3;

#if (SOLVER_GOLD == 1)
		SimpleFemSocket.Init(FEMGoldSolver, 0.5f, 0.01f, Anchor, LastLinearTriple, LastLinearTriple + 1, trIdx[0], trIdx[1], trIdx[2]);
#endif

		SimpleFemSocket.Init(FEMTestSolver, 0.5f, 0.01f, Anchor, LastLinearTriple, LastLinearTriple + 1, trIdx[0], trIdx[1], trIdx[2]);
		ballsockets_fem.push_back(SimpleFemSocket);
#endif

#endif

#elif (DEMO_FEMCUBE == 2)

		FEMTestSolver.SetPrecision(0.0001f);
		FEMTestSolver.SetMaxIterations(1000);

		float mass = 0.1f;

		CMatrix3 InertiaMatrix, ZeroMatrix;
		InertiaMatrix.MakeIdentity();

		// r = 0.5; r*r = 0.25;
		// mr*r * (2/5) = 0.4 * mass * 0.25
		InertiaMatrix.mMatrix[0][0] = 1.0f / (0.1f * mass);
		InertiaMatrix.mMatrix[1][1] = 1.0f / (0.1f * mass);
		InertiaMatrix.mMatrix[2][2] = 1.0f / (0.1f * mass);

		ZeroMatrix.MakeZero();

		int NumLinks = 50;

		// Add rigid bodies
		// 1st
		for (i = 0; i < NumLinks + 1; ++i)
		{
			FEMTestSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(-i*2.0f, 0.0f, 0.0f), NullVec3);
			FEMTestSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
		}



		// 2nd
		for (i = 0; i < NumLinks + 1; ++i)
		{
			FEMTestSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(-i*2.0f, 0.0f, 3.0f), NullVec3);
			FEMTestSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
		}

		// 3rd
		for (i = 0; i < NumLinks + 1; ++i)
		{
			FEMTestSolver.AddTranslationalNode((i == 0 || i == NumLinks+1)?0.0f:(1.0f / mass), Vec3(-i*2.0f, 0.0f, -3.0f), NullVec3);
			FEMTestSolver.AddRotationalNode((i == 0 || i == NumLinks+1)?ZeroMatrix:InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
		}

/*
				// Shifted (zero node is reserved)
				int NodeIdx1 = (i * GridDepth + j)*2 + 1;
				int NodeIdx2 = (i * GridDepth + j - 1)*2 + 1;

				Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
				Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
				Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

				SimpleSocket.Init(FEMTestSolver, JointERP, JointCFM, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

				ballsockets.push_back(SimpleSocket);
*/


		BallSocket SimpleSocket;
		Vec3 Anchor;

		// 1st
		for (i = 0; i < NumLinks; ++i)
		{
			int NodeIdx1 = i * 2 + 1;
			int NodeIdx2 = (i+1) * 2 + 1;

			Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
			Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
			Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

			SimpleSocket.Init(FEMTestSolver, 0.1f, 0.09f, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

			ballsockets.push_back(SimpleSocket);
		}




		// 2nd
		for (i = NumLinks+1; i < (NumLinks)*2; ++i)
		{
			int NodeIdx1 = i * 2 + 1;
			int NodeIdx2 = (i+1) * 2 + 1;

			Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
			Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
			Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

			SimpleSocket.Init(FEMTestSolver, 0.1f, 0.09f, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

			ballsockets.push_back(SimpleSocket);
		}



		// 3rd
		for (i = (NumLinks+1)*2; i < NumLinks*3; ++i)
		{
			int NodeIdx1 = i * 2 + 1;
			int NodeIdx2 = (i+1) * 2 + 1;

			Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
			Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
			Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

			SimpleSocket.Init(FEMTestSolver, 0.1f, 0.09f, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

			ballsockets.push_back(SimpleSocket);
		}

#elif (DEMO_FEMCUBE == 3)

		FEMTestSolver.SetPrecision(0.0001f);
		FEMTestSolver.SetMaxIterations(140);

		float mass = 1.0f;

		CMatrix3 InertiaMatrix, ZeroMatrix;
		InertiaMatrix.MakeIdentity();

		// r = 0.5; r*r = 0.25;
		// mr*r * (2/5) = 0.4 * mass * 0.25
		InertiaMatrix.mMatrix[0][0] = 1.0f / (0.1f * mass);
		InertiaMatrix.mMatrix[1][1] = 1.0f / (0.1f * mass);
		InertiaMatrix.mMatrix[2][2] = 1.0f / (0.1f * mass);

		ZeroMatrix.MakeZero();

		BallSocket SimpleSocket;
		Vec3 Anchor;

		Vec3 TailVelocity(20.0f, 0.0f, 20.0f);

#if 0
		float JointERP = 0.5f, JointCFM = 0.001f;

		int GridWidth = 23;
		int GridDepth = 10;
		int GridHeight = 40;

		// Shift about Y-axis
		float GridShift = 25.0f;

		Vec3 GridSize(100.0f, 40.0f, 100.0f);

		for (int i = 0; i < GridWidth; ++i)
		{
			for (int j = 0; j < GridDepth; ++j)
			{
				for (int k = 0; k < GridHeight; ++k)
				{
					float inv_mass;
					inv_mass = 1.0f / mass;
					if (k == 0) inv_mass = 0.0f;

					CMatrix3 *rot_inv_mass;
					rot_inv_mass = &InertiaMatrix;

					Vec3 Velocity(0.0f, 0.0f, 0.0f);
					if (k == GridHeight - 1)
					{
						Velocity.x = ((rand()%1000) / 500.0f - 1.0f) * TailVelocity.x;
						Velocity.y = ((rand()%1000) / 500.0f - 1.0f) * TailVelocity.y;
						Velocity.z = ((rand()%1000) / 500.0f - 1.0f) * TailVelocity.z;
					}


					FEMTestSolver.AddTranslationalNode(inv_mass,
							Vec3((i - GridWidth / 2) / (float)GridWidth * GridSize.x, GridShift - k / (float)GridHeight * GridSize.y, (j - GridDepth / 2) / (float)GridDepth * GridSize.z),
														Velocity);
					FEMTestSolver.AddRotationalNode(*rot_inv_mass, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
				}
			}
		}

		for (int i = 0; i < GridWidth; ++i)
		{
			for (int j = 0; j < GridDepth; ++j)
			{
				for (int k = 1; k < GridHeight; ++k)
				{
					// Shifted (zero node is reserved)
					int NodeIdx1 = (i * GridDepth * GridHeight + j * GridHeight + k)*2 + 1;
					int NodeIdx2 = (i * GridDepth * GridHeight + j * GridHeight + k - 1)*2 + 1;

					Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
					Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
					Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

					SimpleSocket.Init(FEMTestSolver, JointERP, JointCFM, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

					ballsockets.push_back(SimpleSocket);
				}
			}
		}
#else
		float JointERP = 0.2f, JointCFM = 0.005f;

// 		int GridWidth = 23;
// 		int GridDepth = 10;
// 		int GridHeight = 40;

		int GridWidth = 20;
		int GridDepth = 20;
		int GridHeight = 20;

		// Shift about Y-axis
		float GridShift = 25.0f;

		Vec3 GridSize(100.0f, 80.0f, 100.0f);

		int bodyIdxShift = 0;
		for (int i = 0; i < GridWidth; ++i)
		{
			for (int j = 0; j < GridDepth; ++j)
			{
				float inv_mass;
				inv_mass = 1.0f / mass;

				if (i == 0 || j == 0)
					inv_mass = 0.0f;
				if ((i == GridWidth - 1) || (j == GridDepth - 1))
					inv_mass = 0.0f;

				CMatrix3 *rot_inv_mass;
				rot_inv_mass = &InertiaMatrix;

				Vec3 Velocity(0.0f, 0.0f, 0.0f);

// 				if (k == GridHeight - 1)
// 				{
// 					Velocity.x = ((rand()%1000) / 500.0f - 1.0f) * TailVelocity.x;
// 					Velocity.y = ((rand()%1000) / 500.0f - 1.0f) * TailVelocity.y;
// 					Velocity.z = ((rand()%1000) / 500.0f - 1.0f) * TailVelocity.z;
// 				}


				FEMTestSolver.AddTranslationalNode(inv_mass,
						Vec3((i - GridWidth / 2) / (float)GridWidth * GridSize.x, GridShift, (j - GridDepth / 2) / (float)GridDepth * GridSize.z),
													Velocity);
				FEMTestSolver.AddRotationalNode(*rot_inv_mass, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);

				++bodyIdxShift;
			}
		}

		// NET
		//////////////////////////////////////////////////////////////////////////

		for (int i = 1; i < GridWidth - 1; ++i)
		{
			for (int j = 1; j < GridDepth; ++j)
			{
				// Shifted (zero node is reserved)
				int NodeIdx1 = (i * GridDepth + j)*2 + 1;
				int NodeIdx2 = (i * GridDepth + j - 1)*2 + 1;

				Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
				Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
				Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

				SimpleSocket.Init(FEMTestSolver, JointERP, JointCFM, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

				ballsockets.push_back(SimpleSocket);
			}
		}

		for (int i = 1; i < GridWidth; ++i)
		{
			for (int j = 1; j < GridDepth-1; ++j)
			{
				// Shifted (zero node is reserved)
				int NodeIdx1 = (i * GridDepth + j)*2 + 1;
				int NodeIdx2 = ((i - 1) * GridDepth + j)*2 + 1;

				Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
				Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
				Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

				SimpleSocket.Init(FEMTestSolver, JointERP, JointCFM, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

				ballsockets.push_back(SimpleSocket);
			}
		}

		// CHAINS
		//////////////////////////////////////////////////////////////////////////

		for (int i = 1; i < GridWidth - 1; ++i)
		{
			for (int j = 1; j < GridDepth - 1; ++j)
			{
				for (int k = 0; k < GridHeight; ++k)
				{
					float inv_mass;
					inv_mass = 1.0f / mass;

					Vec3 Velocity(0.0f, 0.0f, 0.0f);

					CMatrix3 *rot_inv_mass;
					rot_inv_mass = &InertiaMatrix;

					FEMTestSolver.AddTranslationalNode(inv_mass,
							Vec3((i - GridWidth / 2) / (float)GridWidth * GridSize.x, GridShift - k / (float)GridHeight * GridSize.y, (j - GridDepth / 2) / (float)GridDepth * GridSize.z),
														Velocity);
					FEMTestSolver.AddRotationalNode(*rot_inv_mass, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);
				}
			}
		}

		for (int i = 1; i < GridWidth-1; ++i)
		{
			for (int j = 1; j < GridDepth-1; ++j)
			{
				for (int k = 0; k < GridHeight; ++k)
				{
					// Shifted (zero node is reserved)
					int NodeIdx1 = ((i-1) * (GridDepth - 2) * GridHeight + (j-1) * GridHeight + k + bodyIdxShift)*2 + 1;
					int NodeIdx2;

					if (k > 0)
						NodeIdx2 = ((i-1) * (GridDepth - 2) * GridHeight + (j-1) * GridHeight + k - 1 + bodyIdxShift)*2 + 1;
					else
						NodeIdx2 = (i * GridDepth + j)*2 + 1;

					Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
					Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
					Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

					SimpleSocket.Init(FEMTestSolver, JointERP, JointCFM, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

					ballsockets.push_back(SimpleSocket);
				}
			}
		}


#endif

#elif (DEMO_FEMCUBE == 4)

/*
	FEMTestSolver.SetPrecision(0.0001f);
	FEMTestSolver.SetMaxIterations(1000);

	float mass = 1.0f;

	CMatrix3 InertiaMatrix, ZeroMatrix;
	InertiaMatrix.MakeIdentity();

	// r = 0.5; r*r = 0.25;
	// mr*r * (2/5) = 0.4 * mass * 0.25
	InertiaMatrix.mMatrix[0][0] = 1.0f / (0.1f * mass);
	InertiaMatrix.mMatrix[1][1] = 1.0f / (0.1f * mass);
	InertiaMatrix.mMatrix[2][2] = 1.0f / (0.1f * mass);

	ZeroMatrix.MakeZero();

	float inv_mass;
	inv_mass = 1.0f / mass;

	CMatrix3 *rot_inv_mass;
	rot_inv_mass = &InertiaMatrix;

	BallSocket SimpleSocket;
	Vec3 Anchor;

	float JointERP = 0.2f, JointCFM = 0.005f;

	// Node 1 - Fixed
	FEMTestSolver.AddTranslationalNode(	0.0f,
										Vec3(0.0f, 0.0f, 0.0f),
										Vec3(0.0f, 0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(ZeroMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);

	FEMTestSolver.AddTranslationalNode(	inv_mass,
										Vec3(0.0f, -1.0f, 0.0f),
										Vec3(0.0f,  0.0f, 0.0f));
	FEMTestSolver.AddRotationalNode(InertiaMatrix, CQuaternion(1.0f, 0.0f, 0.0f, 0.0f), NullVec3);

	int NodeIdx1 = 1;
	int NodeIdx2 = 3;

	Anchor.x = (FEMTestSolver.m_NodePosRot_x[NodeIdx1] + FEMTestSolver.m_NodePosRot_x[NodeIdx2]) * 0.5f;
	Anchor.y = (FEMTestSolver.m_NodePosRot_y[NodeIdx1] + FEMTestSolver.m_NodePosRot_y[NodeIdx2]) * 0.5f;
	Anchor.z = (FEMTestSolver.m_NodePosRot_z[NodeIdx1] + FEMTestSolver.m_NodePosRot_z[NodeIdx2]) * 0.5f;

	SimpleSocket.Init(FEMTestSolver, JointERP, JointCFM, Anchor, NodeIdx1, NodeIdx1+1, NodeIdx2, NodeIdx2+1);

	ballsockets.push_back(SimpleSocket);
*/
	
	InitChandeliere(FEMTestSolver, ballsockets);
	gLog.Print("Number of ballsocket joints: %d;\n", ballsockets.size());

#endif


#if (HARDWARE_SOLVER != HW_SOLVER_NONE)

		FEMTestSolver.SendInitParams_GPU();

#endif
	}
	
	int sort_order;

	void Draw(const GLDrawer &gld, float dt)
	{
		Time += dt;

		Timer OverallTimer;
		OverallTimer.Start();
		gLog.Print("FRAME BEGIN ---------------------------------------");

		if (!gbFreeze || gbStep)
		{
			float solver_dt = 0.016f;

#if (SOLVER_GOLD == 1)
			FEMGoldSolver.m_NumJoints = 0;
#endif

			FEMTestSolver.m_NumJoints = 0;

#if (HARDWARE_SOLVER != HW_SOLVER_NONE)

			FEMTestSolver.m_Num_FEM_Joints = 0;

#endif

			//std::random_shuffle(tetrajoints.begin(), tetrajoints.end());

			Timer PerfTimer;
			PerfTimer.Start();


			std::vector<FEMJoint, GlobalAllocator<FEMJoint>>::iterator itJointEnd = tetrajoints.end(), itJoint;

#if (HARDWARE_SOLVER == HW_SOLVER_NONE)

			for (itJoint = tetrajoints.begin(); itJoint != itJointEnd; ++itJoint)
			{
#if (SOLVER_GOLD == 1)
				itJoint->UpdateCopy(solver_dt, FEMGoldSolver);
#endif

				itJoint->UpdateCopy(solver_dt, FEMTestSolver);
			}

#else

			FEMTestSolver.SendFEMParams_GPU(tetrajoints);
			FEMTestSolver.m_NumJoints = 6 * FEMTestSolver.m_Num_FEM_Joints;

			FEMTestSolver.CalcFEMParams_GPU(solver_dt);

#endif

			float FEM_JacobUpdate = (float)PerfTimer.Time();
 //			gLog.Print("FEM Jacobians Update: %f ms;", FEM_JacobUpdate);


			std::vector<PlaneConstraint, GlobalAllocator<PlaneConstraint>>::iterator itPlaneEnd = planejoints.end(), itPlane;
			for (itPlane = planejoints.begin(); itPlane != itPlaneEnd; ++itPlane)
			{
#if (SOLVER_GOLD == 1)
				itPlane->UpdateCopy(solver_dt, FEMGoldSolver);
#endif

				itPlane->UpdateCopy(solver_dt, FEMTestSolver);
			}


			std::vector<BallSocket, GlobalAllocator<BallSocket>>::iterator itBallSocketEnd = ballsockets.end(), itBallSocket;
			for (itBallSocket = ballsockets.begin(); itBallSocket != itBallSocketEnd; ++itBallSocket)
			{
#if (SOLVER_GOLD == 1)
				itBallSocket->UpdateCopy(solver_dt, FEMGoldSolver);
#endif

				itBallSocket->UpdateCopy(solver_dt, FEMTestSolver);
			}

			std::vector<BallSocket_FEM, GlobalAllocator<BallSocket_FEM>>::iterator itBallSocketFEMEnd = ballsockets_fem.end(), itBallSocketFEM;
			for (itBallSocketFEM = ballsockets_fem.begin(); itBallSocketFEM != itBallSocketFEMEnd; ++itBallSocketFEM)
			{
#if (SOLVER_GOLD == 1)
				itBallSocketFEM->UpdateCopy(solver_dt, FEMGoldSolver);
#endif

				itBallSocketFEM->UpdateCopy(solver_dt, FEMTestSolver);
			}

//			gLog.Print("Num Rows: %d;", FEMTestSolver.m_NumJoints);
			FEMTestSolver.Solve(solver_dt);

#if (SOLVER_GOLD == 1)
			FEMGoldSolver.Solve(solver_dt);

			float error = 0.0f;
			for (int i = 0; i < FEMTestSolver.m_NumJoints; ++i)
			{
				error += math::fastabs(FEMTestSolver.m_lambda[i] - FEMGoldSolver.m_lambda[i]);
			}
			gLog.Print("Error from GOLD: %12.10f;", error);

// 			bool test_succ = true;
// 			for (int i = 0; i < FEMTestSolver.m_NumJoints; ++i)
// 			{
// 				if (math::fastabs(FEMTestSolver.m_lambda[i] - FEMGoldSolver.m_lambda[i]) >  1e-4f)
// 				{
// 					test_succ = false;
// 					gLog.Print("Test: FAILED;");
// 					break;
// 				}
// 			}
#endif


			itJointEnd = tetrajoints.end();
			for (itJoint = tetrajoints.begin(); itJoint != itJointEnd; ++itJoint)
			{
				itJoint->FetchLambdas(FEMTestSolver);
			}

			itPlaneEnd = planejoints.end();
			for (itPlane = planejoints.begin(); itPlane != itPlaneEnd; ++itPlane)
			{
				itPlane->FetchLambdas(FEMTestSolver);
			}


			itBallSocketEnd = ballsockets.end();
			for (itBallSocket = ballsockets.begin(); itBallSocket != itBallSocketEnd; ++itBallSocket)
			{
				itBallSocket->FetchLambdas(FEMTestSolver);
			}

			itBallSocketFEMEnd = ballsockets_fem.end();
			for (itBallSocketFEM = ballsockets_fem.begin(); itBallSocketFEM != itBallSocketFEMEnd; ++itBallSocketFEM)
			{
				itBallSocketFEM->FetchLambdas(FEMTestSolver);
			}

			UpdateFEM = false;

			if (gbStep)
			{
				gbStep = false;
			}
		}

		static float frametimeMax = 0.0f, frametimeSum = 0.0f;

		float PhysicsUpdate = (float)OverallTimer.Time();


		if (PhysicsUpdate > frametimeMax)
			frametimeMax = PhysicsUpdate;

		frametimeSum += PhysicsUpdate;

		gLog.Print("Physics Update: %f ms;", PhysicsUpdate);

		gld.DrawBasis(Vec3(-2.0f, -1.0f, 1.0f));
		gld.SetColor(1.0f, 1.0f, 0.0f, 0.3f);

		int cnt = 0;
		float slowTime = Time * 0.0002f;								// Time "slowed down"

		Vec3 PlaneNormal;
		PlaneNormal.x = sinf(1.5f * slowTime);
		PlaneNormal.y = cosf(1.3f * slowTime);
		PlaneNormal.z = sinf(1.9f * slowTime) * cosf(1.7f * slowTime);

		PlaneNormal.Normalize();

		unsigned jointCnt = 0;

		std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>>::iterator itTetraEnd = FEMCube.FEMTetras.end();

#if 1
		for (std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>>::iterator itTetra = FEMCube.FEMTetras.begin(); itTetra != itTetraEnd; ++itTetra)
		{
#if (RENDER_OBJECT == 1)
			const Vec3	&v1 = FEMCube.FEMNodes[itTetra->ind1],
						&v2 = FEMCube.FEMNodes[itTetra->ind2],
						&v3 = FEMCube.FEMNodes[itTetra->ind3],
						&v4 = FEMCube.FEMNodes[itTetra->ind4];
#else
			Vec3	v1 = FEMTestSolver.GetPosition(itTetra->ind1),
					v2 = FEMTestSolver.GetPosition(itTetra->ind2),
					v3 = FEMTestSolver.GetPosition(itTetra->ind3),
					v4 = FEMTestSolver.GetPosition(itTetra->ind4);

// 			Vec3	v1 = FEMGoldSolver.GetPosition(itTetra->ind1),
// 					v2 = FEMGoldSolver.GetPosition(itTetra->ind2),
// 					v3 = FEMGoldSolver.GetPosition(itTetra->ind3),
// 					v4 = FEMGoldSolver.GetPosition(itTetra->ind4);

#endif

			// Check if all vertices is in positive halfspace of [PlaneNormal, 0]
// 			if ((PlaneNormal.Dot(v1) < 0.0f) || (PlaneNormal.Dot(v2) < 0.0f) ||
// 				(PlaneNormal.Dot(v3) < 0.0f) || (PlaneNormal.Dot(v4) < 0.0f))
// 			{
// 				continue;
// 			}

#if (COLORING_STRAIN == 1)
			float FEstrain = tetrajoints[jointCnt].m_StrainNorm * 100.0f;
			gld.SetColor(1.0f + FEstrain, 1.0f - FEstrain * 0.5f, 1.0f - FEstrain * 0.5f, 1.0f);
#else
			float lambda[6];
			float lambdaRenderScale = 0.001f;
#if 0
			for (int i = 0; i < 6; ++i)
			{
				lambda[i] = 1.0f - lambdaRenderScale * fabsf(tetrajoints[jointCnt].m_lambda0[i]);
			}
#else
			float lambdaNorm = sqrtf(tetrajoints[jointCnt].m_lambda0[0]*tetrajoints[jointCnt].m_lambda0[0] +
									 tetrajoints[jointCnt].m_lambda0[1]*tetrajoints[jointCnt].m_lambda0[1] + 
									 tetrajoints[jointCnt].m_lambda0[2]*tetrajoints[jointCnt].m_lambda0[2]);

			lambda[0] = lambda[1] = lambda[2] = lambdaRenderScale * lambdaNorm;
#endif

			if (!LambdaRender)
			{
				// Linear
				gld.SetColor(lambda[0], lambda[1], lambda[2], 1.0f);
			}
			else
			{
				// Angular
				gld.SetColor(lambda[3], lambda[4], lambda[5], 1.0f);
			}
#endif

			gld.DrawTriangle(v3, v2, v1);
			gld.DrawTriangle(v4, v1, v2);
			gld.DrawTriangle(v4, v2, v3);
			gld.DrawTriangle(v4, v3, v1);

			gld.SetColor(0.0f, 0.0f, 0.0f, 1.0f);
			gld.DrawLine(v1, v2);
			gld.DrawLine(v1, v3);
			gld.DrawLine(v1, v4);
			gld.DrawLine(v2, v3);
			gld.DrawLine(v3, v4);
			gld.DrawLine(v4, v2);

			++jointCnt;
		}
#endif

		//FEMTestSolver.m_NodeF[(8*2)*3 + 2]
//		int idx = FEMTestSolver.NumNodes - 8*2 - 1; +Y
//		int idx = FEMTestSolver.NumNodes - 8*2 - 2; +Z
//		int idx = FEMTestSolver.NumNodes - 8*2 - 5; -Z
//		int idx = FEMTestSolver.NumNodes - 8*2 - 6; -Y

/*
//		int idx = 48;
		int idx = 49;
		gld.SetColor(1.0f, 0.0f, 0.0f, 1.0f);
		Vec3 NodeRenderPos = FEMTestSolver.GetPosition(idx);
		gld.DrawPoint(NodeRenderPos, 50.0f);
//*/

/*
		const int IndexNum = 3;
		int Indices[IndexNum] = { 48, 50, 44 };
//		int Indices[IndexNum] = { 0, 44, 0 };

		gld.SetColor(1.0f, 1.0f, 0.0f, 1.0f);
		for (int i = 0; i < IndexNum; ++i)
		{
			Vec3 NodeRenderPos = FEMTestSolver.GetPosition(Indices[i]);
			gld.DrawPoint(NodeRenderPos, 50.0f);
		}
//*/


		float RenderObject = (float)OverallTimer.Time();
		gLog.Print(" + Render Object: %f ms;", RenderObject);

		gld.SetColor(1.0f, 1.0f, 0.0f, 1.0f);

#if (LIGHT_RENDER == 0)

		std::vector<BallSocket, GlobalAllocator<BallSocket>>::iterator itBallSocketJEnd = ballsockets.end(), itBallSocketJ;
		for (itBallSocketJ = ballsockets.begin(); itBallSocketJ != itBallSocketJEnd; ++itBallSocketJ)
		{
			CMatrix3 rotMatrix;

			int b1_Lidx, b1_Ridx;
			int b2_Lidx, b2_Ridx;

			Vec3 b1_AnchorBS, b2_AnchorBS;

			Vec3 b1_Pos, b2_Pos;
			Vec3 b1_AnchorWS, b2_AnchorWS;

			itBallSocketJ->GetIndices(b1_Lidx, b1_Ridx, b2_Lidx, b2_Ridx);

			itBallSocketJ->GetAnchorPoints_BS(b1_AnchorBS, b2_AnchorBS);

			b1_Pos = FEMTestSolver.GetPosition(b1_Lidx);

			FEMTestSolver.GetOrientation(b1_Ridx).ToMatrix3(rotMatrix);
			b1_AnchorWS = b1_Pos + rotMatrix * b1_AnchorBS;

			gld.DrawLine(b1_Pos, b1_AnchorWS);
			gld.DrawPoint(b1_AnchorWS, 5.0f);

			if (b2_Lidx > 0)
			{
				b2_Pos = FEMTestSolver.GetPosition(b2_Lidx);

				FEMTestSolver.GetOrientation(b2_Ridx).ToMatrix3(rotMatrix);
				b2_AnchorWS = b2_Pos + rotMatrix * b2_AnchorBS;

				gld.DrawLine(b2_Pos, b2_AnchorWS);
				gld.DrawPoint(b2_AnchorWS, 5.0f);
			}
		}

#endif

#if 1
		std::vector<PlaneConstraint, GlobalAllocator<PlaneConstraint>>::iterator itPlaneJEnd = planejoints.end(), itPlaneJ;
		for (itPlaneJ = planejoints.begin(); itPlaneJ != itPlaneJEnd; ++itPlaneJ)
		{
			if (!itPlaneJ->m_bActive) continue;

			Vec3 NodePos = FEMTestSolver.GetPosition(itPlaneJ->m_BodyIdx);
			Vec3 plLambda = itPlaneJ->m_PlaneNormal * itPlaneJ->m_lambda0[0];

			

//			float color = fabsf(itPlaneJ->m_Violation);
			float color = fabsf(FEMTestSolver.m_resid[itPlaneJ->m_StartIdx]) * 1000.0f;

			gld.SetColor(0.7f, color * 0.25f, color * 0.5f, 1.0f);
			gld.DrawLine(NodePos, NodePos+plLambda);
		}
#endif
		gld.SetColor(1.0f, 1.0f, 1.0f, 1.0f);


		gld.SetColor(1.0f, 1.0f, 0.0f, 1.0f);

		Vec3 wsX, wsY, wsZ, Pos;
		for (unsigned int i = 0; i < FEMTestSolver.m_NumNodes - 1; ++i)
		{
			// We need linear nodes, sequenced by rotational ones [full rigid body]
			if ((!FEMTestSolver.m_IsRotational[i+1]) || FEMTestSolver.m_IsRotational[i]) continue;
 			Pos.x = FEMTestSolver.m_NodePosRot_x[i];
 			Pos.y = FEMTestSolver.m_NodePosRot_y[i];
 			Pos.z = FEMTestSolver.m_NodePosRot_z[i];

#if (LIGHT_RENDER == 1)

			gld.DrawPoint(Pos, 5.0f);

#else

			CMatrix3 Rotation;
			CQuaternion(FEMTestSolver.m_NodePosRot_w[i+1], FEMTestSolver.m_NodePosRot_x[i+1],
						FEMTestSolver.m_NodePosRot_y[i+1], FEMTestSolver.m_NodePosRot_z[i+1]).ToMatrix3(Rotation);


			// Get Basis
			wsX = Vec3(Rotation.mMatrix[0][0], Rotation.mMatrix[1][0], Rotation.mMatrix[2][0]);
			wsY = Vec3(Rotation.mMatrix[0][1], Rotation.mMatrix[1][1], Rotation.mMatrix[2][1]);
			wsZ = Vec3(Rotation.mMatrix[0][2], Rotation.mMatrix[1][2], Rotation.mMatrix[2][2]);

			gld.SetColor(1.0f, 0.0f, 0.0f, 1.0f);
			gld.DrawLine(Pos, Pos + wsX);

			gld.SetColor(0.0f, 1.0f, 0.0f, 1.0f);
			gld.DrawLine(Pos, Pos + wsY);

			gld.SetColor(0.0f, 0.0f, 1.0f, 1.0f);
			gld.DrawLine(Pos, Pos + wsZ);

#endif
		}

#if 0
		// DEBUG render nodes to attach to FEM object
		unsigned int i;
			
		i = 1;
		Pos.x = FEMTestSolver.m_NodePosRot_x[i];
		Pos.y = FEMTestSolver.m_NodePosRot_y[i];
		Pos.z = FEMTestSolver.m_NodePosRot_z[i];
		gld.DrawPoint(Pos, 5.0f);

		i = 2;
		Pos.x = FEMTestSolver.m_NodePosRot_x[i];
		Pos.y = FEMTestSolver.m_NodePosRot_y[i];
		Pos.z = FEMTestSolver.m_NodePosRot_z[i];
		gld.DrawPoint(Pos, 5.0f);

		i = 3;
		Pos.x = FEMTestSolver.m_NodePosRot_x[i];
		Pos.y = FEMTestSolver.m_NodePosRot_y[i];
		Pos.z = FEMTestSolver.m_NodePosRot_z[i];
		gld.DrawPoint(Pos, 5.0f);

		gld.SetColor(0.0f, 1.0f, 0.0f, 1.0f);

		i = 49;
		Pos.x = FEMTestSolver.m_NodePosRot_x[i];
		Pos.y = FEMTestSolver.m_NodePosRot_y[i];
		Pos.z = FEMTestSolver.m_NodePosRot_z[i];
		gld.DrawPoint(Pos, 15.0f);

		i = 50;
		Pos.x = FEMTestSolver.m_NodePosRot_x[i];
		Pos.y = FEMTestSolver.m_NodePosRot_y[i];
		Pos.z = FEMTestSolver.m_NodePosRot_z[i];
		gld.DrawPoint(Pos, 15.0f);

		i = 52;
		Pos.x = FEMTestSolver.m_NodePosRot_x[i];
		Pos.y = FEMTestSolver.m_NodePosRot_y[i];
		Pos.z = FEMTestSolver.m_NodePosRot_z[i];
		gld.DrawPoint(Pos, 15.0f);
#endif

		float NodesCycle = (float)OverallTimer.Time();
		gLog.Print(" + Nodes Cycle: %f ms;", NodesCycle);

		gld.PrintText("Cube Subdivision Demo");
		gld.PrintText("w: %f, %f, %f",	FEMTestSolver.m_NodeVel_x[FEMTestSolver.m_NumNodes - 2],
										FEMTestSolver.m_NodeVel_y[FEMTestSolver.m_NumNodes - 2],
										FEMTestSolver.m_NodeVel_z[FEMTestSolver.m_NumNodes - 2] );

		static unsigned int  iterationsMin = 200000, iterationsMax = 0, iterationsSum = 0;
		static unsigned int frameNum = 0;

		unsigned int effIterations = FEMTestSolver.GetEffectiveIterations();

		if (effIterations && (effIterations < iterationsMin))
			iterationsMin = effIterations;

		if (effIterations > iterationsMax)
			iterationsMax = effIterations;

		iterationsSum += effIterations;
		++frameNum;

		float avgIter = iterationsSum / (float)frameNum;
		float avgFrametime = frametimeSum / (float)frameNum;

		static unsigned int frameNumber = 0;

		gld.RenderText(-0.99f, 0.85f, "||Gradient||^2 = %f", FEMTestSolver.GetGradNormSq());
		gld.RenderText(-0.99f, 0.80f, "||Lambda||^2 = %f", FEMTestSolver.GetLambdaNormSq());
		gld.RenderText(-0.99f, 0.75f, "Lambda dot Grad = %f", FEMTestSolver.GetDotLambdaGrad());
		gld.RenderText(-0.99f, 0.70f, "Effective Iterations: %d (%d|%d) [%3.2f]", effIterations, iterationsMin, iterationsMax,avgIter);
		gld.RenderText(-0.99f, 0.65f, "Frametime: %3.2f (%3.2f) [%3.2f]", PhysicsUpdate, frametimeMax, avgFrametime);
		gld.RenderText(-0.99f, 0.60f, "System size: %d", FEMTestSolver.m_NumJoints);
		gld.RenderText(-0.99f, 0.55f, "Frames: %d", frameNumber);
		if (!gbFreeze)
			++frameNumber;

		if (m_FE_NodeNum > 0 && m_FE_TetraNum > 0)
		{
			gld.RenderText(-0.99f, 0.50f, "FE: %d Nodes, %d Tetras", m_FE_NodeNum, m_FE_TetraNum);
		}

		float OverallTime = (float)OverallTimer.Time();
		gLog.Print("Total: %f ms;", OverallTime);
	}

private:

	int m_FE_TetraNum, m_FE_NodeNum;

	int CubeSegWidth, CubeSegHeight, CubeSegDepth;
	FEMObject FEMCube;

	float Time;
} demo;

HWND hMainWindow = NULL;

GLDrawer drawer;
CCamera camera;
bool hasFocus = true;
unsigned char keys[256];				// Array for Key States

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);	// WndProc function prototype

const int windowX = 1280, windowY = 720;
const float znear = 0.01f, zfar = 100.0f;
const float fov = 45.0f;
const char *szTitle = "Physics demo framework ver. 0.01";

bool gbFullScreen;
bool gbWireframe = false;
float gDTime;

// FPS timer ID
#define FPS_TIMER 1
// For correct FPS counting and printing
int debugFPSCount, exFPSCount;

bool ProcessKeys(void)
{
	if(keys[VK_ESCAPE])
		return false;

	if (keys[VK_TAB])
	{
		// To prevent several switches
		keys[VK_TAB] = false;
		gbWireframe = !gbWireframe;
	}

	// Camera Movement
	if(keys['W'])
	{
		if (!keys[VK_SHIFT])
			camera.MoveForward( gDTime / 200.0f);
		else
			camera.MoveForward( gDTime / 20.0f);
	}

	if(keys['S'])
	{
		if (!keys[VK_SHIFT])
			camera.MoveForward(-gDTime / 200.0f);
		else
			camera.MoveForward(-gDTime / 20.0f);
	}

	if(keys['A'])
	{
		if (!keys[VK_SHIFT])
			camera.Strafe(-gDTime / 200.0f);
		else
			camera.Strafe(-gDTime / 20.0f);
	}

	if (keys['D'])
	{
		if (!keys[VK_SHIFT])
			camera.Strafe( gDTime / 200.0f);
		else
			camera.Strafe( gDTime / 20.0f);
	}

	if (keys[VK_PRIOR])
	{
		if (!keys[VK_SHIFT])
			camera.MoveUp( gDTime / 200.0f);
		else
			camera.MoveUp( gDTime / 20.0f);
	}

	if (keys[VK_NEXT] )
	{
		if (!keys[VK_SHIFT])
			camera.MoveUp(-gDTime / 200.0f);
		else
			camera.MoveUp(-gDTime / 20.0f);
	}//*/

	if (keys[VK_SPACE])
	{
		keys[VK_SPACE] = false;
		UpdateFEM = true;
	}

	if (keys['F'])
	{
		keys['F'] = false;

// 		demo.FEMTestSolver.m_NodeF[(demo.FEMTestSolver.NumNodes - 2)*3 + 0] += 100.0f;
//		demo.FEMTestSolver.m_NodeF[(demo.FEMTestSolver.NumNodes - 2)*3 + 2] += 200.0f;


//		demo.FEMTestSolver.m_NodeF_y[demo.FEMTestSolver.m_NumNodes - 2] += 20000.0f;
//		demo.FEMTestSolver.m_NodeF_y[demo.FEMTestSolver.m_NumNodes - 4] += 20000.0f;

//		demo.FEMTestSolver.m_NodeF_x[3] += 20000.0f;	// Rotational
		demo.FEMTestSolver.m_NodeF_y[4] += 20000.0f;	// Linear
	}

	if (keys['G'])
	{
		keys['G'] = false;

// 		demo.FEMTestSolver.m_NodeF[(demo.FEMTestSolver.NumNodes - 2)*3 + 0] -= 100.0f;
//		demo.FEMTestSolver.m_NodeF[(demo.FEMTestSolver.NumNodes - 2)*3 + 2] -= 200.0f;


// 		demo.FEMTestSolver.m_NodeF_y[demo.FEMTestSolver.m_NumNodes - 2] -= 20000.0f;
		demo.FEMTestSolver.m_NodeF_y[demo.FEMTestSolver.m_NumNodes - 4] -= 20000.0f;
	}

	if (keys['T'])
	{
		keys['T'] = false;

		// FAILED TWIST [Need to co-rotate twist vectors for correct twist]
/*
		float TwistMagn = 40000.0f;
		//		int idx = FEMTestSolver.NumNodes - 8*2 - 1; +Y
		//		int idx = FEMTestSolver.NumNodes - 8*2 - 2; +Z
		//		int idx = FEMTestSolver.NumNodes - 8*2 - 5; -Z
		//		int idx = FEMTestSolver.NumNodes - 8*2 - 6; -Y
		demo.FEMTestSolver.m_NodeF_y[demo.FEMTestSolver.m_NumNodes - 8*2 - 1] += TwistMagn;
		demo.FEMTestSolver.m_NodeF_z[demo.FEMTestSolver.m_NumNodes - 8*2 - 2] += TwistMagn;
		demo.FEMTestSolver.m_NodeF_z[demo.FEMTestSolver.m_NumNodes - 8*2 - 5] -= TwistMagn;
		demo.FEMTestSolver.m_NodeF_y[demo.FEMTestSolver.m_NumNodes - 8*2 - 6] -= TwistMagn;
*/
		demo.FEMTestSolver.m_NodeF_x[2] -= 100.0f;
	}


	if (keys['K'])
	{
		keys['K'] = false;

		if (demo.LambdaRender)
			demo.LambdaRender = 0;
		else
			demo.LambdaRender = 1;
	}

	if (keys['O'])
	{
		keys['O'] = false;

		demo.FEMTestSolver.CalcSylvesterCriterion();
	}

	if (keys['Q'])
	{
		keys['Q'] = false;

		if (keys[VK_SHIFT])
			gbFreeze = true;
		else
			gbFreeze = !gbFreeze;
	}
	if (keys['B'])
	{
		keys['B'] = false;

		gbStep = !gbStep;
	}

	if (keys['P'])
	{
		keys['P'] = false;

		std::vector<FEMJoint, GlobalAllocator<FEMJoint>>::iterator itJointEnd = demo.tetrajoints.end(), itJoint;
		for (itJoint = demo.tetrajoints.begin(); itJoint != itJointEnd; ++itJoint)
		{
			float NewYoung = itJoint->GetYoung() + 10000.0f;
			float NewPoisson = itJoint->GetPoisson();// + 0.1f;
			itJoint->SetYoungPoisson(NewYoung, (NewPoisson < 0.45f) ? NewPoisson : 0.45f);
		}
	}

	if (keys['L'])
	{
		keys['L'] = false;

		std::vector<FEMJoint, GlobalAllocator<FEMJoint>>::iterator itJointEnd = demo.tetrajoints.end(), itJoint;
		for (itJoint = demo.tetrajoints.begin(); itJoint != itJointEnd; ++itJoint)
		{
			float NewYoung = itJoint->GetYoung() - 10000.0f;
			float NewPoisson = itJoint->GetPoisson();// - 0.1f;
			itJoint->SetYoungPoisson((NewYoung > 1.0f) ? NewYoung : 1.0f, (NewPoisson > -0.45f) ? NewPoisson : -0.45f);
		}
	}

	return true;
}

// WndProc Body
LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg,			// Handle && Message for Window
						 WPARAM wParam, LPARAM lParam)	// Message parameters
{
	bool bHasFocus = true;
	switch (uMsg)
	{
		case WM_SYSCOMMAND:
		{
			switch (wParam)
			{
				case SC_SCREENSAVE:		// Don't let Screen Saver Start
				case SC_MONITORPOWER:	// Don't let Monitor enter in PowerSafe mode
				return 0;
			}
			break;
		}
		case WM_SIZE:	// Resizing Window
		{
			drawer.SetViewport(LOWORD(lParam), HIWORD(lParam), fov, znear, zfar);
			return 0;
		}
		case WM_CLOSE:	// Closing Window
		{
			PostQuitMessage(0);
			return 0;
		}
		case WM_KEYDOWN:	// Key Pressed
		{
			keys[wParam] = true;
			return 0;
		}
		case WM_KEYUP:		// Key Pressed
		{
			keys[wParam] = false;
			return 0;
		}
		case WM_KILLFOCUS:
		{
			hasFocus = false;
			break;
		}
		case WM_SETFOCUS:
		{
			hasFocus = true;
			break;
		}
	}
	
	// Pass unhandled messages to DefWindowProc
	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

// Function to calculate FPS
void __stdcall TimerProc(struct HWND__ *,unsigned int,unsigned int,unsigned long)
{
	char TempFPS[256];
	sprintf(TempFPS, "%s - [%04d] fps", szTitle, debugFPSCount);
	SetWindowText(hMainWindow, TempFPS);
	exFPSCount = debugFPSCount;
	debugFPSCount = 0;
}

void Render(demoscene &demo, float dt)
{
	drawer.BeginDraw();
	drawer.RenderFPS(exFPSCount);
	drawer.BeginScene(camera, gbWireframe);
	demo.Draw(drawer, dt);
	drawer.EndScene(camera, gbWireframe);
	drawer.EndDraw();
}

// Main Function
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
				   LPSTR lpCmdLine, int nShowCmd)
{
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	srand(GetTickCount());

// 	gLog2.SetLevel(LOGLEV_DEBUG);
// 	gLog2.Initialize("runlog.log");
// 	gLog2.LogDebug("test");

	unsigned int deviceID = strlen(lpCmdLine) ? atoi(lpCmdLine) : 0;

	// Ask user if he want to start fullscreen
	if (0)
	{
		if (MessageBox(NULL, "FullScreen [recommended] ?", "TEngine Question", MB_YESNO | MB_ICONQUESTION) == IDYES)
			gbFullScreen = true;
		else
			gbFullScreen = false;
	}
	else
	{
		gbFullScreen = false;
	}
	
	WNDCLASS wc;
	
	// Setting Window Class
	wc.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Rendraw on move and own DC for Window
	wc.lpfnWndProc   = (WNDPROC)WndProc;					// WndProc function handles messages
	wc.cbClsExtra    = 0;									// 
	wc.cbWndExtra    = 0;									// No extra Window data
	wc.hInstance     = hInstance;							// Instance
	wc.hIcon         = LoadIcon(NULL, IDI_WINLOGO);			// Load the Default Icon for Window
	wc.hCursor       = LoadCursor(NULL, IDC_ARROW);			// Load Arrow Pointer for Mouse
	wc.hbrBackground = NULL;								// No Background
	wc.lpszMenuName	 = NULL;								// No Menu
	wc.lpszClassName = "TEngine";							// Setting Class Name

	if (!RegisterClass(&wc))
	{
		MessageBox(NULL, "Failed To Register The Window Class.", "TEngine Error", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}
	
	RECT		WindowRect;		// Window rectangle

	// Setting Window rectangle
	WindowRect.left		= 0;
	WindowRect.right	= (long) windowX;
	WindowRect.top		= 0;
	WindowRect.bottom	= (long) windowY;

	DWORD dwStyle;
	// Setting Window Style
	if (gbFullScreen)
		dwStyle = WS_POPUP;
	else
		dwStyle = WS_CAPTION;	// We don't need our window to resize

	// Adjust Window to requested size
	AdjustWindowRect(&WindowRect, dwStyle, FALSE);

	long x = gbFullScreen ? 0 : (GetSystemMetrics(SM_CXSCREEN) - (WindowRect.right - WindowRect.left)) / 2;
	long y = gbFullScreen ? 0 : (GetSystemMetrics(SM_CYSCREEN) - (WindowRect.bottom - WindowRect.top)) / 2;

	HWND hwnd;
	// Trying to Create Window
	if (!(hwnd = CreateWindowEx(WS_EX_APPWINDOW,
							   "TEngine",
							   szTitle,
							   dwStyle,
							   x,
							   y,
							   WindowRect.right - WindowRect.left,
							   WindowRect.bottom - WindowRect.top,
							   NULL,
							   NULL,
							   hInstance,
							   NULL)))	
	{
		MessageBox(NULL, "Can't create window.", "TEngine Error", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}

	hMainWindow = hwnd;

	// Toggle Showing Cursor
	if (gbFullScreen)
		ShowCursor(false);
	else 
		ShowCursor(true);
	
	std::string error;

	if(!drawer.Init(hwnd, hInstance, gbFullScreen, windowX, windowY, error))
	{
		MessageBox(NULL, error.c_str(), "TEngine Error", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}
	
	ShowWindow(hwnd, SW_SHOW);	// Show Window
	SetForegroundWindow(hwnd);	// Higher priority
	SetFocus(hwnd);			// Set Focus to Window

	drawer.SetViewport(windowX, windowY, fov, znear, zfar);
	
	int middleX = GetSystemMetrics(SM_CXSCREEN) >> 1;
	int middleY = GetSystemMetrics(SM_CYSCREEN) >> 1;

	camera.SetScreenCenterCoords(middleX, middleY);
	camera.Setup(Vec3(0.0f, 5.0f, 15.0f), Vec3(0.0f, 0.0f, 0.0f));
	
	SetWindowText(hwnd, szTitle);
	
	// Setting Cursor Position in centre Window
	SetCursorPos(middleX, middleY);

	// Setting Timer for FPS counting
	SetTimer(hwnd, FPS_TIMER, 1000, TimerProc);

	DWORD PrgBegin;

	// TimeSync
	PrgBegin = GetTickCount();

	int ElapsedTime = 0;					// Simple Time Syncronization
	
	MSG msg; // Message structure
	bool bLoop = true;

	demo.Init(deviceID);
	gLog.SetEnabled(false);
	gLog2.SetEnabled(false);
	drawer.SetVSyncState(true);
	while (bLoop)
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			// Handles some more messages
			switch (msg.message) 
			{
				case WM_QUIT:
					bLoop = true;
				    break;
			    default:
					TranslateMessage(&msg);
					DispatchMessage(&msg);
			        break;
			};			
		}
		else
		{
			debugFPSCount++;							// For FPS counting
		
			int prevTime = ElapsedTime;
			ElapsedTime = GetTickCount() - PrgBegin;	// Some Time Synchronization

			gDTime = (float) ElapsedTime - prevTime;

			if (hasFocus)
			{
				POINT mousemPos;
				GetCursorPos(&mousemPos);
				camera.OnMouse(mousemPos.x, mousemPos.y);
				camera.GetSavedMouseCoords(mousemPos.x, mousemPos.y);
				SetCursorPos(mousemPos.x, mousemPos.y);
			}

			// We should check Esc key here, not in ProcessKeys function
			if (!ProcessKeys())
				bLoop = false;
			Render(demo, gDTime);								// Render Current Scene
		}
	}
	
	drawer.Destroy();		// Destroy Scene

	return (int)msg.wParam;
}
