#ifndef _DYNAMICS_SOLVER_BASE_H_
#define _DYNAMICS_SOLVER_BASE_H_

#include "all.h"
#include <vector>
#include "FEM/FEM.h"
#include "Math/m_MathFunctions.h"
#include "Math/m_Vector3.h"
#include "math/m_Quaternion.h"
#include <float.h>

template <typename T>
static size_t getAlign(const T & val) { static_assert(false); }

template <> static size_t getAlign(const float & val) { return 4; }
template <> static size_t getAlign(const int & val) { return 4; }
template <> static size_t getAlign(const unsigned int & val) { return 4; }

#define GET_UPPER_ALIGNMENT_BOUND(offset, data) \
			((offset) + getAlign(data) - 1) & ~(getAlign(data) - 1)

const int MAX_NUM_NODES		= 20000;
const int MAX_NUM_JACOBI	= 40000;

const float PLANE_THRESHOLD = 0.00001f;
const float GRAVITY_THRESHOLD = 0.000001f;

// Uncomment to enable cycle multiplication (and enables corresponding arrays)
// !!!!!!!!!!!!!!!!!!!!!!!!!!DO NOT USE WITH SOA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define JOINT_TRIPLES_CYCLE 0
// Number of maximum triples per joint (triple per movement part, e.g. linear/rotational)
#define NUM_JOINT_TRIPLES 6

#if (JOINT_TRIPLES_CYCLE == 1)

#error CAN'T USE CYCLES WITH STRUCTURE OF ARRAYS

#endif

// Additional calculations for testing purposes (gradient/residual vector)
#define SOLVERS_ANALYZE_RESIDUAL 1
// Enable stopping criteria in solvers (w/o it iterations are fixed)
#define SOLVER_STOPPING_CRITERIA 1

// [ DEBUG ] Output objective function to log (CG/LCPCG only)
#define SOLVER_OBJECTIVE_FUNCTION_OUTPUT 0


// Switches Stiffness Warping basis calculation to Polar Decomposition
// (turned off - calculation based on FE edges, Mueller's approach)
#define WARPING_POLAR_DECOMPOSITION 0

#define inc_arr(arr) (*(ptr_##arr)++)
#define val_arr(arr) (*(ptr_##arr))



// Common macros for solvers

//////////////////////////////////////////////////////////////////////////
// INIT A SPARSE
//////////////////////////////////////////////////////////////////////////
#if (NUM_JOINT_TRIPLES == 4)

#define INIT_A_SPARSE(idx)\
		ptr_m_Asp_00 = &m_Asp_00[(idx)];\
		ptr_m_Asp_01 = &m_Asp_01[(idx)];\
		ptr_m_Asp_02 = &m_Asp_02[(idx)];\
		ptr_m_Asp_03 = &m_Asp_03[(idx)];\
		ptr_m_Asp_04 = &m_Asp_04[(idx)];\
		ptr_m_Asp_05 = &m_Asp_05[(idx)];\
		ptr_m_Asp_06 = &m_Asp_06[(idx)];\
		ptr_m_Asp_07 = &m_Asp_07[(idx)];\
		ptr_m_Asp_08 = &m_Asp_08[(idx)];\
		ptr_m_Asp_09 = &m_Asp_09[(idx)];\
		ptr_m_Asp_10 = &m_Asp_10[(idx)];\
		ptr_m_Asp_11 = &m_Asp_11[(idx)];

#elif (NUM_JOINT_TRIPLES == 6)

#define INIT_A_SPARSE(idx)\
		ptr_m_Asp_00 = &m_Asp_00[(idx)];\
		ptr_m_Asp_01 = &m_Asp_01[(idx)];\
		ptr_m_Asp_02 = &m_Asp_02[(idx)];\
		ptr_m_Asp_03 = &m_Asp_03[(idx)];\
		ptr_m_Asp_04 = &m_Asp_04[(idx)];\
		ptr_m_Asp_05 = &m_Asp_05[(idx)];\
		ptr_m_Asp_06 = &m_Asp_06[(idx)];\
		ptr_m_Asp_07 = &m_Asp_07[(idx)];\
		ptr_m_Asp_08 = &m_Asp_08[(idx)];\
		ptr_m_Asp_09 = &m_Asp_09[(idx)];\
		ptr_m_Asp_10 = &m_Asp_10[(idx)];\
		ptr_m_Asp_11 = &m_Asp_11[(idx)];\
		ptr_m_Asp_12 = &m_Asp_12[(idx)];\
		ptr_m_Asp_13 = &m_Asp_13[(idx)];\
		ptr_m_Asp_14 = &m_Asp_14[(idx)];\
		ptr_m_Asp_15 = &m_Asp_15[(idx)];\
		ptr_m_Asp_16 = &m_Asp_16[(idx)];\
		ptr_m_Asp_17 = &m_Asp_17[(idx)];

#endif

//////////////////////////////////////////////////////////////////////////
// INIT J SPARSE
//////////////////////////////////////////////////////////////////////////
#if (NUM_JOINT_TRIPLES == 4)

#define INIT_J_SPARSE(idx)\
	ptr_m_J_00 = &m_J_00[(idx)];\
	ptr_m_J_01 = &m_J_01[(idx)];\
	ptr_m_J_02 = &m_J_02[(idx)];\
	ptr_m_J_03 = &m_J_03[(idx)];\
	ptr_m_J_04 = &m_J_04[(idx)];\
	ptr_m_J_05 = &m_J_05[(idx)];\
	ptr_m_J_06 = &m_J_06[(idx)];\
	ptr_m_J_07 = &m_J_07[(idx)];\
	ptr_m_J_08 = &m_J_08[(idx)];\
	ptr_m_J_09 = &m_J_09[(idx)];\
	ptr_m_J_10 = &m_J_10[(idx)];\
	ptr_m_J_11 = &m_J_11[(idx)];

#elif (NUM_JOINT_TRIPLES == 6)

#define INIT_J_SPARSE(idx)\
	ptr_m_J_00 = &m_J_00[(idx)];\
	ptr_m_J_01 = &m_J_01[(idx)];\
	ptr_m_J_02 = &m_J_02[(idx)];\
	ptr_m_J_03 = &m_J_03[(idx)];\
	ptr_m_J_04 = &m_J_04[(idx)];\
	ptr_m_J_05 = &m_J_05[(idx)];\
	ptr_m_J_06 = &m_J_06[(idx)];\
	ptr_m_J_07 = &m_J_07[(idx)];\
	ptr_m_J_08 = &m_J_08[(idx)];\
	ptr_m_J_09 = &m_J_09[(idx)];\
	ptr_m_J_10 = &m_J_10[(idx)];\
	ptr_m_J_11 = &m_J_11[(idx)];\
	ptr_m_J_12 = &m_J_12[(idx)];\
	ptr_m_J_13 = &m_J_13[(idx)];\
	ptr_m_J_14 = &m_J_14[(idx)];\
	ptr_m_J_15 = &m_J_15[(idx)];\
	ptr_m_J_16 = &m_J_16[(idx)];\
	ptr_m_J_17 = &m_J_17[(idx)];

#endif

//////////////////////////////////////////////////////////////////////////
// INIT JOINTED NODES INDICES
//////////////////////////////////////////////////////////////////////////
#if (NUM_JOINT_TRIPLES == 4)

#define INIT_NODES_IDX(idx)\
	ptr_m_JtdNodes_00 = &m_JtdNodes_00[(idx)];\
	ptr_m_JtdNodes_01 = &m_JtdNodes_01[(idx)];\
	ptr_m_JtdNodes_02 = &m_JtdNodes_02[(idx)];\
	ptr_m_JtdNodes_03 = &m_JtdNodes_03[(idx)];

#elif (NUM_JOINT_TRIPLES == 6)

#define INIT_NODES_IDX(idx)\
	ptr_m_JtdNodes_00 = &m_JtdNodes_00[(idx)];\
	ptr_m_JtdNodes_01 = &m_JtdNodes_01[(idx)];\
	ptr_m_JtdNodes_02 = &m_JtdNodes_02[(idx)];\
	ptr_m_JtdNodes_03 = &m_JtdNodes_03[(idx)];\
	ptr_m_JtdNodes_04 = &m_JtdNodes_04[(idx)];\
	ptr_m_JtdNodes_05 = &m_JtdNodes_05[(idx)];

#endif

//////////////////////////////////////////////////////////////////////////
// CALCULATE VECTOR A
//////////////////////////////////////////////////////////////////////////
#define CALC_VECTOR_A(vec, j, t0, t1, t2)\
			idx = inc_arr( m_JtdNodes_0##j );\
			\
			ptr_m_a_x[idx] += inc_arr(m_Asp_##t0) * val_arr(vec);\
			ptr_m_a_y[idx] += inc_arr(m_Asp_##t1) * val_arr(vec);\
			ptr_m_a_z[idx] += inc_arr(m_Asp_##t2) * val_arr(vec);


//////////////////////////////////////////////////////////////////////////
// CALCULATE VECTOR A -- JACOBIAN PRODUCT
//////////////////////////////////////////////////////////////////////////
#if (NUM_JOINT_TRIPLES == 4)

#define MUL_J_A() \
	inc_arr(m_J_00) * ptr_m_a_x[t0Idx] + inc_arr(m_J_01) * ptr_m_a_y[t0Idx] + inc_arr(m_J_02) * ptr_m_a_z[t0Idx] + \
	inc_arr(m_J_03) * ptr_m_a_x[t1Idx] + inc_arr(m_J_04) * ptr_m_a_y[t1Idx] + inc_arr(m_J_05) * ptr_m_a_z[t1Idx] + \
	inc_arr(m_J_06) * ptr_m_a_x[t2Idx] + inc_arr(m_J_07) * ptr_m_a_y[t2Idx] + inc_arr(m_J_08) * ptr_m_a_z[t2Idx] + \
	inc_arr(m_J_09) * ptr_m_a_x[t3Idx] + inc_arr(m_J_10) * ptr_m_a_y[t3Idx] + inc_arr(m_J_11) * ptr_m_a_z[t3Idx];

#elif (NUM_JOINT_TRIPLES == 6)

#define MUL_J_A() \
	inc_arr(m_J_00) * ptr_m_a_x[t0Idx] + inc_arr(m_J_01) * ptr_m_a_y[t0Idx] + inc_arr(m_J_02) * ptr_m_a_z[t0Idx] + \
	inc_arr(m_J_03) * ptr_m_a_x[t1Idx] + inc_arr(m_J_04) * ptr_m_a_y[t1Idx] + inc_arr(m_J_05) * ptr_m_a_z[t1Idx] + \
	inc_arr(m_J_06) * ptr_m_a_x[t2Idx] + inc_arr(m_J_07) * ptr_m_a_y[t2Idx] + inc_arr(m_J_08) * ptr_m_a_z[t2Idx] + \
	inc_arr(m_J_09) * ptr_m_a_x[t3Idx] + inc_arr(m_J_10) * ptr_m_a_y[t3Idx] + inc_arr(m_J_11) * ptr_m_a_z[t3Idx] + \
	inc_arr(m_J_12) * ptr_m_a_x[t4Idx] + inc_arr(m_J_13) * ptr_m_a_y[t4Idx] + inc_arr(m_J_14) * ptr_m_a_z[t4Idx] + \
	inc_arr(m_J_15) * ptr_m_a_x[t5Idx] + inc_arr(m_J_16) * ptr_m_a_y[t5Idx] + inc_arr(m_J_17) * ptr_m_a_z[t5Idx];

#endif


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
// TEMPORARY CRAP
//////////////////////////////////////////////////////////////////////////
template<typename T, unsigned int S>
class dynarray
{
public:
	dynarray()
	{
		m_pPointer = new T[S];
	}

	dynarray(const dynarray<T, S> &ref)
	{
		m_pPointer = new T[S];
		memcpy(m_pPointer, ref.m_pPointer, S * sizeof(T));
	}

	const dynarray<T, S> & operator =(const dynarray<T, S> &ref)
	{
		if(&ref != this)
		{
			m_pPointer = new T[S];
			memcpy(m_pPointer, ref.m_pPointer, S * sizeof(T));
		}

		return *this;
	}

	~dynarray()
	{
		if (m_pPointer)
			delete[] m_pPointer;
	}

	T& operator[](unsigned int idx)
	{
		return m_pPointer[idx];
	}

	const T& operator[](unsigned int idx) const
	{
		return m_pPointer[idx];
	}

	void Empty()
	{
		memset(m_pPointer, 0, S * sizeof(T));
	}

	void Empty(unsigned int size)
	{
		memset(m_pPointer, 0, size * sizeof(T));
	}

private:
	T *m_pPointer;
};

//////////////////////////////////////////////////////////////////////////



class SolverBase
{
	// list all joints here; in future make some kind of interface for SOlver tompass to joints

	friend class FEMJoint;	
	friend class BallSocket;
	friend class BallSocket_FEM;
	friend class PlaneConstraint;

public:

	/**
	This method adds rotational node to the system

	\param invInertia0 		inversed inertia tensor
	\param orient			quaternion of initial orientation
	\param angvel 			initial angular velocity
	*/
	unsigned int AddRotationalNode(const CMatrix3 &invInertia0, const CQuaternion &orient, const Vec3 &angvel);
	/**
	This method adds linear node to the system

	\param invInertia0 		inversed mass
	\param orient			initial position
	\param angvel 			initial velocity
	*/
	unsigned int AddTranslationalNode(float invMass, const Vec3 &pos, const Vec3 &vel);

	void SetForce(unsigned int node, float x, float y, float z);
	void GetForce(unsigned int node, float &x, float &y, float &z) const;

	float CalcSylvesterCriterion();

	inline Vec3 GetPosition(unsigned int NodeIdx) const
	{
		assert(NodeIdx < m_NumNodes);

		return Vec3(m_NodePosRot_x[NodeIdx], m_NodePosRot_y[NodeIdx], m_NodePosRot_z[NodeIdx]);
	}
	void SetPosition(unsigned int NodeIdx, const Vec3 &Pos);

	inline CQuaternion GetOrientation(unsigned int NodeIdx) const
	{
		assert(NodeIdx < m_NumNodes);

		return CQuaternion(	m_NodePosRot_w[NodeIdx], m_NodePosRot_x[NodeIdx],
							m_NodePosRot_y[NodeIdx], m_NodePosRot_z[NodeIdx] );
	}
	void SetOrientation(unsigned int NodeIdx, const CQuaternion &Quat);

	unsigned int GetEffectiveIterations() const
	{
		return m_EffectiveIterations;
	}

	void SetMaxIterations(unsigned int MaxIterations)
	{
		m_Iterations = MaxIterations;
	}
	unsigned int GetMaxIterations() const
	{
		return m_Iterations;
	}

	void SetPrecision(float Precision)
	{
		m_Precision = Precision;
	}
	float GetPrecision(void) const
	{
		return m_Precision;
	}

	void SetInitialGuessCutFlag(bool Flag)
	{
		m_IG_CutFlag = Flag;
	}
	bool GetInitialGuessCutFlag(void) const
	{
		return m_IG_CutFlag;
	}

	void SetInitialGuessCutThreshold(float Threshold)
	{
		m_IG_CutThreshold = Threshold;
	}
	float GetInitialGuessCutThreshold(void) const
	{
		return m_IG_CutThreshold;
	}

	void SetGravity(const Vec3 &Gravity)
	{
		m_Gravity = Gravity;
	}
	Vec3 GetGravity() const
	{
		return m_Gravity;
	}

#if (SOLVERS_ANALYZE_RESIDUAL == 1)
	// Interface for Analysis
	float GetLambdaNormSq() const
	{
		return m_LambdaNormSq;
	}

	float GetGradNormSq() const
	{
		return m_GradNormSq;
	}

	float GetDotLambdaGrad() const
	{
		return m_DotLambdaGrad;
	}
#else
	// Interface for Analysis [EMPTY]
	float GetLambdaNormSq() const
	{
		return 0.0f;
	}

	float GetGradNormSq() const
	{
		return 0.0f;
	}

	float GetDotLambdaGrad() const
	{
		return 0.0f;
	}
#endif

//protected:
	unsigned int m_NumNodes, m_NumJoints;

#if (NUM_JOINT_TRIPLES == 4)

	// Matrix of derivatives [dC/dx]
	dynarray<float, MAX_NUM_JACOBI> m_J_00, m_J_01, m_J_02, m_J_03, m_J_04, m_J_05,
									m_J_06, m_J_07, m_J_08, m_J_09, m_J_10, m_J_11;

#elif (NUM_JOINT_TRIPLES == 6)

	// Matrix of derivatives [dC/dx]
	dynarray<float, MAX_NUM_JACOBI> m_J_00, m_J_01, m_J_02, m_J_03, m_J_04, m_J_05,
									m_J_06, m_J_07, m_J_08, m_J_09, m_J_10, m_J_11,
									m_J_12, m_J_13, m_J_14, m_J_15, m_J_16, m_J_17;

#endif

	dynarray<float, MAX_NUM_JACOBI> m_J_rhs;								// dC / dt - J * (V + dt * M^-1 * Fext)

#if (NUM_JOINT_TRIPLES == 4)

	// Asp = dt * M^-1 * J^T
	dynarray<float, MAX_NUM_JACOBI> m_Asp_00, m_Asp_01, m_Asp_02, m_Asp_03, m_Asp_04, m_Asp_05,
									m_Asp_06, m_Asp_07, m_Asp_08, m_Asp_09, m_Asp_10, m_Asp_11;

#elif (NUM_JOINT_TRIPLES == 6)

	// Asp = dt * M^-1 * J^T
	dynarray<float, MAX_NUM_JACOBI> m_Asp_00, m_Asp_01, m_Asp_02, m_Asp_03, m_Asp_04, m_Asp_05,
									m_Asp_06, m_Asp_07, m_Asp_08, m_Asp_09, m_Asp_10, m_Asp_11,
									m_Asp_12, m_Asp_13, m_Asp_14, m_Asp_15, m_Asp_16, m_Asp_17;

#endif
	dynarray<float, MAX_NUM_JACOBI> m_CFM;									// Constraint Force Mixing

#if (NUM_JOINT_TRIPLES == 4)

	// Jointed nodes indices
	dynarray<unsigned int, MAX_NUM_JACOBI>	m_JtdNodes_00, m_JtdNodes_01, m_JtdNodes_02, m_JtdNodes_03;

#elif (NUM_JOINT_TRIPLES == 6)
	
	// Jointed nodes indices
	dynarray<unsigned int, MAX_NUM_JACOBI>	m_JtdNodes_00, m_JtdNodes_01, m_JtdNodes_02,
											m_JtdNodes_03, m_JtdNodes_04, m_JtdNodes_05;

#endif

	dynarray<float, MAX_NUM_JACOBI> m_lambda, m_resid;						// Lambda, Initial Guess [lambda from prev solve]
	dynarray<float, MAX_NUM_JACOBI> m_RHS;									// Right Hand Side, Ax = b <- THIS

	dynarray<float, MAX_NUM_JACOBI> m_Lo, m_Hi;								// Lo && Hi limits



	dynarray<bool, MAX_NUM_NODES> m_IsRotational;


	dynarray<float, MAX_NUM_NODES> m_NodePosRot_x, m_NodePosRot_y, m_NodePosRot_z, m_NodePosRot_w;

	dynarray<float, MAX_NUM_NODES>	m_NodeInvMass0_00, m_NodeInvMass0_01, m_NodeInvMass0_02,
									m_NodeInvMass0_10, m_NodeInvMass0_11, m_NodeInvMass0_12,
									m_NodeInvMass0_20, m_NodeInvMass0_21, m_NodeInvMass0_22;

	dynarray<float, MAX_NUM_NODES>	m_NodeInvMass_00, m_NodeInvMass_01, m_NodeInvMass_02,
									m_NodeInvMass_10, m_NodeInvMass_11, m_NodeInvMass_12,
									m_NodeInvMass_20, m_NodeInvMass_21, m_NodeInvMass_22;

	dynarray<float, MAX_NUM_NODES> m_NodeF_x, m_NodeF_y, m_NodeF_z;


	dynarray<float, MAX_NUM_JACOBI> m_invDiag;							// inverse diagonal

	// total force V + dt * M^-1 * Fext - M^-1 * K(x - x0) [???]
	dynarray<float, MAX_NUM_NODES> m_Ftot_x, m_Ftot_y, m_Ftot_z;
	dynarray<float, MAX_NUM_NODES> m_NodeVel_x, m_NodeVel_y, m_NodeVel_z;

	// Odd-job vector
	// PGS: a = Asp * lambda = dt * M^-1 * J^T * lambda
	dynarray<float, MAX_NUM_NODES> m_a_x, m_a_y, m_a_z;

#if (JOINT_TRIPLES_CYCLE == 1)
	dynarray<unsigned int, MAX_NUM_JACOBI> m_JointTriples;
#endif


protected:

	SolverBase(unsigned int MaxIterations, float Precision);

	// Common solver pre-step: calculating Ftotal = V + dt * M^-1 * Fext
	void Prestep(float dt);

	// Symplectic Euler integration, based on solver's output
	void Integrate(float dt);

	// Cut initial guess when needed
	void CutLambdas();

#if (SOLVERS_ANALYZE_RESIDUAL == 1)
	// Solver Analysis
	float m_LambdaNormSq, m_GradNormSq, m_DotLambdaGrad;
#endif

	// Used for statistics measurement
	unsigned int m_EffectiveIterations;

	bool m_IG_CutFlag;
	float m_IG_CutThreshold;

	// Precision of solver's convergence [ if ||residual|| / NumJoints < Precision ==> STOP ]
	float m_Precision;

	// Maximum number of iterations
	unsigned int m_Iterations;

	Vec3 m_Gravity;
};

#endif