#ifndef _FEM_JOINTS_H_
#define _FEM_JOINTS_H_

#include "dynamics/joints.h"

class FEMJoint: public IJoint
{
	friend class SolverLCPCG_CUDA;
	friend class SolverJacobi_CUDA;

public:

	/**
	This method is initializes Joint and calculates internal data

	\param solver 				solver that must be used
	\param NodePos0				vector of initial node positions (needed for calculations)
	\param Node1Idx				Index of the 1st node in Tetrahedral-FE
	\param Node2Idx				Index of the 2nd node in Tetrahedral-FE
	\param Node3Idx				Index of the 3rd node in Tetrahedral-FE
	\param Node4Idx				Index of the 4th node in Tetrahedral-FE
	\param Young				Young's modulus of FE (stiffness)
	\param Poisson				Poisson's ratio (transverse strain to axial strain -- perp. contraction to extension)
	\param Yield				Threshold from where strain goes to plastic creep
	\param Creep				Amount of elastic strain, absorbed by plasticity [0; 1]
	\param MaxPlasticStrain		Maximum amount of strain, absorbed by plasticity
	\param Damping				Oscillation reducing coefficient
	\param Regularization 		constraint force mixing parameter
	*/
	virtual void Init(SolverBase &solver, const std::vector<Vec3, GlobalAllocator<Vec3>> &NodePos0,
		unsigned int Node1Idx, unsigned int Node2Idx, unsigned int Node3Idx, unsigned int Node4Idx,
		float Young, float Poisson, float Yield, float Creep,
		float MaxPlasticStrain, float Damping, float Regularization);

	virtual void UpdateCopy(float dt, SolverBase &solver);

	virtual void FetchLambdas(SolverBase &solver);

	void SetYoungPoisson(float Young, float Poisson);

	float GetYoung() const
	{
		return m_Young;
	}

	float GetPoisson() const
	{
		return m_Poisson;
	}

	float GetYield() const
	{
		return m_Yield;
	}
	void SetYield(float Yield)
	{
		m_Yield = Yield;
	}

	float GetCreep() const
	{
		return m_Creep;
	}
	void SetCreep(float Creep)
	{
		m_Creep = Creep;
	}

	float GetMaxPlasticStrain() const
	{
		return m_MaxPlasticStrain;
	}
	void SetMaxPlasticStrain(float MaxPlasticStrain)
	{
		m_MaxPlasticStrain = MaxPlasticStrain;
	}

	float GetDamping() const
	{
		return m_Damping;
	}
	void SetDamping(float Damping)
	{
		m_Damping = Damping;
	}

	float GetRegularization() const
	{
		return m_Regularization;
	}
	void SetRegularization(float Regularization)
	{
		m_Regularization = Regularization;
	}



	unsigned int m_fetchIndices[6];

	// MEMBERS PUBLIC FOR DEBUG ONLY
	// !!!!!!!!!!DO NOT USE!!!!!!!!!!!!
	float m_lambda0[6];

	// Norm of Elastic Strain for debug purposes [Should be Getter for this?]
	float m_StrainNorm;

private:

	unsigned int m_StartIdx;
	float m_Regularization;
	float m_CFM[6];						// Store precomputed CFM (before multiplication)
	float m_B_loc[6 * 12];

	float m_Damping;					// "beta" parameter

	// !!!! This variable is used ONLY FOR SetYoungPoisson method
	float m_InitialFEVolume;

	// Elasticity parameters
	float m_Young;
	float m_Poisson;

	// Plasticity parameters
	float m_Yield;						// Threshold from where strain goes to plastic creep
	float m_Creep;						// Amount of elastic strain, absorbed by plasticity [0; 1]
	float m_MaxPlasticStrain;			// Maximum amount of "saved" plastic strain for FE

#if (WARPING_POLAR_DECOMPOSITION == 1)
	CMatrix3 m_mp0_inv;					// Matrix of tetrahedron positions;
#else
	CMatrix3 m_N;						// Matrix of original FE basis
#endif

	float m_Jp0[6];						// J * pos0;

	// TRIPLES INDICES
	unsigned int m_Node1Idx, m_Node2Idx, m_Node3Idx, m_Node4Idx;

	float m_E_plastic[6];
};

class BallSocket_FEM: public IJoint
{
public:

	/**
	This method is initializes Joint and calculates internal data

	\param solver 			solver that must be used
	\param ERP				error reduction parameter (stabilization)
	\param CFM 				constraint force mixing parameter (regularization)
	\param AnchorPoint	 	coordinates of the pin-point in world space
	\param Body1LIdx		body1 linear index (should be -1 if attachment to world)
	\param Body1RIdx		body1 rotational index
	\param Node1Idx			FE-point1 linear index
	\param Node2Idx			FE-point2 linear index
	\param Node3Idx			FE-point3 linear index
	*/
	virtual void Init(SolverBase &solver, float ERP, float CFM, const Vec3 &AnchorPoint, int Body1LIdx, int Body1RIdx, int Node1Idx, int Node2Idx, int Node3Idx);
	virtual void UpdateCopy(float dt, SolverBase &solver);

	virtual void FetchLambdas(SolverBase &solver);

private:

	unsigned int m_StartIdx;
	float m_CFM[3];				// Store precomputed CFM (before multiplication)
	float m_ERP;

	float m_lambda0[3];

	// Barycentric coordinates
	float m_Alpha, m_Beta, m_Gamma, m_Delta;

	Vec3 m_BodyAnchorPoints[2];	// BallSocket Anchors, body-space

	// TRIPLES INDICES
	// Shifts of Bodies in solver
	int m_Body1LIdx, m_Body1RIdx;
	// Shift of Nodes in solver
	int m_Node1, m_Node2, m_Node3;
};

#endif