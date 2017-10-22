#ifndef _DYNAMICS_JOINTS_H_
#define _DYNAMICS_JOINTS_H_

#include "solver_base.h"


#define PLANE_FRICTION 1


class IJoint
{
public:

	virtual ~IJoint()
	{
	}

	/**
	This method is called each frame before the joint is used for solving in a particular solver.
	The method must prepare the Jacobians, CFMs, ERPs, low and high bounds, and lambda initial guesses
	for each mathematical constraint of this joint and append this data to the end of the arrays of the solver.

	\param dt 		the length of the time step being processed
	\param solver 	solver that must be used
	*/
	virtual void UpdateCopy(float dt, SolverBase &solver) = 0;

	/**
	This method is called each frame after solving the system. It must store computed lambdas into internal
	storage, so it can restore it later.		

	\param solver 	solver that must used for computations
	*/
	virtual void FetchLambdas(SolverBase &solver) = 0;
};



class BallSocket: public IJoint
{
public:

	/**
	This method is initializes Joint and calculates internal data

	\param solver 			solver that must be used
	\param ERP				error reduction parameter (stabilization)
	\param CFM 				constraint force mixing parameter (regularization)
	\param AnchorPoint	 	coordinates of the pin-point in world space
	\param Body1LIdx		body1 linear index
	\param Body1RIdx		body1 rotational index
	\param Body2LIdx		body2 linear index (should be -1 if attachment to world)
	\param Body2RIdx		body2 rotational index
	*/
	virtual void Init(SolverBase &solver, float ERP, float CFM, const Vec3 &AnchorPoint, int Body1LIdx, int Body1RIdx, int Body2LIdx, int Body2RIdx);
	virtual void UpdateCopy(float dt, SolverBase &solver);

	virtual void FetchLambdas(SolverBase &solver);

	// Temporary function for rendering
	void GetIndices(int &idx1_lin, int&idx1_ang, int &idx2_lin, int&idx2_ang);
	void GetAnchorPoints_BS(Vec3 &bsAnchor1, Vec3 &bsAnchor2);

private:

	unsigned int m_StartIdx;
	float m_CFM[3];			//store precomputed CFM (before multiplication)
	float m_ERP;

	float m_lambda0[3];

	Vec3 m_BodyAnchorPoints[2];	// BallSocket Anchors, body-space

	// TRIPLES INDICES
	// Shifts of Bodies in solver [ shift == -1 for Body2 linear means attachment to world ]
	int m_Body1LIdx, m_Body1RIdx;
	int m_Body2LIdx, m_Body2RIdx;
};


#if (PLANE_FRICTION == 1)

class PlaneConstraint: public IJoint
{
public:

	/**
	This method is initializes Joint and calculates internal data

	\param solver 			solver that must be used
	\param ERP				error reduction parameter (stabilization)
	\param CFM 				constraint force mixing parameter (regularization)
	\param PlanePoint	 	coordinates of any point on the plane
	\param PlaneNormal	 	coordinates of plane's normal vector
	\param BodyIdx			FE-point linear index
	\param Tolerance		acceptable amount of penetration
	*/
	virtual void Init(SolverBase &solver, float ERP, float CFM, const Vec3 &PlanePoint, const Vec3 &PlaneNormal, int BodyIdx, float Tolerance);
	virtual void UpdateCopy(float dt, SolverBase &solver);

	virtual void FetchLambdas(SolverBase &solver);


	// MEMBERS PUBLIC FOR DEBUG ONLY
	// !!!!!!!!!!DO NOT USE!!!!!!!!!!!!
	unsigned int m_StartIdx;
	bool m_bActive;
	float m_lambda0[3];

	// TRIPLES INDICES
	// Shifts of Bodies in solver
	int m_BodyIdx;

	// Plane equation vectors
	Vec3 m_PlaneNormal;
	float m_PlaneD;

private:

	float m_CFM[3];				// Store precomputed CFM (before multiplication)
	float m_ERP[3];
	float m_Tolerance;
};

//////////////////////////////////////////////////////////////////////////
#else
//////////////////////////////////////////////////////////////////////////

class PlaneConstraint: public IJoint
{
public:

	/**
	This method is initializes Joint and calculates internal data

	\param solver 			solver that must be used
	\param ERP				error reduction parameter (stabilization)
	\param CFM 				constraint force mixing parameter (regularization)
	\param PlanePoint	 	coordinates of any point on the plane
	\param PlaneNormal	 	coordinates of plane's normal vector
	\param BodyIdx			FE-point linear index
	\param Tolerance		acceptable amount of penetration
	*/
	virtual void Init(SolverBase &solver, float ERP, float CFM, const Vec3 &PlanePoint, const Vec3 &PlaneNormal, int BodyIdx, float Tolerance);
	virtual void UpdateCopy(float dt, SolverBase &solver);

	virtual void FetchLambdas(SolverBase &solver);


	// MEMBERS PUBLIC FOR DEBUG ONLY
	// !!!!!!!!!!DO NOT USE!!!!!!!!!!!!
	unsigned int m_StartIdx;
	bool m_bActive;
	float m_lambda0[1];

	// TRIPLES INDICES
	// Shifts of Bodies in solver
	int m_BodyIdx;

	// Plane equation vectors
	Vec3 m_PlaneNormal;
	float m_PlaneD;

private:

	float m_CFM[1];				// Store precomputed CFM (before multiplication)
	float m_ERP[1];
	float m_Tolerance;
};

#endif

#endif