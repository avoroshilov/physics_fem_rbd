#include "joints.h"

//////////////////////////////////////////////////////////////////////////
// BallSocket
//////////////////////////////////////////////////////////////////////////

void BallSocket::Init(SolverBase &solver, float ERP, float CFM, const Vec3 &AnchorPoint, int Body1LIdx, int Body1RIdx, int Body2LIdx, int Body2RIdx)
{
	memset(m_lambda0, 0, 3 * sizeof(float));

	m_Body1LIdx = Body1LIdx;
	m_Body1RIdx = Body1RIdx;
	m_Body2LIdx = Body2LIdx;
	m_Body2RIdx = Body2RIdx;

	m_CFM[0] = CFM;
	m_CFM[1] = CFM;
	m_CFM[2] = CFM;

	m_ERP = ERP;

	CMatrix3 Body1Rot, Body2Rot;
	CQuaternion tmpQuaterion;

	// Retrieving Bodies' Orientations
	tmpQuaterion.w = solver.m_NodePosRot_w[m_Body1RIdx];
	tmpQuaterion.x = solver.m_NodePosRot_x[m_Body1RIdx];
	tmpQuaterion.y = solver.m_NodePosRot_y[m_Body1RIdx];
	tmpQuaterion.z = solver.m_NodePosRot_z[m_Body1RIdx];
	tmpQuaterion.ToMatrix3(Body1Rot);

	if (Body2LIdx != -1)
	{
		tmpQuaterion.w = solver.m_NodePosRot_w[m_Body2RIdx];
		tmpQuaterion.x = solver.m_NodePosRot_x[m_Body2RIdx];
		tmpQuaterion.y = solver.m_NodePosRot_y[m_Body2RIdx];
		tmpQuaterion.z = solver.m_NodePosRot_z[m_Body2RIdx];
		tmpQuaterion.ToMatrix3(Body2Rot);
	}

	Vec3 BodyPos;

	// Calculating Anchor Points in Local Space
	BodyPos.x = solver.m_NodePosRot_x[m_Body1LIdx];
	BodyPos.y = solver.m_NodePosRot_y[m_Body1LIdx];
	BodyPos.z = solver.m_NodePosRot_z[m_Body1LIdx];
	m_BodyAnchorPoints[0] = Body1Rot.GetTransposed() * (AnchorPoint - BodyPos);

	if (Body2LIdx != -1)
	{
		BodyPos.x = solver.m_NodePosRot_x[m_Body2LIdx];
		BodyPos.y = solver.m_NodePosRot_y[m_Body2LIdx];
		BodyPos.z = solver.m_NodePosRot_z[m_Body2LIdx];
		m_BodyAnchorPoints[1] = Body2Rot.GetTransposed() * (AnchorPoint - BodyPos);
	}
	else
	{
		Body2Rot.MakeIdentity();
		m_BodyAnchorPoints[1] = AnchorPoint;
	}
}

void BallSocket::UpdateCopy(float dt, SolverBase &solver)
{
	// Receive Number for this joint && Update Solver information about number of joints
	m_StartIdx = solver.m_NumJoints;
	solver.m_NumJoints += 3;

	float idt = 1.0f / dt;

	// Transfer CFM into solver
	float *ptr_m_CFM = &solver.m_CFM[m_StartIdx];

	for (int i = 0; i < 3; ++i)
	{
		*(ptr_m_CFM++) = m_CFM[i];
	}

	CMatrix3 Body1Rot, Body2Rot;
	CQuaternion tmpQuaterion;

	// Body position & orientation indices fuckup
	float	*ptr_m_NodePosRot_w = &solver.m_NodePosRot_w[0],
			*ptr_m_NodePosRot_x = &solver.m_NodePosRot_x[0],
			*ptr_m_NodePosRot_y = &solver.m_NodePosRot_y[0],
			*ptr_m_NodePosRot_z = &solver.m_NodePosRot_z[0];

	// Retrieving Bodies' Orientations
	tmpQuaterion.w = ptr_m_NodePosRot_w[m_Body1RIdx];
	tmpQuaterion.x = ptr_m_NodePosRot_x[m_Body1RIdx];
	tmpQuaterion.y = ptr_m_NodePosRot_y[m_Body1RIdx];
	tmpQuaterion.z = ptr_m_NodePosRot_z[m_Body1RIdx];
	tmpQuaterion.ToMatrix3(Body1Rot);

	if (m_Body2LIdx != -1)
	{
		tmpQuaterion.w = ptr_m_NodePosRot_w[m_Body2RIdx];
		tmpQuaterion.x = ptr_m_NodePosRot_x[m_Body2RIdx];
		tmpQuaterion.y = ptr_m_NodePosRot_y[m_Body2RIdx];
		tmpQuaterion.z = ptr_m_NodePosRot_z[m_Body2RIdx];
		tmpQuaterion.ToMatrix3(Body2Rot);
	}
	else
	{
		Body2Rot.MakeIdentity();
	}

	// Fill Jacobi matrices

	// Jacobian indices fuck-up
	float	*ptr_m_J_00 = &solver.m_J_00[0], *ptr_m_J_01 = &solver.m_J_01[0], *ptr_m_J_02 = &solver.m_J_02[0],
			*ptr_m_J_03 = &solver.m_J_03[0], *ptr_m_J_04 = &solver.m_J_04[0], *ptr_m_J_05 = &solver.m_J_05[0],
			*ptr_m_J_06 = &solver.m_J_06[0], *ptr_m_J_07 = &solver.m_J_07[0], *ptr_m_J_08 = &solver.m_J_08[0],
			*ptr_m_J_09 = &solver.m_J_09[0], *ptr_m_J_10 = &solver.m_J_10[0], *ptr_m_J_11 = &solver.m_J_11[0];

#if (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_J_12 = &solver.m_J_12[0], *ptr_m_J_13 = &solver.m_J_13[0], *ptr_m_J_14 = &solver.m_J_14[0],
			*ptr_m_J_15 = &solver.m_J_15[0], *ptr_m_J_16 = &solver.m_J_16[0], *ptr_m_J_17 = &solver.m_J_17[0];

#endif

	unsigned int Row0Idx, Row1Idx, Row2Idx;

	Row0Idx = m_StartIdx + 0;
	Row1Idx = m_StartIdx + 1;
	Row2Idx = m_StartIdx + 2;

	Vec3 wsAnchorPoints[2];

	// Linear Body1
	ptr_m_J_00[Row0Idx] =  1.0f;
	ptr_m_J_01[Row0Idx] =  0.0f;
	ptr_m_J_02[Row0Idx] =  0.0f;

	ptr_m_J_00[Row1Idx] =  0.0f;
	ptr_m_J_01[Row1Idx] =  1.0f;
	ptr_m_J_02[Row1Idx] =  0.0f;

	ptr_m_J_00[Row2Idx] =  0.0f;
	ptr_m_J_01[Row2Idx] =  0.0f;
	ptr_m_J_02[Row2Idx] =  1.0f;

	// Rotational Body1
	wsAnchorPoints[0] = Body1Rot * m_BodyAnchorPoints[0];

	// -CrossProdMatrix(b1rot * b1anchor_ls)
	ptr_m_J_03[Row0Idx] = 0.0f;
	ptr_m_J_04[Row0Idx] =  wsAnchorPoints[0].z;
	ptr_m_J_05[Row0Idx] = -wsAnchorPoints[0].y;
	ptr_m_J_03[Row1Idx] = -wsAnchorPoints[0].z;
	ptr_m_J_04[Row1Idx] = 0.0f;
	ptr_m_J_05[Row1Idx] =  wsAnchorPoints[0].x;
	ptr_m_J_03[Row2Idx] =  wsAnchorPoints[0].y;
	ptr_m_J_04[Row2Idx] = -wsAnchorPoints[0].x;
	ptr_m_J_05[Row2Idx] = 0.0f;

	wsAnchorPoints[0].x += ptr_m_NodePosRot_x[m_Body1LIdx];
	wsAnchorPoints[0].y += ptr_m_NodePosRot_y[m_Body1LIdx];
	wsAnchorPoints[0].z += ptr_m_NodePosRot_z[m_Body1LIdx];

	if (m_Body2LIdx != -1)
	{
		// Linear Body2
		ptr_m_J_06[Row0Idx] = -1.0f;
		ptr_m_J_07[Row0Idx] =  0.0f;
		ptr_m_J_08[Row0Idx] =  0.0f;

		ptr_m_J_06[Row1Idx] =  0.0f;
		ptr_m_J_07[Row1Idx] = -1.0f;
		ptr_m_J_08[Row1Idx] =  0.0f;

		ptr_m_J_06[Row2Idx] =  0.0f;
		ptr_m_J_07[Row2Idx] =  0.0f;
		ptr_m_J_08[Row2Idx] = -1.0f;

		// Rotational Body2
		wsAnchorPoints[1] = Body2Rot * m_BodyAnchorPoints[1];

		// +CrossProdMatrix(b2rot * b2anchor_ls)
		ptr_m_J_09[Row0Idx] = 0.0f;
		ptr_m_J_10[Row0Idx] = -wsAnchorPoints[1].z;
		ptr_m_J_11[Row0Idx] =  wsAnchorPoints[1].y;
		ptr_m_J_09[Row1Idx] =  wsAnchorPoints[1].z;
		ptr_m_J_10[Row1Idx] = 0.0f;
		ptr_m_J_11[Row1Idx] = -wsAnchorPoints[1].x;
		ptr_m_J_09[Row2Idx] = -wsAnchorPoints[1].y;
		ptr_m_J_10[Row2Idx] =  wsAnchorPoints[1].x;
		ptr_m_J_11[Row2Idx] = 0.0f;

		wsAnchorPoints[1].x += ptr_m_NodePosRot_x[m_Body2LIdx];
		wsAnchorPoints[1].y += ptr_m_NodePosRot_y[m_Body2LIdx];
		wsAnchorPoints[1].z += ptr_m_NodePosRot_z[m_Body2LIdx];
	}
	else
	{
		// Null out Jacobians (if attached to world)
		ptr_m_J_06[Row0Idx] = 0.0f;
		ptr_m_J_06[Row1Idx] = 0.0f;
		ptr_m_J_06[Row2Idx] = 0.0f;
		ptr_m_J_07[Row0Idx] = 0.0f;
		ptr_m_J_07[Row1Idx] = 0.0f;
		ptr_m_J_07[Row2Idx] = 0.0f;
		ptr_m_J_08[Row0Idx] = 0.0f;
		ptr_m_J_08[Row1Idx] = 0.0f;
		ptr_m_J_08[Row2Idx] = 0.0f;
		ptr_m_J_09[Row0Idx] = 0.0f;
		ptr_m_J_09[Row1Idx] = 0.0f;
		ptr_m_J_09[Row2Idx] = 0.0f;
		ptr_m_J_10[Row0Idx] = 0.0f;
		ptr_m_J_10[Row1Idx] = 0.0f;
		ptr_m_J_10[Row2Idx] = 0.0f;
		ptr_m_J_11[Row0Idx] = 0.0f;
		ptr_m_J_11[Row1Idx] = 0.0f;
		ptr_m_J_11[Row2Idx] = 0.0f;

		// Anchor Point is in World Space already
		wsAnchorPoints[1] = m_BodyAnchorPoints[1];
	}

	// Null remaining Jacobian
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

		ptr_m_J_12[Row0Idx] = 0.0f;
		ptr_m_J_12[Row1Idx] = 0.0f;
		ptr_m_J_12[Row2Idx] = 0.0f;
		ptr_m_J_13[Row0Idx] = 0.0f;
		ptr_m_J_13[Row1Idx] = 0.0f;
		ptr_m_J_13[Row2Idx] = 0.0f;
		ptr_m_J_14[Row0Idx] = 0.0f;
		ptr_m_J_14[Row1Idx] = 0.0f;
		ptr_m_J_14[Row2Idx] = 0.0f;
		ptr_m_J_15[Row0Idx] = 0.0f;
		ptr_m_J_15[Row1Idx] = 0.0f;
		ptr_m_J_15[Row2Idx] = 0.0f;
		ptr_m_J_16[Row0Idx] = 0.0f;
		ptr_m_J_16[Row1Idx] = 0.0f;
		ptr_m_J_16[Row2Idx] = 0.0f;
		ptr_m_J_17[Row0Idx] = 0.0f;
		ptr_m_J_17[Row1Idx] = 0.0f;
		ptr_m_J_17[Row2Idx] = 0.0f;

#endif

	float *ptr_m_J_rhs = &solver.m_J_rhs[m_StartIdx];

	// Baumgarte Stabilization
	*(ptr_m_J_rhs++) = idt * m_ERP * (wsAnchorPoints[1].x - wsAnchorPoints[0].x);
	*(ptr_m_J_rhs++) = idt * m_ERP * (wsAnchorPoints[1].y - wsAnchorPoints[0].y);
	*(ptr_m_J_rhs++) = idt * m_ERP * (wsAnchorPoints[1].z - wsAnchorPoints[0].z);

	float	*ptr_m_Lo = &solver.m_Lo[m_StartIdx],
			*ptr_m_Hi = &solver.m_Hi[m_StartIdx];

	*(ptr_m_Lo++) = -FLT_MAX;
	*(ptr_m_Hi++) =  FLT_MAX;
	*(ptr_m_Lo++) = -FLT_MAX;
	*(ptr_m_Hi++) =  FLT_MAX;
	*(ptr_m_Lo++) = -FLT_MAX;
	*(ptr_m_Hi++) =  FLT_MAX;

	// Calculate Jointed Nodes Indices
	unsigned int	*ptr_m_JtdNodes_00 = &solver.m_JtdNodes_00[Row0Idx],
					*ptr_m_JtdNodes_01 = &solver.m_JtdNodes_01[Row0Idx],
					*ptr_m_JtdNodes_02 = &solver.m_JtdNodes_02[Row0Idx],
					*ptr_m_JtdNodes_03 = &solver.m_JtdNodes_03[Row0Idx];

	// First Triple Index
	*(ptr_m_JtdNodes_00++) = m_Body1LIdx;
	*(ptr_m_JtdNodes_00++) = m_Body1LIdx;
	*(ptr_m_JtdNodes_00++) = m_Body1LIdx;

	// Second Triple Index
	*(ptr_m_JtdNodes_01++) = m_Body1RIdx;
	*(ptr_m_JtdNodes_01++) = m_Body1RIdx;
	*(ptr_m_JtdNodes_01++) = m_Body1RIdx;

	if (m_Body2LIdx != -1)
	{
		// Third Triple Index
		*(ptr_m_JtdNodes_02++) = m_Body2LIdx;
		*(ptr_m_JtdNodes_02++) = m_Body2LIdx;
		*(ptr_m_JtdNodes_02++) = m_Body2LIdx;

		// Fourth Triple Index
		*(ptr_m_JtdNodes_03++) = m_Body2RIdx;
		*(ptr_m_JtdNodes_03++) = m_Body2RIdx;
		*(ptr_m_JtdNodes_03++) = m_Body2RIdx;
	} else
	{
		// No triple -- "world" anchor

		// Third Triple Index
		*(ptr_m_JtdNodes_02++) = 0;
		*(ptr_m_JtdNodes_02++) = 0;
		*(ptr_m_JtdNodes_02++) = 0;

		// Fourth Triple Index
		*(ptr_m_JtdNodes_03++) = 0;
		*(ptr_m_JtdNodes_03++) = 0;
		*(ptr_m_JtdNodes_03++) = 0;
	}

	// Null out remaining indices
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

	unsigned int	*ptr_m_JtdNodes_04 = &solver.m_JtdNodes_04[Row0Idx],
					*ptr_m_JtdNodes_05 = &solver.m_JtdNodes_05[Row0Idx];

	// Fifth Triple Index (should be empty)
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;

	// Sixth Triple Index (should be empty)
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;

#endif


	float *ptr_m_lambda = &solver.m_lambda[m_StartIdx];

#if (JOINT_TRIPLES_CYCLE == 1)
	unsigned int *ptr_m_JointTriples = &solver.m_JointTriples[m_StartIdx];
#endif

	for (unsigned int i = 0; i < 3; ++i)
	{
		// Initial Guess
		*(ptr_m_lambda++) = m_lambda0[i];

#if (JOINT_TRIPLES_CYCLE == 1)
		if (m_Body2LIdx != -1)
		{
			*(ptr_m_JointTriples++) = 4;
		}
		else
		{
			// We're joined with World, so only 2 triples used
			*(ptr_m_JointTriples++) = 2;
		}
#endif
	}
}

void BallSocket::FetchLambdas(SolverBase &solver)
{
	m_lambda0[0] = solver.m_lambda[m_StartIdx    ];
	m_lambda0[1] = solver.m_lambda[m_StartIdx + 1];
	m_lambda0[2] = solver.m_lambda[m_StartIdx + 2];
}

void BallSocket::GetIndices(int &idx1_lin, int&idx1_ang, int &idx2_lin, int&idx2_ang)
{
	idx1_lin = m_Body1LIdx;
	idx1_ang = m_Body1RIdx;
	idx2_lin = m_Body2LIdx;
	idx2_ang = m_Body2RIdx;
}

void BallSocket::GetAnchorPoints_BS(Vec3 &bsAnchor1, Vec3 &bsAnchor2)
{
	bsAnchor1 = m_BodyAnchorPoints[0];
	bsAnchor2 = m_BodyAnchorPoints[1];
}


//////////////////////////////////////////////////////////////////////////
// PlaneConstraint
//////////////////////////////////////////////////////////////////////////

#if (PLANE_FRICTION == 1)

void PlaneConstraint::Init(SolverBase &solver, float ERP, float CFM, const Vec3 &PlanePoint, const Vec3 &PlaneNormal, int BodyIdx, float Tolerance)
{
	m_Tolerance = Tolerance;

	m_bActive = false;

	memset(m_lambda0, 0, 3 * sizeof(float));

	m_BodyIdx = BodyIdx;

	m_CFM[0] = m_CFM[1] = m_CFM[2] = CFM;
	m_ERP[0] = m_ERP[1] = m_ERP[2] = ERP;

	m_PlaneNormal = PlaneNormal.GetNormalized();
	m_PlaneD = PlanePoint.Dot(m_PlaneNormal);
}

void PlaneConstraint::UpdateCopy(float dt, SolverBase &solver)
{
	// Determine, if plane is in active
	Vec3 NodePos;

	// Body position & orientation indices fuck-up
	float	*ptr_m_NodePosRot_w = &solver.m_NodePosRot_w[0],
			*ptr_m_NodePosRot_x = &solver.m_NodePosRot_x[0],
			*ptr_m_NodePosRot_y = &solver.m_NodePosRot_y[0],
			*ptr_m_NodePosRot_z = &solver.m_NodePosRot_z[0];

	NodePos.x = ptr_m_NodePosRot_x[m_BodyIdx];
	NodePos.y = ptr_m_NodePosRot_y[m_BodyIdx];
	NodePos.z = ptr_m_NodePosRot_z[m_BodyIdx];

	float Violation = (NodePos).Dot(m_PlaneNormal) - m_PlaneD;
	if (Violation > PLANE_THRESHOLD)
	{
		m_bActive = false;

		m_lambda0[0] = 0.0f;
		m_lambda0[1] = 0.0f;
		m_lambda0[2] = 0.0f;

		return;
	}

	m_bActive = true;

	//////////////////////////////////////////////////////////////////////////
	// Standard Jacobi-fill and stuff
	//////////////////////////////////////////////////////////////////////////

	// Receive Number for this joint && Update Solver information about number of joints
	m_StartIdx = solver.m_NumJoints;
	solver.m_NumJoints += 3;

	float idt = 1.0f / dt;

	// Transfer CFM into solver
	solver.m_CFM[m_StartIdx  ] = m_CFM[0];
	solver.m_CFM[m_StartIdx+1] = m_CFM[1];
	solver.m_CFM[m_StartIdx+2] = m_CFM[2];

	// Pointer to Jacobi row
	float	*ptr_m_J_00 = &solver.m_J_00[m_StartIdx], *ptr_m_J_01 = &solver.m_J_01[m_StartIdx], *ptr_m_J_02 = &solver.m_J_02[m_StartIdx],
			*ptr_m_J_03 = &solver.m_J_03[m_StartIdx], *ptr_m_J_04 = &solver.m_J_04[m_StartIdx], *ptr_m_J_05 = &solver.m_J_05[m_StartIdx],
			*ptr_m_J_06 = &solver.m_J_06[m_StartIdx], *ptr_m_J_07 = &solver.m_J_07[m_StartIdx], *ptr_m_J_08 = &solver.m_J_08[m_StartIdx],
			*ptr_m_J_09 = &solver.m_J_09[m_StartIdx], *ptr_m_J_10 = &solver.m_J_10[m_StartIdx], *ptr_m_J_11 = &solver.m_J_11[m_StartIdx];

#if (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_J_12 = &solver.m_J_12[m_StartIdx], *ptr_m_J_13 = &solver.m_J_13[m_StartIdx], *ptr_m_J_14 = &solver.m_J_14[m_StartIdx],
			*ptr_m_J_15 = &solver.m_J_15[m_StartIdx], *ptr_m_J_16 = &solver.m_J_16[m_StartIdx], *ptr_m_J_17 = &solver.m_J_17[m_StartIdx];

#endif


	// Fill Jacobi matrix
	*ptr_m_J_00++ = m_PlaneNormal.x;
	*ptr_m_J_01++ = m_PlaneNormal.y;
	*ptr_m_J_02++ = m_PlaneNormal.z;

	Vec3 Tan1, Tan2;
	m_PlaneNormal.TangentSpace(Tan1, Tan2);

	*ptr_m_J_00++ = Tan1.x;
	*ptr_m_J_01++ = Tan1.y;
	*ptr_m_J_02++ = Tan1.z;

	*ptr_m_J_00++ = Tan2.x;
	*ptr_m_J_01++ = Tan2.y;
	*ptr_m_J_02++ = Tan2.z;

	for (int i = 0; i < 3; ++i)
	{
		*ptr_m_J_03++ = 0.0f;
		*ptr_m_J_04++ = 0.0f;
		*ptr_m_J_05++ = 0.0f;
		*ptr_m_J_06++ = 0.0f;
		*ptr_m_J_07++ = 0.0f;
		*ptr_m_J_08++ = 0.0f;
		*ptr_m_J_09++ = 0.0f;
		*ptr_m_J_10++ = 0.0f;
		*ptr_m_J_11++ = 0.0f;

	// Null remaining Jacobian
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

		*ptr_m_J_12++ = 0.0f;
		*ptr_m_J_13++ = 0.0f;
		*ptr_m_J_14++ = 0.0f;
		*ptr_m_J_15++ = 0.0f;
		*ptr_m_J_16++ = 0.0f;
		*ptr_m_J_17++ = 0.0f;

#endif

	}

	// Allow some penetration for stability
	Violation += m_Tolerance;

	if (Violation > 0.0f)
		Violation = 0.0f;

	// Baumgarte Stabilization
	solver.m_J_rhs[m_StartIdx  ] = (Violation < 0.0f) ? (-idt * m_ERP[0] * Violation) : 0.0f;
	solver.m_J_rhs[m_StartIdx+1] = 0.0f;
	solver.m_J_rhs[m_StartIdx+2] = 0.0f;

	// Limits
	solver.m_Lo[m_StartIdx  ] = 0.0f;
	solver.m_Hi[m_StartIdx  ] = FLT_MAX;

	float friction_limit = 10.1f;

	solver.m_Lo[m_StartIdx+1] = -friction_limit;
	solver.m_Hi[m_StartIdx+1] =  friction_limit;
	solver.m_Lo[m_StartIdx+2] = -friction_limit;
	solver.m_Hi[m_StartIdx+2] =  friction_limit;

// 	solver.m_Lo[m_StartIdx+1] = -FLT_MAX;
// 	solver.m_Hi[m_StartIdx+1] =  FLT_MAX;
// 	solver.m_Lo[m_StartIdx+2] = -FLT_MAX;
// 	solver.m_Hi[m_StartIdx+2] =  FLT_MAX;


	// Calculate Jointed Nodes Indices
	unsigned int	*ptr_m_JtdNodes_00 = &solver.m_JtdNodes_00[m_StartIdx],
					*ptr_m_JtdNodes_01 = &solver.m_JtdNodes_01[m_StartIdx],
					*ptr_m_JtdNodes_02 = &solver.m_JtdNodes_02[m_StartIdx],
					*ptr_m_JtdNodes_03 = &solver.m_JtdNodes_03[m_StartIdx];

#if (NUM_JOINT_TRIPLES == 6)

	unsigned int	*ptr_m_JtdNodes_04 = &solver.m_JtdNodes_04[m_StartIdx],
					*ptr_m_JtdNodes_05 = &solver.m_JtdNodes_05[m_StartIdx];

#endif

	for (int i = 0; i < 3; ++i)
	{
		*ptr_m_JtdNodes_00++ = m_BodyIdx;

		*ptr_m_JtdNodes_01++ = 0;
		*ptr_m_JtdNodes_02++ = 0;
		*ptr_m_JtdNodes_03++ = 0;


	// Null out remaining indices
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

		*ptr_m_JtdNodes_04++ = 0;
		*ptr_m_JtdNodes_05++ = 0;


#endif
	}


	// Lambda Initial Guess
	solver.m_lambda[m_StartIdx  ] = m_lambda0[0];
	solver.m_lambda[m_StartIdx+1] = m_lambda0[1];
	solver.m_lambda[m_StartIdx+2] = m_lambda0[2];

#if (JOINT_TRIPLES_CYCLE == 1)
	// If we're using cycles - number of triples used in joint
	solver.m_JointTriples[m_StartIdx  ] = 1;
	solver.m_JointTriples[m_StartIdx+1] = 1;
	solver.m_JointTriples[m_StartIdx+2] = 1;
#endif
}

void PlaneConstraint::FetchLambdas(SolverBase &solver)
{
	if (!m_bActive)
		return;

	m_lambda0[0] = solver.m_lambda[m_StartIdx  ];
	m_lambda0[1] = solver.m_lambda[m_StartIdx+1];
	m_lambda0[2] = solver.m_lambda[m_StartIdx+2];
}

//////////////////////////////////////////////////////////////////////////
#else
//////////////////////////////////////////////////////////////////////////

void PlaneConstraint::Init(SolverBase &solver, float ERP, float CFM, const Vec3 &PlanePoint, const Vec3 &PlaneNormal, int BodyIdx, float Tolerance)
{
	m_Tolerance = Tolerance;

	m_bActive = false;

	memset(m_lambda0, 0, 1 * sizeof(float));

	m_BodyIdx = BodyIdx;

	m_CFM[0] = CFM;
	m_ERP[0] = ERP;

	m_PlaneNormal = PlaneNormal.GetNormalized();
	m_PlaneD = PlanePoint.Dot(m_PlaneNormal);
}

void PlaneConstraint::UpdateCopy(float dt, SolverBase &solver)
{
	// Determine, if plane is in active
	Vec3 NodePos;

	// Body position & orientation indices fuck-up
	float	*ptr_m_NodePosRot_w = &solver.m_NodePosRot_w[0],
			*ptr_m_NodePosRot_x = &solver.m_NodePosRot_x[0],
			*ptr_m_NodePosRot_y = &solver.m_NodePosRot_y[0],
			*ptr_m_NodePosRot_z = &solver.m_NodePosRot_z[0];

	NodePos.x = ptr_m_NodePosRot_x[m_BodyIdx];
	NodePos.y = ptr_m_NodePosRot_y[m_BodyIdx];
	NodePos.z = ptr_m_NodePosRot_z[m_BodyIdx];

	float Violation = (NodePos).Dot(m_PlaneNormal) - m_PlaneD;
	if (Violation > PLANE_THRESHOLD)
	{
		m_bActive = false;

		m_lambda0[0] = 0.0f;

		return;
	}

	m_bActive = true;

	//////////////////////////////////////////////////////////////////////////
	// Standard Jacobi-fill and stuff
	//////////////////////////////////////////////////////////////////////////

	// Receive Number for this joint && Update Solver information about number of joints
	m_StartIdx = solver.m_NumJoints;
	solver.m_NumJoints += 1;

	float idt = 1.0f / dt;

	// Transfer CFM into solver
	solver.m_CFM[m_StartIdx  ] = m_CFM[0];
	solver.m_CFM[m_StartIdx+1] = m_CFM[1];
	solver.m_CFM[m_StartIdx+2] = m_CFM[2];

	// Pointer to Jacobi row
	float	*ptr_m_J_00 = &solver.m_J_00[m_StartIdx], *ptr_m_J_01 = &solver.m_J_01[m_StartIdx], *ptr_m_J_02 = &solver.m_J_02[m_StartIdx],
			*ptr_m_J_03 = &solver.m_J_03[m_StartIdx], *ptr_m_J_04 = &solver.m_J_04[m_StartIdx], *ptr_m_J_05 = &solver.m_J_05[m_StartIdx],
			*ptr_m_J_06 = &solver.m_J_06[m_StartIdx], *ptr_m_J_07 = &solver.m_J_07[m_StartIdx], *ptr_m_J_08 = &solver.m_J_08[m_StartIdx],
			*ptr_m_J_09 = &solver.m_J_09[m_StartIdx], *ptr_m_J_10 = &solver.m_J_10[m_StartIdx], *ptr_m_J_11 = &solver.m_J_11[m_StartIdx];

#if (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_J_12 = &solver.m_J_12[m_StartIdx], *ptr_m_J_13 = &solver.m_J_13[m_StartIdx], *ptr_m_J_14 = &solver.m_J_14[m_StartIdx],
			*ptr_m_J_15 = &solver.m_J_15[m_StartIdx], *ptr_m_J_16 = &solver.m_J_16[m_StartIdx], *ptr_m_J_17 = &solver.m_J_17[m_StartIdx];

#endif


	// Fill Jacobi matrix
	*ptr_m_J_00++ = m_PlaneNormal.x;
	*ptr_m_J_01++ = m_PlaneNormal.y;
	*ptr_m_J_02++ = m_PlaneNormal.z;

	*ptr_m_J_03++ = 0.0f;
	*ptr_m_J_04++ = 0.0f;
	*ptr_m_J_05++ = 0.0f;
	*ptr_m_J_06++ = 0.0f;
	*ptr_m_J_07++ = 0.0f;
	*ptr_m_J_08++ = 0.0f;
	*ptr_m_J_09++ = 0.0f;
	*ptr_m_J_10++ = 0.0f;
	*ptr_m_J_11++ = 0.0f;

	// Null remaining Jacobian
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

	*ptr_m_J_12++ = 0.0f;
	*ptr_m_J_13++ = 0.0f;
	*ptr_m_J_14++ = 0.0f;
	*ptr_m_J_15++ = 0.0f;
	*ptr_m_J_16++ = 0.0f;
	*ptr_m_J_17++ = 0.0f;

#endif

	// Allow some penetration for stability
	Violation += m_Tolerance;

	if (Violation > 0.0f)
		Violation = 0.0f;

	// Baumgarte Stabilization
	solver.m_J_rhs[m_StartIdx  ] = (Violation < 0.0f) ? (-idt * m_ERP[0] * Violation) : 0.0f;

	// Limits
	solver.m_Lo[m_StartIdx  ] = 0.0f;
	solver.m_Hi[m_StartIdx  ] = FLT_MAX;

	// Calculate Jointed Nodes Indices
	unsigned int	*ptr_m_JtdNodes_00 = &solver.m_JtdNodes_00[m_StartIdx],
					*ptr_m_JtdNodes_01 = &solver.m_JtdNodes_01[m_StartIdx],
					*ptr_m_JtdNodes_02 = &solver.m_JtdNodes_02[m_StartIdx],
					*ptr_m_JtdNodes_03 = &solver.m_JtdNodes_03[m_StartIdx];

#if (NUM_JOINT_TRIPLES == 6)

	unsigned int	*ptr_m_JtdNodes_04 = &solver.m_JtdNodes_04[m_StartIdx],
					*ptr_m_JtdNodes_05 = &solver.m_JtdNodes_05[m_StartIdx];

#endif

	*ptr_m_JtdNodes_00++ = m_BodyIdx;

	*ptr_m_JtdNodes_01++ = 0;
	*ptr_m_JtdNodes_02++ = 0;
	*ptr_m_JtdNodes_03++ = 0;


	// Null out remaining indices
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

	*ptr_m_JtdNodes_04++ = 0;
	*ptr_m_JtdNodes_05++ = 0;


#endif

	// Lambda Initial Guess
	solver.m_lambda[m_StartIdx  ] = m_lambda0[0];

#if (JOINT_TRIPLES_CYCLE == 1)
	// If we're using cycles - number of triples used in joint
	solver.m_JointTriples[m_StartIdx  ] = 1;
#endif
}

void PlaneConstraint::FetchLambdas(SolverBase &solver)
{
	if (!m_bActive)
		return;

	m_lambda0[0] = solver.m_lambda[m_StartIdx  ];
}

#endif