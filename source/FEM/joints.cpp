#include "joints.h"

//////////////////////////////////////////////////////////////////////////
// FEMJoint
//////////////////////////////////////////////////////////////////////////

void FEMJoint::SetYoungPoisson(float Young, float Poisson)
{
	m_Young = Young;
	m_Poisson = Poisson;

	float D_loc[6 * 6];
	memset(D_loc, 0, 6 * 6 * sizeof(float));

	// D
	float sigma1 = m_Poisson / (1 - m_Poisson);
	float premult = m_Young / ((1 + m_Poisson) * ( 1 - sigma1));
	float coef1 = premult * sigma1;
	float coef2 = 0.5f * premult * (1 - sigma1);

	D_loc[0 * 6 + 0] = premult * (1.0f - sigma1);
	D_loc[1 * 6 + 1] = premult * (1.0f + 2.0f * sigma1);
	D_loc[2 * 6 + 2] = premult * (1.0f - sigma1);

	D_loc[3 * 6 + 3] = coef2;
	D_loc[4 * 6 + 4] = coef2;
	D_loc[5 * 6 + 5] = coef2;

	float CFM_mult = 1.0f / m_InitialFEVolume;

	m_CFM[0] = CFM_mult / D_loc[0 * 6 + 0];
	m_CFM[1] = CFM_mult / D_loc[1 * 6 + 1];
	m_CFM[2] = CFM_mult / D_loc[2 * 6 + 2];
	m_CFM[3] = CFM_mult / D_loc[3 * 6 + 3];
	m_CFM[4] = CFM_mult / D_loc[4 * 6 + 4];
	m_CFM[5] = CFM_mult / D_loc[5 * 6 + 5];
}

void FEMJoint::Init(SolverBase &solver, const std::vector<Vec3, GlobalAllocator<Vec3>> &NodePos0,
					unsigned int Node1Idx, unsigned int Node2Idx, unsigned int Node3Idx, unsigned int Node4Idx,
					float Young, float Poisson, float Yield, float Creep,
					float MaxPlasticStrain, float Damping, float Regularization)
{
	m_Young = Young;
	m_Poisson = Poisson;

	m_Yield = Yield;
	m_Creep = Creep;
	m_MaxPlasticStrain = MaxPlasticStrain;
	m_Damping = Damping;

	memset(m_lambda0, 0, 6 * sizeof(float));

	float D_loc[6 * 6];
	memset(D_loc, 0, 6 * 6 * sizeof(float));
	memset(m_B_loc, 0, 6 * 12 * sizeof(float));
	memset(m_E_plastic, 0, 6 * sizeof(float));

	m_Node1Idx = Node1Idx;
	m_Node2Idx = Node2Idx;
	m_Node3Idx = Node3Idx;
	m_Node4Idx = Node4Idx;

	// Calculate B for Tetra [Old Pos]
	float det =	(NodePos0[m_Node2Idx].x - NodePos0[m_Node1Idx].x) * ((NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z) - (NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z)) +
				(NodePos0[m_Node3Idx].x - NodePos0[m_Node1Idx].x) * ((NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z) - (NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z)) +
				(NodePos0[m_Node4Idx].x - NodePos0[m_Node1Idx].x) * ((NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z) - (NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z));

	float invDet = 1.0f / det;

	float dNBdx = invDet * (  (NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z)
							- (NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z) );
	float dNBdy = invDet * (  (NodePos0[m_Node4Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z)
							- (NodePos0[m_Node3Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z) );
	float dNBdz = invDet * (  (NodePos0[m_Node3Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y)
							- (NodePos0[m_Node4Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y) );

	float dNCdx = invDet * (  (NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z)
							- (NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z) );
	float dNCdy = invDet * (  (NodePos0[m_Node2Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z)
							- (NodePos0[m_Node4Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z) );
	float dNCdz = invDet * (  (NodePos0[m_Node4Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y)
							- (NodePos0[m_Node2Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y) );

	float dNDdx = invDet * (  (NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z)
							- (NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y) * (NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z) );
	float dNDdy = invDet * (  (NodePos0[m_Node3Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z)
							- (NodePos0[m_Node2Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z) );
	float dNDdz = invDet * (  (NodePos0[m_Node2Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y)
							- (NodePos0[m_Node3Idx].x - NodePos0[m_Node1Idx].x) * (NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y) );

	float dNAdx = -(dNBdx + dNCdx + dNDdx);
	float dNAdy = -(dNBdy + dNCdy + dNDdy);
	float dNAdz = -(dNBdz + dNCdz + dNDdz);

	float sqrt_2 = 0.5f * sqrtf(2.0f);
	float sqrt_3 = sqrtf(3.0f) / 3.0f;
	float sqrt_6 = sqrtf(6.0f) / 6.0f;

	// B1
	m_B_loc[0*12+ 0] =  sqrt_2 * dNAdx;
	m_B_loc[0*12+ 2] = -sqrt_2 * dNAdz;

	m_B_loc[1*12+ 0] = sqrt_3 * dNAdx;
	m_B_loc[1*12+ 1] = sqrt_3 * dNAdy;
	m_B_loc[1*12+ 2] = sqrt_3 * dNAdz;

	m_B_loc[2*12+ 0] = sqrt_6 * dNAdx;
	m_B_loc[2*12+ 1] = -2.0f * sqrt_6 * dNAdy;
	m_B_loc[2*12+ 2] = sqrt_6 * dNAdz;

	m_B_loc[3*12+ 1] = dNAdz;
	m_B_loc[3*12+ 2] = dNAdy;
	m_B_loc[4*12+ 0] = dNAdz;
	m_B_loc[4*12+ 2] = dNAdx;
	m_B_loc[5*12+ 0] = dNAdy;
	m_B_loc[5*12+ 1] = dNAdx;

	// B2
	m_B_loc[0*12+ 3] =  sqrt_2 * dNBdx;
	m_B_loc[0*12+ 5] = -sqrt_2 * dNBdz;

	m_B_loc[1*12+ 3] = sqrt_3 * dNBdx;
	m_B_loc[1*12+ 4] = sqrt_3 * dNBdy;
	m_B_loc[1*12+ 5] = sqrt_3 * dNBdz;

	m_B_loc[2*12+ 3] = sqrt_6 * dNBdx;
	m_B_loc[2*12+ 4] = -2.0f * sqrt_6 * dNBdy;
	m_B_loc[2*12+ 5] = sqrt_6 * dNBdz;

	m_B_loc[3*12+ 4] = dNBdz;
	m_B_loc[3*12+ 5] = dNBdy;
	m_B_loc[4*12+ 3] = dNBdz;
	m_B_loc[4*12+ 5] = dNBdx;
	m_B_loc[5*12+ 3] = dNBdy;
	m_B_loc[5*12+ 4] = dNBdx;

	// B3
	m_B_loc[0*12+ 6] =  sqrt_2 * dNCdx;
	m_B_loc[0*12+ 8] = -sqrt_2 * dNCdz;

	m_B_loc[1*12+ 6] = sqrt_3 * dNCdx;
	m_B_loc[1*12+ 7] = sqrt_3 * dNCdy;
	m_B_loc[1*12+ 8] = sqrt_3 * dNCdz;

	m_B_loc[2*12+ 6] = sqrt_6 * dNCdx;
	m_B_loc[2*12+ 7] = -2.0f * sqrt_6 * dNCdy;
	m_B_loc[2*12+ 8] = sqrt_6 * dNCdz;

	m_B_loc[3*12+ 7] = dNCdz;
	m_B_loc[3*12+ 8] = dNCdy;
	m_B_loc[4*12+ 6] = dNCdz;
	m_B_loc[4*12+ 8] = dNCdx;
	m_B_loc[5*12+ 6] = dNCdy;
	m_B_loc[5*12+ 7] = dNCdx;

	// B4
	m_B_loc[0*12+ 9] =  sqrt_2 * dNDdx;
	m_B_loc[0*12+11] = -sqrt_2 * dNDdz;

	m_B_loc[1*12+ 9] = sqrt_3 * dNDdx;
	m_B_loc[1*12+10] = sqrt_3 * dNDdy;
	m_B_loc[1*12+11] = sqrt_3 * dNDdz;

	m_B_loc[2*12+ 9] = sqrt_6 * dNDdx;
	m_B_loc[2*12+10] = -2.0f * sqrt_6 * dNDdy;
	m_B_loc[2*12+11] = sqrt_6 * dNDdz;

	m_B_loc[3*12+10] = dNDdz;
	m_B_loc[3*12+11] = dNDdy;
	m_B_loc[4*12+ 9] = dNDdz;
	m_B_loc[4*12+11] = dNDdx;
	m_B_loc[5*12+ 9] = dNDdy;
	m_B_loc[5*12+10] = dNDdx;

	// D
	float sigma1 = m_Poisson / (1 - m_Poisson);
	float premult = m_Young / ((1 + m_Poisson) * (1 - sigma1));
	float coef1 = premult * sigma1;
	float coef2 = 0.5f * premult * (1 - sigma1);

	D_loc[0 * 6 + 0] = premult * (1.0f - sigma1);
	D_loc[1 * 6 + 1] = premult * (1.0f + 2.0f * sigma1);
	D_loc[2 * 6 + 2] = premult * (1.0f - sigma1);

	D_loc[3 * 6 + 3] = coef2;
	D_loc[4 * 6 + 4] = coef2;
	D_loc[5 * 6 + 5] = coef2;

	m_InitialFEVolume = (1.0f / 6.0f) * (NodePos0[m_Node4Idx] - NodePos0[m_Node1Idx]).Dot((NodePos0[m_Node3Idx] - NodePos0[m_Node1Idx]).Cross(NodePos0[m_Node2Idx] - NodePos0[m_Node1Idx]));

	float CFM_mult = 1.0f / m_InitialFEVolume;

	m_CFM[0] = CFM_mult / D_loc[0 * 6 + 0];
	m_CFM[1] = CFM_mult / D_loc[1 * 6 + 1];
	m_CFM[2] = CFM_mult / D_loc[2 * 6 + 2];
	m_CFM[3] = CFM_mult / D_loc[3 * 6 + 3];
	m_CFM[4] = CFM_mult / D_loc[4 * 6 + 4];
	m_CFM[5] = CFM_mult / D_loc[5 * 6 + 5];

	m_Regularization = Regularization;

#if (WARPING_POLAR_DECOMPOSITION == 1)

	m_mp0_inv.mMatrix[0][0] = NodePos0[m_Node2Idx].x - NodePos0[m_Node1Idx].x; 
	m_mp0_inv.mMatrix[1][0] = NodePos0[m_Node2Idx].y - NodePos0[m_Node1Idx].y; 
	m_mp0_inv.mMatrix[2][0] = NodePos0[m_Node2Idx].z - NodePos0[m_Node1Idx].z; 

	m_mp0_inv.mMatrix[0][1] = NodePos0[m_Node3Idx].x - NodePos0[m_Node1Idx].x; 
	m_mp0_inv.mMatrix[1][1] = NodePos0[m_Node3Idx].y - NodePos0[m_Node1Idx].y; 
	m_mp0_inv.mMatrix[2][1] = NodePos0[m_Node3Idx].z - NodePos0[m_Node1Idx].z; 

	m_mp0_inv.mMatrix[0][2] = NodePos0[m_Node4Idx].x - NodePos0[m_Node1Idx].x; 
	m_mp0_inv.mMatrix[1][2] = NodePos0[m_Node4Idx].y - NodePos0[m_Node1Idx].y; 
	m_mp0_inv.mMatrix[2][2] = NodePos0[m_Node4Idx].z - NodePos0[m_Node1Idx].z; 

	m_mp0_inv = m_mp0_inv.GetInverse();

#else

	// Calculating matrix N
	Vec3 edge0 = NodePos0[m_Node2Idx] - NodePos0[m_Node1Idx];

	Vec3 N1 = edge0 + (NodePos0[m_Node3Idx] - NodePos0[m_Node1Idx]) + (NodePos0[m_Node4Idx] - NodePos0[m_Node1Idx]);
	N1.Normalize();

	Vec3 N2 = N1.Cross(edge0);
	N2.Normalize();

	Vec3 N3 = N2.Cross(N1);
	N3.Normalize();

	m_N.mMatrix[0][0] = N1.x;
	m_N.mMatrix[1][0] = N1.y;
	m_N.mMatrix[2][0] = N1.z;

	m_N.mMatrix[0][1] = N2.x;
	m_N.mMatrix[1][1] = N2.y;
	m_N.mMatrix[2][1] = N2.z;

	m_N.mMatrix[0][2] = N3.x;
	m_N.mMatrix[1][2] = N3.y;
	m_N.mMatrix[2][2] = N3.z;

#endif

	for (unsigned int i = 0; i < 6; ++i)
	{
		m_Jp0[i] = 
			NodePos0[m_Node1Idx].x * m_B_loc[i * 12 + 0] +
			NodePos0[m_Node1Idx].y * m_B_loc[i * 12 + 1] +
			NodePos0[m_Node1Idx].z * m_B_loc[i * 12 + 2] +

			NodePos0[m_Node2Idx].x * m_B_loc[i * 12 + 3] +
			NodePos0[m_Node2Idx].y * m_B_loc[i * 12 + 4] +
			NodePos0[m_Node2Idx].z * m_B_loc[i * 12 + 5] +

			NodePos0[m_Node3Idx].x * m_B_loc[i * 12 + 6] +
			NodePos0[m_Node3Idx].y * m_B_loc[i * 12 + 7] +
			NodePos0[m_Node3Idx].z * m_B_loc[i * 12 + 8] +

			NodePos0[m_Node4Idx].x * m_B_loc[i * 12 + 9] +
			NodePos0[m_Node4Idx].y * m_B_loc[i * 12 + 10] +
			NodePos0[m_Node4Idx].z * m_B_loc[i * 12 + 11];

	}
}

void FEMJoint::UpdateCopy(float dt, SolverBase &solver)
{
	// Receive Number for this joint && Update Solver information about number of joints
	m_StartIdx = solver.m_NumJoints;
	solver.m_NumJoints += 6;

	m_fetchIndices[0] = m_StartIdx;
	m_fetchIndices[1] = m_StartIdx + 1;
	m_fetchIndices[2] = m_StartIdx + 2;
	m_fetchIndices[3] = m_StartIdx + 3;
	m_fetchIndices[4] = m_StartIdx + 4;
	m_fetchIndices[5] = m_StartIdx + 5;

	float idt = 1.0f / dt;

	float gamma = 1.0f / (1.0f + idt * m_Damping);

	// Transfer CFM into solver
	float *ptr_m_CFM = &solver.m_CFM[m_StartIdx];

	for (int i = 0; i < 6; ++i)
		*(ptr_m_CFM++) = (idt * gamma * m_CFM[i] + m_Regularization);


	// FE-nodes indices fuck up
	float	*ptr_m_NodePosRot_w = &solver.m_NodePosRot_w[0],
			*ptr_m_NodePosRot_x = &solver.m_NodePosRot_x[0],
			*ptr_m_NodePosRot_y = &solver.m_NodePosRot_y[0],
			*ptr_m_NodePosRot_z = &solver.m_NodePosRot_z[0];

	Vec3 v1, v2, v3;

	v1.x = ptr_m_NodePosRot_x[m_Node4Idx] - ptr_m_NodePosRot_x[m_Node1Idx];
	v1.y = ptr_m_NodePosRot_y[m_Node4Idx] - ptr_m_NodePosRot_y[m_Node1Idx];
	v1.z = ptr_m_NodePosRot_z[m_Node4Idx] - ptr_m_NodePosRot_z[m_Node1Idx];

	v2.x = ptr_m_NodePosRot_x[m_Node3Idx] - ptr_m_NodePosRot_x[m_Node1Idx];
	v2.y = ptr_m_NodePosRot_y[m_Node3Idx] - ptr_m_NodePosRot_y[m_Node1Idx];
	v2.z = ptr_m_NodePosRot_z[m_Node3Idx] - ptr_m_NodePosRot_z[m_Node1Idx];

	v3.x = ptr_m_NodePosRot_x[m_Node2Idx] - ptr_m_NodePosRot_x[m_Node1Idx];
	v3.y = ptr_m_NodePosRot_y[m_Node2Idx] - ptr_m_NodePosRot_y[m_Node1Idx];
	v3.z = ptr_m_NodePosRot_z[m_Node2Idx] - ptr_m_NodePosRot_z[m_Node1Idx];

#if (WARPING_POLAR_DECOMPOSITION == 1)

	CMatrix3 mpt, a,  p, u;

	mpt.mMatrix[0][0] = (ptr_m_NodePosRot_x[m_Node2Idx] - ptr_m_NodePosRot_x[m_Node1Idx]);
	mpt.mMatrix[0][1] = (ptr_m_NodePosRot_x[m_Node3Idx] - ptr_m_NodePosRot_x[m_Node1Idx]);
	mpt.mMatrix[0][2] = (ptr_m_NodePosRot_x[m_Node4Idx] - ptr_m_NodePosRot_x[m_Node1Idx]);

	mpt.mMatrix[1][0] = (ptr_m_NodePosRot_y[m_Node2Idx] - ptr_m_NodePosRot_y[m_Node1Idx]);
	mpt.mMatrix[1][1] = (ptr_m_NodePosRot_y[m_Node3Idx] - ptr_m_NodePosRot_y[m_Node1Idx]);
	mpt.mMatrix[1][2] = (ptr_m_NodePosRot_y[m_Node4Idx] - ptr_m_NodePosRot_y[m_Node1Idx]);

	mpt.mMatrix[2][0] = (ptr_m_NodePosRot_z[m_Node2Idx] - ptr_m_NodePosRot_z[m_Node1Idx]);
	mpt.mMatrix[2][1] = (ptr_m_NodePosRot_z[m_Node3Idx] - ptr_m_NodePosRot_z[m_Node1Idx]);
	mpt.mMatrix[2][2] = (ptr_m_NodePosRot_z[m_Node4Idx] - ptr_m_NodePosRot_z[m_Node1Idx]);

	a = mpt * m_mp0_inv;
	a.PolarDecompose(p, u);

#else

	// Calculating new matrix N

	CMatrix3 N_new;

	Vec3 edge0(	ptr_m_NodePosRot_x[m_Node2Idx] - ptr_m_NodePosRot_x[m_Node1Idx],
				ptr_m_NodePosRot_y[m_Node2Idx] - ptr_m_NodePosRot_y[m_Node1Idx],
				ptr_m_NodePosRot_z[m_Node2Idx] - ptr_m_NodePosRot_z[m_Node1Idx] );

	Vec3 N1;
	N1.x = edge0.x	+ (ptr_m_NodePosRot_x[m_Node3Idx] - ptr_m_NodePosRot_x[m_Node1Idx])
					+ (ptr_m_NodePosRot_x[m_Node4Idx] - ptr_m_NodePosRot_x[m_Node1Idx]);
	N1.y = edge0.y	+ (ptr_m_NodePosRot_y[m_Node3Idx] - ptr_m_NodePosRot_y[m_Node1Idx])
					+ (ptr_m_NodePosRot_y[m_Node4Idx] - ptr_m_NodePosRot_y[m_Node1Idx]);
	N1.z = edge0.z	+ (ptr_m_NodePosRot_z[m_Node3Idx] - ptr_m_NodePosRot_z[m_Node1Idx])
					+ (ptr_m_NodePosRot_z[m_Node4Idx] - ptr_m_NodePosRot_z[m_Node1Idx]);
	N1.Normalize();

	Vec3 N2 = N1.Cross(edge0);
	N2.Normalize();

	Vec3 N3 = N2.Cross(N1);
	N3.Normalize();

	N_new.mMatrix[0][0] = N1.x;
	N_new.mMatrix[1][0] = N1.y;
	N_new.mMatrix[2][0] = N1.z;

	N_new.mMatrix[0][1] = N2.x;
	N_new.mMatrix[1][1] = N2.y;
	N_new.mMatrix[2][1] = N2.z;

	N_new.mMatrix[0][2] = N3.x;
	N_new.mMatrix[1][2] = N3.y;
	N_new.mMatrix[2][2] = N3.z;

	CMatrix3 u = N_new * m_N.GetTransposed();

#endif	

	if (u.Det() < 0.0f)
		u = u * -1.0f;

	// Jacobian indices fuck up
	float	*ptr_m_J_00 = &solver.m_J_00[0], *ptr_m_J_01 = &solver.m_J_01[0], *ptr_m_J_02 = &solver.m_J_02[0],
			*ptr_m_J_03 = &solver.m_J_03[0], *ptr_m_J_04 = &solver.m_J_04[0], *ptr_m_J_05 = &solver.m_J_05[0],
			*ptr_m_J_06 = &solver.m_J_06[0], *ptr_m_J_07 = &solver.m_J_07[0], *ptr_m_J_08 = &solver.m_J_08[0],
			*ptr_m_J_09 = &solver.m_J_09[0], *ptr_m_J_10 = &solver.m_J_10[0], *ptr_m_J_11 = &solver.m_J_11[0];

#if (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_J_12 = &solver.m_J_12[0], *ptr_m_J_13 = &solver.m_J_13[0], *ptr_m_J_14 = &solver.m_J_14[0],
			*ptr_m_J_15 = &solver.m_J_15[0], *ptr_m_J_16 = &solver.m_J_16[0], *ptr_m_J_17 = &solver.m_J_17[0];

#endif

	unsigned int Jrow0, Jrow1, Jrow2;
	unsigned int Brow0, Brow1, Brow2;

#define TRANSFER_B_SUBMATRIX(i, j, t0, t1, t2) \
			Jrow0 = (m_StartIdx + i * 3 + 0);\
			Jrow1 = (m_StartIdx + i * 3 + 1);\
			Jrow2 = (m_StartIdx + i * 3 + 2);\
			\
			Brow0 = (i * 3 + 0) * 12 + j * 3;\
			Brow1 = (i * 3 + 1) * 12 + j * 3;\
			Brow2 = (i * 3 + 2) * 12 + j * 3;\
			\
			ptr_m_J_##t0[Jrow0] =	m_B_loc[Brow0 + 0] * u.mMatrix[0][0] + m_B_loc[Brow0 + 1] * u.mMatrix[0][1] + \
									m_B_loc[Brow0 + 2] * u.mMatrix[0][2];\
			\
			ptr_m_J_##t1[Jrow0] =	m_B_loc[Brow0 + 0] * u.mMatrix[1][0] + m_B_loc[Brow0 + 1] * u.mMatrix[1][1] + \
									m_B_loc[Brow0 + 2] * u.mMatrix[1][2];\
			\
			ptr_m_J_##t2[Jrow0] =	m_B_loc[Brow0 + 0] * u.mMatrix[2][0] + m_B_loc[Brow0 + 1] * u.mMatrix[2][1] + \
									m_B_loc[Brow0 + 2] * u.mMatrix[2][2];\
			\
			ptr_m_J_##t0[Jrow1] =	m_B_loc[Brow1 + 0] * u.mMatrix[0][0] + m_B_loc[Brow1 + 1] * u.mMatrix[0][1] + \
									m_B_loc[Brow1 + 2] * u.mMatrix[0][2];\
			\
			ptr_m_J_##t1[Jrow1] =	m_B_loc[Brow1 + 0] * u.mMatrix[1][0] + m_B_loc[Brow1 + 1] * u.mMatrix[1][1] + \
									m_B_loc[Brow1 + 2] * u.mMatrix[1][2];\
			\
			ptr_m_J_##t2[Jrow1] =	m_B_loc[Brow1 + 0] * u.mMatrix[2][0] + m_B_loc[Brow1 + 1] * u.mMatrix[2][1] + \
									m_B_loc[Brow1 + 2] * u.mMatrix[2][2];\
			\
			ptr_m_J_##t0[Jrow2] =	m_B_loc[Brow2 + 0] * u.mMatrix[0][0] + m_B_loc[Brow2 + 1] * u.mMatrix[0][1] + \
									m_B_loc[Brow2 + 2] * u.mMatrix[0][2];\
			\
			ptr_m_J_##t1[Jrow2] =	m_B_loc[Brow2 + 0] * u.mMatrix[1][0] + m_B_loc[Brow2 + 1] * u.mMatrix[1][1] + \
									m_B_loc[Brow2 + 2] * u.mMatrix[1][2];\
			\
			ptr_m_J_##t2[Jrow2] =	m_B_loc[Brow2 + 0] * u.mMatrix[2][0] + m_B_loc[Brow2 + 1] * u.mMatrix[2][1] + \
									m_B_loc[Brow2 + 2] * u.mMatrix[2][2];

	TRANSFER_B_SUBMATRIX(0, 0, 00, 01, 02);
	TRANSFER_B_SUBMATRIX(0, 1, 03, 04, 05);
	TRANSFER_B_SUBMATRIX(0, 2, 06, 07, 08);
	TRANSFER_B_SUBMATRIX(0, 3, 09, 10, 11);
	TRANSFER_B_SUBMATRIX(1, 0, 00, 01, 02);
	TRANSFER_B_SUBMATRIX(1, 1, 03, 04, 05);
	TRANSFER_B_SUBMATRIX(1, 2, 06, 07, 08);
	TRANSFER_B_SUBMATRIX(1, 3, 09, 10, 11);

#undef TRANSFER_B_SUBMATRIX

	// Null remaining Jacobian
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

	for (int i = 0; i < 6; ++i)
	{
		ptr_m_J_12[m_StartIdx + i] = 0.0f;
		ptr_m_J_13[m_StartIdx + i] = 0.0f;
		ptr_m_J_14[m_StartIdx + i] = 0.0f;
		ptr_m_J_15[m_StartIdx + i] = 0.0f;
		ptr_m_J_16[m_StartIdx + i] = 0.0f;
		ptr_m_J_17[m_StartIdx + i] = 0.0f;
	}

#endif

	unsigned int i;

	float e_elastic[6];

	// Calculate Elastic Strain = Total Strain - Plastic Strain [ Bq - e_plastic ]
	for (i = 0; i < 6; ++i)
	{
		unsigned int RowIdx = m_StartIdx + i;

		e_elastic[i] = (
					ptr_m_NodePosRot_x[m_Node1Idx] * ptr_m_J_00[RowIdx] +
					ptr_m_NodePosRot_y[m_Node1Idx] * ptr_m_J_01[RowIdx] +
					ptr_m_NodePosRot_z[m_Node1Idx] * ptr_m_J_02[RowIdx] +

					ptr_m_NodePosRot_x[m_Node2Idx] * ptr_m_J_03[RowIdx] +
					ptr_m_NodePosRot_y[m_Node2Idx] * ptr_m_J_04[RowIdx] +
					ptr_m_NodePosRot_z[m_Node2Idx] * ptr_m_J_05[RowIdx] +

					ptr_m_NodePosRot_x[m_Node3Idx] * ptr_m_J_06[RowIdx] +
					ptr_m_NodePosRot_y[m_Node3Idx] * ptr_m_J_07[RowIdx] +
					ptr_m_NodePosRot_z[m_Node3Idx] * ptr_m_J_08[RowIdx] +

					ptr_m_NodePosRot_x[m_Node4Idx] * ptr_m_J_09[RowIdx] +
					ptr_m_NodePosRot_y[m_Node4Idx] * ptr_m_J_10[RowIdx] +
					ptr_m_NodePosRot_z[m_Node4Idx] * ptr_m_J_11[RowIdx]
			- m_Jp0[i] - m_E_plastic[i]);
	}

	float el_norm = sqrtf(
		e_elastic[0] * e_elastic[0] + e_elastic[1] * e_elastic[1] + e_elastic[2] * e_elastic[2] +
		e_elastic[3] * e_elastic[3] + e_elastic[4] * e_elastic[4] + e_elastic[5] * e_elastic[5]
	);

	m_StrainNorm = el_norm;

	// Creep is in [0; 1/dt]
	float Creep = m_Creep * idt;

	if (el_norm > m_Yield)
	{
		for (i = 0; i < 6; ++i)
			m_E_plastic[i] += dt * Creep * e_elastic[i];
	}

	float pl_norm = sqrtf(
		m_E_plastic[0] * m_E_plastic[0] + m_E_plastic[1] * m_E_plastic[1] + m_E_plastic[2] * m_E_plastic[2] +
		m_E_plastic[3] * m_E_plastic[3] + m_E_plastic[4] * m_E_plastic[4] + m_E_plastic[5] * m_E_plastic[5]
	);

	if (pl_norm > m_MaxPlasticStrain)
	{
		for (i = 0; i < 6; ++i)
		{
			m_E_plastic[i] *= m_MaxPlasticStrain / pl_norm;
		}
	}

	float *ptr_m_J_rhs = &solver.m_J_rhs[m_StartIdx];

	for (i = 0; i < 6; ++i)
	{
		*(ptr_m_J_rhs++) = -idt * gamma * e_elastic[i];
	}

	// Limits
	float liForces = FLT_MAX;

	float	*ptr_m_Lo = &solver.m_Lo[m_StartIdx],
			*ptr_m_Hi = &solver.m_Hi[m_StartIdx];

	for (i = 0; i < 6; ++i)
	{
		*(ptr_m_Lo++) = -liForces;
		*(ptr_m_Hi++) =  liForces;
	}

	// Calculate Jointed Nodes Indices
	unsigned int	*ptr_m_JtdNodes_00 = &solver.m_JtdNodes_00[m_StartIdx],
					*ptr_m_JtdNodes_01 = &solver.m_JtdNodes_01[m_StartIdx],
					*ptr_m_JtdNodes_02 = &solver.m_JtdNodes_02[m_StartIdx],
					*ptr_m_JtdNodes_03 = &solver.m_JtdNodes_03[m_StartIdx];

	// First Triple Index
	*(ptr_m_JtdNodes_00++) = m_Node1Idx;
	*(ptr_m_JtdNodes_00++) = m_Node1Idx;
	*(ptr_m_JtdNodes_00++) = m_Node1Idx;
	*(ptr_m_JtdNodes_00++) = m_Node1Idx;
	*(ptr_m_JtdNodes_00++) = m_Node1Idx;
	*(ptr_m_JtdNodes_00++) = m_Node1Idx;

	// Second Triple Index
	*(ptr_m_JtdNodes_01++) = m_Node2Idx;
	*(ptr_m_JtdNodes_01++) = m_Node2Idx;
	*(ptr_m_JtdNodes_01++) = m_Node2Idx;
	*(ptr_m_JtdNodes_01++) = m_Node2Idx;
	*(ptr_m_JtdNodes_01++) = m_Node2Idx;
	*(ptr_m_JtdNodes_01++) = m_Node2Idx;

	// Third Triple Index
	*(ptr_m_JtdNodes_02++) = m_Node3Idx;
	*(ptr_m_JtdNodes_02++) = m_Node3Idx;
	*(ptr_m_JtdNodes_02++) = m_Node3Idx;
	*(ptr_m_JtdNodes_02++) = m_Node3Idx;
	*(ptr_m_JtdNodes_02++) = m_Node3Idx;
	*(ptr_m_JtdNodes_02++) = m_Node3Idx;

	// Fourth Triple Index
	*(ptr_m_JtdNodes_03++) = m_Node4Idx;
	*(ptr_m_JtdNodes_03++) = m_Node4Idx;
	*(ptr_m_JtdNodes_03++) = m_Node4Idx;
	*(ptr_m_JtdNodes_03++) = m_Node4Idx;
	*(ptr_m_JtdNodes_03++) = m_Node4Idx;
	*(ptr_m_JtdNodes_03++) = m_Node4Idx;


	// Null out remaining indices
#if (NUM_JOINT_TRIPLES == 4)

#elif (NUM_JOINT_TRIPLES == 6)

	unsigned int	*ptr_m_JtdNodes_04 = &solver.m_JtdNodes_04[m_StartIdx],
					*ptr_m_JtdNodes_05 = &solver.m_JtdNodes_05[m_StartIdx];

	// Fifth Triple Index (should be empty)
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;
	*(ptr_m_JtdNodes_04++) = 0;

	// Sixth Triple Index (should be empty)
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;
	*(ptr_m_JtdNodes_05++) = 0;

#endif


	float *ptr_m_lambda = &solver.m_lambda[m_StartIdx];

#if (JOINT_TRIPLES_CYCLE == 1)
	unsigned int *ptr_m_JointTriples = &solver.m_JointTriples[m_StartIdx];
#endif

	for (unsigned int i = 0; i < 6; ++i)
	{
		// Initial Guess
		*(ptr_m_lambda++) = m_lambda0[i];

#if (JOINT_TRIPLES_CYCLE == 1)
		*(ptr_m_JointTriples++) = 4;
#endif
	}
}

void FEMJoint::FetchLambdas(SolverBase &solver)
{
	/*
#if (HARDWARE_SOLVER == 0)
	m_lambda0[0] = solver.m_lambda[m_StartIdx    ];
	m_lambda0[1] = solver.m_lambda[m_StartIdx + 1];
	m_lambda0[2] = solver.m_lambda[m_StartIdx + 2];
	m_lambda0[3] = solver.m_lambda[m_StartIdx + 3];
	m_lambda0[4] = solver.m_lambda[m_StartIdx + 4];
	m_lambda0[5] = solver.m_lambda[m_StartIdx + 5];
#else
	*/
	m_lambda0[0] = solver.m_lambda[m_fetchIndices[0]];
	m_lambda0[1] = solver.m_lambda[m_fetchIndices[1]];
	m_lambda0[2] = solver.m_lambda[m_fetchIndices[2]];
	m_lambda0[3] = solver.m_lambda[m_fetchIndices[3]];
	m_lambda0[4] = solver.m_lambda[m_fetchIndices[4]];
	m_lambda0[5] = solver.m_lambda[m_fetchIndices[5]];
}



//////////////////////////////////////////////////////////////////////////
// BallSocket_FEM
//////////////////////////////////////////////////////////////////////////


void BallSocket_FEM::Init(SolverBase &solver, float ERP, float CFM, const Vec3 &AnchorPoint, int Body1LIdx, int Body1RIdx, int Node1Idx, int Node2Idx, int Node3Idx)
{
	memset(m_lambda0, 0, 3 * sizeof(float));

	m_Body1LIdx = Body1LIdx;
	m_Body1RIdx = Body1RIdx;
	m_Node1 = Node1Idx;
	m_Node2 = Node2Idx;
	m_Node3 = Node3Idx;

	m_CFM[0] = CFM;
	m_CFM[1] = CFM;
	m_CFM[2] = CFM;

	m_ERP = ERP;

	// FE setup

	float	*ptr_m_NodePosRot_w = &solver.m_NodePosRot_w[0],
			*ptr_m_NodePosRot_x = &solver.m_NodePosRot_x[0],
			*ptr_m_NodePosRot_y = &solver.m_NodePosRot_y[0],
			*ptr_m_NodePosRot_z = &solver.m_NodePosRot_z[0];

	// Calculating Barycentric for
	Vec3 v_p1p0, v_p2p0, p0;
	p0.x = ptr_m_NodePosRot_x[m_Node1];
	p0.y = ptr_m_NodePosRot_y[m_Node1];
	p0.z = ptr_m_NodePosRot_z[m_Node1];

	v_p1p0.x = ptr_m_NodePosRot_x[m_Node2] - ptr_m_NodePosRot_x[m_Node1];
	v_p1p0.y = ptr_m_NodePosRot_y[m_Node2] - ptr_m_NodePosRot_y[m_Node1];
	v_p1p0.z = ptr_m_NodePosRot_z[m_Node2] - ptr_m_NodePosRot_z[m_Node1];

	v_p2p0.x = ptr_m_NodePosRot_x[m_Node3] - ptr_m_NodePosRot_x[m_Node1];
	v_p2p0.y = ptr_m_NodePosRot_y[m_Node3] - ptr_m_NodePosRot_y[m_Node1];
	v_p2p0.z = ptr_m_NodePosRot_z[m_Node3] - ptr_m_NodePosRot_z[m_Node1];

	Vec3 Normal = v_p1p0.Cross(v_p2p0);
	Normal.Normalize();

	m_Delta = (AnchorPoint - p0).Dot(Normal);

	// Calculate tmpVectors
	Vec3 vPA, vPB, vPC, vPPproj;

	vPA.x = ptr_m_NodePosRot_x[m_Node1] - AnchorPoint.x;
	vPA.y = ptr_m_NodePosRot_y[m_Node1] - AnchorPoint.y;
	vPA.z = ptr_m_NodePosRot_z[m_Node1] - AnchorPoint.z;
	vPB.x = ptr_m_NodePosRot_x[m_Node2] - AnchorPoint.x;
	vPB.y = ptr_m_NodePosRot_y[m_Node2] - AnchorPoint.y;
	vPB.z = ptr_m_NodePosRot_z[m_Node2] - AnchorPoint.z;
	vPC.x = ptr_m_NodePosRot_x[m_Node3] - AnchorPoint.x;
	vPC.y = ptr_m_NodePosRot_y[m_Node3] - AnchorPoint.y;
	vPC.z = ptr_m_NodePosRot_z[m_Node3] - AnchorPoint.z;

	vPPproj.x = -m_Delta * Normal.x;
	vPPproj.y = -m_Delta * Normal.y;
	vPPproj.z = -m_Delta * Normal.z;

	// Volumes of tetrahedrons for barycentric coordinates calculation
	float tABPProj, tBCPProj, tCAPProj, tABPC;

	tABPProj = vPA.Triple(vPB, vPPproj);
	tBCPProj = vPB.Triple(vPC, vPPproj);
	tCAPProj = vPC.Triple(vPA, vPPproj);
	tABPC = vPA.Triple(vPB, vPC);

	m_Alpha = tBCPProj / tABPC;
	m_Beta  = tCAPProj / tABPC;
	m_Gamma = tABPProj / tABPC;

	if (m_Body1LIdx != -1)
	{
		CMatrix3 Body1Rot;
		CQuaternion tmpQuaterion;

		// Retrieving Bodies' Orientations
		tmpQuaterion.w = ptr_m_NodePosRot_w[m_Body1RIdx];
		tmpQuaterion.x = ptr_m_NodePosRot_x[m_Body1RIdx];
		tmpQuaterion.y = ptr_m_NodePosRot_y[m_Body1RIdx];
		tmpQuaterion.z = ptr_m_NodePosRot_z[m_Body1RIdx];
		tmpQuaterion.ToMatrix3(Body1Rot);

		Vec3 BodyPos;

		// Calculating Anchor Points in Local Space
		BodyPos.x = ptr_m_NodePosRot_x[m_Body1LIdx];
		BodyPos.y = ptr_m_NodePosRot_y[m_Body1LIdx];
		BodyPos.z = ptr_m_NodePosRot_z[m_Body1LIdx];
		m_BodyAnchorPoints[0] = Body1Rot.GetTransposed() * (AnchorPoint - BodyPos);
	}
	else
	{
		m_BodyAnchorPoints[0] = AnchorPoint;
	}
}

void BallSocket_FEM::UpdateCopy(float dt, SolverBase &solver)
{
	// Receive Number for this joint && Update Solver information about number of joints
	m_StartIdx = solver.m_NumJoints;
	solver.m_NumJoints += 3;

	float idt = 1.0f / dt;

	// Transfer CFM into solver
	float *ptr_m_CFM = &solver.m_CFM[m_StartIdx];

	for (int i = 0; i < 3; ++i)
		*(ptr_m_CFM++) = m_CFM[i];

	// Jacobian indices fuck-up
#if (NUM_JOINT_TRIPLES == 4)

#error Ball-socket FEM-version should be compiled ONLY with NUM_JOINT_TRIPLES more than 4

#elif (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_J_00 = &solver.m_J_00[0], *ptr_m_J_01 = &solver.m_J_01[0], *ptr_m_J_02 = &solver.m_J_02[0],
			*ptr_m_J_03 = &solver.m_J_03[0], *ptr_m_J_04 = &solver.m_J_04[0], *ptr_m_J_05 = &solver.m_J_05[0],
			*ptr_m_J_06 = &solver.m_J_06[0], *ptr_m_J_07 = &solver.m_J_07[0], *ptr_m_J_08 = &solver.m_J_08[0],
			*ptr_m_J_09 = &solver.m_J_09[0], *ptr_m_J_10 = &solver.m_J_10[0], *ptr_m_J_11 = &solver.m_J_11[0],
			*ptr_m_J_12 = &solver.m_J_12[0], *ptr_m_J_13 = &solver.m_J_13[0], *ptr_m_J_14 = &solver.m_J_14[0],
			*ptr_m_J_15 = &solver.m_J_15[0], *ptr_m_J_16 = &solver.m_J_16[0], *ptr_m_J_17 = &solver.m_J_17[0];

#endif

	unsigned int Row0Idx, Row1Idx, Row2Idx;

	Row0Idx = m_StartIdx + 0;
	Row1Idx = m_StartIdx + 1;
	Row2Idx = m_StartIdx + 2;

	// Body position & orientation indices fuck-up
	float *ptr_m_NodePosRot;

	// FE-nodes indices fuck up
	float	*ptr_m_NodePosRot_w = &solver.m_NodePosRot_w[0],
			*ptr_m_NodePosRot_x = &solver.m_NodePosRot_x[0],
			*ptr_m_NodePosRot_y = &solver.m_NodePosRot_y[0],
			*ptr_m_NodePosRot_z = &solver.m_NodePosRot_z[0];

	// Filling Jacobians for FE-nodes

	// (p1 - p0), (p2 - p0)
	Vec3 v_p1p0, v_p2p0, v_p2p1;
	Vec3 vec1, vec2, vec3;

	v_p1p0.x = ptr_m_NodePosRot_x[m_Node2] - ptr_m_NodePosRot_x[m_Node1];
	v_p1p0.y = ptr_m_NodePosRot_y[m_Node2] - ptr_m_NodePosRot_y[m_Node1];
	v_p1p0.z = ptr_m_NodePosRot_z[m_Node2] - ptr_m_NodePosRot_z[m_Node1];

	v_p2p0.x = ptr_m_NodePosRot_x[m_Node3] - ptr_m_NodePosRot_x[m_Node1];
	v_p2p0.y = ptr_m_NodePosRot_y[m_Node3] - ptr_m_NodePosRot_y[m_Node1];
	v_p2p0.z = ptr_m_NodePosRot_z[m_Node3] - ptr_m_NodePosRot_z[m_Node1];

	Vec3 Normal = v_p1p0.Cross(v_p2p0);

	// due to internal formulation, MOD = |n|
	float mod = Normal.GetLength();

	vec3 = -(1.0f / (mod * mod * mod)) * Normal;

	// temporary storage
	float vecTmp1[3], vecTmp2[3];

	// J_p0 [ starts from 6 ]
	v_p2p1 = v_p2p0 - v_p1p0;
	vecTmp1[0] = (v_p2p1.x) / mod;
	vecTmp1[1] = (v_p2p1.y) / mod;
	vecTmp1[2] = (v_p2p1.z) / mod;

	vec1 = Normal.Cross(v_p2p1);
	vecTmp2[0] = vec1.x;
	vecTmp2[1] = vec1.y;
	vecTmp2[2] = vec1.z;

	// row1
	ptr_m_J_00[Row0Idx] = -m_Delta *  vec3.x * vecTmp2[0] - m_Alpha;
	ptr_m_J_01[Row0Idx] = -m_Delta * (vec3.x * vecTmp2[1] - vecTmp1[2]);
	ptr_m_J_02[Row0Idx] = -m_Delta * (vec3.x * vecTmp2[2] + vecTmp1[1]);

	// row2
	ptr_m_J_00[Row1Idx] = -m_Delta * (vec3.y * vecTmp2[0] + vecTmp1[2]);
	ptr_m_J_01[Row1Idx] = -m_Delta *  vec3.y * vecTmp2[1] - m_Alpha;
	ptr_m_J_02[Row1Idx] = -m_Delta * (vec3.y * vecTmp2[2] - vecTmp1[0]);

	// row3
	ptr_m_J_00[Row2Idx] = -m_Delta * (vec3.z * vecTmp2[0] - vecTmp1[1]);
	ptr_m_J_01[Row2Idx] = -m_Delta * (vec3.z * vecTmp2[1] + vecTmp1[0]);
	ptr_m_J_02[Row2Idx] = -m_Delta *  vec3.z * vecTmp2[2] - m_Alpha;


	// J_p1 [ starts from 9 ]
	vecTmp1[0] = -(v_p2p0.x) / mod;
	vecTmp1[1] = -(v_p2p0.y) / mod;
	vecTmp1[2] = -(v_p2p0.z) / mod;

	vec1 = Normal.Cross(-v_p2p0);
	vecTmp2[0] = vec1.x;
	vecTmp2[1] = vec1.y;
	vecTmp2[2] = vec1.z;

	// row1
	ptr_m_J_03[Row0Idx] = -m_Delta *  vec3.x * vecTmp2[0] - m_Beta;
	ptr_m_J_04[Row0Idx] = -m_Delta * (vec3.x * vecTmp2[1] - vecTmp1[2]);
	ptr_m_J_05[Row0Idx] = -m_Delta * (vec3.x * vecTmp2[2] + vecTmp1[1]);

	// row2
	ptr_m_J_03[Row1Idx] = -m_Delta * (vec3.y * vecTmp2[0] + vecTmp1[2]);
	ptr_m_J_04[Row1Idx] = -m_Delta *  vec3.y * vecTmp2[1] - m_Beta;
	ptr_m_J_05[Row1Idx] = -m_Delta * (vec3.y * vecTmp2[2] - vecTmp1[0]);

	// row3
	ptr_m_J_03[Row2Idx] = -m_Delta * (vec3.z * vecTmp2[0] - vecTmp1[1]);
	ptr_m_J_04[Row2Idx] = -m_Delta * (vec3.z * vecTmp2[1] + vecTmp1[0]);
	ptr_m_J_05[Row2Idx] = -m_Delta *  vec3.z * vecTmp2[2] - m_Beta;

	// J_p2 [ starts from 12 ]
	vecTmp1[0] = (v_p1p0.x) / mod;
	vecTmp1[1] = (v_p1p0.y) / mod;
	vecTmp1[2] = (v_p1p0.z) / mod;

	vec1 = Normal.Cross(v_p1p0);
	vecTmp2[0] = vec1.x;
	vecTmp2[1] = vec1.y;
	vecTmp2[2] = vec1.z;

	// row1
	ptr_m_J_06[Row0Idx] = -m_Delta *  vec3.x * vecTmp2[0] - m_Gamma;
	ptr_m_J_07[Row0Idx] = -m_Delta * (vec3.x * vecTmp2[1] - vecTmp1[2]);
	ptr_m_J_08[Row0Idx] = -m_Delta * (vec3.x * vecTmp2[2] + vecTmp1[1]);

	// row2
	ptr_m_J_06[Row1Idx] = -m_Delta * (vec3.y * vecTmp2[0] + vecTmp1[2]);
	ptr_m_J_07[Row1Idx] = -m_Delta *  vec3.y * vecTmp2[1] - m_Gamma;
	ptr_m_J_08[Row1Idx] = -m_Delta * (vec3.y * vecTmp2[2] - vecTmp1[0]);

	// row3
	ptr_m_J_06[Row2Idx] = -m_Delta * (vec3.z * vecTmp2[0] - vecTmp1[1]);
	ptr_m_J_07[Row2Idx] = -m_Delta * (vec3.z * vecTmp2[1] + vecTmp1[0]);
	ptr_m_J_08[Row2Idx] = -m_Delta *  vec3.z * vecTmp2[2] - m_Gamma;

	// Filling Jacobians for rigid body (if World-attachment not used)
	Vec3 wsAnchorPoint;

	if (m_Body1LIdx != -1)
	{
		CMatrix3 Body1Rot;
		CQuaternion tmpQuaterion;

		// Retrieving Bodies' Orientations
		tmpQuaterion.w = ptr_m_NodePosRot_w[m_Body1RIdx];
		tmpQuaterion.x = ptr_m_NodePosRot_x[m_Body1RIdx];
		tmpQuaterion.y = ptr_m_NodePosRot_y[m_Body1RIdx];
		tmpQuaterion.z = ptr_m_NodePosRot_z[m_Body1RIdx];
		tmpQuaterion.ToMatrix3(Body1Rot);

		// Linear Body1
		ptr_m_J_09[Row0Idx] =  1.0f;
		ptr_m_J_10[Row0Idx] =  0.0f;
		ptr_m_J_11[Row0Idx] =  0.0f;

		ptr_m_J_09[Row1Idx] =  0.0f;
		ptr_m_J_10[Row1Idx] =  1.0f;
		ptr_m_J_11[Row1Idx] =  0.0f;

		ptr_m_J_09[Row2Idx] =  0.0f;
		ptr_m_J_10[Row2Idx] =  0.0f;
		ptr_m_J_11[Row2Idx] =  1.0f;

		// Rotational Body1
		wsAnchorPoint = Body1Rot * m_BodyAnchorPoints[0];

		// -CrossProdMatrix(b1rot * b1anchor_ls)
		ptr_m_J_12[Row0Idx] = 0.0f;
		ptr_m_J_13[Row0Idx] =  wsAnchorPoint.z;
		ptr_m_J_14[Row0Idx] = -wsAnchorPoint.y;
		ptr_m_J_12[Row1Idx] = -wsAnchorPoint.z;
		ptr_m_J_13[Row1Idx] = 0.0f;
		ptr_m_J_14[Row1Idx] =  wsAnchorPoint.x;
		ptr_m_J_12[Row2Idx] =  wsAnchorPoint.y;
		ptr_m_J_13[Row2Idx] = -wsAnchorPoint.x;
		ptr_m_J_14[Row2Idx] = 0.0f;

		wsAnchorPoint.x += ptr_m_NodePosRot_x[m_Body1LIdx];
		wsAnchorPoint.y += ptr_m_NodePosRot_y[m_Body1LIdx];
		wsAnchorPoint.z += ptr_m_NodePosRot_z[m_Body1LIdx];
	}
	else
	{
		ptr_m_J_09[Row0Idx] = 0.0f;
		ptr_m_J_09[Row1Idx] = 0.0f;
		ptr_m_J_09[Row2Idx] = 0.0f;
		ptr_m_J_10[Row0Idx] = 0.0f;
		ptr_m_J_10[Row1Idx] = 0.0f;
		ptr_m_J_10[Row2Idx] = 0.0f;
		ptr_m_J_11[Row0Idx] = 0.0f;
		ptr_m_J_11[Row1Idx] = 0.0f;
		ptr_m_J_11[Row2Idx] = 0.0f;
		ptr_m_J_12[Row0Idx] = 0.0f;
		ptr_m_J_12[Row1Idx] = 0.0f;
		ptr_m_J_12[Row2Idx] = 0.0f;
		ptr_m_J_13[Row0Idx] = 0.0f;
		ptr_m_J_13[Row1Idx] = 0.0f;
		ptr_m_J_13[Row2Idx] = 0.0f;
		ptr_m_J_14[Row0Idx] = 0.0f;
		ptr_m_J_14[Row1Idx] = 0.0f;
		ptr_m_J_14[Row2Idx] = 0.0f;

		wsAnchorPoint = m_BodyAnchorPoints[0];
	}


	// Clear the rest of the Jacobi
#if (NUM_JOINT_TRIPLES == 6)

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


	// World-space coordinates of socket, attached to FE (where rigid body's ball should be)
	Vec3 FE_Socket;

	FE_Socket.x = m_Alpha * ptr_m_NodePosRot_x[m_Node1] + m_Beta * ptr_m_NodePosRot_x[m_Node2]
				+ m_Gamma * ptr_m_NodePosRot_x[m_Node3] + m_Delta * Normal.x / mod;

	FE_Socket.y = m_Alpha * ptr_m_NodePosRot_y[m_Node1] + m_Beta * ptr_m_NodePosRot_y[m_Node2]
				+ m_Gamma * ptr_m_NodePosRot_y[m_Node3] + m_Delta * Normal.y / mod;

	FE_Socket.z = m_Alpha * ptr_m_NodePosRot_z[m_Node1] + m_Beta * ptr_m_NodePosRot_z[m_Node2]
				+ m_Gamma * ptr_m_NodePosRot_z[m_Node3] + m_Delta * Normal.z / mod;

	// Baumgarte Stabilization
	float *ptr_m_J_rhs = &solver.m_J_rhs[m_StartIdx];

	*(ptr_m_J_rhs++) = idt * m_ERP * (FE_Socket.x - wsAnchorPoint.x);
	*(ptr_m_J_rhs++) = idt * m_ERP * (FE_Socket.y - wsAnchorPoint.y);
	*(ptr_m_J_rhs++) = idt * m_ERP * (FE_Socket.z - wsAnchorPoint.z);

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
					*ptr_m_JtdNodes_03 = &solver.m_JtdNodes_03[Row0Idx],
					*ptr_m_JtdNodes_04 = &solver.m_JtdNodes_04[Row0Idx];

	// First Triple Index
	*(ptr_m_JtdNodes_00++) = m_Node1;
	*(ptr_m_JtdNodes_00++) = m_Node1;
	*(ptr_m_JtdNodes_00++) = m_Node1;

	// Second Triple Index
	*(ptr_m_JtdNodes_01++) = m_Node2;
	*(ptr_m_JtdNodes_01++) = m_Node2;
	*(ptr_m_JtdNodes_01++) = m_Node2;

	// Third Triple Index
	*(ptr_m_JtdNodes_02++) = m_Node3;
	*(ptr_m_JtdNodes_02++) = m_Node3;
	*(ptr_m_JtdNodes_02++) = m_Node3;

	if (m_Body1LIdx != -1)
	{
		// Fourth Triple Index
		*(ptr_m_JtdNodes_03++) = m_Body1LIdx;
		*(ptr_m_JtdNodes_03++) = m_Body1LIdx;
		*(ptr_m_JtdNodes_03++) = m_Body1LIdx;

		// Fifth Triple Index
		*(ptr_m_JtdNodes_04++) = m_Body1RIdx;
		*(ptr_m_JtdNodes_04++) = m_Body1RIdx;
		*(ptr_m_JtdNodes_04++) = m_Body1RIdx;
	} else
	{
		// No triple -- "world" anchor

		// Fourth Triple Index
		*(ptr_m_JtdNodes_03++) = 0;
		*(ptr_m_JtdNodes_03++) = 0;
		*(ptr_m_JtdNodes_03++) = 0;

		// Fifth Triple Index
		*(ptr_m_JtdNodes_04++) = 0;
		*(ptr_m_JtdNodes_04++) = 0;
		*(ptr_m_JtdNodes_04++) = 0;
	}

#if (NUM_JOINT_TRIPLES == 6)

	unsigned int *ptr_m_JtdNodes_05 = &solver.m_JtdNodes_05[Row0Idx];

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
		if (m_Body1LIdx != -1)
		{
			*(ptr_m_JointTriples++) = 5;
		}
		else
		{
			// We're joined with World, so only 2 triples used
			*(ptr_m_JointTriples++) = 3;
		}
#endif
	}
}

void BallSocket_FEM::FetchLambdas(SolverBase &solver)
{
	m_lambda0[0] = solver.m_lambda[m_StartIdx  ];
	m_lambda0[1] = solver.m_lambda[m_StartIdx+1];
	m_lambda0[2] = solver.m_lambda[m_StartIdx+2];
}