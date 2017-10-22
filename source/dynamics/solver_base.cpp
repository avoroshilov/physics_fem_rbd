#include "solver_base.h"

SolverBase::SolverBase(unsigned int MaxIterations, float Precision) : m_Gravity(0.0f, -9.81f, 0.0f)
{
	m_IG_CutFlag = false;
	m_IG_CutThreshold = 10000.0f;

	m_Precision = Precision;
	m_Iterations = MaxIterations;

	m_LambdaNormSq = m_GradNormSq = m_DotLambdaGrad = 0.0f;

	m_NumJoints = 0;
	m_NumNodes = 1;

	// Initialize Node 0 (world anchor)
	m_IsRotational[0] = false;

	m_NodePosRot_x[0] = m_NodePosRot_y[0] = m_NodePosRot_z[0] = m_NodePosRot_w[0] = 0.0f;

	m_NodeInvMass0_00[0] = m_NodeInvMass0_01[0] = m_NodeInvMass0_02[0] = 0.0f;
	m_NodeInvMass0_10[0] = m_NodeInvMass0_11[0] = m_NodeInvMass0_12[0] = 0.0f;
	m_NodeInvMass0_20[0] = m_NodeInvMass0_21[0] = m_NodeInvMass0_22[0] = 0.0f;

	m_NodeInvMass_00[0] = m_NodeInvMass_01[0] = m_NodeInvMass_02[0] = 0.0f;
	m_NodeInvMass_10[0] = m_NodeInvMass_11[0] = m_NodeInvMass_12[0] = 0.0f;
	m_NodeInvMass_20[0] = m_NodeInvMass_21[0] = m_NodeInvMass_22[0] = 0.0f;

	m_NodeF_x[0] =  m_NodeF_y[0] =  m_NodeF_z[0] = 0.0f;

	m_Ftot_x[0] =  m_Ftot_y[0] =  m_Ftot_z[0] = 0.0f;
	m_NodeVel_x[0] = m_NodeVel_y[0] = m_NodeVel_z[0] = 0.0f;

	m_a_x[0] = m_a_y[0] =  m_a_z[0] = 0.0f;


	for (int i = 0; i < MAX_NUM_JACOBI; ++i)
	{
		m_lambda[i] = 0.0f;

#if (NUM_JOINT_TRIPLES == 4)

		for (int j = 0; j < 6; ++j)
		{
			m_JtdNodes_00[i] = -1;
			m_JtdNodes_01[i] = -1;
			m_JtdNodes_02[i] = -1;
			m_JtdNodes_03[i] = -1;
		}

#elif (NUM_JOINT_TRIPLES == 6)

		for (int j = 0; j < 6; ++j)
		{
			m_JtdNodes_00[i] = -1;
			m_JtdNodes_01[i] = -1;
			m_JtdNodes_02[i] = -1;
			m_JtdNodes_03[i] = -1;
			m_JtdNodes_04[i] = -1;
			m_JtdNodes_05[i] = -1;
		}

#endif
	}
}

void SolverBase::Prestep(float dt)
{
	CutLambdas();

	float *ptr_m_a_x = &m_a_x[0];
	float *ptr_m_a_y = &m_a_y[0];
	float *ptr_m_a_z = &m_a_z[0];

	float *ptr_m_Ftot_x = &m_Ftot_x[0];
	float *ptr_m_Ftot_y = &m_Ftot_y[0];
	float *ptr_m_Ftot_z = &m_Ftot_z[0];

	float *ptr_m_NodeF_x = &m_NodeF_x[0];
	float *ptr_m_NodeF_y = &m_NodeF_y[0];
	float *ptr_m_NodeF_z = &m_NodeF_z[0];

	float *ptr_m_NodeVel_x = &m_NodeVel_x[0];
	float *ptr_m_NodeVel_y = &m_NodeVel_y[0];
	float *ptr_m_NodeVel_z = &m_NodeVel_z[0];

	float *ptr_m_NodeInvMass_00 = &m_NodeInvMass_00[0];
	float *ptr_m_NodeInvMass_01 = &m_NodeInvMass_01[0];
	float *ptr_m_NodeInvMass_02 = &m_NodeInvMass_02[0];
	float *ptr_m_NodeInvMass_10 = &m_NodeInvMass_10[0];
	float *ptr_m_NodeInvMass_11 = &m_NodeInvMass_11[0];
	float *ptr_m_NodeInvMass_12 = &m_NodeInvMass_12[0];
	float *ptr_m_NodeInvMass_20 = &m_NodeInvMass_20[0];
	float *ptr_m_NodeInvMass_21 = &m_NodeInvMass_21[0];
	float *ptr_m_NodeInvMass_22 = &m_NodeInvMass_22[0];

	bool *ptr_m_IsRotational = &m_IsRotational[0];

	float tmpMass;

	for (unsigned int i = 0; i < m_NumNodes; ++i)
	{
		inc_arr(m_a_x) = 0.0f;
		inc_arr(m_a_y) = 0.0f;
		inc_arr(m_a_z) = 0.0f;

		tmpMass = val_arr(m_NodeInvMass_00);

		*ptr_m_Ftot_x = inc_arr(m_NodeVel_x) + dt *
			( inc_arr(m_NodeInvMass_00) * (*ptr_m_NodeF_x)
			+ inc_arr(m_NodeInvMass_01) * (*ptr_m_NodeF_y) 
			+ inc_arr(m_NodeInvMass_02) * (*ptr_m_NodeF_z) );

		*ptr_m_Ftot_y = inc_arr(m_NodeVel_y) + dt *
			( inc_arr(m_NodeInvMass_10) * (*ptr_m_NodeF_x)
			+ inc_arr(m_NodeInvMass_11) * (*ptr_m_NodeF_y) 
			+ inc_arr(m_NodeInvMass_12) * (*ptr_m_NodeF_z) );

		*ptr_m_Ftot_z = inc_arr(m_NodeVel_z) + dt *
			( inc_arr(m_NodeInvMass_20) * (*ptr_m_NodeF_x)
			+ inc_arr(m_NodeInvMass_21) * (*ptr_m_NodeF_y) 
			+ inc_arr(m_NodeInvMass_22) * (*ptr_m_NodeF_z) );

		if ((!(inc_arr(m_IsRotational))) && (tmpMass > GRAVITY_THRESHOLD))
		{
			*ptr_m_Ftot_x += dt * m_Gravity.x;
			*ptr_m_Ftot_y += dt * m_Gravity.y;
			*ptr_m_Ftot_z += dt * m_Gravity.z;
		}

		++ptr_m_Ftot_x;
		++ptr_m_Ftot_y;
		++ptr_m_Ftot_z;

		++ptr_m_NodeF_x;
		++ptr_m_NodeF_y;
		++ptr_m_NodeF_z;
	}
}

void SolverBase::Integrate(float dt)
{
	float *ptr_m_a_x = &m_a_x[1];
	float *ptr_m_a_y = &m_a_y[1];
	float *ptr_m_a_z = &m_a_z[1];

	float *ptr_m_Ftot_x = &m_Ftot_x[1];
	float *ptr_m_Ftot_y = &m_Ftot_y[1];
	float *ptr_m_Ftot_z = &m_Ftot_z[1];

	float *ptr_m_NodeF_x = &m_NodeF_x[1];
	float *ptr_m_NodeF_y = &m_NodeF_y[1];
	float *ptr_m_NodeF_z = &m_NodeF_z[1];

	float *ptr_m_NodeVel_x = &m_NodeVel_x[1];
	float *ptr_m_NodeVel_y = &m_NodeVel_y[1];
	float *ptr_m_NodeVel_z = &m_NodeVel_z[1];

	float	*ptr_m_NodePosRot_w = &m_NodePosRot_w[1],
			*ptr_m_NodePosRot_x = &m_NodePosRot_x[1],
			*ptr_m_NodePosRot_y = &m_NodePosRot_y[1],
			*ptr_m_NodePosRot_z = &m_NodePosRot_z[1];

	float *ptr_m_NodeInvMass_00 = &m_NodeInvMass_00[1];
	float *ptr_m_NodeInvMass_01 = &m_NodeInvMass_01[1];
	float *ptr_m_NodeInvMass_02 = &m_NodeInvMass_02[1];
	float *ptr_m_NodeInvMass_10 = &m_NodeInvMass_10[1];
	float *ptr_m_NodeInvMass_11 = &m_NodeInvMass_11[1];
	float *ptr_m_NodeInvMass_12 = &m_NodeInvMass_12[1];
	float *ptr_m_NodeInvMass_20 = &m_NodeInvMass_20[1];
	float *ptr_m_NodeInvMass_21 = &m_NodeInvMass_21[1];
	float *ptr_m_NodeInvMass_22 = &m_NodeInvMass_22[1];

	float *ptr_m_NodeInvMass0_00 = &m_NodeInvMass0_00[1];
	float *ptr_m_NodeInvMass0_01 = &m_NodeInvMass0_01[1];
	float *ptr_m_NodeInvMass0_02 = &m_NodeInvMass0_02[1];
	float *ptr_m_NodeInvMass0_10 = &m_NodeInvMass0_10[1];
	float *ptr_m_NodeInvMass0_11 = &m_NodeInvMass0_11[1];
	float *ptr_m_NodeInvMass0_12 = &m_NodeInvMass0_12[1];
	float *ptr_m_NodeInvMass0_20 = &m_NodeInvMass0_20[1];
	float *ptr_m_NodeInvMass0_21 = &m_NodeInvMass0_21[1];
	float *ptr_m_NodeInvMass0_22 = &m_NodeInvMass0_22[1];

	bool *ptr_m_IsRotational = &m_IsRotational[1];

	// Integrate [semi-implicit Euler]
	for (unsigned int i = 1; i < m_NumNodes; ++i)
	{
		val_arr(m_NodeVel_x) = inc_arr(m_a_x) + inc_arr(m_Ftot_x);
		val_arr(m_NodeVel_y) = inc_arr(m_a_y) + inc_arr(m_Ftot_y);
		val_arr(m_NodeVel_z) = inc_arr(m_a_z) + inc_arr(m_Ftot_z);

		if (inc_arr(m_IsRotational))
		{
			CQuaternion orient(	val_arr(m_NodePosRot_w), val_arr(m_NodePosRot_x),
								val_arr(m_NodePosRot_y), val_arr(m_NodePosRot_z) ),
						velocity(0.0f, inc_arr(m_NodeVel_x), inc_arr(m_NodeVel_y), inc_arr(m_NodeVel_z));

			orient += dt * 0.5f * (velocity * orient);
			orient.Normalize();

			inc_arr(m_NodePosRot_w) = orient.w;
			inc_arr(m_NodePosRot_x) = orient.x;
			inc_arr(m_NodePosRot_y) = orient.y;
			inc_arr(m_NodePosRot_z) = orient.z;

			CMatrix3 rot, tensor;
			orient.ToMatrix3(rot);

			tensor.mMatrix[0][0] = inc_arr(m_NodeInvMass0_00);
			tensor.mMatrix[0][1] = inc_arr(m_NodeInvMass0_01);
			tensor.mMatrix[0][2] = inc_arr(m_NodeInvMass0_02);
			tensor.mMatrix[1][0] = inc_arr(m_NodeInvMass0_10);
			tensor.mMatrix[1][1] = inc_arr(m_NodeInvMass0_11);
			tensor.mMatrix[1][2] = inc_arr(m_NodeInvMass0_12);
			tensor.mMatrix[2][0] = inc_arr(m_NodeInvMass0_20);
			tensor.mMatrix[2][1] = inc_arr(m_NodeInvMass0_21);
			tensor.mMatrix[2][2] = inc_arr(m_NodeInvMass0_22);

			tensor = rot * tensor * rot.GetTransposed();

			inc_arr(m_NodeInvMass_00) = tensor.mMatrix[0][0];
			inc_arr(m_NodeInvMass_01) = tensor.mMatrix[0][1];
			inc_arr(m_NodeInvMass_02) = tensor.mMatrix[0][2];
			inc_arr(m_NodeInvMass_10) = tensor.mMatrix[1][0];
			inc_arr(m_NodeInvMass_11) = tensor.mMatrix[1][1];
			inc_arr(m_NodeInvMass_12) = tensor.mMatrix[1][2];
			inc_arr(m_NodeInvMass_20) = tensor.mMatrix[2][0];
			inc_arr(m_NodeInvMass_21) = tensor.mMatrix[2][1];
			inc_arr(m_NodeInvMass_22) = tensor.mMatrix[2][2];
		}
		else
		{
			++ptr_m_NodePosRot_w;
			inc_arr(m_NodePosRot_x) += inc_arr(m_NodeVel_x) * dt;
			inc_arr(m_NodePosRot_y) += inc_arr(m_NodeVel_y) * dt;
			inc_arr(m_NodePosRot_z) += inc_arr(m_NodeVel_z) * dt;
		}

		inc_arr(m_NodeF_x) = 0.0f;
		inc_arr(m_NodeF_y) = 0.0f;
		inc_arr(m_NodeF_z) = 0.0f;
	}
}

void SolverBase::CutLambdas()
{
	unsigned int i;

	if (m_IG_CutFlag)
	{
		if (m_LambdaNormSq / m_NumJoints > m_IG_CutThreshold)
		{
			for (i = 0; i < m_NumJoints; ++i)
			{
				m_lambda[i] = 0.0f;
			}
		}
	}
}

float SolverBase::CalcSylvesterCriterion()
{
	// J[CDoFs x 3*n] * M^-1 [3*n x 3*n] * J^T [3*n x CDoFs]
	float *A_Matrix = new float[m_NumJoints * m_NumJoints];

	float	*L = new float [m_NumJoints * m_NumJoints],
			*D = new float [m_NumJoints];

#define AMat(i, j) (A_Matrix[(i) * m_NumJoints + (j)])
#define LMat(i, j) (L[(i) * m_NumJoints + (j)])

	unsigned int i, j, k;

	// Set A_Matrix to zero
	memset(A_Matrix, 0, m_NumJoints * m_NumJoints * sizeof(float));

	// Asp already calculated on "Solve" step
	unsigned int	*ptr_JtdNodes_00 = &m_JtdNodes_00[0],
					*ptr_JtdNodes_01 = &m_JtdNodes_01[0],
					*ptr_JtdNodes_02 = &m_JtdNodes_02[0],
					*ptr_JtdNodes_03 = &m_JtdNodes_03[0];

#if (NUM_JOINT_TRIPLES == 6)

	unsigned int	*ptr_JtdNodes_04 = &m_JtdNodes_04[0],
					*ptr_JtdNodes_05 = &m_JtdNodes_05[0];

#endif

	float	*ptr_m_Asp_00, *ptr_m_Asp_01, *ptr_m_Asp_02, *ptr_m_Asp_03, *ptr_m_Asp_04, *ptr_m_Asp_05,
			*ptr_m_Asp_06, *ptr_m_Asp_07, *ptr_m_Asp_08, *ptr_m_Asp_09, *ptr_m_Asp_10, *ptr_m_Asp_11;
	
#if (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_Asp_12, *ptr_m_Asp_13, *ptr_m_Asp_14, *ptr_m_Asp_15, *ptr_m_Asp_16, *ptr_m_Asp_17;

#endif

	float	*ptr_m_J_00, *ptr_m_J_01, *ptr_m_J_02, *ptr_m_J_03, *ptr_m_J_04, *ptr_m_J_05,
			*ptr_m_J_06, *ptr_m_J_07, *ptr_m_J_08, *ptr_m_J_09, *ptr_m_J_10, *ptr_m_J_11;
	
#if (NUM_JOINT_TRIPLES == 6)

	float	*ptr_m_J_12, *ptr_m_J_13, *ptr_m_J_14, *ptr_m_J_15, *ptr_m_J_16, *ptr_m_J_17;

#endif


	ptr_m_Asp_00 = &m_Asp_00[0];
	ptr_m_Asp_01 = &m_Asp_01[0];
	ptr_m_Asp_02 = &m_Asp_02[0];
	ptr_m_Asp_03 = &m_Asp_03[0];
	ptr_m_Asp_04 = &m_Asp_04[0];
	ptr_m_Asp_05 = &m_Asp_05[0];
	ptr_m_Asp_06 = &m_Asp_06[0];
	ptr_m_Asp_07 = &m_Asp_07[0];
	ptr_m_Asp_08 = &m_Asp_08[0];
	ptr_m_Asp_09 = &m_Asp_09[0];
	ptr_m_Asp_10 = &m_Asp_10[0];
	ptr_m_Asp_11 = &m_Asp_11[0];

#if (NUM_JOINT_TRIPLES == 6)

	ptr_m_Asp_12 = &m_Asp_12[0];
	ptr_m_Asp_13 = &m_Asp_13[0];
	ptr_m_Asp_14 = &m_Asp_14[0];
	ptr_m_Asp_15 = &m_Asp_15[0];
	ptr_m_Asp_16 = &m_Asp_16[0];
	ptr_m_Asp_17 = &m_Asp_17[0];

#endif


	ptr_m_J_00 = &m_J_00[0];
	ptr_m_J_01 = &m_J_01[0];
	ptr_m_J_02 = &m_J_02[0];
	ptr_m_J_03 = &m_J_03[0];
	ptr_m_J_04 = &m_J_04[0];
	ptr_m_J_05 = &m_J_05[0];
	ptr_m_J_06 = &m_J_06[0];
	ptr_m_J_07 = &m_J_07[0];
	ptr_m_J_08 = &m_J_08[0];
	ptr_m_J_09 = &m_J_09[0];
	ptr_m_J_10 = &m_J_10[0];
	ptr_m_J_11 = &m_J_11[0];

#if (NUM_JOINT_TRIPLES == 6)

	ptr_m_J_12 = &m_J_12[0];
	ptr_m_J_13 = &m_J_13[0];
	ptr_m_J_14 = &m_J_14[0];
	ptr_m_J_15 = &m_J_15[0];
	ptr_m_J_16 = &m_J_16[0];
	ptr_m_J_17 = &m_J_17[0];

#endif


	for (i = 0; i < m_NumJoints; ++i)
	{
		for (j = 0; j < m_NumJoints; ++j)
		{
			AMat(i, j) = 0.0f;

#define CALCULATE_JWJ(k, l, kt0, kt1, kt2, lt0, lt1, lt2)\
					if (ptr_JtdNodes_0##k[i] == ptr_JtdNodes_0##l[j])\
					{\
						int node_idx = ptr_JtdNodes_0##k[i];\
						\
						AMat(i, j) += ptr_m_J_##kt0[i] * ptr_m_Asp_##lt0[j]\
									+ ptr_m_J_##kt1[i] * ptr_m_Asp_##lt1[j]\
									+ ptr_m_J_##kt2[i] * ptr_m_Asp_##lt2[j];\
					}

//////////////////////////////////////////////////////////////////////////
// NUM_JOINT_TRIPPLES == 4
//////////////////////////////////////////////////////////////////////////
#if (NUM_JOINT_TRIPLES == 4)

#define CALCULATE_JWJ_L(k, kt0, kt1, kt2)\
					CALCULATE_JWJ(k, 0, kt0, kt1, kt2, 00, 01, 02);\
					CALCULATE_JWJ(k, 1, kt0, kt1, kt2, 03, 04, 05);\
					CALCULATE_JWJ(k, 2, kt0, kt1, kt2, 06, 07, 08);\
					CALCULATE_JWJ(k, 3, kt0, kt1, kt2, 09, 10, 11);

			CALCULATE_JWJ_L(0, 00, 01, 02);
			CALCULATE_JWJ_L(1, 03, 04, 05);
			CALCULATE_JWJ_L(2, 06, 07, 08);
			CALCULATE_JWJ_L(3, 09, 10, 11);

#undef CALCULATE_JWJ_L

//////////////////////////////////////////////////////////////////////////
// NUM_JOINT_TRIPPLES == 6
//////////////////////////////////////////////////////////////////////////
#elif (NUM_JOINT_TRIPLES == 6)

#define CALCULATE_JWJ_L(k, kt0, kt1, kt2)\
					CALCULATE_JWJ(k, 0, kt0, kt1, kt2, 00, 01, 02);\
					CALCULATE_JWJ(k, 1, kt0, kt1, kt2, 03, 04, 05);\
					CALCULATE_JWJ(k, 2, kt0, kt1, kt2, 06, 07, 08);\
					CALCULATE_JWJ(k, 3, kt0, kt1, kt2, 09, 10, 11);\
					CALCULATE_JWJ(k, 4, kt0, kt1, kt2, 12, 13, 14);\
					CALCULATE_JWJ(k, 5, kt0, kt1, kt2, 15, 16, 17);

			CALCULATE_JWJ_L(0, 00, 01, 02);
			CALCULATE_JWJ_L(1, 03, 04, 05);
			CALCULATE_JWJ_L(2, 06, 07, 08);
			CALCULATE_JWJ_L(3, 09, 10, 11);
			CALCULATE_JWJ_L(4, 12, 13, 14);
			CALCULATE_JWJ_L(5, 15, 16, 17);

#undef CALCULATE_JWJ_L

#endif

#undef CALCULATE_JWJ

		}
	}

	// CFM fuck up
	float *CFM_ptr = &m_CFM[0];
	for (i = 0; i < m_NumJoints; ++i)
	{
		AMat(i, i) += CFM_ptr[i];
	}

	for (i = 0; i < m_NumJoints; ++i)
	{
		for (j = 0; j < m_NumJoints; ++j)
		{
			AMat(i, j) *= m_invDiag[i];
		}
	}

	// Symmetry testing
	/* DBG: PRINT */
	int Sym = 1;

#define DBG_PRINT_MATRIX	0

#if (DBG_PRINT_MATRIX == 1)
	FILE *fp_matrix_log;
	fp_matrix_log = fopen("matrix.log", "w");
	fprintf(fp_matrix_log, "{");
#endif
	for (i = 0; i < m_NumJoints; ++i)
	{
#if (DBG_PRINT_MATRIX == 1)
		fprintf(fp_matrix_log, "{");
		for (j = 0; j < m_NumJoints; ++j)
		{
			if (j != m_NumJoints - 1)
				fprintf(fp_matrix_log, "%f,", AMat(i, j));
			else
				fprintf(fp_matrix_log, "%f", AMat(i, j));
		}
		if (i != m_NumJoints - 1)
			fprintf(fp_matrix_log, "},");
		else
			fprintf(fp_matrix_log, "}");
#endif

		for (j = i + 1; j < m_NumJoints; ++j)
		{
			float tmp1 = AMat(i, j);
			float tmp2 = AMat(j, i);
			if (math::fastabs(AMat(i, j) - AMat(j, i)) > 0.01f)
			{
				assert(false);
				Sym = 0;
				break;
			}
		}
	}
#if (DBG_PRINT_MATRIX == 1)
	fprintf(fp_matrix_log, "}\n\n");
	fclose(fp_matrix_log);
#endif

	if (Sym == 0)
	{
		gLog.Print("Matrix is not self-conjugate!");
	}
	else
	{
		gLog.Print("Matrix IS self-conjugate!");
	}

	int PD = 1;
	int CurDetDimension;
	float lastDet = 0.0f;

	// Determinant, Sylvester Criterion
	// We must find NumJoints determinants for the upper left corner of matrix
	// [ Determinant is found via Cholesky Decomposition and U-matrix diagonal mult ]

	// Cholesky Decomposition [ WIKI version ]

	gLog.Print("Diagonal in Cholesky Decomposition:");

	lastDet = 1.0f;

	// General case: i = 0 .. n
	for (i = 0; i < m_NumJoints; ++i)
	{
		float temp_var;

		for (j = 0; j < i; ++j)
		{
			temp_var = AMat(i, j);
			for (k = 0; k < j; ++k)
			{
				temp_var -= LMat(i, k) * LMat(j, k) * D[k];
			}

			LMat(i, j) = temp_var / D[j];
			LMat(j, i) = 0.0f;
		}

		temp_var = AMat(i, i);
		for (k = 0; k < i; ++k)
		{
			temp_var -= LMat(i, k) * LMat(i, k) * D[k];
		}

		D[i] = temp_var;
		lastDet *= temp_var;
		gLog.Print("%d: %f; [%f]", i, temp_var, lastDet);
		LMat(i, i) = 1.0f;
	}

	gLog.Print("\n\n");

	// Determinant: det(L * D * L^T) = det L * det D * det L = (MUL Lii) ^2 * (MUL Dii)
	lastDet = 1.0f;
	for (i = 0; i < m_NumJoints; ++i)
	{
		lastDet *= LMat(i, i) * D[i] * LMat(i, i);

		if (D[i] < 0.0f)
		{
			CurDetDimension = i;
			PD = 0;
			break;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Test if A == L * D * L^T
	float *LD = new float [m_NumJoints * m_NumJoints];

#define LDMat(i, j) (LD[(i) * m_NumJoints + (j)])

	for (i = 0; i < m_NumJoints; ++i)
	{
		for (j = 0; j < m_NumJoints; ++j)
		{
			LDMat(i, j) = LMat(i, j) * D[j];
		}
	}

	int DecompEq = 1;
	for (i = 0; i < m_NumJoints; ++i)
	{
		for (j = 0; j < m_NumJoints; ++j)
		{
			float Res = 0.0f;
			for (k = 0; k < m_NumJoints; ++k)
			{
				// LD * L^T
				Res += LDMat(i, k) * LMat(j, k);
			}

			float tmp = AMat(i, j);
			if (math::fastabs(Res - tmp) > 0.1f)
			{
				assert(false);
				DecompEq = 0;
				break;
			}
		}
	}

	if (DecompEq == 0)
	{
		gLog.Print("Decomposition is INCORRECT!");
	}
	else
	{
		gLog.Print("Decomposition is correct!");
	}

#undef LDMat

	delete [] LD;
	//////////////////////////////////////////////////////////////////////////

#undef RMat
#undef AMat

	delete [] A_Matrix;

	delete [] L;
	delete [] D;

	if (PD == 0)
	{
		gLog.Print("Sylvester's criterion failed (step %d, det * 1000.0 = %.16f)!", CurDetDimension, lastDet);
	}
	else
	{
		gLog.Print("Sylvester's criterion OK, det A * 1000.0 = %.16f;", lastDet);
	}

	return lastDet;
}


unsigned int SolverBase::AddRotationalNode(const CMatrix3 &invInertia0, const CQuaternion &orient, const Vec3 &angvel)
{
	unsigned int i = m_NumNodes++;

	m_IsRotational[i] = true;

	m_NodePosRot_w[i] = orient.w;
	m_NodePosRot_x[i] = orient.x;
	m_NodePosRot_y[i] = orient.y;
	m_NodePosRot_z[i] = orient.z;

	m_NodeVel_x[i] = angvel.x;
	m_NodeVel_y[i] = angvel.y;
	m_NodeVel_z[i] = angvel.z;

	m_NodeF_x[i] = 0.0f;
	m_NodeF_y[i] = 0.0f;
	m_NodeF_z[i] = 0.0f;

	m_NodeInvMass0_00[i] = invInertia0.mMatrix[0][0];
	m_NodeInvMass0_01[i] = invInertia0.mMatrix[0][1];
	m_NodeInvMass0_02[i] = invInertia0.mMatrix[0][2];
	m_NodeInvMass0_10[i] = invInertia0.mMatrix[1][0];
	m_NodeInvMass0_11[i] = invInertia0.mMatrix[1][1];
	m_NodeInvMass0_12[i] = invInertia0.mMatrix[1][2];
	m_NodeInvMass0_20[i] = invInertia0.mMatrix[2][0];
	m_NodeInvMass0_21[i] = invInertia0.mMatrix[2][1];
	m_NodeInvMass0_22[i] = invInertia0.mMatrix[2][2];

	CMatrix3 rot, tensor;
	orient.ToMatrix3(rot);

	tensor = rot * invInertia0 * rot.GetTransposed();

	m_NodeInvMass_00[i] = tensor.mMatrix[0][0];
	m_NodeInvMass_01[i] = tensor.mMatrix[0][1];
	m_NodeInvMass_02[i] = tensor.mMatrix[0][2];
	m_NodeInvMass_10[i] = tensor.mMatrix[1][0];
	m_NodeInvMass_11[i] = tensor.mMatrix[1][1];
	m_NodeInvMass_12[i] = tensor.mMatrix[1][2];
	m_NodeInvMass_20[i] = tensor.mMatrix[2][0];
	m_NodeInvMass_21[i] = tensor.mMatrix[2][1];
	m_NodeInvMass_22[i] = tensor.mMatrix[2][2];

	return i;
}

unsigned int SolverBase::AddTranslationalNode(float invMass, const Vec3 &pos, const Vec3 &vel)
{
	unsigned int i = m_NumNodes++;

	m_IsRotational[i] = false;

	m_NodePosRot_w[i] = 0.0f;
	m_NodePosRot_x[i] = pos.x;
	m_NodePosRot_y[i] = pos.y;
	m_NodePosRot_z[i] = pos.z;

	m_NodeVel_x[i] = vel.x;
	m_NodeVel_y[i] = vel.y;
	m_NodeVel_z[i] = vel.z;

	m_NodeF_x[i] = 0.0f;
	m_NodeF_y[i] = 0.0f;
	m_NodeF_z[i] = 0.0f;

	m_NodeInvMass_00[i] = m_NodeInvMass0_00[i] = invMass;
	m_NodeInvMass_01[i] = m_NodeInvMass0_01[i] = 0.0f;
	m_NodeInvMass_02[i] = m_NodeInvMass0_02[i] = 0.0f;
	m_NodeInvMass_10[i] = m_NodeInvMass0_10[i] = 0.0f;
	m_NodeInvMass_11[i] = m_NodeInvMass0_11[i] = invMass;
	m_NodeInvMass_12[i] = m_NodeInvMass0_12[i] = 0.0f;
	m_NodeInvMass_20[i] = m_NodeInvMass0_20[i] = 0.0f;
	m_NodeInvMass_21[i] = m_NodeInvMass0_21[i] = 0.0f;
	m_NodeInvMass_22[i] = m_NodeInvMass0_22[i] = invMass;

	return i;
}

void SolverBase::SetForce(unsigned int node, float x, float y, float z)
{
	assert(node < m_NumNodes);

	m_NodeF_x[node] = x;
	m_NodeF_y[node] = y;
	m_NodeF_z[node] = z;
}

void SolverBase::GetForce(unsigned int node, float &x, float &y, float &z) const
{
	assert(node < m_NumNodes);

	x = m_NodeF_x[node];
	y = m_NodeF_y[node];
	z = m_NodeF_z[node];
}

//inline Vec3 SolverBase::GetPosition(unsigned int NodeIdx) const
void SolverBase::SetPosition(unsigned int NodeIdx, const Vec3 &Pos)
{
	assert(NodeIdx < m_NumNodes);

	m_NodePosRot_w[NodeIdx] = 0.0f;
	m_NodePosRot_x[NodeIdx] = Pos.x;
	m_NodePosRot_y[NodeIdx] = Pos.y;
	m_NodePosRot_z[NodeIdx] = Pos.z;
}

//inline CQuaternion SolverBase::GetOrientation(unsigned int NodeIdx) const
void SolverBase::SetOrientation(unsigned int NodeIdx, const CQuaternion &Quat)
{
	assert(NodeIdx < m_NumNodes);

	m_NodePosRot_w[NodeIdx] = Quat.w;
	m_NodePosRot_x[NodeIdx] = Quat.x;
	m_NodePosRot_y[NodeIdx] = Quat.y;
	m_NodePosRot_z[NodeIdx] = Quat.z;
}
