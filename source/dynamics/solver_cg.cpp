#include "solver_cg.h"

#define FORM_A_SPARSE(j, t0, t1, t2)\
			idx = ptr_m_JtdNodes_0##j[i];\
			\
			val_arr(m_RHS) -= (	ptr_m_J_##t0[i] * ptr_m_Ftot_x[idx] + \
								ptr_m_J_##t1[i] * ptr_m_Ftot_y[idx] +\
								ptr_m_J_##t2[i] * ptr_m_Ftot_z[idx] );\
			\
			inc_arr( m_Asp_##t0 ) = dt * \
				( ptr_m_NodeInvMass_00[idx] * ptr_m_J_##t0[i] \
				+ ptr_m_NodeInvMass_01[idx] * ptr_m_J_##t1[i] \
				+ ptr_m_NodeInvMass_02[idx] * ptr_m_J_##t2[i] );\
			\
			inc_arr( m_Asp_##t1 ) = dt * \
				( ptr_m_NodeInvMass_10[idx] * ptr_m_J_##t0[i] \
				+ ptr_m_NodeInvMass_11[idx] * ptr_m_J_##t1[i] \
				+ ptr_m_NodeInvMass_12[idx] * ptr_m_J_##t2[i] );\
			\
			inc_arr( m_Asp_##t2 ) = dt * \
				( ptr_m_NodeInvMass_20[idx] * ptr_m_J_##t0[i] \
				+ ptr_m_NodeInvMass_21[idx] * ptr_m_J_##t1[i] \
				+ ptr_m_NodeInvMass_22[idx] * ptr_m_J_##t2[i] );



#define SOLVER_OBJECTIVE_FUNCTION_OUTPUT 0



void SolverCG::Solve(float dt)
{
	if (dt < 0.00001f)
	{
		return;
	}

	Prestep(dt);

	unsigned int i, idx;

	// PTRs for quicker array traverse
	float *ptr_m_a_x;
	float *ptr_m_a_y;
	float *ptr_m_a_z;

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

	float *ptr_m_lambda;
	float *ptr_m_d, *ptr_m_r;

	unsigned int *ptr_m_JtdNodes_00, *ptr_m_JtdNodes_01, *ptr_m_JtdNodes_02, *ptr_m_JtdNodes_03;

#if (NUM_JOINT_TRIPLES == 6)

	unsigned int *ptr_m_JtdNodes_04, *ptr_m_JtdNodes_05;

#endif

	float *ptr_m_J_rhs;
	float *ptr_m_ad_vec;
	float *ptr_m_CFM;
	float *ptr_m_RHS;

	float idt = 1.0f / dt;

	INIT_A_SPARSE(0);
	
	ptr_m_J_rhs = &m_J_rhs[0];

	INIT_NODES_IDX(0);

	ptr_m_RHS = &m_RHS[0];


	float	*ptr_m_NodeInvMass_00 = &m_NodeInvMass_00[0],
			*ptr_m_NodeInvMass_01 = &m_NodeInvMass_01[0],
			*ptr_m_NodeInvMass_02 = &m_NodeInvMass_02[0],
			*ptr_m_NodeInvMass_10 = &m_NodeInvMass_10[0],
			*ptr_m_NodeInvMass_11 = &m_NodeInvMass_11[0],
			*ptr_m_NodeInvMass_12 = &m_NodeInvMass_12[0],
			*ptr_m_NodeInvMass_20 = &m_NodeInvMass_20[0],
			*ptr_m_NodeInvMass_21 = &m_NodeInvMass_21[0],
			*ptr_m_NodeInvMass_22 = &m_NodeInvMass_22[0];

	float	*ptr_m_Ftot_x = &m_Ftot_x[0],
			*ptr_m_Ftot_y = &m_Ftot_y[0],
			*ptr_m_Ftot_z = &m_Ftot_z[0];

	INIT_J_SPARSE(0);

	for (i = 0; i < m_NumJoints; ++i)
	{
		val_arr(m_RHS) = inc_arr(m_J_rhs);

		FORM_A_SPARSE(0, 00, 01, 02);
		FORM_A_SPARSE(1, 03, 04, 05);
		FORM_A_SPARSE(2, 06, 07, 08);
		FORM_A_SPARSE(3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

		FORM_A_SPARSE(4, 12, 13, 14);
		FORM_A_SPARSE(5, 15, 16, 17);

#endif

		++ptr_m_RHS;
	}

	INIT_A_SPARSE(0);

	ptr_m_lambda = &m_lambda[0];

	INIT_NODES_IDX(0);

	ptr_m_a_x = &m_a_x[0];
	ptr_m_a_y = &m_a_y[0];
	ptr_m_a_z = &m_a_z[0];

	for (i = 0; i < m_NumJoints; ++i)
	{
		CALC_VECTOR_A(m_lambda, 0, 00, 01, 02);
		CALC_VECTOR_A(m_lambda, 1, 03, 04, 05);
		CALC_VECTOR_A(m_lambda, 2, 06, 07, 08);
		CALC_VECTOR_A(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

		CALC_VECTOR_A(m_lambda, 4, 12, 13, 14);
		CALC_VECTOR_A(m_lambda, 5, 15, 16, 17);

#endif

		++ptr_m_lambda;
	}

	float res_scalar = 0.0f, res_scalar2;

	// Pre-step: calculating r0, d0, res_scalar = r0^T * r0;

	ptr_m_d = &m_d[0];
	ptr_m_r = &m_r[0];

	ptr_m_a_x = &m_a_x[0];
	ptr_m_a_y = &m_a_y[0];
	ptr_m_a_z = &m_a_z[0];

	INIT_J_SPARSE(0);

	ptr_m_RHS = &m_RHS[0];

	INIT_NODES_IDX(0);

	for (i = 0; i < m_NumJoints; ++i)
	{
		unsigned int base = i * 3 * NUM_JOINT_TRIPLES;

		// Calculating residual, direction, d0 = r0 = dx = b - Ax

		val_arr(m_r) = inc_arr(m_RHS);

		// dLam -= <(J*M^-1*J^T)*lambda>
		// dLam -= <J*A*lambda>
		// dLam -= <J*a>

		unsigned int t0Idx = inc_arr(m_JtdNodes_00);
		unsigned int t1Idx = inc_arr(m_JtdNodes_01);
		unsigned int t2Idx = inc_arr(m_JtdNodes_02);
		unsigned int t3Idx = inc_arr(m_JtdNodes_03);

#if (NUM_JOINT_TRIPLES == 6)

		unsigned int t4Idx = inc_arr(m_JtdNodes_04);
		unsigned int t5Idx = inc_arr(m_JtdNodes_05);

#endif

		// J * a
		val_arr(m_r) -=	MUL_J_A();

		// CFM fuck up
		val_arr(m_r) -= m_CFM[i] * m_lambda[i];

		// Direction
		inc_arr(m_d) = val_arr(m_r);

		// Inner product
		res_scalar += val_arr(m_r) * val_arr(m_r);
		++ptr_m_r;
	}

	float dot_dAd;
	float alpha, beta;

	m_EffectiveIterations = m_Iterations;

	for (unsigned k = 0; k < m_Iterations; ++k)
	{
		dot_dAd = 0.0f;

		m_a_x.Empty();
		m_a_y.Empty();
		m_a_z.Empty();

		INIT_A_SPARSE(0);

		ptr_m_d = &m_d[0];

		INIT_NODES_IDX(0);

		for (i = 0; i < m_NumJoints; ++i)
		{
			CALC_VECTOR_A(m_d, 0, 00, 01, 02);
			CALC_VECTOR_A(m_d, 1, 03, 04, 05);
			CALC_VECTOR_A(m_d, 2, 06, 07, 08);
			CALC_VECTOR_A(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

			CALC_VECTOR_A(m_d, 4, 12, 13, 14);
			CALC_VECTOR_A(m_d, 5, 15, 16, 17);

#endif

			++ptr_m_d;
		}

		ptr_m_d = &m_d[0];

		INIT_J_SPARSE(0);

		ptr_m_CFM = &m_CFM[0];
		ptr_m_ad_vec = &m_ad_vec[0];

		INIT_NODES_IDX(0);

		// calculate dot_dAd, ad_vec
		for (i = 0; i < m_NumJoints; ++i)
		{
			unsigned int t0Idx = inc_arr(m_JtdNodes_00);
			unsigned int t1Idx = inc_arr(m_JtdNodes_01);
			unsigned int t2Idx = inc_arr(m_JtdNodes_02);
			unsigned int t3Idx = inc_arr(m_JtdNodes_03);

#if (NUM_JOINT_TRIPLES == 6)

			unsigned int t4Idx = inc_arr(m_JtdNodes_04);
			unsigned int t5Idx = inc_arr(m_JtdNodes_05);

#endif

			// J * a
			val_arr(m_ad_vec) =	MUL_J_A();

			// CFM fuck up
			val_arr(m_ad_vec) += ptr_m_CFM[i] * val_arr(m_d);

			// Calculating inner product <d, Ad>
			dot_dAd += inc_arr(m_d) * inc_arr(m_ad_vec);
		}

		// Calculate alpha
		alpha = res_scalar / dot_dAd;

		res_scalar2 = 0.0f;

		ptr_m_d = &m_d[0];
		ptr_m_r = &m_r[0];
		ptr_m_ad_vec = &m_ad_vec[0];
		ptr_m_lambda = &m_lambda[0];

		// Calculate x', r'
		for (i = 0; i < m_NumJoints; ++i)
		{
			inc_arr(m_lambda) += alpha * inc_arr(m_d);
			val_arr(m_r) -= alpha * inc_arr(m_ad_vec);

			res_scalar2 += val_arr(m_r) * val_arr(m_r);
			++ptr_m_r;
		}

		// Calculate beta
		beta = res_scalar2 / res_scalar;

		// Calculate new direction d
		ptr_m_d = &m_d[0];
		ptr_m_r = &m_r[0];

#if (SOLVER_STOPPING_CRITERIA == 1)
		float grad_norm = 0.0f;
#endif

		for (i = 0; i < m_NumJoints; ++i)
		{
#if (SOLVER_STOPPING_CRITERIA == 1)
			grad_norm += val_arr(m_r)*val_arr(m_r);
//			grad_norm += val_arr(m_d)*val_arr(m_d);
#endif

			val_arr(m_d) = inc_arr(m_r) + beta * val_arr(m_d);
			++ptr_m_d;
		}

		res_scalar = res_scalar2;


#if (SOLVER_OBJECTIVE_FUNCTION_OUTPUT == 1)

		ptr_m_r = &m_r[0];
		ptr_m_RHS = &m_RHS[0];
		ptr_m_lambda = &m_lambda[0];

		float lambda_norm = 0.0f;
		float interm_lambda_dot_grad = 0.0f;
		for (i = 0; i < m_NumJoints; ++i)
		{
			lambda_norm += val_arr(m_lambda) * val_arr(m_lambda);

			interm_lambda_dot_grad += -inc_arr(m_r) * val_arr(m_lambda);
			
			//  - x dot b
			interm_lambda_dot_grad -= inc_arr(m_RHS) * inc_arr(m_lambda);
		}

		gLog.Print("gradient: %f; alpha_cg: %12.10f;", grad_norm, alpha);
		gLog.Print("[%12.10f]---------------------------------------> 0.5 * (x * (Ax - b) - x * b) = %f;", lambda_norm, 0.5f * interm_lambda_dot_grad);

#endif

#if (SOLVER_STOPPING_CRITERIA == 1)
		// Stopping criteria
		if (grad_norm / m_NumJoints < m_Precision)
		{
			m_EffectiveIterations = k;
			break;
		}
#endif
	}

	m_a_x.Empty();
	m_a_y.Empty();
	m_a_z.Empty();

	// Calculating a = W * J^T * lambda
	INIT_A_SPARSE(0);
	
	ptr_m_lambda = &m_lambda[0];

	INIT_NODES_IDX(0);

#if (SOLVERS_ANALYZE_RESIDUAL == 1)
	ptr_m_r = &m_r[0];

	m_GradNormSq = 0.0f;
	m_LambdaNormSq = 0.0f;
	m_DotLambdaGrad = 0.0f;
#endif

	for (i = 0; i < m_NumJoints; ++i)
	{
#if (SOLVERS_ANALYZE_RESIDUAL == 1)
		// grad = -residual [resid = b - Ax]
		m_GradNormSq += val_arr(m_r) * val_arr(m_r);
		m_LambdaNormSq += val_arr(m_lambda) * val_arr(m_lambda);
		m_DotLambdaGrad += val_arr(m_lambda) * val_arr(m_r);

		++ptr_m_r;
#endif

			CALC_VECTOR_A(m_lambda, 0, 00, 01, 02);
			CALC_VECTOR_A(m_lambda, 1, 03, 04, 05);
			CALC_VECTOR_A(m_lambda, 2, 06, 07, 08);
			CALC_VECTOR_A(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

			CALC_VECTOR_A(m_lambda, 4, 12, 13, 14);
			CALC_VECTOR_A(m_lambda, 5, 15, 16, 17);

#endif

		++ptr_m_lambda;
	}

	Integrate(dt);
}

#undef FORM_A_SPARSE