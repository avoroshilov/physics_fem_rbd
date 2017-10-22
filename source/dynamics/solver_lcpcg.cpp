#include "solver_lcpcg.h"

#define SW_LCPCG_PRECONDITIONER_LEFTJACOBI		0
#define SW_LCPCG_PRECONDITIONER_FACEJACOBI		0
#define SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI	1

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

#define CALC_VECTOR_A_PRECOND(vec, j, t0, t1, t2)\
			idx = inc_arr( m_JtdNodes_0##j );\
			\
			ptr_m_a_x[idx] += inc_arr(m_Asp_##t0) * (ptr_m_invDiag[i] * (val_arr(vec)));\
			ptr_m_a_y[idx] += inc_arr(m_Asp_##t1) * (ptr_m_invDiag[i] * (val_arr(vec)));\
			ptr_m_a_z[idx] += inc_arr(m_Asp_##t2) * (ptr_m_invDiag[i] * (val_arr(vec)));

#endif


#define FORM_A_SPARSE(j, t0, t1, t2)\
			idx = ptr_m_JtdNodes_0##j[i];\
			\
			val_arr(m_RHS) -= (	ptr_m_J_##t0[i] * ptr_m_Ftot_x[idx] + \
								ptr_m_J_##t1[i] * ptr_m_Ftot_y[idx] + \
								ptr_m_J_##t2[i] * ptr_m_Ftot_z[idx] );\
			\
			val_arr( m_Asp_##t0 ) = dt * \
				( ptr_m_NodeInvMass_00[idx] * ptr_m_J_##t0[i] \
				+ ptr_m_NodeInvMass_01[idx] * ptr_m_J_##t1[i] \
				+ ptr_m_NodeInvMass_02[idx] * ptr_m_J_##t2[i] );\
			ptr_m_invDiag[i] += ptr_m_J_##t0[i] * inc_arr( m_Asp_##t0 );\
			\
			val_arr( m_Asp_##t1 ) = dt * \
				( ptr_m_NodeInvMass_10[idx] * ptr_m_J_##t0[i] \
				+ ptr_m_NodeInvMass_11[idx] * ptr_m_J_##t1[i] \
				+ ptr_m_NodeInvMass_12[idx] * ptr_m_J_##t2[i] );\
			ptr_m_invDiag[i] += ptr_m_J_##t1[i] * inc_arr( m_Asp_##t1 );\
			\
			val_arr( m_Asp_##t2 ) = dt * \
				( ptr_m_NodeInvMass_20[idx] * ptr_m_J_##t0[i] \
				+ ptr_m_NodeInvMass_21[idx] * ptr_m_J_##t1[i] \
				+ ptr_m_NodeInvMass_22[idx] * ptr_m_J_##t2[i] );\
			ptr_m_invDiag[i] += ptr_m_J_##t2[i] * inc_arr( m_Asp_##t2 );

#if (NUM_JOINT_TRIPLES == 4)

#define ADVANCE_A_SPARSE_NODES_IDX()\
		++ptr_m_Asp_00;\
		++ptr_m_Asp_01;\
		++ptr_m_Asp_02;\
		++ptr_m_Asp_03;\
		++ptr_m_Asp_04;\
		++ptr_m_Asp_05;\
		++ptr_m_Asp_06;\
		++ptr_m_Asp_07;\
		++ptr_m_Asp_08;\
		++ptr_m_Asp_09;\
		++ptr_m_Asp_10;\
		++ptr_m_Asp_11;\
		++ptr_m_JtdNodes_00;\
		++ptr_m_JtdNodes_01;\
		++ptr_m_JtdNodes_02;\
		++ptr_m_JtdNodes_03;

#elif (NUM_JOINT_TRIPLES == 6)

#define ADVANCE_A_SPARSE_NODES_IDX()\
		++ptr_m_Asp_00;\
		++ptr_m_Asp_01;\
		++ptr_m_Asp_02;\
		++ptr_m_Asp_03;\
		++ptr_m_Asp_04;\
		++ptr_m_Asp_05;\
		++ptr_m_Asp_06;\
		++ptr_m_Asp_07;\
		++ptr_m_Asp_08;\
		++ptr_m_Asp_09;\
		++ptr_m_Asp_10;\
		++ptr_m_Asp_11;\
		++ptr_m_Asp_12;\
		++ptr_m_Asp_13;\
		++ptr_m_Asp_14;\
		++ptr_m_Asp_15;\
		++ptr_m_Asp_16;\
		++ptr_m_Asp_17;\
		++ptr_m_JtdNodes_00;\
		++ptr_m_JtdNodes_01;\
		++ptr_m_JtdNodes_02;\
		++ptr_m_JtdNodes_03;\
		++ptr_m_JtdNodes_04;\
		++ptr_m_JtdNodes_05;

#endif



#define PROPORTIONING_DOUBLE_STEP 1


void SolverLCPCG::Solve(float dt)
{
	unsigned i, idx;

	if (dt < 0.00001f)
		return;

	Prestep(dt);

	float	ACT_SET_THRESHOLD = 0.001f,
			ZERO_THRESHOLD = 0.0f;

	float BetaNorm, Phi_Phi;

	// Additional variables for calculations
	float tmp_float, tmp_float2, tmp_float3;

	// Alpha_fucked should be in (0.0f, 2 * ||A||^-1];
	float alpha_fucked = 0.02f, alpha_cg, beta;

	float GAMMA = 1.0f;
	float GAMMA_Sq = GAMMA * GAMMA;

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
	float *ptr_m_g, *ptr_m_p, *ptr_m_d;

	unsigned int *ptr_m_JtdNodes_00, *ptr_m_JtdNodes_01, *ptr_m_JtdNodes_02, *ptr_m_JtdNodes_03;

#if (NUM_JOINT_TRIPLES == 6)

	unsigned int *ptr_m_JtdNodes_04, *ptr_m_JtdNodes_05;

#endif

	float *ptr_m_J_rhs;
	float *ptr_m_ad_vec, *ptr_m_ap_vec;
	float *ptr_m_CFM;
	float *ptr_m_RHS;
	float *ptr_m_Lo, *ptr_m_Hi;

	float idt = 1.0f / dt;

	INIT_J_SPARSE(0);

	ptr_m_d = &m_d[0];

	INIT_A_SPARSE(0);

	ptr_m_RHS = &m_RHS[0];
	ptr_m_J_rhs = &m_J_rhs[0];

	INIT_NODES_IDX(0);

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

	float	*ptr_m_invDiag = &m_invDiag[0];

	float diagNorm = 0.0f;

	ptr_m_CFM = &m_CFM[0];

	// y = [1 1 1 1 <...> 1]^T -- for power iteration step
	for (i = 0; i < m_NumJoints; ++i)
	{
		val_arr(m_RHS) = inc_arr(m_J_rhs);

		inc_arr(m_d) = 1.0f;

		ptr_m_invDiag[i] = 0.0f;

		FORM_A_SPARSE(0, 00, 01, 02);
		FORM_A_SPARSE(1, 03, 04, 05);
		FORM_A_SPARSE(2, 06, 07, 08);
		FORM_A_SPARSE(3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

		FORM_A_SPARSE(4, 12, 13, 14);
		FORM_A_SPARSE(5, 15, 16, 17);

#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
		ptr_m_invDiag[i] = 1.0f / sqrtf(ptr_m_invDiag[i] + ptr_m_CFM[i]);
#else
		ptr_m_invDiag[i] = 1.0f / (ptr_m_invDiag[i] + ptr_m_CFM[i]);
#endif

		//assert(ptr_m_invDiag[i] >= 0.0f);

		diagNorm += ptr_m_invDiag[i] * ptr_m_invDiag[i];

		++ptr_m_RHS;
	}

	// First norm value should be 1.0
	float y_norm_inv = 1.0f, y_norm_sq = 0.0001f;

	// d - temporary storage, having same size
	float interm_y;

	// Power iteration pre-step
	for (unsigned int pow_iters = 0; pow_iters < 15; ++pow_iters)
	{
		ptr_m_d = &m_d[0];

		ptr_m_a_x = &m_a_x[0];
		ptr_m_a_y = &m_a_y[0];
		ptr_m_a_z = &m_a_z[0];

		ptr_m_Lo = &m_Lo[0];
		ptr_m_Hi = &m_Hi[0];
		
		INIT_A_SPARSE(0);
		INIT_NODES_IDX(0);

		m_a_x.Empty();
		m_a_y.Empty();
		m_a_z.Empty();

		// Calculate a = W * J^T * y;
		for (i = 0; i < m_NumJoints; ++i)
		{
			val_arr(m_d) *= y_norm_inv;

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

			CALC_VECTOR_A_PRECOND(m_d, 0, 00, 01, 02);
			CALC_VECTOR_A_PRECOND(m_d, 1, 03, 04, 05);
			CALC_VECTOR_A_PRECOND(m_d, 2, 06, 07, 08);
			CALC_VECTOR_A_PRECOND(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

			CALC_VECTOR_A_PRECOND(m_d, 4, 12, 13, 14);
			CALC_VECTOR_A_PRECOND(m_d, 5, 15, 16, 17);

#endif

#else

			CALC_VECTOR_A(m_d, 0, 00, 01, 02);
			CALC_VECTOR_A(m_d, 1, 03, 04, 05);
			CALC_VECTOR_A(m_d, 2, 06, 07, 08);
			CALC_VECTOR_A(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

			CALC_VECTOR_A(m_d, 4, 12, 13, 14);
			CALC_VECTOR_A(m_d, 5, 15, 16, 17);

#endif

#endif

			++ptr_m_d;
		}

		y_norm_sq = 0.0f;

		ptr_m_d = &m_d[0];

		ptr_m_a_x = &m_a_x[0];
		ptr_m_a_y = &m_a_y[0];
		ptr_m_a_z = &m_a_z[0];

		ptr_m_CFM = &m_CFM[0];
		ptr_m_RHS = &m_RHS[0];
		ptr_m_lambda = &m_lambda[0];

		INIT_J_SPARSE(0);
		INIT_NODES_IDX(0);

		for (i = 0; i < m_NumJoints; ++i)
		{
			interm_y = val_arr(m_d);

			unsigned int t0Idx = inc_arr(m_JtdNodes_00);
			unsigned int t1Idx = inc_arr(m_JtdNodes_01);
			unsigned int t2Idx = inc_arr(m_JtdNodes_02);
			unsigned int t3Idx = inc_arr(m_JtdNodes_03);

#if (NUM_JOINT_TRIPLES == 6)

			unsigned int t4Idx = inc_arr(m_JtdNodes_04);
			unsigned int t5Idx = inc_arr(m_JtdNodes_05);

#endif

			// J * a
			val_arr(m_d) = MUL_J_A();

			// CFM fuck up
			val_arr(m_d) += ptr_m_CFM[i] * interm_y;

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
			// Preconditioning
			val_arr(m_d) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
			// Double Jacobi:
			// D^-1/2 * (A * vec') = D^-1/2 * A * D^-1/2 * vec
			val_arr(m_d) *= ptr_m_invDiag[i];
#endif

// 			assert(val_arr(m_d) == val_arr(m_d));
// 			assert(val_arr(m_d) < FLT_MAX);
// 			assert(val_arr(m_d) > -FLT_MAX);

			y_norm_sq += val_arr(m_d) * val_arr(m_d);

			++ptr_m_d;
		}

		// Calculate ||y||^-1
		y_norm_inv = 1.0f / sqrtf(y_norm_sq);
	}
	// We have already y = Ax, so ||A|| == ||y|| (which is already calculated too)
	alpha_fucked = 1.6f / sqrtf(y_norm_sq);

	float prev_grad_norm = FLT_MAX;

	//////////////////////////////////////////////////////////////////////////
	// CG-based solver itself
	//////////////////////////////////////////////////////////////////////////
	ptr_m_a_x = &m_a_x[0];
	ptr_m_a_y = &m_a_y[0];
	ptr_m_a_z = &m_a_z[0];

	ptr_m_Lo = &m_Lo[0];
	ptr_m_Hi = &m_Hi[0];

	ptr_m_lambda = &m_lambda[0];

	INIT_A_SPARSE(0);
	INIT_NODES_IDX(0);

	// Erase data from Power Iteration in a
	m_a_x.Empty();
	m_a_y.Empty();
	m_a_z.Empty();

	// Calculate a = W * J^T * lambda0;
	for (i = 0; i < m_NumJoints; ++i)
	{
#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
		// Double Jacobi:
		// multiply IG and LoHi by D^1/2 to scale to the preconditioned unknown

		val_arr(m_Lo) /= ptr_m_invDiag[i];
		val_arr(m_Hi) /= ptr_m_invDiag[i];

		val_arr(m_lambda) /= ptr_m_invDiag[i];
#endif

		// Project initial guess so that it fit criteria
		if (val_arr(m_lambda) < val_arr(m_Lo))
		{
			val_arr(m_lambda) = val_arr(m_Lo);
		}
		else if (val_arr(m_lambda) > val_arr(m_Hi))
		{
			val_arr(m_lambda) = val_arr(m_Hi);
		}

		++ptr_m_Lo;
		++ptr_m_Hi;

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

		CALC_VECTOR_A_PRECOND(m_lambda, 0, 00, 01, 02);
		CALC_VECTOR_A_PRECOND(m_lambda, 1, 03, 04, 05);
		CALC_VECTOR_A_PRECOND(m_lambda, 2, 06, 07, 08);
		CALC_VECTOR_A_PRECOND(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

		CALC_VECTOR_A_PRECOND(m_lambda, 4, 12, 13, 14);
		CALC_VECTOR_A_PRECOND(m_lambda, 5, 15, 16, 17);

#endif

#else

		CALC_VECTOR_A(m_lambda, 0, 00, 01, 02);
		CALC_VECTOR_A(m_lambda, 1, 03, 04, 05);
		CALC_VECTOR_A(m_lambda, 2, 06, 07, 08);
		CALC_VECTOR_A(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

		CALC_VECTOR_A(m_lambda, 4, 12, 13, 14);
		CALC_VECTOR_A(m_lambda, 5, 15, 16, 17);

#endif

#endif

		++ptr_m_lambda;
	}

	ptr_m_p = &m_p[0];
	ptr_m_g = &m_g[0];

	ptr_m_a_x = &m_a_x[0];
	ptr_m_a_y = &m_a_y[0];
	ptr_m_a_z = &m_a_z[0];

	ptr_m_Lo = &m_Lo[0];
	ptr_m_Hi = &m_Hi[0];
	ptr_m_RHS = &m_RHS[0];
	ptr_m_lambda = &m_lambda[0];

	INIT_J_SPARSE(0);
	INIT_NODES_IDX(0);

	BetaNorm = 0.0f;
	Phi_Phi = 0.0f;

	// Calculate g = J * a - b; p = phi(lambda0);
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
		val_arr(m_g) = MUL_J_A();

		// CFM fuck up
		val_arr(m_g) += m_CFM[i] * ptr_m_lambda[i];
		val_arr(m_g) -= inc_arr(m_RHS);

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
		// Preconditioning
		val_arr(m_g) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
		val_arr(m_g) *= ptr_m_invDiag[i];
#endif

// 		assert(val_arr(m_g) < 1e2);
// 
// 		assert(val_arr(m_g) == val_arr(m_g));
// 		assert(val_arr(m_g) < FLT_MAX);
// 		assert(val_arr(m_g) > -FLT_MAX);

		// Calculate stepping direction p (CG)

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
		// Preconditioner in face:
		// p = z
		val_arr(m_p) = ptr_m_invDiag[i] * val_arr(m_g);
#else
		val_arr(m_p) = val_arr(m_g);
#endif

		// UPD: LoHi
		tmp_float2 = ptr_m_lambda[i] - ptr_m_Lo[i];
		tmp_float3 = ptr_m_lambda[i] - ptr_m_Hi[i];

		// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
		if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
		{
			val_arr(m_p) = 0.0f;

			// Active set
			// Here we calculate chopped gradient Beta, hence BetaNorm updated
			// beta = g- = min{g, 0.0f};

			tmp_float = inc_arr(m_g);

			if (tmp_float < 0.0f)
			{
#if (PROPORTIONING_DOUBLE_STEP == 1)

				if (tmp_float3 / alpha_fucked > tmp_float)
				{
					BetaNorm += tmp_float * tmp_float3 / alpha_fucked;
				}
				else
				{
					BetaNorm += tmp_float * tmp_float;
				}

#else

				BetaNorm += tmp_float * tmp_float;

#endif

			}
		}
		// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
		else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
		{
			val_arr(m_p) = 0.0f;

			// Active set
			// Here we calculate chopped gradient Beta, hence BetaNorm updated
			// beta = g+ = max{g, 0.0f};

			tmp_float = inc_arr(m_g);

			if (tmp_float > 0.0f)
			{
#if (PROPORTIONING_DOUBLE_STEP == 1)

				if (tmp_float2 / alpha_fucked < tmp_float)
				{
					BetaNorm += tmp_float * tmp_float2 / alpha_fucked;
				}
				else
				{
					BetaNorm += tmp_float * tmp_float;
				}

#else

						BetaNorm += tmp_float * tmp_float;

#endif
			}
		}
		else
		{
			// Free set
			// Here we calculate free gradient phi and reduced free gradient phi~
			// phi~(alpha) = min{(x_i - lo_i), alpha}

			tmp_float = inc_arr(m_g);

// 			if (!M_EQUAL(ptr_m_Lo[i], -FLT_MAX) && tmp_float2 / alpha_fucked < tmp_float)
// 			{
// 				Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
// 			}
// 			else
// 			{
// 				Phi_Phi += tmp_float * tmp_float;
// 			}

			// We're going to Lo
			if (tmp_float > 0.0f)
			{
				if (tmp_float2 / alpha_fucked < tmp_float)
				{
	 				Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
				} else
				{
	 				Phi_Phi += tmp_float * tmp_float;
				}
			}
			// We're going to Hi
			else
			{
				if (tmp_float3 / alpha_fucked > tmp_float)
				{
	 				Phi_Phi += tmp_float * tmp_float3 / alpha_fucked;
				} else
				{
	 				Phi_Phi += tmp_float * tmp_float;
				}
			}
		}

		++ptr_m_p;
	}

	// MAIN ITERATION CYCLE
	m_EffectiveIterations = m_Iterations;

	for (unsigned k = 0; k < m_Iterations; ++k)
	{
		// grad_norm = ||beta + phi||
		float grad_norm = 0.0f;

		if (BetaNorm <= GAMMA_Sq * Phi_Phi)
		{
			// CG prestep

			float dot_gp = 0.0f, dot_pAp = 0.0f;

			ptr_m_p = &m_p[0];
			ptr_m_g = &m_g[0];

			ptr_m_a_x = &m_a_x[0];
			ptr_m_a_y = &m_a_y[0];
			ptr_m_a_z = &m_a_z[0];

			INIT_A_SPARSE(0);
			INIT_NODES_IDX(0);
			
			m_a_x.Empty();
			m_a_y.Empty();
			m_a_z.Empty();

			ptr_m_Lo = &m_Lo[0];
			ptr_m_Hi = &m_Hi[0];
			ptr_m_lambda = &m_lambda[0];

			// Calculate:
			//		temp_vec = W * J^T * p;
			//		g^T * p;
			for (i = 0; i < m_NumJoints; ++i)
			{
#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

				CALC_VECTOR_A_PRECOND(m_p, 0, 00, 01, 02);
				CALC_VECTOR_A_PRECOND(m_p, 1, 03, 04, 05);
				CALC_VECTOR_A_PRECOND(m_p, 2, 06, 07, 08);
				CALC_VECTOR_A_PRECOND(m_p, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

				CALC_VECTOR_A_PRECOND(m_p, 4, 12, 13, 14);
				CALC_VECTOR_A_PRECOND(m_p, 5, 15, 16, 17);

#endif

#else

				CALC_VECTOR_A(m_p, 0, 00, 01, 02);
				CALC_VECTOR_A(m_p, 1, 03, 04, 05);
				CALC_VECTOR_A(m_p, 2, 06, 07, 08);
				CALC_VECTOR_A(m_p, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

				CALC_VECTOR_A(m_p, 4, 12, 13, 14);
				CALC_VECTOR_A(m_p, 5, 15, 16, 17);

#endif

#endif

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
				// Preconditioner in face:
				// z dot g = D^-1 * phi dot g
				if ( (math::fastabs(ptr_m_lambda[i] - ptr_m_Lo[i]) < ACT_SET_THRESHOLD) ||
					 (math::fastabs(ptr_m_lambda[i] - ptr_m_Hi[i]) < ACT_SET_THRESHOLD) )
				{
					// Active set
				}
				else
				{
					// Free set
					dot_gp += val_arr(m_g) * (ptr_m_invDiag[i] * val_arr(m_g));
				}
#else
				dot_gp += val_arr(m_g) * val_arr(m_p);
#endif

				++ptr_m_g;
				++ptr_m_p;
			}

			ptr_m_p = &m_p[0];
			
			ptr_m_a_x = &m_a_x[0];
			ptr_m_a_y = &m_a_y[0];
			ptr_m_a_z = &m_a_z[0];

			ptr_m_Lo = &m_Lo[0];
			ptr_m_Hi = &m_Hi[0];
			ptr_m_RHS = &m_RHS[0];
			ptr_m_ap_vec = &m_ap_vec[0];
			ptr_m_lambda = &m_lambda[0];

			ptr_m_CFM = &m_CFM[0];

			INIT_J_SPARSE(0);
			INIT_NODES_IDX(0);

			float alpha_f = FLT_MAX;

			// Calculate:
			//		p^T * (J * temp_vec);
			for (i = 0; i < m_NumJoints; ++i)
			{
				// Calculating residual, direction, d0 = r0 = dx = b - Ax
				unsigned int t0Idx = inc_arr(m_JtdNodes_00);
				unsigned int t1Idx = inc_arr(m_JtdNodes_01);
				unsigned int t2Idx = inc_arr(m_JtdNodes_02);
				unsigned int t3Idx = inc_arr(m_JtdNodes_03);

#if (NUM_JOINT_TRIPLES == 6)

				unsigned int t4Idx = inc_arr(m_JtdNodes_04);
				unsigned int t5Idx = inc_arr(m_JtdNodes_05);

#endif

				// J * a
				val_arr(m_ap_vec) =	MUL_J_A();

				// CFM fuck up
				val_arr(m_ap_vec) += ptr_m_CFM[i] * ptr_m_p[i];

				float precond_track = val_arr(m_ap_vec);

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
				// Preconditioning
				val_arr(m_ap_vec) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
				val_arr(m_ap_vec) *= ptr_m_invDiag[i];
#endif
// 				assert(val_arr(m_ap_vec) == val_arr(m_ap_vec));
// 				assert(val_arr(m_ap_vec) < FLT_MAX);
// 				assert(val_arr(m_ap_vec) > -FLT_MAX);

				dot_pAp += ptr_m_p[i] * inc_arr(m_ap_vec);

				//assert((dot_pAp) > -0.01f);

// 				assert((dot_pAp) == (dot_pAp));
// 				assert((dot_pAp) < FLT_MAX);
// 				assert((dot_pAp) > -FLT_MAX);

				// Calculate alpha_f
				// UPD: LoHi
				if (ptr_m_p[i] > 0.0f)
				{
					tmp_float = (ptr_m_lambda[i] - ptr_m_Lo[i]) / ptr_m_p[i];

					if (tmp_float < alpha_f)
					{
						alpha_f = tmp_float;
					}
				}
				else if (ptr_m_p[i] < 0.0f)
				{
					tmp_float = (ptr_m_lambda[i] - ptr_m_Hi[i]) / ptr_m_p[i];

					if (tmp_float < alpha_f)
					{
						alpha_f = tmp_float;
					}
				}
			}

			alpha_cg = dot_gp / dot_pAp;
			if (alpha_cg < alpha_f)
			{
				//////////////////////////////////////////////////////////////////////////
				// Conjugate Gradient Step
				//////////////////////////////////////////////////////////////////////////
// 				gLog.Print("CG   Step [%d], B[%f]; P[%f]", k, BetaNorm, Phi_Phi);

				ptr_m_p = &m_p[0];
				ptr_m_g = &m_g[0];
				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];
				ptr_m_lambda = &m_lambda[0];
				ptr_m_ap_vec = &m_ap_vec[0];

				BetaNorm = 0.0f;
				Phi_Phi = 0.0f;

				float dot_phiAp = 0.0f;

				// Calculate:
				//		x_k+1 = x_k - alpha_cg * p;
				//		g = g - alpha_cg * A * p;
				//		recalculate BetaNorm, Phi_phi
				for (i = 0; i < m_NumJoints; ++i)
				{
					val_arr(m_lambda) -= alpha_cg * val_arr(m_p);
					val_arr(m_g) -= alpha_cg * val_arr(m_ap_vec);

					// UPD: LoHi
					tmp_float2 = val_arr(m_lambda) - ptr_m_Lo[i];
					tmp_float3 = val_arr(m_lambda) - ptr_m_Hi[i];

					++ptr_m_p;

					// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
					if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
					{
						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g- = min{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float < 0.0f)
						{
#if (PROPORTIONING_DOUBLE_STEP == 1)

							if (tmp_float3 / alpha_fucked > tmp_float)
							{
								BetaNorm += tmp_float * tmp_float3 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}
#else

							BetaNorm += tmp_float * tmp_float;

#endif

							grad_norm += tmp_float * tmp_float;
						}
					}
					// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
					else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
					{
						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g+ = max{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float > 0.0f)
						{
#if (PROPORTIONING_DOUBLE_STEP == 1)

							if (tmp_float2 / alpha_fucked < tmp_float)
							{
								BetaNorm += tmp_float * tmp_float2 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}

#else

							BetaNorm += tmp_float * tmp_float;

#endif

							grad_norm += tmp_float * tmp_float;
						}
					}
					else
					{
						// Free set
						// Here we calculate free gradient phi and reduced free gradient phi~
						// phi~(alpha) = min{(x_i - lo_i), alpha}

						tmp_float = inc_arr(m_g);

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
						// Preconditioner in face:
						// (phi dot A * p) replaced with (z dot A * p)
						dot_phiAp += (ptr_m_invDiag[i] * tmp_float) * val_arr(m_ap_vec);
#else
						dot_phiAp += tmp_float * val_arr(m_ap_vec);
#endif

						grad_norm += tmp_float * tmp_float;

						// We're going to Lo
						if (tmp_float > 0.0f)
						{
							if (tmp_float2 / alpha_fucked < tmp_float)
							{
	 							Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
							} else
							{
	 							Phi_Phi += tmp_float * tmp_float;
							}
						}
						// We're going to Hi
						else
						{
							if (tmp_float3 / alpha_fucked > tmp_float)
							{
	 							Phi_Phi += tmp_float * tmp_float3 / alpha_fucked;
							} else
							{
	 							Phi_Phi += tmp_float * tmp_float;
							}
						}
					}

// 					assert(Phi_Phi >= 0.0f);
// 
// 					assert((dot_phiAp) == (dot_phiAp));
// 					assert((dot_phiAp) < FLT_MAX);
// 					assert((dot_phiAp) > -FLT_MAX);
// 
// 					assert((dot_pAp) == (dot_pAp));
// 					assert((dot_pAp) < FLT_MAX);
// 					assert((dot_pAp) > -FLT_MAX);

					++ptr_m_ap_vec;
					++ptr_m_lambda;
				}

// 				assert(Phi_Phi >= 0.0f);

				beta = dot_phiAp / dot_pAp;

				ptr_m_p = &m_p[0];
				ptr_m_g = &m_g[0];
				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];
				ptr_m_lambda = &m_lambda[0];

				// Calculate p = phi(x_k+1) - beta * p;
				for (i = 0; i < m_NumJoints; ++i)
				{
					// Calculate phi from gradient
					tmp_float = inc_arr(m_g);

					if ( (math::fastabs(val_arr(m_lambda) - ptr_m_Lo[i]) < ACT_SET_THRESHOLD) ||
						 (math::fastabs(val_arr(m_lambda) - ptr_m_Hi[i]) < ACT_SET_THRESHOLD) )
					{
						tmp_float = 0.0f;
					}

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
					// Preconditioner in face:
					// p = z - beta * p;
					val_arr(m_p) = (ptr_m_invDiag[i] * tmp_float) - beta * val_arr(m_p);
#else
					val_arr(m_p) = tmp_float - beta * val_arr(m_p);
#endif

// 					assert(val_arr(m_p) == val_arr(m_p));
// 					assert(val_arr(m_p) < FLT_MAX);
// 					assert(val_arr(m_p) > -FLT_MAX);

					++ptr_m_p;
					++ptr_m_lambda;
				}

// 				gLog.Print("CG   Step Gradient[%f]: Phi[%f], Beta[%f], alpha: %12.10f;", grad_norm, Phi_Phi, BetaNorm, alpha_cg);
			}
			else
			{
				//////////////////////////////////////////////////////////////////////////
				// Expansion Step
				//////////////////////////////////////////////////////////////////////////
// 				gLog.Print("Exp  Step [%d], B[%f]; P[%f]", k, BetaNorm, Phi_Phi);

				ptr_m_p = &m_p[0];
				ptr_m_g = &m_g[0];
				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];

				ptr_m_lambda = &m_lambda[0];
				ptr_m_ap_vec = &m_ap_vec[0];
				
				INIT_A_SPARSE(0);
				INIT_NODES_IDX(0);

				ptr_m_a_x = &m_a_x[0];
				ptr_m_a_y = &m_a_y[0];
				ptr_m_a_z = &m_a_z[0];

				m_a_x.Empty();
				m_a_y.Empty();
				m_a_z.Empty();

				float lambda_sq = 0.0f;

				// Calculate:
				//		x_k+1/2 = x_k - alpha_f * p;
				//		g = g - alpha_f * A * p;
				//		x_k+1 = project(x_k+1/2 - alpha_f * phi(x_k+1/2))
				//		a = W * J^T * x_k+1
				for (i = 0; i < m_NumJoints; ++i)
				{
					val_arr(m_g) -= alpha_f * inc_arr(m_ap_vec);

					// calculate x_k+1/2 {needed for phi} [TODO: optimize]
					val_arr(m_lambda) -= alpha_f * inc_arr(m_p);

					// Calculate phi(x_k+1/2)
					if ( (math::fastabs(val_arr(m_lambda) - ptr_m_Lo[i]) > ACT_SET_THRESHOLD) &&
						 (math::fastabs(val_arr(m_lambda) - ptr_m_Hi[i]) > ACT_SET_THRESHOLD) )
					{
						// Free set
						val_arr(m_lambda) -= alpha_fucked * val_arr(m_g);
					}

					// Project
					if (val_arr(m_lambda) < ptr_m_Lo[i])
					{
						val_arr(m_lambda) = ptr_m_Lo[i];
					}
					else if (val_arr(m_lambda) > ptr_m_Hi[i])
					{
						val_arr(m_lambda) = ptr_m_Hi[i];
					}

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

					CALC_VECTOR_A_PRECOND(m_lambda, 0, 00, 01, 02);
					CALC_VECTOR_A_PRECOND(m_lambda, 1, 03, 04, 05);
					CALC_VECTOR_A_PRECOND(m_lambda, 2, 06, 07, 08);
					CALC_VECTOR_A_PRECOND(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

					CALC_VECTOR_A_PRECOND(m_lambda, 4, 12, 13, 14);
					CALC_VECTOR_A_PRECOND(m_lambda, 5, 15, 16, 17);

#endif

#else

					CALC_VECTOR_A(m_lambda, 0, 00, 01, 02);
					CALC_VECTOR_A(m_lambda, 1, 03, 04, 05);
					CALC_VECTOR_A(m_lambda, 2, 06, 07, 08);
					CALC_VECTOR_A(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

					CALC_VECTOR_A(m_lambda, 4, 12, 13, 14);
					CALC_VECTOR_A(m_lambda, 5, 15, 16, 17);

#endif

#endif

					lambda_sq += val_arr(m_lambda) * val_arr(m_lambda);

					++ptr_m_lambda;
					++ptr_m_g;
				}

				ptr_m_p = &m_p[0];
				ptr_m_g = &m_g[0];

				ptr_m_a_x = &m_a_x[0];
				ptr_m_a_y = &m_a_y[0];
				ptr_m_a_z = &m_a_z[0];

				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];
				ptr_m_RHS = &m_RHS[0];
				ptr_m_lambda = &m_lambda[0];
				
				INIT_J_SPARSE(0);
				INIT_NODES_IDX(0);

				BetaNorm = 0.0f;
				Phi_Phi = 0.0f;

				// Calculate:
				//		g = A * x_k+1 - b;
				//		p = phi(x_k+1)
				//		recalculate BetaNorm, Phi_phi
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
					val_arr(m_g) = MUL_J_A();

					// CFM fuck up
					val_arr(m_g) += m_CFM[i] * ptr_m_lambda[i];

					val_arr(m_g) -= inc_arr(m_RHS);

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
					// Preconditioning
					val_arr(m_g) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
					val_arr(m_g) *= ptr_m_invDiag[i];
#endif

// 					assert(val_arr(m_g) == val_arr(m_g));
// 					assert(val_arr(m_g) < FLT_MAX);
// 					assert(val_arr(m_g) > -FLT_MAX);

					// Calculate phi from gradient, and p = phi(x_k+1)
					// zero later if x[i] is not in free set

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
					// Preconditioner in face:
					// p = z
					val_arr(m_p) = ptr_m_invDiag[i] * val_arr(m_g);
#else
					val_arr(m_p) = val_arr(m_g);
#endif

					// UPD: LoHi
					tmp_float2 = ptr_m_lambda[i] - ptr_m_Lo[i];
					tmp_float3 = ptr_m_lambda[i] - ptr_m_Hi[i];

					// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
					if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
					{
						val_arr(m_p) = 0.0f;

						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g- = min{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float < 0.0f)
						{
#if (PROPORTIONING_DOUBLE_STEP == 1)

							if (tmp_float3 / alpha_fucked > tmp_float)
							{
								BetaNorm += tmp_float * tmp_float3 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}

#else

						BetaNorm += tmp_float * tmp_float;

#endif

							grad_norm += tmp_float * tmp_float;
						}
					}
					// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
					else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
					{
						val_arr(m_p) = 0.0f;

						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g+ = max{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float > 0.0f)
						{

#if (PROPORTIONING_DOUBLE_STEP == 1)

							if (tmp_float2 / alpha_fucked < tmp_float)
							{
								BetaNorm += tmp_float * tmp_float2 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}

#else

						BetaNorm += tmp_float * tmp_float;

#endif
							grad_norm += tmp_float * tmp_float;
						}
					}
					else
					{
						// Free set
						// Here we calculate free gradient phi and reduced free gradient phi~
						// phi~(alpha) = min{(x_i - lo_i), alpha}

						tmp_float = inc_arr(m_g);

						grad_norm += tmp_float * tmp_float;

						// We're going to Lo
						if (tmp_float > 0.0f)
						{
							if (tmp_float2 / alpha_fucked < tmp_float)
							{
	 							Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
							} else
							{
	 							Phi_Phi += tmp_float * tmp_float;
							}
						}
						// We're going to Hi
						else
						{
							if (tmp_float3 / alpha_fucked > tmp_float)
							{
	 							Phi_Phi += tmp_float * tmp_float3 / alpha_fucked;
							} else
							{
	 							Phi_Phi += tmp_float * tmp_float;
							}
						}
					}

// 					assert(val_arr(m_p) == val_arr(m_p));
// 					assert(val_arr(m_p) < FLT_MAX);
// 					assert(val_arr(m_p) > -FLT_MAX);

					++ptr_m_p;
				}

// 				assert(Phi_Phi >= 0.0f);
//  				gLog.Print("Exp  Step [%f]: Phi[%f], Beta[%f]", grad_norm, Phi_Phi, BetaNorm);
			}
		}
		else
		{
			//////////////////////////////////////////////////////////////////////////
			// Proportioning step
			//////////////////////////////////////////////////////////////////////////
// 			gLog.Print("Prop Step [%d], B[%f]; P[%f]", k, BetaNorm, Phi_Phi);

#if (PROPORTIONING_DOUBLE_STEP == 1)

			float dot_gd = 0.0f, dot_dAd = 0.0f;

			ptr_m_p = &m_p[0];
			ptr_m_g = &m_g[0];
			ptr_m_d = &m_d[0];

			ptr_m_a_x = &m_a_x[0];
			ptr_m_a_y = &m_a_y[0];
			ptr_m_a_z = &m_a_z[0];

			ptr_m_Lo = &m_Lo[0];
			ptr_m_Hi = &m_Hi[0];

			ptr_m_lambda = &m_lambda[0];

			INIT_A_SPARSE(0);
			INIT_NODES_IDX(0);

			m_a_x.Empty();
			m_a_y.Empty();
			m_a_z.Empty();


			float alpha_f = FLT_MAX;

			// Calculate:
			//		d = Beta(x_k);
			//		temp_vec = W * J^T * d;
			//		g^T * d;

			unsigned int temp_cnt = 0;

			for (i = 0; i < m_NumJoints; ++i)
			{
				val_arr(m_d) = 0.0f;

				// UPD: LoHi
				tmp_float2 = ptr_m_lambda[i] - ptr_m_Lo[i];
				tmp_float3 = ptr_m_lambda[i] - ptr_m_Hi[i];


				// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
				if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
				{
					// Active set [Lo]
					// Here we calculate chopped gradient Beta, hence BetaNorm updated
					// beta = g- = min{g, 0.0f};

					tmp_float = val_arr(m_g);
					if (tmp_float > 0.0f)
					{
						// Beta = 0.0
						tmp_float = 0.0f;

						// Advance arrays with no action
						ADVANCE_A_SPARSE_NODES_IDX();
					}
					else
					{
						float tmp_alpha_f = tmp_float3 / tmp_float;
						
						if (tmp_alpha_f < alpha_f)
						{
							alpha_f = tmp_alpha_f;
						}

						// Beta = g_i
						val_arr(m_d) = tmp_float;

						++temp_cnt;

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

						CALC_VECTOR_A_PRECOND(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A_PRECOND(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A_PRECOND(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A_PRECOND(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A_PRECOND(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A_PRECOND(m_d, 5, 15, 16, 17);

#endif

#else

						CALC_VECTOR_A(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A(m_d, 5, 15, 16, 17);

#endif

#endif

						dot_gd += val_arr(m_g) * val_arr(m_d);
					}
				}
				// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
				else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
				{
					// Active set [Hi]
					// Here we calculate chopped gradient Beta, hence BetaNorm updated
					// beta = g+ = max{g, 0.0f};

					tmp_float = val_arr(m_g);
					if (tmp_float < 0.0f)
					{
						// Beta = 0.0
						tmp_float = 0.0f;

						// Advance arrays with no action
						ADVANCE_A_SPARSE_NODES_IDX();
					}
					else
					{
						float tmp_alpha_f = tmp_float2 / tmp_float;

						if (tmp_alpha_f < alpha_f)
						{
							alpha_f = tmp_alpha_f;
						}

						// Beta = g_i
						val_arr(m_d) = tmp_float;
						++temp_cnt;

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

						CALC_VECTOR_A_PRECOND(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A_PRECOND(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A_PRECOND(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A_PRECOND(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A_PRECOND(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A_PRECOND(m_d, 5, 15, 16, 17);

#endif

#else

						CALC_VECTOR_A(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A(m_d, 5, 15, 16, 17);

#endif

#endif

						dot_gd += val_arr(m_g) * val_arr(m_d);
					}
				}
				else
				{
					// Beta = 0.0

					// Advance arrays with no action
					ADVANCE_A_SPARSE_NODES_IDX();
				}

				++ptr_m_g;
				++ptr_m_d;
			}

// 			assert(temp_cnt > 0);

			ptr_m_d = &m_d[0];

			ptr_m_a_x = &m_a_x[0];
			ptr_m_a_y = &m_a_y[0];
			ptr_m_a_z = &m_a_z[0];

			ptr_m_RHS = &m_RHS[0];
			ptr_m_ad_vec = &m_ad_vec[0];
			ptr_m_lambda = &m_lambda[0];

			INIT_J_SPARSE(0);
			INIT_NODES_IDX(0);


			// Calculate:
			//		A * d
			//		d^T * A * d
			for (i = 0; i < m_NumJoints; ++i)
			{
				// Calculating residual, direction, d0 = r0 = dx = b - Ax
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
				val_arr(m_ad_vec) += m_CFM[i] * ptr_m_d[i];

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
				// Preconditioning
				val_arr(m_ad_vec) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
				val_arr(m_ad_vec) *= ptr_m_invDiag[i];
#endif

// 				assert(val_arr(m_ad_vec) == val_arr(m_ad_vec));
// 				assert(val_arr(m_ad_vec) < FLT_MAX);
// 				assert(val_arr(m_ad_vec) > -FLT_MAX);

				dot_dAd += ptr_m_d[i] * inc_arr(m_ad_vec);
			}

			alpha_cg = dot_gd / dot_dAd;

// 			assert(dot_gd > 1e-4f);
// 			assert(dot_dAd > 1e-4f);

			if (alpha_f < alpha_cg)
			{
				alpha_cg = alpha_f;


				ptr_m_a_x = &m_a_x[0];
				ptr_m_a_y = &m_a_y[0];
				ptr_m_a_z = &m_a_z[0];

				m_a_x.Empty();
				m_a_y.Empty();
				m_a_z.Empty();

				ptr_m_g = &m_g[0];
				ptr_m_p = &m_p[0];
				ptr_m_d = &m_d[0];
				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];
				ptr_m_ad_vec = &m_ad_vec[0];
				ptr_m_lambda = &m_lambda[0];


				INIT_A_SPARSE(0);
				INIT_NODES_IDX(0);

				// Calculate:
				//		x_k+1/2 = x_k - alpha_cg * d;
				//		g = g - alpha_cg * A * d;
				//		x_k+1 = Proj(x_k+1/2 - alpha * Beta(x_k+1/2);
				for (i = 0; i < m_NumJoints; ++i)
				{
					val_arr(m_lambda) -= alpha_cg * inc_arr(m_d);
					val_arr(m_g) -= alpha_cg * inc_arr(m_ad_vec);


					tmp_float2 = val_arr(m_lambda) - ptr_m_Lo[i];
					tmp_float3 = val_arr(m_lambda) - ptr_m_Hi[i];

					// Calculate x_k+1
					if ( math::fastabs(tmp_float2) < ACT_SET_THRESHOLD )
					{
						// Active set [Lo]
						// Beta = g- = min(g, 0.0);
						if (val_arr(m_g) < 0.0f)
							val_arr(m_lambda) -= alpha_fucked * val_arr(m_g);
					}
					else if ( math::fastabs(tmp_float3) < ACT_SET_THRESHOLD )
					{
						// Active set [Hi]
						// Beta = g+ = max(g, 0.0);
						if (val_arr(m_g) > 0.0f)
							val_arr(m_lambda) -= alpha_fucked * val_arr(m_g);
					}

					// Project
					if (val_arr(m_lambda) < ptr_m_Lo[i])
					{
						val_arr(m_lambda) = ptr_m_Lo[i];
					}
					else if (val_arr(m_lambda) > ptr_m_Hi[i])
					{
						val_arr(m_lambda) = ptr_m_Hi[i];
					}

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

					CALC_VECTOR_A_PRECOND(m_lambda, 0, 00, 01, 02);
					CALC_VECTOR_A_PRECOND(m_lambda, 1, 03, 04, 05);
					CALC_VECTOR_A_PRECOND(m_lambda, 2, 06, 07, 08);
					CALC_VECTOR_A_PRECOND(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

					CALC_VECTOR_A_PRECOND(m_lambda, 4, 12, 13, 14);
					CALC_VECTOR_A_PRECOND(m_lambda, 5, 15, 16, 17);

#endif

#else

					CALC_VECTOR_A(m_lambda, 0, 00, 01, 02);
					CALC_VECTOR_A(m_lambda, 1, 03, 04, 05);
					CALC_VECTOR_A(m_lambda, 2, 06, 07, 08);
					CALC_VECTOR_A(m_lambda, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

					CALC_VECTOR_A(m_lambda, 4, 12, 13, 14);
					CALC_VECTOR_A(m_lambda, 5, 15, 16, 17);

#endif

#endif

					++ptr_m_lambda;
					++ptr_m_g;
				}


				ptr_m_a_x = &m_a_x[0];
				ptr_m_a_y = &m_a_y[0];
				ptr_m_a_z = &m_a_z[0];

				ptr_m_g = &m_g[0];
				ptr_m_p = &m_p[0];
				ptr_m_d = &m_d[0];
				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];
				ptr_m_RHS = &m_RHS[0];
				ptr_m_lambda = &m_lambda[0];

				INIT_J_SPARSE(0);
				INIT_NODES_IDX(0);

				BetaNorm = 0.0f;
				Phi_Phi = 0.0f;

				// Calculate:
				//		g = A * x_k+1 - b;
				//		p = phi(x_k+1);
				//		recalculate BetaNorm, Phi_phi
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
					val_arr(m_g) = MUL_J_A();

					// CFM fuck up
					val_arr(m_g) += m_CFM[i] * ptr_m_lambda[i];

					val_arr(m_g) -= inc_arr(m_RHS);

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
					// Preconditioning
					val_arr(m_g) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
					val_arr(m_g) *= ptr_m_invDiag[i];
#endif

// 					assert(val_arr(m_g) == val_arr(m_g));
// 					assert(val_arr(m_g) < FLT_MAX);
// 					assert(val_arr(m_g) > -FLT_MAX);

					// Calculate phi from gradient, and p = phi(x_k+1)
					// zero later if x[i] is not in free set

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
					// Preconditioner in face:
					// p = z
					val_arr(m_p) = ptr_m_invDiag[i] * val_arr(m_g);
#else
					val_arr(m_p) = val_arr(m_g);
#endif

					// UPD: LoHi
					tmp_float2 = ptr_m_lambda[i] - ptr_m_Lo[i];
					tmp_float3 = ptr_m_lambda[i] - ptr_m_Hi[i];

					// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
					if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
					{
						val_arr(m_p) = 0.0f;

						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g- = min{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float < 0.0f)
						{
							if (tmp_float3 / alpha_fucked > tmp_float)
							{
								BetaNorm += tmp_float * tmp_float3 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}
							grad_norm += tmp_float * tmp_float;
						}
					}
					// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
					else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
					{
						val_arr(m_p) = 0.0f;

						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g+ = max{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float > 0.0f)
						{
							if (tmp_float2 / alpha_fucked < tmp_float)
							{
								BetaNorm += tmp_float * tmp_float2 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}
							grad_norm += tmp_float * tmp_float;
						}
					}
					else
					{
						// Free set
						// Here we calculate free gradient phi and reduced free gradient phi~
						// phi~(alpha) = min{(x_i - lo_i) / alpha, phi}

						tmp_float = inc_arr(m_g);

						grad_norm += tmp_float*tmp_float;

						// We're going to Lo
						if (tmp_float > 0.0f)
						{
							if (tmp_float2 / alpha_fucked < tmp_float)
							{
 								Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
							} else
							{
 								Phi_Phi += tmp_float * tmp_float;
							}
						}
						// We're going to Hi
						else
						{
							if (tmp_float3 / alpha_fucked > tmp_float)
							{
 								Phi_Phi += tmp_float * tmp_float3 / alpha_fucked;
							} else
							{
 								Phi_Phi += tmp_float * tmp_float;
							}
						}
					}

// 					assert(val_arr(m_p) == val_arr(m_p));
// 					assert(val_arr(m_p) < FLT_MAX);
// 					assert(val_arr(m_p) > -FLT_MAX);

					++ptr_m_p;
				}

// 				assert(Phi_Phi >= 0.0f);
			}
			else
			{
				ptr_m_g = &m_g[0];
				ptr_m_p = &m_p[0];
				ptr_m_d = &m_d[0];
				ptr_m_Lo = &m_Lo[0];
				ptr_m_Hi = &m_Hi[0];
				ptr_m_ad_vec = &m_ad_vec[0];
				ptr_m_lambda = &m_lambda[0];

				BetaNorm = 0.0f;
				Phi_Phi = 0.0f;

				// Calculate:
				//		x_k+1 = x_k - alpha_cg * d;
				//		g = g - alpha_cg * A * d;
				//		p = phi(x_k+1);
				//		recalculate BetaNorm, Phi_phi
				for (i = 0; i < m_NumJoints; ++i)
				{
					val_arr(m_lambda) -= alpha_cg * inc_arr(m_d);
					val_arr(m_g) -= alpha_cg * inc_arr(m_ad_vec);

					// Calculate phi from gradient, and p = phi(x_k+1)

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
					// Preconditioner in face:
					// p = z
					val_arr(m_p) = ptr_m_invDiag[i] * val_arr(m_g);
#else
					val_arr(m_p) = val_arr(m_g);
#endif

					// UPD: LoHi
					tmp_float2 = val_arr(m_lambda) - ptr_m_Lo[i];
					tmp_float3 = val_arr(m_lambda) - ptr_m_Hi[i];

					++ptr_m_lambda;

					// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
					if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
					{
						val_arr(m_p) = 0.0f;

						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g- = min{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float < 0.0f)
						{
							if (tmp_float3 / alpha_fucked > tmp_float)
							{
								BetaNorm += tmp_float * tmp_float3 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}
							grad_norm += tmp_float * tmp_float;
						}
					}
					// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
					else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
					{
						val_arr(m_p) = 0.0f;

						// Active set
						// Here we calculate chopped gradient Beta, hence BetaNorm updated
						// beta = g+ = max{g, 0.0f};

						tmp_float = inc_arr(m_g);
						if (tmp_float > 0.0f)
						{
							if (tmp_float2 / alpha_fucked < tmp_float)
							{
								BetaNorm += tmp_float * tmp_float2 / alpha_fucked;
							}
							else
							{
								BetaNorm += tmp_float * tmp_float;
							}
							grad_norm += tmp_float * tmp_float;
						}
					}
					else
					{
						// Free set
						// Here we calculate free gradient phi and reduced free gradient phi~
						// phi~(alpha) = min{(x_i - lo_i) / alpha, phi}

						tmp_float = inc_arr(m_g);

						grad_norm += tmp_float*tmp_float;

						// We're going to Lo
						if (tmp_float > 0.0f)
						{
							if (tmp_float2 / alpha_fucked < tmp_float)
							{
 								Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
							} else
							{
 								Phi_Phi += tmp_float * tmp_float;
							}
						}
						// We're going to Hi
						else
						{
							if (tmp_float3 / alpha_fucked > tmp_float)
							{
 								Phi_Phi += tmp_float * tmp_float3 / alpha_fucked;
							} else
							{
 								Phi_Phi += tmp_float * tmp_float;
							}
						}
					}

// 					assert(val_arr(m_p) == val_arr(m_p));
// 					assert(val_arr(m_p) < FLT_MAX);
// 					assert(val_arr(m_p) > -FLT_MAX);

					++ptr_m_p;
				}

// 				assert(Phi_Phi >= 0.0f);
			}

#else


			float dot_gd = 0.0f, dot_dAd = 0.0f;

			ptr_m_p = &m_p[0];
			ptr_m_g = &m_g[0];
			ptr_m_d = &m_d[0];

			ptr_m_a_x = &m_a_x[0];
			ptr_m_a_y = &m_a_y[0];
			ptr_m_a_z = &m_a_z[0];

			ptr_m_Lo = &m_Lo[0];
			ptr_m_Hi = &m_Hi[0];

			ptr_m_lambda = &m_lambda[0];

			INIT_A_SPARSE(0);
			INIT_NODES_IDX(0);

			m_a_x.Empty();
			m_a_y.Empty();
			m_a_z.Empty();


			float alpha_f = FLT_MAX;

			// Calculate:
			//		d = Beta(x_k);
			//		temp_vec = W * J^T * d;
			//		g^T * d;
			for (i = 0; i < m_NumJoints; ++i)
			{
				val_arr(m_d) = 0.0f;

				// UPD: LoHi
				tmp_float2 = ptr_m_lambda[i] - ptr_m_Lo[i];
				tmp_float3 = ptr_m_lambda[i] - ptr_m_Hi[i];


				// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
				if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
				{
					// Active set [Lo]
					// Here we calculate chopped gradient Beta, hence BetaNorm updated
					// beta = g- = min{g, 0.0f};

					tmp_float = val_arr(m_g);
					if (tmp_float > 0.0f)
					{
						// Beta = 0.0
						tmp_float = 0.0f;

						// Advance arrays with no action
						ADVANCE_A_SPARSE_NODES_IDX();
					}
					else
					{
						float tmp_alpha_f = tmp_float3 / tmp_float;
						
						if (tmp_alpha_f < alpha_f)
						{
							alpha_f = tmp_alpha_f;
						}

						// Beta = g_i
						val_arr(m_d) = tmp_float;

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

						CALC_VECTOR_A_PRECOND(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A_PRECOND(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A_PRECOND(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A_PRECOND(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A_PRECOND(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A_PRECOND(m_d, 5, 15, 16, 17);

#endif

#else

						CALC_VECTOR_A(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A(m_d, 5, 15, 16, 17);

#endif

#endif

						dot_gd += val_arr(m_g) * val_arr(m_d);
					}
				}
				// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
				else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
				{
					// Active set [Hi]
					// Here we calculate chopped gradient Beta, hence BetaNorm updated
					// beta = g+ = max{g, 0.0f};

					tmp_float = val_arr(m_g);
					if (tmp_float < 0.0f)
					{
						// Beta = 0.0
						tmp_float = 0.0f;

						// Advance arrays with no action
						ADVANCE_A_SPARSE_NODES_IDX();
					}
					else
					{
						float tmp_alpha_f = tmp_float2 / tmp_float;

						if (tmp_alpha_f < alpha_f)
						{
							alpha_f = tmp_alpha_f;
						}

						// Beta = g_i
						val_arr(m_d) = tmp_float;

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)

						CALC_VECTOR_A_PRECOND(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A_PRECOND(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A_PRECOND(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A_PRECOND(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A_PRECOND(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A_PRECOND(m_d, 5, 15, 16, 17);

#endif

#else

						CALC_VECTOR_A(m_d, 0, 00, 01, 02);
						CALC_VECTOR_A(m_d, 1, 03, 04, 05);
						CALC_VECTOR_A(m_d, 2, 06, 07, 08);
						CALC_VECTOR_A(m_d, 3, 09, 10, 11);

#if (NUM_JOINT_TRIPLES == 6)

						CALC_VECTOR_A(m_d, 4, 12, 13, 14);
						CALC_VECTOR_A(m_d, 5, 15, 16, 17);

#endif

#endif

						dot_gd += val_arr(m_g) * val_arr(m_d);
					}
				}
				else
				{
					// Beta = 0.0

					// Advance arrays with no action
					ADVANCE_A_SPARSE_NODES_IDX();
				}

				++ptr_m_g;
				++ptr_m_d;
			}

			ptr_m_d = &m_d[0];

			ptr_m_a_x = &m_a_x[0];
			ptr_m_a_y = &m_a_y[0];
			ptr_m_a_z = &m_a_z[0];

			ptr_m_RHS = &m_RHS[0];
			ptr_m_ad_vec = &m_ad_vec[0];
			ptr_m_lambda = &m_lambda[0];

			INIT_J_SPARSE(0);
			INIT_NODES_IDX(0);

			// Calculate:
			//		A * d
			//		d^T * A * d
			for (i = 0; i < m_NumJoints; ++i)
			{
				// Calculating residual, direction, d0 = r0 = dx = b - Ax
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
				val_arr(m_ad_vec) += m_CFM[i] * ptr_m_d[i];

#if (SW_LCPCG_PRECONDITIONER_LEFTJACOBI == 1)
				// Preconditioning
				val_arr(m_ad_vec) *= ptr_m_invDiag[i];
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
				val_arr(m_ad_vec) *= ptr_m_invDiag[i];
#endif

// 				assert(val_arr(m_ad_vec) == val_arr(m_ad_vec));
// 				assert(val_arr(m_ad_vec) < FLT_MAX);
// 				assert(val_arr(m_ad_vec) > -FLT_MAX);

				dot_dAd += ptr_m_d[i] * inc_arr(m_ad_vec);
			}

			alpha_cg = dot_gd / dot_dAd;

			if (alpha_f < alpha_cg)
			{
				alpha_cg = alpha_f;
			}

			ptr_m_g = &m_g[0];
			ptr_m_p = &m_p[0];
			ptr_m_d = &m_d[0];
			ptr_m_Lo = &m_Lo[0];
			ptr_m_Hi = &m_Hi[0];
			ptr_m_ad_vec = &m_ad_vec[0];
			ptr_m_lambda = &m_lambda[0];

			BetaNorm = 0.0f;
			Phi_Phi = 0.0f;

			// Calculate:
			//		x_k+1 = x_k - alpha_cg * d;
			//		g = g - alpha_cg * A * d;
			//		p = phi(x_k+1);
			//		recalculate BetaNorm, Phi_phi
			for (i = 0; i < m_NumJoints; ++i)
			{
				val_arr(m_lambda) -= alpha_cg * inc_arr(m_d);
				val_arr(m_g) -= alpha_cg * inc_arr(m_ad_vec);

				// Calculate phi from gradient, and p = phi(x_k+1)

#if (SW_LCPCG_PRECONDITIONER_FACEJACOBI == 1)
				// Preconditioner in face:
				// p = z
				val_arr(m_p) = ptr_m_invDiag[i] * val_arr(m_g);
#else
				val_arr(m_p) = val_arr(m_g);
#endif

				// UPD: LoHi
				tmp_float2 = val_arr(m_lambda) - ptr_m_Lo[i];
				tmp_float3 = val_arr(m_lambda) - ptr_m_Hi[i];

				++ptr_m_lambda;

				// Check if we're at Active set (--> Lo <--)[TODO: Optimize]
				if (math::fastabs(tmp_float2) < ACT_SET_THRESHOLD)
				{
					val_arr(m_p) = 0.0f;

					// Active set
					// Here we calculate chopped gradient Beta, hence BetaNorm updated
					// beta = g- = min{g, 0.0f};

					tmp_float = inc_arr(m_g);
					if (tmp_float < 0.0f)
					{
						BetaNorm += tmp_float * tmp_float;
						grad_norm += tmp_float * tmp_float;
					}
				}
				// Check if we're at Active set (--> Hi <--)[TODO: Optimize]
				else if (math::fastabs(tmp_float3) < ACT_SET_THRESHOLD)
				{
					val_arr(m_p) = 0.0f;

					// Active set
					// Here we calculate chopped gradient Beta, hence BetaNorm updated
					// beta = g+ = max{g, 0.0f};

					tmp_float = inc_arr(m_g);
					if (tmp_float > 0.0f)
					{
						BetaNorm += tmp_float * tmp_float;
						grad_norm += tmp_float * tmp_float;
					}
				}
				else
				{
					// Free set
					// Here we calculate free gradient phi and reduced free gradient phi~
					// phi~(alpha) = min{(x_i - lo_i) / alpha, phi}

					tmp_float = inc_arr(m_g);

					grad_norm += tmp_float*tmp_float;

// 					if (!M_EQUAL(ptr_m_Lo[i], -FLT_MAX) && tmp_float2 / alpha_fucked < tmp_float)
// 					{
// 						Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
// 					}
// 					else
// 					{
// 						Phi_Phi += tmp_float * tmp_float;
// 					}

					// We're going to Lo
					if (tmp_float > 0.0f)
					{
						if (tmp_float2 / alpha_fucked < tmp_float)
						{
 							Phi_Phi += tmp_float * tmp_float2 / alpha_fucked;
						} else
						{
 							Phi_Phi += tmp_float * tmp_float;
						}
					}
					// We're going to Hi
					else
					{
						if (tmp_float3 / alpha_fucked > tmp_float)
						{
 							Phi_Phi += tmp_float * tmp_float3 / alpha_fucked;
						} else
						{
 							Phi_Phi += tmp_float * tmp_float;
						}
					}
				}

				++ptr_m_p;
			}

#endif
//  			gLog.Print("Prop Step [%f]: Phi[%f], Beta[%f]", grad_norm, Phi_Phi, BetaNorm);
		}

		
#if (SOLVER_OBJECTIVE_FUNCTION_OUTPUT == 1)

		ptr_m_g = &m_g[0];
		ptr_m_RHS = &m_RHS[0];
		ptr_m_lambda = &m_lambda[0];

		float lambda_norm = 0.0f;
		float interm_lambda_dot_grad = 0.0f;
		for (i = 0; i < m_NumJoints; ++i)
		{
			lambda_norm += val_arr(m_lambda) * val_arr(m_lambda);

			interm_lambda_dot_grad += inc_arr(m_g) * val_arr(m_lambda);
			
			//  - x dot b
			interm_lambda_dot_grad -= inc_arr(m_RHS) * inc_arr(m_lambda);
		}

		gLog.Print("[%12.10f]----------------------------------------> 0.5 * (x * (Ax - b) - x * b) = %f;", lambda_norm, 0.5f * interm_lambda_dot_grad);

#endif

		//assert(grad_norm < prev_grad_norm);
		//prev_grad_norm = grad_norm;

#if (SOLVER_STOPPING_CRITERIA == 1)
		// Stopping criteria
		if (grad_norm / m_NumJoints < m_Precision)
		{
			m_EffectiveIterations = k;
			break;
		}
#endif
	}

	ptr_m_a_x = &m_a_x[0];
	ptr_m_a_y = &m_a_y[0];
	ptr_m_a_z = &m_a_z[0];

	ptr_m_lambda = &m_lambda[0];
	
	INIT_A_SPARSE(0);
	INIT_NODES_IDX(0);

	m_a_x.Empty();
	m_a_y.Empty();
	m_a_z.Empty();

#if (SOLVERS_ANALYZE_RESIDUAL == 1)
	ptr_m_g = &m_g[0];

	m_GradNormSq = 0.0f;
	m_LambdaNormSq = 0.0f;
	m_DotLambdaGrad = 0.0f;
#endif

	for (i = 0; i < m_NumJoints; ++i)
	{
#if (SOLVERS_ANALYZE_RESIDUAL == 1)
		// g = Ax - b;
		m_resid[i] = val_arr(m_g);

		m_GradNormSq += val_arr(m_g) * val_arr(m_g);
		m_LambdaNormSq += val_arr(m_lambda) * val_arr(m_lambda);
		m_DotLambdaGrad += val_arr(m_lambda) * val_arr(m_g); 

		++ptr_m_g;
#endif

#if (SW_LCPCG_PRECONDITIONER_DOUBLEJACOBI == 1)
		// Calculate original lambda:
		// lambda_or = D^1/2 * lambda
		// and get a vector with it
		val_arr(m_lambda) *= ptr_m_invDiag[i];
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