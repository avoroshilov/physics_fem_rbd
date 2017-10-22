#include "helpers/log.h"
#include "FEM.h"
#include <float.h>

#define TEST_TETRA_SUITABLE 0

extern CLogChannel gLog;

inline bool isTetraSuitable(Vec3 &v1, Vec3 &v2, Vec3 &v3, Vec3 &v4, float Volume,
							float *oblongMin, float *oblongMax, float *volumeMin, float *volumeMax,
							float oblongLo, float oblongHi, float volumeLo, float volumeHi)
{
	float volAbs = fabsf(Volume);

	// Oblongness

	// Area of faces
	float faceAreas[4];

	faceAreas[0] = 0.5f * (v2 - v1).Cross(v3 - v1).GetLength();
	faceAreas[1] = 0.5f * (v4 - v1).Cross(v2 - v1).GetLength();
	faceAreas[2] = 0.5f * (v4 - v1).Cross(v3 - v1).GetLength();
	faceAreas[3] = 0.5f * (v4 - v3).Cross(v2 - v3).GetLength();

	// Heights to respective faces
	float heights[4];

	for (int i = 0; i < 4; ++i)
		heights[i] = 3.0f * volAbs / faceAreas[i];

	float oblMinLoc = FLT_MAX, oblMaxLoc = -FLT_MAX;

	float tmp_val;
	
	for (int i = 0; i < 4; ++i)
	{
		tmp_val = heights[i] / faceAreas[i];

		if (tmp_val < oblMinLoc)
			oblMinLoc = tmp_val;
		if (tmp_val > oblMaxLoc)
			oblMaxLoc = tmp_val;

	}


	if (Volume < volumeLo || Volume > volumeHi)
		return false;

	if (oblMinLoc < oblongLo || oblMaxLoc > oblongHi)
		return false;


	//////////////////////////////////////////////////////////////////////////
	if (oblMinLoc < *oblongMin)
		*oblongMin = oblMinLoc;
	if (oblMaxLoc > *oblongMax)
		*oblongMax = oblMaxLoc;

	if (volAbs < *volumeMin)
		*volumeMin = volAbs;
	if (volAbs > *volumeMax)
		*volumeMax = volAbs;

	return true;
}

// pyr_ind - pyramid's base indices, CW (from the vertex)
INLINE void FormTetraFromPyramid(FEMObject &Object, unsigned pyr_vertex,
								 unsigned pyr_ind1, unsigned pyr_ind2, unsigned pyr_ind3, unsigned pyr_ind4)
{
	FEMTetrahedron tmpTetra;

	// Tetra1
	tmpTetra.ind1 = pyr_ind1;
	tmpTetra.ind2 = pyr_ind2;
	tmpTetra.ind3 = pyr_ind3;
	tmpTetra.ind4 = pyr_vertex;
	Object.FEMTetras.push_back(tmpTetra);

	// Tetra2
	tmpTetra.ind1 = pyr_ind3;
	tmpTetra.ind2 = pyr_ind4;
	tmpTetra.ind3 = pyr_ind1;
	tmpTetra.ind4 = pyr_vertex;
	Object.FEMTetras.push_back(tmpTetra);
}

void BuildCube(FEMObject &Object, const Vec3 &Pos, const Vec3 &Extents, const CMatrix3 &Rot,
			   unsigned SegWidth, unsigned SegHeight, unsigned SegDepth)
{
	unsigned i, j, k;

	float dX = 2.0f * Extents.x / (float)(SegWidth),
		  dY = 2.0f * Extents.y / (float)(SegHeight),
		  dZ = 2.0f * Extents.z / (float)(SegDepth);

	Vec3 CurNode;

	unsigned VertexNum = Object.FEMNodes.size();
	Object.FEMNodes.reserve(VertexNum + (SegWidth + 1) * (SegHeight + 1) * (SegDepth + 1) + SegWidth * SegHeight * SegDepth);

	// Place cubes' border
	CurNode.x = -Extents.x;
	for (i = 0; i < SegWidth + 1; ++i, CurNode.x += dX)
	{
		CurNode.y = -Extents.y;
		for (j = 0; j < SegHeight + 1; ++j, CurNode.y += dY)
		{
			CurNode.z = -Extents.z;
			for (k = 0; k < SegDepth + 1; ++k, CurNode.z += dZ)
			{
				Object.FEMNodes.push_back(Rot * CurNode + Pos);
			}
		}
	}

	// Place cubes' middle
	CurNode.x = -Extents.x + dX * 0.5f;
	for (i = 0; i < SegWidth; ++i, CurNode.x += dX)
	{
		CurNode.y = -Extents.y + dY * 0.5f;
		for (j = 0; j < SegHeight; ++j, CurNode.y += dY)
		{
			CurNode.z = -Extents.z + dZ * 0.5f;
			for (k = 0; k < SegDepth; ++k, CurNode.z += dZ)
			{
				Object.FEMNodes.push_back(Rot * CurNode + Pos);
			}
		}
	}

	unsigned MidShift = VertexNum + (SegWidth + 1) * (SegHeight + 1) * (SegDepth + 1);

	// 12 - number of tetrahedrons per cube (6 pyramids, 2 tetra per pyramid)
	Object.FEMTetras.reserve(Object.FEMTetras.size() + SegWidth * SegHeight * SegDepth * 12);

#define idx(i, j, k)	(VertexNum + (i) * (SegHeight + 1) * (SegDepth + 1) + (j) * (SegDepth + 1) + (k))

	// Mark up tetrahedrons
	for (i = 0; i < SegWidth; ++i)
	{
		for (j = 0; j < SegHeight; ++j)
		{
			for (k = 0; k < SegDepth; ++k)
			{
				unsigned MidIndex = (i * (SegHeight * SegDepth) + j * SegDepth + k) + MidShift;

				// PYRAMID TOP
				FormTetraFromPyramid(Object, MidIndex,	idx( i , j+1,  k ), idx( i , j+1, k+1),
														idx(i+1, j+1, k+1), idx(i+1, j+1,  k ));

				// PYRAMID BOTTOM
				FormTetraFromPyramid(Object, MidIndex,  idx( i ,  j , k+1), idx( i ,  j ,  k ),
														idx(i+1,  j ,  k ), idx(i+1,  j , k+1));

				// PYRAMID RIGHT
				FormTetraFromPyramid(Object, MidIndex,  idx(i+1, j+1, k+1), idx(i+1,  j , k+1),
														idx(i+1,  j ,  k ), idx(i+1, j+1,  k ));

				// PYRAMID LEFT
				FormTetraFromPyramid(Object, MidIndex,  idx( i , j+1,  k ), idx( i ,  j ,  k ),
														idx( i ,  j , k+1), idx( i , j+1, k+1));

				// PYRAMID BACK
				FormTetraFromPyramid(Object, MidIndex,  idx(i+1, j+1,  k ), idx(i+1,  j ,  k ),
														idx( i ,  j ,  k ), idx( i , j+1,  k ));

				// PYRAMID FRONT
				FormTetraFromPyramid(Object, MidIndex,  idx( i , j+1, k+1), idx( i ,  j , k+1),
														idx(i+1,  j , k+1), idx(i+1, j+1, k+1));
			}
		}
	}

#undef idx
}

void ReadLineNextSection(FILE *fp, char *lineBuffer, unsigned int bufSize)
{
	int readComments = 1;

	do
	{
		fgets(lineBuffer, 512, fp);

		// If comment symbol or empty string
		if ((lineBuffer[0] != '#') && (strlen(lineBuffer) > 3))
			readComments = 0;

	} while (readComments != 0);
}

int BuildFromFile(char *fileName, FEMObject &Object, const Vec3 &Pos, const Vec3 &Scale, const CMatrix3 &Rot)
{
	FILE *fp = fopen(fileName, "r");

	if (fp == NULL)
		return 0;

	const int bufSize = 512;
	char fileLineBuf[bufSize];

	// Vertices
	//////////////////////////////////////////////////////////////////////////

	ReadLineNextSection(fp, fileLineBuf, bufSize);

	// We already have 1 vertex in buffer
	unsigned int vertCount = 1;

	// READ VERTICES
	Vec3 CurNode;
	while (true)
	{
		// Analyze && append vertex into FEMObject
		sscanf(fileLineBuf, "v %f %f %f\n", &CurNode.x, &CurNode.y, &CurNode.z);

		CurNode.x *= Scale.x;
		CurNode.y *= Scale.y;
		CurNode.z *= Scale.z;

		Object.FEMNodes.push_back(Rot * CurNode + Pos);
		//////////////////////////////////////////////////////////////////////////

		char *res_ptr = fgets(fileLineBuf, 512, fp);

		if ((fileLineBuf[0] != 'v') || (res_ptr == NULL))
			break;

		++vertCount;
	}

//	printf("> Vertex count: %u;\n\n", vertCount);



	// Tetrahedra
	//////////////////////////////////////////////////////////////////////////

	ReadLineNextSection(fp, fileLineBuf, bufSize);

	// We already have 1 vertex in buffer
	unsigned int tetCount = 1;

	float volMin = FLT_MAX, volMax = -FLT_MAX;
	float oblMin = FLT_MAX, oblMax = -FLT_MAX;

	// READ VERTICES
	FEMTetrahedron CurTetra;
	while (true)
	{
		// Analyze && append tetrahedra into FEMObject
		sscanf(fileLineBuf, "t %u %u %u %u\n", &CurTetra.ind1, &CurTetra.ind2, &CurTetra.ind3, &CurTetra.ind4);

		// Shift due to first zero-vertex
		++CurTetra.ind1;
		++CurTetra.ind2;
		++CurTetra.ind3;
		++CurTetra.ind4;

		Vec3 & v1 = Object.FEMNodes[CurTetra.ind1];
		Vec3 & v2 = Object.FEMNodes[CurTetra.ind2];
		Vec3 & v3 = Object.FEMNodes[CurTetra.ind3];
		Vec3 & v4 = Object.FEMNodes[CurTetra.ind4];

		float Volume = (v2 - v1).Triple(v3 - v1, v4 - v1);
		if (Volume > 0.0f)
		{
			// Swap
			unsigned int tmp = CurTetra.ind1;
			CurTetra.ind1 = CurTetra.ind3;
			CurTetra.ind3 = tmp;
		}

		//////////////////////////////////////////////////////////////////////////

		char *res_ptr = fgets(fileLineBuf, 512, fp);

		if ((fileLineBuf[0] != 't') || (res_ptr == NULL))
			break;

#if (TEST_TETRA_SUITABLE == 1)

		/* DBG */
		Vec3 planePoint = Vec3(-4.5f, 2.5f, 0.0f);
		Vec3 planeNormal = Vec3(-1.0f, 1.0f, 0.0f);
		planeNormal.Normalize();

/*
		if ( (v1 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
		if ( (v2 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
		if ( (v3 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
		if ( (v4 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
*/

		float TetraSuitable = isTetraSuitable(v1, v2, v3, v4, Volume, &oblMin, &oblMax, &volMin, &volMax,
//												10.0f, 150.0f, 0.005f, FLT_MAX);
												15.0f, 150.0f, 0.005f, 0.05f);

		if (!TetraSuitable)
			continue;

#endif

		//////////////////////////////////////////////////////////////////////////

		Object.FEMTetras.push_back(CurTetra);
		++tetCount;

/*
		if (tetCount > 2000)
			break;
*/
	}

	gLog.Print("[TET] Mesh: %s;", fileName);
	gLog.Print("Volume min/max: %f/%f", volMin, volMax);
	gLog.Print("Oblongness min/max: %f/%f;", oblMin, oblMax);

//	printf("> Tetrahedra count: %u;\n", tetCount);

	fclose(fp);	
}

int BuildFromFile_Tetgen(char *fileNameNode, char *fileNameEle, FEMObject &Object, const Vec3 &Pos, const Vec3 &Scale, const CMatrix3 &Rot)
{
	FILE *fpNode = fopen(fileNameNode, "r");

	if (fpNode == NULL)
		return 0;

	FILE *fpEle = fopen(fileNameEle, "r");

	if (fpEle == NULL)
		return 0;

	const int bufSize = 512;
	char fileLineBuf[bufSize];

	// Vertices
	//////////////////////////////////////////////////////////////////////////

	// First line is not needed
	fgets(fileLineBuf, 512, fpNode);

	ReadLineNextSection(fpNode, fileLineBuf, bufSize);

	// We already have 1 vertex in buffer
	unsigned int vertCount = 1;

	// READ VERTICES
	Vec3 CurNode;
	while (true)
	{
		unsigned int vertIdx;

		// Analyze && append vertex into FEMObject
		sscanf(fileLineBuf, "%u %f %f %f\n", &vertIdx, &CurNode.x, &CurNode.y, &CurNode.z);

		CurNode.x *= Scale.x;
		CurNode.y *= Scale.y;
		CurNode.z *= Scale.z;

		Object.FEMNodes.push_back(Rot * CurNode + Pos);
		//////////////////////////////////////////////////////////////////////////

		char *res_ptr = fgets(fileLineBuf, 512, fpNode);

		if ((fileLineBuf[0] == '#') || (res_ptr == NULL))
			break;

		++vertCount;
	}

//	printf("> Vertex count: %u;\n\n", vertCount);

	fclose(fpNode);	


	// Tetrahedra
	//////////////////////////////////////////////////////////////////////////

	// First line is not needed
	fgets(fileLineBuf, 512, fpEle);

	ReadLineNextSection(fpEle, fileLineBuf, bufSize);

	// We already have 1 vertex in buffer
	unsigned int tetCount = 1;

	float volMin = FLT_MAX, volMax = -FLT_MAX;
	float oblMin = FLT_MAX, oblMax = -FLT_MAX;

	// READ VERTICES
	FEMTetrahedron CurTetra;
	while (true)
	{
		unsigned int tetIdx;

		// Analyze && append tetrahedra into FEMObject
		sscanf(fileLineBuf, "%u %u %u %u %u\n", &tetIdx, &CurTetra.ind1, &CurTetra.ind2, &CurTetra.ind3, &CurTetra.ind4);

		// Shift due to first zero-vertex
		++CurTetra.ind1;
		++CurTetra.ind2;
		++CurTetra.ind3;
		++CurTetra.ind4;

		Vec3 & v1 = Object.FEMNodes[CurTetra.ind1];
		Vec3 & v2 = Object.FEMNodes[CurTetra.ind2];
		Vec3 & v3 = Object.FEMNodes[CurTetra.ind3];
		Vec3 & v4 = Object.FEMNodes[CurTetra.ind4];

		float Volume = (v2 - v1).Triple(v3 - v1, v4 - v1);
		if (Volume > 0.0f)
		{
			// Swap
			unsigned int tmp = CurTetra.ind1;
			CurTetra.ind1 = CurTetra.ind3;
			CurTetra.ind3 = tmp;
		}

		char *res_ptr = fgets(fileLineBuf, 512, fpEle);

		if ((fileLineBuf[0] == '#') || (res_ptr == NULL))
			break;

#if 0

		float TetraSuitable = isTetraSuitable(v1, v2, v3, v4, Volume, &oblMin, &oblMax, &volMin, &volMax,
												10.0f, 150.0f, 0.05f, FLT_MAX);
//												0.01f, FLT_MAX, 0.01f, FLT_MAX);

		if (!TetraSuitable)
			continue;
		//////////////////////////////////////////////////////////////////////////

		/* DBG */
		Vec3 planePoint = Vec3(-4.0f, 2.5f, 0.0f);
		Vec3 planeNormal = Vec3(-1.0f, 1.0f, 0.0f);
		planeNormal.Normalize();

		if ( (v1 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
		if ( (v2 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
		if ( (v3 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
		if ( (v4 - planePoint).Dot(planeNormal) < 0.0f )
			continue;
#endif

		Object.FEMTetras.push_back(CurTetra);
		++tetCount;

/*
		if (tetCount > 2000)
			break;
*/
	}

//	printf("> Tetrahedra count: %u;\n", tetCount);
	gLog.Print("[NODE/ELE] Mesh: %s / %s;", fileNameNode, fileNameEle);
	gLog.Print("Volume min/max: %f/%f", volMin, volMax);
	gLog.Print("Oblongness min/max: %f/%f;", oblMin, oblMax);

	fclose(fpEle);	
}

void FuncAABBMin(Vec3 &out, Vec3 **tetraVerts)
{
	for (int i = 0; i < 3; ++i)
	{
		out.v[i] = tetraVerts[0]->v[i];
		if (out.v[i] > tetraVerts[1]->v[i])
			out.v[i] = tetraVerts[1]->v[i];
		if (out.v[i] > tetraVerts[2]->v[i])
			out.v[i] = tetraVerts[2]->v[i];
		if (out.v[i] > tetraVerts[3]->v[i])
			out.v[i] = tetraVerts[3]->v[i];
	}
}

void FuncAABBMax(Vec3 &out, Vec3 **tetraVerts)
{
	for (int i = 0; i < 3; ++i)
	{
		out.v[i] = tetraVerts[0]->v[i];
		if (out.v[i] < tetraVerts[1]->v[i])
			out.v[i] = tetraVerts[1]->v[i];
		if (out.v[i] < tetraVerts[2]->v[i])
			out.v[i] = tetraVerts[2]->v[i];
		if (out.v[i] < tetraVerts[3]->v[i])
			out.v[i] = tetraVerts[3]->v[i];
	}
}

bool CheckAxis(const Vec3 &Axis, Vec3 **tetra0Verts, Vec3 **tetra1Verts)
{
	float t0Min, t0Max;
	float t1Min, t1Max;
	float tmp;

	// Project && find minima
	t0Min = tetra0Verts[0]->Dot(Axis);
	t0Max = t0Min;

	for (int i = 1; i < 4; ++i)
	{
		tmp = tetra0Verts[i]->Dot(Axis);

		if (tmp < t0Min)
			t0Min = tmp;

		if (tmp > t0Max)
			t0Max = tmp;
	}

	t1Min = tetra1Verts[0]->Dot(Axis);
	t1Max = t1Min;

	for (int i = 1; i < 4; ++i)
	{
		tmp = tetra1Verts[i]->Dot(Axis);

		if (tmp < t1Min)
			t1Min = tmp;

		if (tmp > t1Max)
			t1Max = tmp;
	}

	float SAT_EPS = 1e-5f;

	t0Min += SAT_EPS;
	t0Max -= SAT_EPS;
	t1Min += SAT_EPS;
	t1Max -= SAT_EPS;

	if (t0Max < t1Min)
		return true;

	if (t0Min > t1Max)
		return true;

	return false;
}

// TRUE if separation
bool CheckFaceNormal(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, Vec3 **tetra0Verts, Vec3 **tetra1Verts)
{
	Vec3 Normal = ( (v2 - v1).Cross(v0 - v1) ).GetNormalized();

	return CheckAxis(Normal, tetra0Verts, tetra1Verts);
}

bool CheckEdges(Vec3 **tetra0Verts, Vec3 **tetra1Verts)
{
	for (unsigned t0v0 = 0; t0v0 < 4; ++t0v0)
	{
		for (unsigned t0v1 = t0v0 + 1; t0v1 < 4; ++t0v1)
		{
			Vec3 edge0 = *(tetra0Verts[t0v1]) - *(tetra0Verts[t0v0]);

			for (unsigned t1v0 = 0; t1v0 < 4; ++t1v0)
			{
				for (unsigned t1v1 = t1v0 + 1; t1v1 < 4; ++t1v1)
				{
					Vec3 edge1 = *(tetra1Verts[t1v1]) - *(tetra1Verts[t1v0]);

					Vec3 cross = edge0.Cross(edge1);

					// Edges share points
					if (cross.GetSquareLength() < 1e-4)
						continue;

					cross.Normalize();

					bool isSeparation = CheckAxis(cross, tetra0Verts, tetra1Verts);

					if (isSeparation)
						return true;
				}
			}
		}
	}

	return false;
}

void NoSeparation(unsigned tetra0, unsigned tetra1, const FEMObject &Object)
{
	gLog.Print("Self-collision! %u vs %u;", tetra0, tetra1);
}

#define CHECK_AXIS(axis) \
	if (t0_aabbMax.axis < t1_aabbMin.axis) \
		continue; \
	if (t0_aabbMin.axis > t1_aabbMax.axis) \
		continue;
	
#define CHECK_CALL(call) \
	isSeparation = call; \
	if (isSeparation) \
		continue;


void CheckSelfIntersections(FEMObject &Object)
{
	unsigned int tetrasNum = Object.FEMTetras.size();
	for (unsigned int i = 0; i < tetrasNum; ++i)
	{
		Vec3 *t0[4] = 
			{
				&Object.FEMNodes[Object.FEMTetras[i].ind1],
				&Object.FEMNodes[Object.FEMTetras[i].ind2],
				&Object.FEMNodes[Object.FEMTetras[i].ind3],
				&Object.FEMNodes[Object.FEMTetras[i].ind4]
			};

		Vec3 t0_aabbMin, t0_aabbMax;

		FuncAABBMin(t0_aabbMin, t0);
		FuncAABBMax(t0_aabbMax, t0);

		for (unsigned int j = i + 1; j < tetrasNum; ++j)
		{
			Vec3 *t1[4] = 
				{
					&Object.FEMNodes[Object.FEMTetras[j].ind1],
					&Object.FEMNodes[Object.FEMTetras[j].ind2],
					&Object.FEMNodes[Object.FEMTetras[j].ind3],
					&Object.FEMNodes[Object.FEMTetras[j].ind4]
				};

			Vec3 t1_aabbMin, t1_aabbMax;

			FuncAABBMin(t1_aabbMin, t1);
			FuncAABBMax(t1_aabbMax, t1);

			// Brute-force O(n*n) broadphase
			CHECK_AXIS(x);
			CHECK_AXIS(y);
			CHECK_AXIS(z);

			// Run SAT
			bool isSeparation;

			// Face Normals
			CHECK_CALL( CheckFaceNormal(*(t0[0]), *(t0[3]), *(t0[2]), t0, t1) );
			CHECK_CALL( CheckFaceNormal(*(t0[2]), *(t0[3]), *(t0[1]), t0, t1) );
			CHECK_CALL( CheckFaceNormal(*(t0[3]), *(t0[1]), *(t0[0]), t0, t1) );
			CHECK_CALL( CheckFaceNormal(*(t0[0]), *(t0[2]), *(t0[1]), t0, t1) );

			CHECK_CALL( CheckFaceNormal(*(t1[0]), *(t1[3]), *(t1[2]), t0, t1) );
			CHECK_CALL( CheckFaceNormal(*(t1[2]), *(t1[3]), *(t1[1]), t0, t1) );
			CHECK_CALL( CheckFaceNormal(*(t1[3]), *(t1[1]), *(t1[0]), t0, t1) );
			CHECK_CALL( CheckFaceNormal(*(t1[0]), *(t1[2]), *(t1[1]), t0, t1) );

			// Edge Crosses (36 tests)
			CHECK_CALL( CheckEdges(t0, t1) );

			NoSeparation(i, j, Object);
		}
	}


}

#undef CHECK_AXIS