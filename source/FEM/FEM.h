#ifndef _FEM_H_
#define _FEM_H_

#include <vector>
#include "helpers/globalallocator.h"
#include "math/m_Vector3.h"
#include "math/m_Matrix3.h"

struct FEMTetrahedron
{
	unsigned ind1, ind2, ind3, ind4;
};

struct FEMObject
{
	std::vector<Vec3, GlobalAllocator<Vec3>> FEMNodes;
	std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>> FEMTetras;
};

// Extents - Cube's half-size
void BuildCube(FEMObject &Object, const Vec3 &Pos, const Vec3 &Extents, const CMatrix3 &Rot,
			   unsigned SegWidth, unsigned SegHeight, unsigned SegDepth);

int BuildFromFile(char *fileName, FEMObject &Object, const Vec3 &Pos, const Vec3 &Scale, const CMatrix3 &Rot);

int BuildFromFile_Tetgen(char *fileNameNode, char *fileNameEle, FEMObject &Object, const Vec3 &Pos, const Vec3 &Scale, const CMatrix3 &Rot);

void CheckSelfIntersections(FEMObject &Object);

#endif