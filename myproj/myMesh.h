#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();
	bool readFile(std::string filename);
	void computeNormals();
	void normalize();

	void triangulate();
	bool triangulate(myFace *);

	void clear();

	myMesh(void);
	~myMesh(void);

	void inflateMesh(double dist);
	void smoothenMesh(double dist);

	bool subdivisionCatmullClark();

	void splitFace(myFace *f, myPoint3D *p);
	void splitEdge(myHalfedge *, myPoint3D *);

	bool simplification();
};
