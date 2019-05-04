#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(0.0, 0.0, 0.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
	myHalfedge *this_edge = adjacent_halfedge;
	myHalfedge *next_edge = this_edge->next;
	myHalfedge *prev_edge = this_edge->prev;

	myPoint3D *prev_point = prev_edge->source->point;
	myPoint3D *this_point = this_edge->source->point;
	myPoint3D *next_point = next_edge->source->point;

	myVector3D result = myVector3D(0.0, 0.0, 0.0);
	result.setNormal(prev_point, this_point, next_point);

	normal->dX = result.dX;
	normal->dY = result.dY;
	normal->dZ = result.dZ;
}
