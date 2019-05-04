#pragma once
#include "mypoint3d.h"

class myHalfedge;
class myVector3D;

class myVertex
{
public:
	myPoint3D *point;
	myHalfedge *originof;

	int index;  //use as you wish.
	// It is originally the reading order from file.

	myVector3D *normal;

	void computeNormal();
	myVertex(void);
	~myVertex(void);
};
