#include <iostream>
#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,0.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
	/**** TODO ****/

	myHalfedge* sample = originof; int count = 0;
	myVector3D result = myVector3D(0, 0 ,0);

	do {//sample->adjacent_face->computeNormal();
		myVector3D *faceNormal = sample->adjacent_face->normal;
		result += *faceNormal; count += 1;
		if (sample->prev->twin == NULL)
			std::cout << "NULL twin" << std::endl;

		sample = sample->prev->twin;

	} while (sample != originof);

	result = result/count;

	normal->dX = result.dX;
	normal->dY = result.dY;
	normal->dZ = result.dZ;
}
