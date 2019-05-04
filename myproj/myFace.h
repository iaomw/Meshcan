#pragma once

class myHalfedge;
class myVector3D;

class myFace
{
public:
	myHalfedge *adjacent_halfedge;

	myVector3D *normal;

	int index; //use this variable as you wish.

	void computeNormal();
	myFace(void);
	~myFace(void);
};