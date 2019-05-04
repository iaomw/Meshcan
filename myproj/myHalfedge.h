#pragma once

#include <stdio.h>
#include <tchar.h>

class myVertex;
class myFace;
class myPoint3D;

class myHalfedge
{
public:
	myVertex *source; 
	myFace *adjacent_face; 
	myHalfedge *next;  
	myHalfedge *prev;  
	myHalfedge *twin;  

	int index; //use as you wish.

	myHalfedge(void);
	void copy(myHalfedge *);
	~myHalfedge(void);
};