#include "myHalfedge.h"

myHalfedge::myHalfedge(void)
{
	source = NULL; 
	adjacent_face = NULL; 
	next = NULL;  
	prev = NULL;  
	twin = NULL;  
}

void myHalfedge::copy(myHalfedge *ie)
{
/**** TODO ****/
}

myHalfedge::~myHalfedge(void)
{
}
