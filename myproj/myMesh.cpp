#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	vector<myHalfedge *>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL)
			break;
	}
	if (it != halfedges.end())
		std:: cout << "Error! Not all edges have their twins!\n";
	else std::cout << "Each edge has a twin!\n";
}



bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;
	myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v")
		{
			
			double x, y, z; 
			myline >> u; x = atof(u.c_str()); 
			myline >> u; y = atof(u.c_str()); 
			myline >> u; z = atof(u.c_str());

			myVertex *vector = new myVertex();
			vector->point = new myPoint3D(x, y, z);
			vector->index = vertices.size();
			vertices.push_back(vector);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			myFace *face = new myFace();
			face->index = faces.size();
			faces.push_back(face);
			
			myHalfedge *prev_edge = NULL, *first_edge = NULL;
			int index_first = halfedges.size();
			while (myline >> u) {
				int index = atoi((u.substr(0, u.find("/"))).c_str()); index -= 1;

				myHalfedge* this_edge = new myHalfedge();
				this_edge->adjacent_face = face;
				this_edge->source = vertices[index];
				if (vertices[index] == NULL) {
					cout << index << " is NULL" << endl;
				}
				if (prev_edge != NULL) {
					this_edge->prev = prev_edge;
					prev_edge->next = this_edge;
				}
				if (vertices[index]->originof == NULL) {
					vertices[index]->originof = this_edge;
				}

				this_edge->index = halfedges.size();
				halfedges.push_back(this_edge);
				prev_edge = this_edge;
			}

			first_edge = halfedges[index_first];
			first_edge->prev = prev_edge;
			prev_edge->next = first_edge;

			face->adjacent_halfedge = first_edge; 

			myHalfedge *twin_edge = NULL, *this_edge = first_edge; 

			do {
				int source = this_edge->source->index;
				int dest = this_edge->next->source->index;

				pair<int, int> this_key = make_pair(source, dest);
				pair<int, int> twin_key = make_pair(dest, source);

				twin_map[twin_key] = this_edge;

				map<pair<int, int>, myHalfedge*>::const_iterator it = twin_map.find(this_key);
				if (it == twin_map.end()) {
					//handle the error
				} else {

					twin_edge = it->second; //twin_map[this_key];
					twin_edge->twin = this_edge;
					this_edge->twin = twin_edge;
				}

				this_edge = this_edge->next;
			} while (this_edge != first_edge);
		}
	}

	checkMesh();
	normalize();

	return true;
}


void myMesh::computeNormals()
{
	/**** TODO ****/

	for (auto element : faces) {
		element->computeNormal();
	}

	for (auto element : vertices) {
		element->computeNormal();
	}
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}

bool myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
	//triangulate();

	std::vector<myPoint3D*> facePoints(faces.size()); 
	std::vector<myPoint3D*> middPoints(halfedges.size());
	std::vector<myVertex*> edgeVertexes(halfedges.size()); 

	//int i = 0; for (auto v : vertices) { v->index = i; i++; }

	int fvi = vertices.size(); int fvc = 0;
	for (auto f : faces) {
		myPoint3D *r = new myPoint3D(0,0,0);
		myHalfedge *e = f->adjacent_halfedge;
		myHalfedge *sample = e;
		int count = 0; do {
			*r += *(sample->source->point);
			count++; sample = sample->next;
		} while (sample != e);
		*r /= count;
		facePoints[f->index] = r;

		auto fv = new myVertex();
		fv->index = fvi+fvc; fvc++;
		vertices.push_back(fv);
		fv->point = r; 
	}

	int evi = vertices.size(); int evc = 0;
	for (auto e : halfedges) {
		myHalfedge* t = e->twin;
		myPoint3D* s = e->source->point;
		myPoint3D* z = t->source->point;

		if (middPoints[e->index] == NULL && middPoints[t->index] == NULL) {
			myPoint3D *midd = new myPoint3D(0,0,0); *midd += (*s + *z)/2;
			middPoints[e->index] = midd;
			middPoints[t->index] = midd;
		} // midd point

		if (edgeVertexes[e->index] == NULL && edgeVertexes[t->index] == NULL) {

			int index_b = e->adjacent_face->index;
			myPoint3D* u = facePoints[index_b];
			int index_d = t->adjacent_face->index;
			myPoint3D* v = facePoints[index_d];

			myPoint3D *r = new myPoint3D(0,0,0);
			*r += (*s + *z + *u + *v)/4;
			myVertex *ev = new myVertex();
			ev->point = r;

			edgeVertexes[e->index] = ev;
			edgeVertexes[t->index] = ev;

			ev->index = evi+evc; evc++;
			vertices.push_back(ev);
		} // edge point
	}

	for (int i = 0; i < fvi; i++) {
		auto v = vertices[i];
		auto q = new myPoint3D(0,0,0);
		auto r = new myPoint3D(0,0,0);
		auto o = new myPoint3D(0,0,0);
		auto e = v->originof;

		auto n = 0; do {
			auto index = e->adjacent_face->index;
			*r += *middPoints[e->index];
			*q += *facePoints[index];
			n++; e = e->twin->next;
		} while (e != v->originof);

		*q = *q/n; *r = *r/n;
		*o += *q/n + (*r/n) * 2;
		*o += (*(v->point)/n)*(n-3);

		v->point->X = o->X; 
		v->point->Y = o->Y; 
		v->point->Z = o->Z; 
		delete q, r, o;
	}

	int nfc = 0, nec = 0;
	int fc = faces.size();
	int ec = halfedges.size();
	for (int i = 0; i < fc; i++) {
		auto f = faces[i];
		auto fv = vertices[fvi+i];
		auto e = f->adjacent_halfedge;
		auto ev = edgeVertexes[e->index];

		do {
			auto a = new myHalfedge(); a->source = fv;
			auto b = new myHalfedge(); b->source = ev;
			fv->originof = a; ev->originof = b;

			e = e->next; ev = edgeVertexes[e->index];

			auto vv = e->source;
			auto c = new myHalfedge(); c->source = vv;
			auto d = new myHalfedge(); d->source = ev;

			vv->originof = c; ev->originof = d;

			a->next = b; b->next = c; c->next = d; d->next = a;
			a->prev = d; b->prev = a; c->prev = b; d->prev = c;

			auto face = new myFace(); face->index = nfc;
			face->adjacent_halfedge = a; nfc++;
			faces.push_back(face);

			a->adjacent_face = face; b->adjacent_face = face;
			c->adjacent_face = face; d->adjacent_face = face;

			a->index = nec; halfedges.push_back(a); nec++;
			b->index = nec; halfedges.push_back(b); nec++;
			c->index = nec; halfedges.push_back(c); nec++;
			d->index = nec; halfedges.push_back(d); nec++;

		} while (e != f->adjacent_halfedge);
	}

	faces.erase(faces.begin(), faces.begin()+fc);
	halfedges.erase(halfedges.begin(), halfedges.begin()+ec);

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	for (auto f : faces) {
		myHalfedge *first_edge = f->adjacent_halfedge;
		myHalfedge *twin_edge = NULL, *this_edge = first_edge;

		do {
			int sour = this_edge->source->index; 
			int dest = this_edge->next->source->index;
			pair<int, int> this_key = make_pair(sour, dest);
			pair<int, int> twin_key = make_pair(dest, sour);

			twin_map[twin_key] = this_edge;
			map<pair<int, int>, myHalfedge*>::const_iterator it = twin_map.find(this_key);
			if (it == twin_map.end()) {
				//handle the error
			} else {
				twin_edge = it->second; //twin_map[this_key];
				twin_edge->twin = this_edge;
				this_edge->twin = twin_edge;
			}
			this_edge = this_edge->next;
		} while (this_edge != first_edge);
	}

	checkMesh();
	return true;
}

void myMesh::triangulate()
{
	/**** TODO ****/
	int count = faces.size();
	for (int i = 0; i < count; i++) {
		triangulate(faces[i]);
	}
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	/**** TODO ****/

	myHalfedge *first = f->adjacent_halfedge;

	if (first->next->next->next == first) {
		return false;
	} 

	myHalfedge *sample = first;
	
		myHalfedge *prev = sample->prev;
		myHalfedge *next = sample->next->next;

		myHalfedge *z = new myHalfedge(); 
		z->index = halfedges.size(); halfedges.push_back(z);
		myHalfedge *s = new myHalfedge(); 
		s->index = halfedges.size(); halfedges.push_back(s);

		z->source = sample->source; z->twin = s;
		s->source = next->source; s->twin = z;

		prev->next = z; z->prev = prev;
		next->prev = z; z->next = next;

		z->adjacent_face = f;
		f->adjacent_halfedge = z;

		triangulate(f);

		myFace *t = new myFace(); t->normal = f->normal;
		t->index = faces.size(); faces.push_back(t);
		t->adjacent_halfedge = s;
		s->adjacent_face = t;
		next = sample->next;

		sample->adjacent_face = t;
		next->adjacent_face = t;

		s->next = sample; s->prev = next; 
		sample->prev = s; next->next = s;

		checkMesh();

		return true;
}

// TP2

void myMesh::inflateMesh(double dist) {

	// Given the input parameter dist, this function moves each myVertex *v in the direction of its normal by distance dist :
	// *(v->point) = *(v->point) + (*(v->normal))*dist.
	//	Calling it gives the appearance of inﬂating the mesh.
	
}

void myMesh::smoothenMesh(double dist) {
	// Moves each vertex point in the direction of the average of it’s neighbors in the mesh.
	// If neighbors average is myPoint *X, your new point for myVertex *v should be :
	// *(v->point) = (*(v->point))*(1 - dist) + (*X)*dist.
}

void myMesh::splitEdge(myHalfedge *e, myPoint3D *p) {

	/**** TODO ****/

	myVertex *sour = e->source;
	myVertex *dest = e->next->source;

	myVertex *midd = new myVertex(); midd->point = new myPoint3D(p->X, p->Y, p->Z);
	midd->index = vertices.size(); vertices.push_back(midd);
	
	myHalfedge *s = new myHalfedge(); s->source = sour; 
	myHalfedge *z = new myHalfedge(); z->source = midd;
	myHalfedge *u = new myHalfedge(); u->source = midd;
	myHalfedge *v = new myHalfedge(); v->source = dest;

	s->index = halfedges.size(); halfedges.push_back(s);
	z->index = halfedges.size(); halfedges.push_back(z);
	u->index = halfedges.size(); halfedges.push_back(u);
	v->index = halfedges.size(); halfedges.push_back(v);

	s->twin = z; z->twin = s;
	u->twin = v; v->twin = u;

	s->next = u; u->prev = s;
	z->prev = v; v->next = z;

	s->prev = e->prev; u->next = e->next;
	s->prev->next = s; u->next->prev = u;

	s->adjacent_face = e->adjacent_face;
	u->adjacent_face = e->adjacent_face;
	e->adjacent_face->adjacent_halfedge = s;

	myFace *newFace = new myFace();
	newFace->adjacent_halfedge = e;
	newFace->index = faces.size();
	e->adjacent_face = newFace;
	v->adjacent_face = newFace;
	z->adjacent_face = newFace;
	e->next = v; e->prev = z;
	v->prev = e; z->next = e;
	faces.push_back(newFace);

	midd->originof = z;

	checkMesh();
}

void myMesh::splitFace(myFace *f, myPoint3D *p) {

	myVertex *origin = new myVertex();
	origin->point = new myPoint3D(p->X, p->Y, p->Z);

	myHalfedge *this_edge = f->adjacent_halfedge;
	myHalfedge *next_edge = this_edge->next;
	myHalfedge *prev_edge = this_edge->prev;

	myFace *t = new myFace(); t->index = faces.size(); faces.push_back(t);
	t->adjacent_halfedge = this_edge; this_edge->adjacent_face = t;

	myHalfedge *s = new myHalfedge(); s->source = origin;
	myHalfedge *z = new myHalfedge(); z->source = this_edge->source;
	myHalfedge *u = new myHalfedge(); u->source = next_edge->source;
	myHalfedge *v = new myHalfedge(); v->source = origin;

	s->index = halfedges.size(); halfedges.push_back(s);
	z->index = halfedges.size(); halfedges.push_back(z);
	u->index = halfedges.size(); halfedges.push_back(u);
	v->index = halfedges.size(); halfedges.push_back(v);

	s->twin = z; z->twin = s; u->twin = v; v->twin = u;
	z->next = v; v->prev = z; s->prev = u; u->next = s;

	z->adjacent_face = f; v->adjacent_face = f;
	s->adjacent_face = t; u->adjacent_face = t;
	
	this_edge->next = u; u->prev = this_edge;
	this_edge->prev = s; s->next = this_edge;

	z->prev = prev_edge; prev_edge->next = z;
	v->next = next_edge; next_edge->prev = v;

	f->adjacent_halfedge = v;
	triangulate(f);

	checkMesh();
}

bool myMesh::simplification() {

	if (faces.size() < 5) return false;
	int vertex_count = vertices.size();
	std::vector<int> rmEdge; 
	std::vector<int> rmVect; 
	std::vector<int> rmFace;

	int **index_map = new int*[vertex_count];
	for (int i = 0; i < vertex_count; ++i) {
		index_map[i] = new int[vertex_count];
		for (int j = 0; j < vertex_count; ++j) {
			index_map[i][j] = -1;
		}
	}

	int *vertex_mark = new int[vertex_count];
	for (int i = 0; i < vertex_count; ++i) {
		vertex_mark[i] = 0;
	}

	typedef pair<double, pair<int, int>> DPType;
	std::vector<DPType> distance_pair;
	int edge_count = halfedges.size();

	// tail is the source, head is the arrow
	for (int i = 0; i < edge_count; i++ ) {
		auto he = halfedges[i];
		int tail = he->source->index;
		int head = he->next->source->index;
		auto head_p = vertices[head]->point;
		auto tail_p = vertices[tail]->point;

		auto delta = *head_p-*tail_p;
		auto distance = delta.length();

		if (index_map[head][tail] != -1 && index_map[tail][head] != -1) { continue; }
		index_map[tail][head] = he->index; index_map[head][tail]= he->twin->index;
		auto i2d = make_pair(head, tail); auto d2p = make_pair(distance, i2d);

		distance_pair.push_back(d2p);
	}

	std::sort(distance_pair.begin(), distance_pair.end(), [](DPType a, DPType b) {
		return b.first > a.first;
	}); 

	int size_e = distance_pair.size();
	for (int i = 0; i < size_e; i++) {

		auto dp = distance_pair[i].second;
		auto h = dp.first, t = dp.second;
		auto edge_index = index_map[t][h];

		if (vertex_mark[h] == 1 || vertex_mark[t] == 1) {
			continue; // give up
		} else {
			vertex_mark[h] = 1;
			vertex_mark[t] = 1;
			// remove vertex
			rmVect.push_back(h);
		}

		auto edge = halfedges[edge_index];
		auto twin = edge->twin;

		// remove edge
		rmEdge.push_back(edge->index);
		rmEdge.push_back(twin->index);

		edge->adjacent_face->adjacent_halfedge = edge->prev;
		twin->adjacent_face->adjacent_halfedge = twin->next;

		auto prev = edge->prev; auto next = edge->next;
		prev->next = next; next->prev = prev;
		next->source = edge->source;

		auto shouldRemoveFace = (prev->next->next == prev);

		if (shouldRemoveFace) {
			//remove face
			rmEdge.push_back(prev->index);
			rmEdge.push_back(next->index);
			prev->twin->twin = prev->next->twin;
			prev->next->twin->twin = prev->twin;
			rmFace.push_back(prev->adjacent_face->index);
		}

		auto v_h = vertices[h]; auto v_t = vertices[t];
		auto v_n = (*v_h->point + *v_t->point) / 2;

		v_t->point->X = v_n.X;
		v_t->point->Y = v_n.Y;
		v_t->point->Z = v_n.Z;

		auto thisSample = shouldRemoveFace ? edge->next->twin->next : edge->next;
		v_t->originof = thisSample;

		auto sample = thisSample;
		do {
			sample->source = v_t;
			sample = sample->twin->next;
		} while (sample != thisSample);

		prev = twin->prev; next = twin->next;
		prev->next = next; next->prev = prev;
		
		shouldRemoveFace = (prev->next->next == prev);

		if (shouldRemoveFace) {
			//remove face
			rmEdge.push_back(prev->index);
			rmEdge.push_back(next->index);
			prev->twin->twin = prev->next->twin;
			prev->next->twin->twin = prev->twin;
			rmFace.push_back(prev->adjacent_face->index);
		}
	}

	delete[] vertex_mark;
	for (int i = 0; i < vertex_count; ++i) {
		delete[] index_map[i];
	}
	delete[] index_map;

	std::sort(rmVect.begin(), rmVect.end(), [](int a, int b) { return a > b; });
	for (auto v : rmVect) {
		auto dv = vertices[v];
		dv->index = -1;
		dv->originof = NULL;
		vertices.erase(vertices.begin()+v); delete dv;
		std::for_each(vertices.begin()+v, vertices.end(), [](auto myV) { myV->index -= 1;});
	}

	std::sort(rmEdge.begin(), rmEdge.end(), [](int a, int b) { return a > b; });
	for (auto e : rmEdge) {
		auto de = halfedges[e];
		de->index = -1;
		de->adjacent_face = NULL;
		de->source = NULL;
		de->next = NULL;
		de->prev = NULL;
		//auto he = halfedges[e];
		//if (he->source->originof == he) {
			//The edge will be removed, but the source vect still has a reference to this edge.
		//}
		halfedges.erase(halfedges.begin() + e); delete de;
		std::for_each(halfedges.begin()+e, halfedges.end(), [](auto myE) { myE->index -= 1; });
	}

	for (auto e : halfedges) { // reset any invalide reference.
		e->source->originof = e;
	}

	std::sort(rmFace.begin(), rmFace.end(), [](int a, int b) { return a > b; });
	for (auto f : rmFace) {
		faces.erase(faces.begin()+f);
		std::for_each(faces.begin()+f, faces.end(), [](auto myF) { myF->index -= 1; });
	}

	checkMesh();
	return true;
}