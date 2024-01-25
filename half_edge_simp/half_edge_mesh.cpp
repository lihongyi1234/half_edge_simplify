#pragma once
#include "half_edge_mesh.h"

#define DEBUG1

#define CV_PI   3.1415926535897932384626433832795

namespace half_edge {

	void Mesh::loadVF(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
	{
		if (m_V.size() > 0)m_V.clear();
		if (m_F.size() > 0)m_F.clear();
		m_V.resize(V.rows());
		loopi(0, V.rows()) {
			Vertex v;
			v.p = V.row(i);
			v.id = i;
			m_V[i] = v;
		}

		m_F.resize(F.rows());
		loopi(0, F.rows()) {
			Face f;
			f.a = F(i, 0); f.b = F(i, 1); f.c = F(i, 2);
			m_F[i] = f;
		}
	}
		
	void Mesh::writeMF(std::string filename) {
		std::ofstream outfile;
		outfile.open(filename, std::ios::out);
		std::string strbuff;
		for (unsigned int i = 0; i < m_V.size(); i++) {
			outfile << "v " << m_V[i].p[0] << " " << m_V[i].p[1] << " " << m_V[i].p[2] << std::endl;
		}
		for (unsigned int i = 0; i < m_F.size(); i++) {
			outfile << "f " << m_F[i].a << " " << m_F[i].b << " " << m_F[i].c << std::endl;
		}
		outfile.close();
	}

	std::vector<HEVertex*> Mesh::getNeighbourVertices(HEVertex* v) {

		std::vector<HEVertex*> tmp;

		HEEdge* curr = v->edge;
		HEEdge* first = v->edge;
		if (!first || !first->exists) {
			//std::cout << "Invalid edge, no neighbours" << std::endl;
			return tmp;
		}

		bool border = false;
		int count = 0;
		do {
			tmp.emplace_back(curr->hevertex);
			if (curr->opposite == nullptr) {
				border = true;
				break;
			}
			curr = curr->opposite->next;
		} while (*curr != *first);
		if (!border)
			return tmp;

		//v->edge is border
		curr = v->edge;
		do
		{
			tmp.emplace_back(curr->next->hevertex);
			curr = curr->next->next->opposite;
		} while (curr != nullptr);
		return tmp;
	}

	std::vector<HEVertex*> Mesh::getNeighbourVerticesV2(HEVertex* v) {

		std::vector<HEVertex*> tmp;

		HEEdge* curr = v->edge;
		HEEdge* first = v->edge;
		if (!first || !first->exists) {
			//std::cout << "Invalid edge, no neighbours" << std::endl;
			return tmp;
		}

		bool border = false;
		int count = 0;
		do {
			tmp.emplace_back(curr->hevertex);
			if (curr->opposite == nullptr) {
				border = true;
				break;
			}
			curr = curr->opposite->next;
		} while (*curr != *first);
		if (!border)
			return tmp;

		//v->edge is border
		curr = v->edge;
		do
		{
			tmp.emplace_back(curr->next->hevertex);
			curr = curr->next->next->opposite;
		} while (curr != nullptr);
		return tmp;
	}

	void Mesh::printNeighborVertex(HEVertex*v)
	{
		std::vector<HEVertex*> neighbor_vertices = getNeighbourVertices(v);
		printf("[Mesh][printNeighborVertex] vertex (%f,%f,%f) neigbor is:\n", v->vertex.p[0], v->vertex.p[1], v->vertex.p[2]);
		for (int i = 0; i < neighbor_vertices.size(); i++) {
			Vertex hev= neighbor_vertices[i]->vertex;
			printf("(%f,%f,%f)  ", hev.p[0], hev.p[1], hev.p[2]);
		}
		printf("\n");
	}

	int Mesh::countCommonVertices(HEVertex* v1, HEVertex* v2) {

		int count = 0;
		std::vector<HEVertex*> v1_neighbours = getNeighbourVertices(v1);
		std::vector<HEVertex*> v2_neighbours = getNeighbourVertices(v2);

		for (int i = 0; i < v1_neighbours.size(); i++) {

			for (int j = 0; j < v2_neighbours.size(); j++) {

				if (*v1_neighbours[i] == *v2_neighbours[j])
					count++;
			}
		}
		return count;
	}

	void Mesh::collapseVertices(HEVertex* v1, HEVertex* v2) {

		std::cout << "Deleted " << std::endl;
		//delete &v1;
		//v1 = v2;

	}

	bool Mesh::canEdgeCollapse(HEEdge* heedge) {

		if (heedge && heedge->exists)
			return true;
		return false;
	}

	void Mesh::debugAdjacentEdges() {

		for (auto v : m_HEV) {

			//HEEdge* curr = v->edge->opposite->next;
			HEEdge* curr = v->edge;
			HEEdge* first = v->edge;

			if (!first->exists)
				continue;

			do {

				curr = curr->opposite->next;

			} while (*curr != *first);
		}
	}

	void Mesh::collapseInnerEdge(HEEdge* edgeToRemove) {
		clock_t start = clock();
		HEEdge* edgeToRemoveTwin = edgeToRemove->opposite;
		HEVertex* hevertex1 = edgeToRemove->hevertex;
		HEVertex* hevertex2 = edgeToRemoveTwin->hevertex;

		// Get faces that share the end point of this edge
		// Get faces that borders this edge
		HEEdge* e_next = edgeToRemove->next;
		HEEdge* e_prev = edgeToRemove->next->next;
		HEEdge* e_next_twin = e_next->opposite;
		HEEdge* e_prev_twin = e_prev->opposite;

		HEEdge* et_next = edgeToRemoveTwin->next;
		HEEdge* et_prev = edgeToRemoveTwin->next->next;
		HEEdge* et_next_twin = et_next->opposite;
		HEEdge* et_prev_twin = et_prev->opposite;

		// Fix surrounding edges that points to vertex that's removed
		HEEdge* curr = hevertex1->edge->opposite;
		HEEdge* first = hevertex1->edge->opposite;

		// Reassign all incoming vertices from removed one to next
		int count = 0;
		do {
			curr->hevertex = hevertex2;
			curr = curr->next->opposite;
			if (count++ > 100) {
					
				bool ret = checkEdgeValid(first);

				getNeighbourVerticesV2(hevertex1);

				int count1 = 0;
				for (auto iter = m_HEE.begin(); iter != m_HEE.end(); iter++) {
					if ((*iter)->id == first->id)
						std::cout << " " << std::endl;
					if (!checkEdgeValid(*iter))
						count1++;
				}
				std::cout << "count:" << count1 << std::endl;

			}	
		} while (*curr != *first);

		if (*hevertex2->edge == *edgeToRemove || *hevertex2->edge == *et_next) {
			hevertex2->edge = et_prev->opposite;
			//std::cout << "** Reassignment of remaining vertex outedge" << std::endl;

			if (hevertex2 == hevertex2->edge->hevertex) {
				printf("[Mesh][collapse edge] error 000\n");
			}
		}

		HEVertex* left_hev = e_next->hevertex;
		if (*left_hev->edge == *e_prev) {
			left_hev->edge = e_next_twin;
			if (hevertex2 == hevertex2->edge->hevertex) {
				printf("[Mesh][collapse edge] error 111\n");
			}
		}
		HEVertex* right_hev = et_next->hevertex;
		if (*right_hev->edge == *et_prev) {
			right_hev->edge = et_next_twin;
			if (hevertex2 == hevertex2->edge->hevertex) {
				printf("[Mesh][collapse edge] error 222\n");
			}
		}

		// Bind together these
		e_next_twin->opposite = e_prev_twin;
		e_prev_twin->opposite = e_next_twin;

		et_next_twin->opposite = et_prev_twin;
		et_prev_twin->opposite = et_next_twin;
			
		if (hevertex2 == hevertex2->edge->hevertex) {
			printf("[Mesh][collapse edge] error 333\n");
		}

		//for (auto edge : m_HEE) {
		//	if (*edge->hevertex == *hevertex1) {
		//		printf("[Mesh][collapse edge] error *edge->hevertex == *hevertex1\n");
		//	}
		//}

		// Run some tests
		if (*et_prev_twin->opposite != *et_next_twin) {
			printf("[Mesh][collapse edge] error *et_prev_twin->opposite != *et_next_twin\n");
		}
		if(*et_next_twin->opposite != *et_prev_twin) {
			printf("[Mesh][collapse edge] error *et_next_twin->opposite != *et_prev_twin\n");
		}
		if (*e_next_twin->opposite != *e_prev_twin) {
			printf("[Mesh][collapse edge] error *e_next_twin->opposite != *e_prev_twin\n");
		}
		if (*e_prev_twin->opposite != *e_next_twin) {
			printf("[Mesh][collapse edge] error *e_prev_twin->opposite != *e_next_twin\n");
		}
		if (*hevertex2->edge == *edgeToRemove) {
			printf("[Mesh][collapse edge] error *hevertex2->edge == *edgeToRemove\n");
		}
		if (*hevertex2->edge == *et_next) {
			printf("[Mesh][collapse edge] error *hevertex2->edge == *et_next\n");
		}
		if (*e_next_twin->hevertex == *hevertex1) {
			printf("[Mesh][collapse edge] error *e_next_twin->hevertex == *hevertex1\n");
		}

		if (hevertex2 == hevertex2->edge->hevertex) {
			printf("[Mesh][collapse edge] error hevertex2 == hevertex2->edge->hevertex\n");
		}

		m_HEV.erase(hevertex1);

		m_HEF.erase(edgeToRemove->face);
		m_HEF.erase(edgeToRemoveTwin->face);

		m_HEE.erase(edgeToRemove);
		m_HEE.erase(e_next);
		m_HEE.erase(e_prev);
		m_HEE.erase(edgeToRemoveTwin);
		m_HEE.erase(et_next);
		m_HEE.erase(et_prev);

		edgeToRemove->exists = false;
		e_next->exists = false;
		e_prev->exists = false;
		edgeToRemoveTwin->exists = false;
		et_next->exists = false;
		et_prev->exists = false;

		//update 
		hevertex2->vertex.p = edgeToRemove->collapse_point;
		hevertex2->vertex.q = hevertex2->vertex.q + hevertex1->vertex.q;
		
		updateQuadricGaussError(hevertex2);
	}

	void Mesh::cleanUp() {

		// TODO. Leave it to garbage collector for now

	}

	void Mesh::simplifyMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int target_count, std::string out)
	{
		loadVF(V, F);
		convertMesh();
		calcBorder();

		calcSymetricMatrix();

		printf("[Mesh][simplifyMesh] simplify from %d to %d.\n", int(m_HEF.size()), target_count);
		std::vector<HEEdge*> edges_collapse_rank;

		double agressiveness = 4.0;
		int num_iterations = 100;

		for (int iter = 0; iter < num_iterations; iter++)
		{
			printf("################## Iter:%d, %d ##################\n", iter, m_HEF.size());
			int cur_count = m_HEF.size();
			if (cur_count <= target_count)break;
			edges_collapse_rank.clear();
			if (iter == 0) {
				allGaussCurvature();
				allEdgeQuadricGaussError();
				calcQuadricGaussWeight();
				//m_weight = 0.f;
				allEdgeCost();
			}

			double threshold = 0.000000001 * pow(double(iter + 3), agressiveness);
			int start_face_count = m_HEF.size();
			//collapse edge
			rankEdgeByCost(edges_collapse_rank);
			for (int i = 0; i < edges_collapse_rank.size(); i++) {
				HEEdge* edge_remove = edges_collapse_rank[i];
				if (!canEdgeCollapse(edge_remove))
					continue;
				if (edge_remove->Cost > threshold)
					continue;
				if (edge_remove->dirty)
					continue;
				Eigen::Vector3d p = edge_remove->collapse_point;
				HEVertex* hev0 = edge_remove->hevertex;
				HEVertex* hev1 = edge_remove->next->next->hevertex;

				if (flipped(p, hev0, hev1))
					continue;
				if (flipped(p, hev1, hev0))
					continue;

				if (!checkVertexsNeighbor(hev0, hev1))
					continue;

				int neighbor_count = calcCommonNeighborCount(hev0, hev1);
				if (neighbor_count > 3) {
					continue;
				}
				if (neighbor_count == 3) {
					HEEdge* early_deleted_edge = toBedeletedNeighborEdge(hev0, hev1);
					if (!checkEdgeEarlyDeleted(early_deleted_edge))
						continue;

					collapseInnerEdge(early_deleted_edge);

					if (!checkVertexsNeighbor(hev0, hev1)) {
						checkVertexsNeighbor(hev0, hev1);
					}

					if (calcCommonNeighborCount(hev0, hev1) > 2) {
						calcCommonNeighborCount(hev0, hev1);
					}
				}

				int border0 = hev0->vertex.border;
				int border1 = hev1->vertex.border;
				if (!border0 && !border1) {
					collapseInnerEdge(edge_remove);
				}
				else
					continue;
					
				if (m_HEF.size() <= target_count)
					break;
			}
		}
		
		revertMesh();
		revertEigenVF();

		writeMF(out);
		cleanUp();
	}

	/* Convert the mesh into the half-edge data structure*/
	void Mesh::convertMesh() {

		std::map<std::pair<unsigned int, unsigned int>, HEEdge*> edges; // Keep track of opposite edges
		std::map<std::pair<unsigned int, unsigned int>, int> edges_number;
		std::vector<HEVertex*> added_hevertices;
		int edge_index = 0;

		// Add vertices
		for (int i = 0; i < m_V.size(); i++) {
			HEVertex* _hevertex = new HEVertex();
			_hevertex->vertex = m_V[i];

			//HEVToIndex[_hevertex] = i;
			added_hevertices.emplace_back(_hevertex);
			m_HEV.insert(_hevertex);

		}

		// Go through all faces
		for (int i = 0; i < m_F.size(); i++) {
			// Declare the variables,
			// Assume that the vertices are declared in a clockwise order
			HEFace* _heface = new HEFace();

			HEEdge* _halfedge1 = new HEEdge();
			HEEdge* _halfedge2 = new HEEdge();
			HEEdge* _halfedge3 = new HEEdge();

			_halfedge1->exists = true;
			_halfedge2->exists = true;
			_halfedge3->exists = true;

			_halfedge1->is_front = true;
			_halfedge2->is_front = true;
			_halfedge3->is_front = true;

			HEVertex* _hevertex1 = added_hevertices[m_F[i].a];
			HEVertex* _hevertex2 = added_hevertices[m_F[i].b];
			HEVertex* _hevertex3 = added_hevertices[m_F[i].c];

			// HEFace
			_heface->edge = _halfedge1;

			// Update pointers from halfedge to vertex 
			_halfedge3->hevertex = _hevertex1;
			_hevertex1->edge = _halfedge1;

			//_hevertex2->edge = _halfedge2;
			_hevertex2->edge = _halfedge3;
			_halfedge2->hevertex = _hevertex2;

			//_hevertex3->edge = _halfedge3;
			_hevertex3->edge = _halfedge2;
			_halfedge1->hevertex = _hevertex3;

			// Next
			_halfedge1->next = _halfedge2;
			_halfedge2->next = _halfedge3;
			_halfedge3->next = _halfedge1;

			// Face
			_halfedge1->face = _heface;
			_halfedge2->face = _heface;
			_halfedge3->face = _heface;


			edges[std::make_pair(m_F[i].b, m_F[i].a)] = _halfedge3;
			edges[std::make_pair(m_F[i].c, m_F[i].b)] = _halfedge2;
			edges[std::make_pair(m_F[i].a, m_F[i].c)] = _halfedge1;

			// If opposite edge has been added
			if (edges.find(std::make_pair(m_F[i].a, m_F[i].b)) != edges.end()) {

				edges[std::make_pair(m_F[i].b, m_F[i].a)]->opposite = edges[std::make_pair(m_F[i].a, m_F[i].b)];
				edges[std::make_pair(m_F[i].a, m_F[i].b)]->opposite = edges[std::make_pair(m_F[i].b, m_F[i].a)];

			}

			if (edges.find(std::make_pair(m_F[i].b, m_F[i].c)) != edges.end()) {

				edges[std::make_pair(m_F[i].c, m_F[i].b)]->opposite = edges[std::make_pair(m_F[i].b, m_F[i].c)];
				edges[std::make_pair(m_F[i].b, m_F[i].c)]->opposite = edges[std::make_pair(m_F[i].c, m_F[i].b)];

			}

			if (edges.find(std::make_pair(m_F[i].c, m_F[i].a)) != edges.end()) {

				edges[std::make_pair(m_F[i].a, m_F[i].c)]->opposite = edges[std::make_pair(m_F[i].c, m_F[i].a)];
				edges[std::make_pair(m_F[i].c, m_F[i].a)]->opposite = edges[std::make_pair(m_F[i].a, m_F[i].c)];

			}

			if (edges_number.find(std::make_pair(m_F[i].b, m_F[i].a)) != edges_number.end()) {
				printf(" \n");
			}
			edges_number[std::make_pair(m_F[i].b, m_F[i].a)] = 1;
				
			if (edges_number.find(std::make_pair(m_F[i].c, m_F[i].b)) != edges_number.end()) {
				printf(" \n");
			}
			edges_number[std::make_pair(m_F[i].c, m_F[i].b)] = 1;

			if (edges_number.find(std::make_pair(m_F[i].a, m_F[i].c)) != edges_number.end()) {
				printf(" \n");
			}
			edges_number[std::make_pair(m_F[i].a, m_F[i].c)] = 1;



			//m_HEF.push_back(_heface);
			m_HEF.insert(_heface);


			_halfedge1->id = edge_index++;
			_halfedge2->id = edge_index++;
			_halfedge3->id = edge_index++;

			m_HEE.insert(_halfedge1);
			m_HEE.insert(_halfedge2);
			m_HEE.insert(_halfedge3);

			if (*_halfedge1->next != *_halfedge2) {
				printf("[Mesh][convertMesh] error *_halfedge1->next != *_halfedge2\n");
			}
			if (*_halfedge2->next != *_halfedge3) {
				printf("[Mesh][convertMesh] error *_halfedge2->next != *_halfedge3\n");
			}
			if (*_halfedge3->next != *_halfedge1) {
				printf("[Mesh][convertMesh] error *_halfedge3->next != *_halfedge1\n");
			}

			if (*_halfedge1->next->next != *_halfedge3) {
				printf("[Mesh][convertMesh] error *_halfedge1->next->next != *_halfedge3\n");
			}
			if (*_halfedge2->next->next != *_halfedge1) {
				printf("[Mesh][convertMesh] error *_halfedge2->next->next != *_halfedge1\n");
			}
			if (*_halfedge3->next->next != *_halfedge2) {
				printf("[Mesh][convertMesh] error *_halfedge3->next->next != *_halfedge2\n");
			}
		}

		for (auto edge : m_HEE) {
			if (*edge->opposite != *(*m_HEE.find(edge))->opposite) {
				printf("[Mesh][convertMesh] error *edge->opposite != *(*m_HEE.find(edge))->opposite\n");
			}
			//assert(*edge->opposite == *(*it)->opposite);
		}

		// Everything ok?
		if (m_V.size() != m_HEV.size()) {
			printf("[Mesh][convertMesh] m_V.size() != m_HEV.size()\n");
		}
		if (m_F.size() != m_HEF.size()) {
			printf("[Mesh][convertMesh] m_F.size() != m_HEF.size()\n");
		}
		if (m_HEE.size() != m_HEF.size() * 3)
		{
			printf("[Mesh][convertMesh] m_HEE.size() != m_HEF.size() * 3\n");
		}

		m_V.clear();
		m_F.clear();
	}

	bool Mesh::isBorderVertex(HEVertex* v)
	{
		HEEdge* curr = v->edge;
		HEEdge* first = v->edge;

		if (!first || !first->exists) {
			//std::cout << "Invalid edge, no neighbours" << std::endl;
			return true;
		}
		if (curr->opposite == nullptr)
			return true;
		do {
			if (curr->opposite == nullptr) {
				return true;
			}
			curr = curr->opposite->next;
		} while (*curr != *first);
			
		return false;
	}

	void Mesh::calcBorder()
	{			
		for (auto iter = m_HEV.begin(); iter != m_HEV.end(); iter++){
			int border = isBorderVertex(*iter);
			(*iter)->vertex.border = border;
		}
	}

	void Mesh::calcSymetricMatrix()
	{
		for (auto iter = m_HEF.begin(); iter != m_HEF.end(); iter++) {
			calcFaceSymetricMatrix(*iter);
		}
	}

	void Mesh::calcFaceSymetricMatrix(HEFace* hef)
	{
		HEEdge* edge0 = hef->edge; Vertex& vtx0 = edge0->hevertex->vertex;
		HEEdge* edge1 = edge0->next; Vertex& vtx1 = edge1->hevertex->vertex;
		HEEdge* edge2 = edge1->next; Vertex& vtx2 = edge2->hevertex->vertex;

		Eigen::Vector3d n = calcFaceNormal(vtx0.p, vtx2.p, vtx1.p);
		vtx0.q += SymetricMatrix(n[0], n[1], n[2], -n.dot(vtx0.p));
		vtx1.q += SymetricMatrix(n[0], n[1], n[2], -n.dot(vtx0.p));
		vtx2.q += SymetricMatrix(n[0], n[1], n[2], -n.dot(vtx0.p));
	}

	void Mesh::vertexGaussCurvature(HEVertex* hev)
	{
		int border = hev->vertex.border;
		std::vector<HEEdge*> hevs = getNeighborFaceEdges(hev);
		std::vector<float> angles;
		std::vector<float> areas;
		int num_face = hevs.size();
		loopi(0, num_face) {
			HEEdge* edge0 = hevs[i]; Vertex vtx0 = edge0->hevertex->vertex;
			HEEdge* edge1 = edge0->next; Vertex vtx1 = edge1->hevertex->vertex;
				
			HEEdge* edge2 = edge1->next; Vertex vtx2 = edge2->hevertex->vertex;
			assert(hev == vtx2);

			float angle = calcAngle(vtx2.p, vtx0.p, vtx1.p); angles.emplace_back(angle);
			float area = calcTriArea(vtx2.p, vtx0.p, vtx1.p); areas.emplace_back(area);
		}

		float sum_angles = accumulate(angles.begin(), angles.end(), 0.f);
		float sum_areas = accumulate(areas.begin(), areas.end(), 0.f);
		hev->vertex.sum_areas = sum_areas;
		float c_v;
		if (!border)
			c_v = (2 * CV_PI - sum_angles) / (1.0 / 3 * sum_areas);
		else
			c_v = (CV_PI - sum_angles) / (1.0 / 3 * sum_areas);
		hev->vertex.c_v = c_v;
	}

	void Mesh::allGaussCurvature()
	{
		for (auto iter = m_HEV.begin(); iter != m_HEV.end(); iter++) {
			vertexGaussCurvature(*iter);
		} 
	}

	void Mesh::revertMesh() {

		std::map<HEVertex, int> indices;

		int tempindex = 1;

		std::set<HEVertex*>::iterator it;

		for (it = m_HEV.begin(); it != m_HEV.end(); ++it) {

			HEVertex* he_vertex = *it;
			Vertex v;
			v = (*it)->vertex;
			m_V.emplace_back(v);
			indices[**it] = tempindex++;

		}

		assert(m_V.size() == m_HEV.size());

		for (auto heface : m_HEF) {

			Face f;

			HEVertex* _hevertex3 = heface->edge->hevertex;
			HEVertex* _hevertex2 = heface->edge->next->hevertex;
			HEVertex* _hevertex1 = heface->edge->next->next->hevertex;

			f.a = indices[*_hevertex1];
			f.b = indices[*_hevertex2];
			f.c = indices[*_hevertex3];

			m_F.emplace_back(f);
		}
		assert(m_HEF.size() == m_F.size());
	}

	std::vector<HEEdge*> Mesh::getNeighborFaceEdges(HEVertex* v) {
		std::vector<HEEdge*> tmp;

		HEEdge* curr = v->edge;
		HEEdge* first = v->edge;

		//tmp.emplace_back(v->edge);

		if (!first || !first->exists) {
			return tmp;
		}

		bool border = false;
		bool is_first = true;
		while (*curr != *first || is_first) {
			is_first = false;
			if (!curr->opposite) {
				border = true;
				break;
			}
			tmp.emplace_back(curr);
			curr = curr->opposite->next;
		}
		if (!border)return tmp;
		while (curr && curr->next->next->opposite) {
			curr = curr->next->next->opposite;
			tmp.emplace_back(curr);

		}
		return tmp;

	}

	std::vector<HEVertex*> Mesh::adjacentVertices(HEFace* f) {
		return std::vector<HEVertex*>();
	}

	Eigen::Vector3d Mesh::calcFaceNormal(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2)
	{
		Eigen::Vector3d a = p1 - p0;
		Eigen::Vector3d b = p2 - p0;

		Eigen::Vector3d n(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
		n.normalize();
		return n;
	}

	float Mesh::calcTriArea(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
	{
		float l0 = sqrtf((v1 - v0).squaredNorm());
		float l1 = sqrtf((v2 - v1).squaredNorm());
		float l2 = sqrtf((v0 - v2).squaredNorm());
		//check l0 l1 l2 is a triangle
		float p = (l0 + l1 + l2) / 2;
		float s = sqrtf(p * (p - l0) * (p - l1) * (p - l2));
		return s;
	}

	float Mesh::calcAngle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
	{
		Eigen::Vector3d v0_v1 = v1 - v0; v0_v1.normalize();
		Eigen::Vector3d v0_v2 = v2 - v0; v0_v2.normalize();
		return acosf(v0_v1.dot(v0_v2));
	}

	double Mesh::vertexQuadricError(SymetricMatrix q, double x, double y, double z)
	{
		return q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x + q[4] * y * y
			+ 2 * q[5] * y * z + 2 * q[6] * y + q[7] * z * z + 2 * q[8] * z + q[9];
	}

	// Error for one edge

	double Mesh::edgeQuadricError(Vertex v1, Vertex v2, Eigen::Vector3d& p_result)
	{
		// compute interpolated vertex

		SymetricMatrix q = v1.q + v2.q;
		bool border = v1.border & v2.border;
		double error = 0;
		double det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);
		if (det != 0 && !border)
		{

			// q_delta is invertible
			p_result[0] = -1 / det * (q.det(1, 2, 3, 4, 5, 6, 5, 7, 8));	// vx = A41/det(q_delta)
			p_result[1] = 1 / det * (q.det(0, 2, 3, 1, 5, 6, 2, 7, 8));	// vy = A42/det(q_delta)
			p_result[2] = -1 / det * (q.det(0, 1, 3, 1, 4, 6, 2, 5, 8));	// vz = A43/det(q_delta)

			error = vertexQuadricError(q, p_result[0], p_result[1], p_result[2]);
		}
		else
		{
			// det = 0 -> try to find best result
			//Eigen::Vector3d p1 = Eigen::Vector3d(v1.x, v1.y, v1.z);
			//Eigen::Vector3d p2 = Eigen::Vector3d(v2.x, v2.y, v2.z);
			Eigen::Vector3d p1 = v1.p;
			Eigen::Vector3d p2 = v2.p;

			Eigen::Vector3d p3 = (p1 + p2) / 2;
			double error1 = vertexQuadricError(q, p1[0], p1[1], p1[2]);
			double error2 = vertexQuadricError(q, p2[0], p2[1], p2[2]);
			double error3 = vertexQuadricError(q, p3[0], p3[1], p3[2]);
			error = std::min(error1, std::min(error2, error3));
			if (error1 == error) p_result = p1;
			if (error2 == error) p_result = p2;
			if (error3 == error) p_result = p3;
		}
		return error;
	}

	void Mesh::calcEdgeQuadricAndGaussError(HEEdge* edgeToSplit)
	{
		if (!edgeToSplit->is_front)
			return;
		HEVertex* hev0 = edgeToSplit->hevertex;
		Vertex vtx0 = hev0->vertex;
		HEVertex* hev1 = edgeToSplit->next->next->hevertex;
		Vertex vtx1 = hev1->vertex;

		float w0 = vtx0.sum_areas / (vtx0.sum_areas + vtx1.sum_areas);
		float w1 = vtx1.sum_areas / (vtx0.sum_areas + vtx1.sum_areas);

		edgeToSplit->C_error = fabsf(vtx0.c_v) * w0 + fabsf(vtx1.c_v) * w1;
		//Eigen::Vector3d d_point;
		edgeToSplit->Q_error = edgeQuadricError(vtx0, vtx1, edgeToSplit->collapse_point);

		edgeToSplit->Cost = edgeToSplit->C_error * m_weight + edgeToSplit->Q_error;

		if (edgeToSplit->opposite) {
			auto oppo = edgeToSplit->opposite;
			oppo->is_front = false;
			oppo->C_error = edgeToSplit->C_error;
			oppo->Q_error = edgeToSplit->Q_error;
			oppo->collapse_point = edgeToSplit->collapse_point;
			oppo->Cost = edgeToSplit->Cost;
		}
	}

	void Mesh::allEdgeQuadricGaussError()
	{
		for (auto iter = m_HEE.begin(); iter != m_HEE.end(); iter++) {
			calcEdgeQuadricAndGaussError(*iter);
		}
	}

	void Mesh::calcQuadricGaussWeight()
	{
		float Q_sum = 0;
		float C_sum = 0;
		int count = 0;
		for (auto iter = m_HEE.begin(); iter != m_HEE.end(); iter++) {
			if ((*iter)->is_front) {
				Q_sum += (*iter)->Q_error;
				C_sum += (*iter)->C_error;
				count++;
			}
		}
			
		float Q_avg = Q_sum / count;
		float C_avg = C_sum / count;
		m_weight = m_w_coeff * Q_avg / C_avg;
		printf("[Mesh][calcQuadricGaussWeight] Gauss and Quadric weight %f\n", m_weight);
	}

	void Mesh::allEdgeCost()
	{
		for (auto iter = m_HEE.begin(); iter != m_HEE.end(); iter++) {
			(*iter)->Cost = (*iter)->C_error * m_weight + (*iter)->Q_error;
		}
	}
		
	double Mesh::cot(Eigen::Vector3d v, Eigen::Vector3d w)
	{
		return v.dot(w) / (v.cross(w).norm());
	}

	void Mesh::rankEdgeByCost(std::vector<HEEdge*>& edges_rank)
	{
		for (auto iter = m_HEE.begin(); iter != m_HEE.end(); iter++) {
			(*iter)->dirty = false;
			if ((*iter)->is_front && (*iter)->exists) {
				edges_rank.emplace_back(*iter);
			}
		}
		sort(edges_rank.begin(), edges_rank.end(), [](HEEdge* e0, HEEdge* e1) {return e0->Cost < e1->Cost; });
	}

	bool Mesh::flipped(Eigen::Vector3d p, HEVertex* v0, HEVertex* v1)
	{
		std::vector<HEEdge*> heedges = getNeighborFaceEdges(v0);
		int length = heedges.size();
		for (size_t i = 0; i < length; i++) {
			HEEdge* edge = heedges[i];
			HEVertex* hev0 = edge->hevertex;
			HEVertex* hev1 = edge->next->hevertex;
			HEVertex* hev2 = edge->next->next->hevertex;
			assert(hev2 == v0);
			if (hev0 == v1 || hev1 == v0)
				continue;
			Eigen::Vector3d n0 = calcFaceNormal(hev2->vertex.p, hev1->vertex.p, hev0->vertex.p);

			Eigen::Vector3d d1 = hev0->vertex.p - p; d1.normalize();
			Eigen::Vector3d d2 = hev1->vertex.p - p; d2.normalize();
			if (fabs(d1.dot(d2)) > 0.999)
				return true;
			Eigen::Vector3d n1 = calcFaceNormal(p, hev1->vertex.p, hev0->vertex.p);
			if (n0.dot(n1) < 0.2)
				return true;
		}
		return false;
	}

	void Mesh::updateQuadricGaussError(HEVertex* hevertex)
	{
		vertexGaussCurvature(hevertex);
		std::vector<HEVertex*> hevs = getNeighbourVertices(hevertex);
		for (int i = 0; i < hevs.size(); i++) {
			HEVertex* hev = hevs[i];
			vertexGaussCurvature(hev);
		}

		for (int i = 0; i < hevs.size(); i++) {
			HEVertex* hev = hevs[i];
			std::vector<HEEdge*> neighbor_hees = getNeighborFaceEdges(hev);

			for (int j = 0; j < neighbor_hees.size(); j++) {
				HEEdge* hee = neighbor_hees[j];
				calcEdgeQuadricAndGaussError(hee);
				hee->dirty = true;
				if (hee->opposite)
					hee->opposite->dirty = true;
			}
		}
	}

	void Mesh::writeVertexObj(std::string fn, HEVertex* hev)
	{
		Eigen::MatrixXd V(1, 3);
		V.row(0) = hev->vertex.p;

		std::ofstream ofs(fn, std::ios::out);
		for (int i = 0; i < V.rows(); i++) {
			Eigen::Vector3d v = V.row(i);
			std::string str = "v " + std::to_string(v[0]) + " " + std::to_string(v[1]) + " " + std::to_string(v[2]);
			ofs << str << std::endl;
		}
		ofs.close();
	}

	bool Mesh::checkEdgeValid(HEEdge* edge)
	{
		HEEdge* curr = edge;
		int count = 0;
		bool valid = true;
		do {
			curr = curr->next->opposite;
			if (count++ > 100) {
				valid = false;
				break;
			}
		} while (*curr != *edge);
		return valid;
	}

	bool Mesh::checkVertexValid(HEVertex* hev)
	{
		std::vector<HEVertex*> hevs = getNeighbourVertices(hev);

		std::vector<int> vecInt1;
		for (size_t i = 0; i < hevs.size(); i++)
			vecInt1.emplace_back(hevs[i]->vertex.id);

		for (int i = 0; i < hevs.size(); i++) {
			std::vector<HEVertex*> hev2s = getNeighbourVertices(hevs[i]);
			std::vector<int> vecInt2;
			for (size_t i = 0; i < hev2s.size(); i++) {
				vecInt2.emplace_back(hev2s[i]->vertex.id);
			}

			int count = 0;
			for (int i = 0; i < vecInt1.size(); i++)
			{
				if (find(vecInt2.begin(), vecInt2.end(), vecInt1[i]) != vecInt2.end())
					count++;
			}
			if (count != 2) {
				std::cout << "Error. " << std::endl;
				return false;
			}

		}
		return true;
	}

	bool Mesh::checkVertexsNeighbor(HEVertex* hev0, HEVertex* hev1)
	{
		std::vector<HEVertex*> hev0_neighbors = getNeighbourVertices(hev0);
		std::vector<HEVertex*> hev1_neighbors = getNeighbourVertices(hev1);
		if (hev0_neighbors.size() < 3 || hev1_neighbors.size() < 3)
			return false;
		return true;
	}

	int Mesh::calcCommonNeighborCount(HEVertex* hev0, HEVertex* hev1)
	{
		std::vector<HEVertex*> hev0_neighbors = getNeighbourVertices(hev0);
		std::vector<HEVertex*> hev1_neighbors = getNeighbourVertices(hev1);
		int count = 0;
		for (size_t i = 0; i < hev0_neighbors.size(); i++) {
			HEVertex* hev0_i = hev0_neighbors[i];
			for (size_t j = 0; j < hev1_neighbors.size(); j++)
			{
				HEVertex* hev1_j = hev1_neighbors[j];
				if (hev0_i == hev1_j) {
					count++;
					break;
				}
			}
		}
		return count;
	}

	HEEdge* Mesh::toBedeletedNeighborEdge(HEVertex* hev0, HEVertex* hev1)
	{
		std::vector<HEEdge*> hev0_neighbor_edges = getNeighborFaceEdges(hev0);
		std::vector<HEVertex*> hev1_neighbors_vtxs = getNeighbourVertices(hev1);
		HEEdge* edge = NULL;
		for (size_t i = 0; i < hev0_neighbor_edges.size(); i++) {
			HEVertex* hev0_i0 = hev0_neighbor_edges[i]->hevertex;
			HEVertex* hev0_i1 = hev0_neighbor_edges[i]->next->hevertex;
			bool finded = false;
			for (size_t j = 0; j < hev1_neighbors_vtxs.size(); j++)
			{
				HEVertex* hev1_j = hev1_neighbors_vtxs[j];
				if (hev0_i0 == hev1_j) {
					finded = true;
					break;
				}
			}
			if (!finded)continue;
			finded = false;
			for (size_t j = 0; j < hev1_neighbors_vtxs.size(); j++)
			{
				HEVertex* hev1_j = hev1_neighbors_vtxs[j];
				if (hev0_i1 == hev1_j) {
					finded = true;
					break;
				}
			}
			if (!finded)continue;
			edge = hev0_neighbor_edges[i]->next;
			if (!edge->is_front)
				edge = hev0_neighbor_edges[i]->next->opposite;
			break;
		}
		return edge;
	}


	bool Mesh::checkEdgeEarlyDeleted(HEEdge* edge_remove)
	{
		if (!canEdgeCollapse(edge_remove))
			return false;

		if (edge_remove->dirty)
			return false;
		Eigen::Vector3d p = edge_remove->collapse_point;
		HEVertex* hev0 = edge_remove->hevertex;
		HEVertex* hev1 = edge_remove->next->next->hevertex;

		if (flipped(p, hev0, hev1))
			return false;
		if (flipped(p, hev1, hev0))
			return false;

		if (!checkVertexsNeighbor(hev0, hev1))
			return false;

		if (calcCommonNeighborCount(hev0, hev1) > 2)
			return false;

		int border0 = hev0->vertex.border;
		int border1 = hev1->vertex.border;
		if (border0 || border1)
			return false;
		return true;
	}

	void Mesh::revertEigenVF()
	{
		m_simp_V.resize(m_V.size(), 3);
		m_simp_F.resize(m_F.size(), 3);

		for (unsigned int i = 0; i < m_V.size(); i++) {
			//outfile << "v " << m_V[i].p[0] << " " << m_V[i].p[1] << " " << m_V[i].p[2] << std::endl;
			m_simp_V(i, 0) = m_V[i].p[0];
			m_simp_V(i, 1) = m_V[i].p[1];
			m_simp_V(i, 2) = m_V[i].p[2];
		}
		for (unsigned int i = 0; i < m_F.size(); i++) {
			//outfile << "f " << m_F[i].a << " " << m_F[i].b << " " << m_F[i].c << std::endl;
			m_simp_F(i, 0) = m_F[i].a - 1;
			m_simp_F(i, 1) = m_F[i].b - 1;
			m_simp_F(i, 2) = m_F[i].c - 1;
		}
		//outfile.close();
	}

	Eigen::MatrixXd& Mesh::getSimplifiedV()
	{
		return m_simp_V;
	}

	Eigen::MatrixXi& Mesh::getSimplifiedF()
	{
		return m_simp_F;
	}

}