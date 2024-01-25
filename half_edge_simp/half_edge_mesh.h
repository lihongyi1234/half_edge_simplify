#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <assert.h>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <numeric>

namespace half_edge
{
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
#define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

	class SymetricMatrix {

	public:

		// Constructor

		SymetricMatrix(double c = 0) { loopi(0, 10) m[i] = c; }

		SymetricMatrix(double m11, double m12, double m13, double m14,
			double m22, double m23, double m24,
			double m33, double m34,
			double m44) {
			m[0] = m11;  m[1] = m12;  m[2] = m13;  m[3] = m14;
			m[4] = m22;  m[5] = m23;  m[6] = m24;
			m[7] = m33;  m[8] = m34;
			m[9] = m44;
		}

		// Make plane

		SymetricMatrix(double a, double b, double c, double d)
		{
			m[0] = a * a;  m[1] = a * b;  m[2] = a * c;  m[3] = a * d;
			m[4] = b * b;  m[5] = b * c;  m[6] = b * d;
			m[7] = c * c; m[8] = c * d;
			m[9] = d * d;
		}

		double operator[](int c) const { return m[c]; }

		// Determinant

		double det(int a11, int a12, int a13,
			int a21, int a22, int a23,
			int a31, int a32, int a33)
		{
			double det = m[a11] * m[a22] * m[a33] + m[a13] * m[a21] * m[a32] + m[a12] * m[a23] * m[a31]
				- m[a13] * m[a22] * m[a31] - m[a11] * m[a23] * m[a32] - m[a12] * m[a21] * m[a33];
			return det;
		}

		const SymetricMatrix operator+(const SymetricMatrix& n) const
		{
			return SymetricMatrix(m[0] + n[0], m[1] + n[1], m[2] + n[2], m[3] + n[3],
				m[4] + n[4], m[5] + n[5], m[6] + n[6],
				m[7] + n[7], m[8] + n[8],
				m[9] + n[9]);
		}

		SymetricMatrix& operator+=(const SymetricMatrix& n)
		{
			m[0] += n[0];   m[1] += n[1];   m[2] += n[2];   m[3] += n[3];
			m[4] += n[4];   m[5] += n[5];   m[6] += n[6];   m[7] += n[7];
			m[8] += n[8];   m[9] += n[9];
			return *this;
		}

		double m[10];
	};

	// Forward declare these
	typedef struct HEEdge HEEdge;
	typedef struct HEFace HEFace;

	//data structure of indexed face set
	struct Vertex {
		//3d coordinates
		int id;
		Eigen::Vector3d p;

		float c_v;
		float sum_areas;
		SymetricMatrix q= SymetricMatrix(0.0);
		int border = 0;

		bool operator==(Vertex v) const {
			return p[0] == v.p[0] && p[1] == v.p[1] && p[2] == v.p[2];
		}

		friend std::ostream& operator<<(std::ostream& stream, const Vertex& v) {
			return stream << "Vertex \n"
				<< " x : " << v.p[0] << "\n"
				<< " y : " << v.p[1] << "\n"
				<< " z : " << v.p[2] << "\n";
		}
	};

	struct Face {
		//three vertex ids
		unsigned int a, b, c;
		int id;
		Eigen::Vector3d n;
	};

	// Half-edge data structure
	struct HEVertex {
		//add members here
		Vertex vertex;
		HEEdge* edge; // One outgoing half-edge

		bool operator<(const HEVertex& hev) const {

			if (vertex.p[0] < hev.vertex.p[0]) return true;
			if (vertex.p[0] > hev.vertex.p[0]) return false;
			// Otherwise equal
			if (vertex.p[1] < hev.vertex.p[1]) return true;
			if (vertex.p[1] > hev.vertex.p[1]) return false;

			// Otherwise equal
			if (vertex.p[2] < hev.vertex.p[2]) return true;
			if (vertex.p[2] > hev.vertex.p[2]) return false;

			return false;

		}

		bool operator==(const HEVertex& other) const {

			return vertex.p[0] == other.vertex.p[0] && vertex.p[1] == other.vertex.p[1]
				&& vertex.p[2] == other.vertex.p[2];

		}

		bool operator!=(const HEVertex& other) const {

			return !(vertex.p[0] == other.vertex.p[0] && vertex.p[1] == other.vertex.p[1]
				&& vertex.p[2] == other.vertex.p[2]);

		}

		friend std::ostream& operator<<(std::ostream& stream, const HEVertex& v) {

			return stream << "\n" << v.vertex;

		}
	};


	struct HEEdge {
		//add members here

		HEFace* face;
		HEVertex* hevertex;
		HEEdge* opposite;
		HEEdge* next;
		int id;
		bool exists = true;
		float Q_error;
		float C_error;
		float Cost;
		Eigen::Vector3d collapse_point;
			
		//I(e) = cot(alpha1) + cot(alpha2)
		//more irregular, I_e more smaller.
		float I_e;  
		float quality0, quality1;
		//one edge consist of two half edge. 
		bool is_front = true;
		bool dirty = false;

		bool operator==(HEEdge e) const {
			return id == e.id;
		}

		bool operator!=(HEEdge e) const {
			return !(id == e.id);
		}

		bool operator<(const HEEdge& other) const {

			if (hevertex->vertex.p[0] < other.hevertex->vertex.p[0]) return true;
			if (hevertex->vertex.p[0] > other.hevertex->vertex.p[0]) return false;
			// Otherwise equal
			if (hevertex->vertex.p[1] < other.hevertex->vertex.p[1]) return true;
			if (hevertex->vertex.p[1] > other.hevertex->vertex.p[1]) return false;

			// Otherwise equal
			if (hevertex->vertex.p[2] < other.hevertex->vertex.p[2]) return true;
			if (hevertex->vertex.p[2] > other.hevertex->vertex.p[2]) return false;

			return false;

		}

	};

	struct HEFace {

		HEEdge* edge;
		int id;

		bool operator<(const HEFace& other) const {

			if (edge->hevertex->vertex.p[0] < other.edge->hevertex->vertex.p[0]) return true;
			if (edge->hevertex->vertex.p[0] > other.edge->hevertex->vertex.p[0]) return false;
			// Otherwise equal
			if (edge->hevertex->vertex.p[1] < other.edge->hevertex->vertex.p[1]) return true;
			if (edge->hevertex->vertex.p[1] > other.edge->hevertex->vertex.p[1]) return false;

			// Otherwise equal
			if (edge->hevertex->vertex.p[2] < other.edge->hevertex->vertex.p[2]) return true;
			if (edge->hevertex->vertex.p[2] > other.edge->hevertex->vertex.p[2]) return false;

			return false;
		}
	};

	class Mesh {
	public:
		Mesh() {};
		//Mesh(const char*);
		//load a Mesh from .mesh file
		//void loadMF(const char*);
		void loadVF(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
		//write a Mesh to .mesh file (no header)
		void writeMF(std::string);
		//simplify a mesh
		//void simplifyMesh(const char* input, const char* output, int faceCnt);
		void simplifyMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int target_count, std::string out);
		//turn indexed face set to halfedge
		void convertMesh();
		void calcBorder();
		bool isBorderVertex(HEVertex*);
		void calcFaceSymetricMatrix(HEFace*);
		void calcSymetricMatrix();
		void allGaussCurvature();
		void vertexGaussCurvature(HEVertex*);
		//turn halfedge to indexed face set
		void revertMesh();
		//collapse random
		void collapseVertices(HEVertex* v1, HEVertex* v2);
		void debugAdjacentEdges();
		bool canEdgeCollapse(HEEdge* edge);
		void collapseInnerEdge(HEEdge* edgeToRemove);
		int countCommonVertices(HEVertex* v1, HEVertex* v2);
		void cleanUp();
		//helper methods
		void printNeighborVertex(HEVertex*);
		std::vector<HEVertex*> getNeighbourVertices(HEVertex* v);
		std::vector<HEVertex*> getNeighbourVerticesV2(HEVertex* v);
		std::vector<HEEdge*> getNeighborFaceEdges(HEVertex* v);
		std::vector<HEVertex*> adjacentVertices(HEFace* f);

		Eigen::Vector3d calcFaceNormal(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);
		float calcTriArea(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);
		float calcAngle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

		double vertexQuadricError(SymetricMatrix q, double x, double y, double z);
		double edgeQuadricError(Vertex id_v1, Vertex id_v2, Eigen::Vector3d& p_result);
		void calcEdgeQuadricAndGaussError(HEEdge* edgeToSplit);
		void allEdgeQuadricGaussError();
		void calcQuadricGaussWeight();
		void allEdgeCost();
		double cot(Eigen::Vector3d v, Eigen::Vector3d w);

		void rankEdgeByCost(std::vector<HEEdge*>&);
		bool flipped(Eigen::Vector3d p, HEVertex* v0, HEVertex* v1);
		void updateQuadricGaussError(HEVertex*);
			
		//debug
		void writeVertexObj(std::string fn, HEVertex* hev);
		bool checkEdgeValid(HEEdge* edge);
		bool checkVertexValid(HEVertex* hev);
		bool checkVertexsNeighbor(HEVertex* hev0, HEVertex* hev1);
		int calcCommonNeighborCount(HEVertex* hev0, HEVertex* hev1);
		HEEdge* toBedeletedNeighborEdge(HEVertex* hev0, HEVertex* hev1);
		bool checkEdgeEarlyDeleted(HEEdge* edge);
		void revertEigenVF();

		Eigen::MatrixXd& getSimplifiedV();
		Eigen::MatrixXi& getSimplifiedF();


	private:
		std::vector<Vertex> m_V;
		std::vector<Face> m_F;
		std::set<HEVertex*> m_HEV;
		std::set<HEEdge*> m_HEE;
		std::set<HEFace*> m_HEF;

		float m_weight = -1;
		float m_w_coeff = 0.01;

		Eigen::MatrixXd m_simp_V;
		Eigen::MatrixXi m_simp_F;
	};

}