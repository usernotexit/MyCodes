#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

MinSurf::MinSurf(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMesh);
}

void MinSurf::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool MinSurf::Init(Ptr<TriMesh> triMesh) {
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;// 三角形顶点序号
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) { // 边界检测 已完成
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}

bool MinSurf::Run() {
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Minimize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	triMesh->Init(indice, positions);

	return true;
}

void MinSurf::Minimize() {
	// TODO
	/*
	cout << "WARNING::MinSurf::Minimize:" << endl
		<< "\t" << "not implemented" << endl;
	*/
	// 草，heMesh的意思是“half edge Mesh”
	// half edge mesh is ready
	// 两层vector，第一层表示每个亏格，第二层表示亏格中的每条半边
	size_t nV = heMesh->NumVertices();

	MatrixX3f B(nV, 3); B.setZero();
	vector<Triplet<float>> ijv;
	for (auto& b : heMesh->Boundaries())
		for (auto& he : b)
		{
			size_t i = heMesh->Index(he->End());
			cout << i << endl;
			for (size_t k = 0; k < 3; k++)
				B(i, k) = he->End()->pos[k];
			ijv.push_back(Triplet<float>(i, i, 1));
		}

	for (auto& v : heMesh->Vertices())
	{
		if (v->IsBoundary()) continue;
		else
		{
			auto i = heMesh->Index(v);
			ijv.push_back(Triplet<float>(i, i, 1));
			for (auto& vv : v->AdjVertices())
			{
				auto j = heMesh->Index(vv);
				ijv.push_back(Triplet<float>(i, j, -1./v->Degree()));
			}
		}
	}
	
	SparseMatrix<float, ColMajor> A(nV, nV);
	A.setFromTriplets(ijv.cbegin(), ijv.cend());

	SparseLU<SparseMatrix<float>> solver(A);
	auto x = solver.solve(B);

	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = vecf3(x(i, 0), x(i, 1), x(i, 2));
	}

}
