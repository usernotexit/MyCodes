#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

using namespace Ubpa;

using namespace std;

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>()){
	// TODO 2024/1/23
	/*
	cout << "Paramaterize::Paramaterize:" << endl
		<< "\t" << "not implemented" << endl;*/
	Init(triMesh);
}

void Paramaterize::Clear() {
	// TODO 2024/1/23
	/*
	cout << "Paramaterize::Clear:" << endl
		<< "\t" << "not implemented" << endl;*/
	heMesh->Clear();
	triMesh = nullptr;
}

bool Paramaterize::Init(Ptr<TriMesh> triMesh) {
	// TODO 2024/1/23
	/*
	cout << "Paramaterize::Init:" << endl
		<< "\t" << "not implemented" << endl;*/

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

bool Paramaterize::Run() {
	// TODO
	cout << "Paramaterize::Init:" << endl
		<< "\t" << "not implemented" << endl;

	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	auto texcoord = param();
	triMesh->Update(texcoord);

	return false;
}

vector<pointf2> Ubpa::Paramaterize::param()
{
	// TBD: kernal of the algorithm
	size_t nV = heMesh->NumVertices();
	vector<pointf2> coord;

	return coord;
}
