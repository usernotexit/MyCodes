#pragma once
#include <QImage>
#include <QtGui>
#include <queue>
#include <Eigen/Sparse>

struct ET //ALE
{
	int ymax;
	double x;	// �¶˵��x����
	double k;	// б�ʵĵ���
	ET* next;
public:
	ET():ymax(-1), x(-1), k(0), next(nullptr) {};
	ET(int _y, double _x, double _k) : ymax(_y), x(_x), k(_k), next(nullptr) {};
	ET(const ET& e) :ymax(e.ymax), x(e.x), k(e.k), next(nullptr) {};
	ET(const QPoint& p1, const QPoint& p2)
	{
		ymax = p1.y() > p2.y() ? p1.y() : p2.y();
		x = p1.y() < p2.y() ? double(p1.x()) : double(p2.x());
		k = p1.y() == p2.y() ? DBL_MAX : (static_cast<double>(p1.x()) - p2.x()) / (p1.y() - p2.y());
		next = nullptr;
	}
	~ET() {}

	bool operator<(const ET& right);

	void uplevel();
};

struct SET // Sorted Edge Table
{
private:
	int h;
	std::vector<ET> list; // һ����ͷ��������
public:
	SET() :h(0) {};
	SET(int _h) : h(_h)
	{
		for (int i = 0; i < _h; i++) 
		{
			list.push_back(ET()); list[i].next = nullptr;
		}
	};
	~SET() {}// { for (int i = 0; i < h;i++)if (list[i])list[i]->~ET(); }

	void add(ET& e, int pos);

	void add(QPoint p1, QPoint p2);

	void read(ET& acti_list, int pos); // ����pos����ӵ�acti_list������ͷ��㣩��ֱ�Ӳ���

	bool isEmpty();
};

struct Mask
{
	Mask(QPoint pos)
		: lefttop(pos), rightbottom(pos)
	{
		kpoints.clear(); kpoints.push_back(pos);
	};
	Mask() :lefttop(0, 0), rightbottom(100000, 100000) {};
	~Mask() {};
	Mask& operator= (const Mask& rhs);

	void restore();

	void add(QPoint pos_new);

	void getET(SET& list);
	void getMask();
	void get_laplace_precomposed();

	void get_boundary(Eigen::MatrixXd& boun, const QImage& src, QRect zone) const;
	void get_boundaryAndInner(Eigen::MatrixXd& pix, const QImage& src, QRect zone) const;
	void fill(QImage& dst, const Eigen::MatrixXd& color, QRect zone) const;// fill the masked pixels (including boundaries)


	void poission_editing(QImage& dst, const QImage& src, const QRect& zone_dst, const QRect& zone_src);

	int getWidth();
	int getHeight();
	QPoint getLefttop() const;
	int getVertexNumber() { return kpoints.size(); }

	void Draw(QPainter&); // ����mask������
private:
	QPoint lefttop;
	QPoint rightbottom;
	std::vector<QPoint> kpoints;
	std::vector<std::vector<int>> map;
	int n_bound;
	int n_inner;
	Eigen::SparseMatrix<double> A; // laplace matrix
	Eigen::SparseLU < Eigen::SparseMatrix < double >, Eigen::AMDOrdering < int > > qr; // precomposed matrix
};