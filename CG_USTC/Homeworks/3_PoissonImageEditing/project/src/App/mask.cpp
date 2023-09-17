#include "mask.h"


bool ET::operator<(const ET& right)
{
	return this->x < right.x;
}

void ET::uplevel()
{
	x += k;
}

void SET::add(ET& e, int pos)
{
	if (pos >= h || pos < 0) return;
	ET* end = &(list[pos]);
	while (end && end->next) end = end->next;

	end->next = &e; e.next = nullptr; //if (end->next->next)printf("ERROR!\n");
}

void SET::add(QPoint p1, QPoint p2)
{
	if (p1.y() == p2.y())return;
	int ymin = p1.y() < p2.y() ? p1.y() : p2.y();
	ET* e = new ET(p1, p2);
	add(*e, ymin);
}

void SET::read(ET& acti_list, int pos)
{
	if (pos >= h || pos < 0) return;

	ET* e = &(list[pos]);
	ET* acti_tmp = &acti_list;

	while (acti_tmp->next)acti_tmp = acti_tmp->next;
	while (e->next)
	{
		e = e->next;
		acti_tmp->next = new ET[1];
		*(acti_tmp->next) = ET(*e);
		acti_tmp = acti_tmp->next;
	}
}

bool SET::isEmpty()
{
	bool isempty = true;
	for (int i = 0;i < h;i++) if (list[i].next) isempty = false;
	return isempty;
}

void Mask::add(QPoint pos_new)
{
	int x_lefttop = lefttop.x() < pos_new.x() ? lefttop.x() : pos_new.x();
	int y_lefttop = lefttop.y() < pos_new.y() ? lefttop.y() : pos_new.y();
	int x_rightbottom = rightbottom.x() > pos_new.x() ? rightbottom.x() : pos_new.x();
	int y_rightbottom = rightbottom.y() > pos_new.y() ? rightbottom.y() : pos_new.y();
	lefttop.setX(x_lefttop); lefttop.setY(y_lefttop);
	rightbottom.setX(x_rightbottom); rightbottom.setY(y_rightbottom);
	kpoints.push_back(pos_new);
}

void Mask::getET(SET& list)
{
	int h = rightbottom.y() - lefttop.y(), w = rightbottom.x() - lefttop.x();
	size_t n_kpoint = kpoints.size();

	list = SET(h);
	for (int k = 0; k < n_kpoint; k++)
	{
		int kk = (k == n_kpoint - 1) ? 0 : k + 1;
		auto p1 = kpoints[k] - lefttop, p2 = kpoints[kk] - lefttop;
		if (p1.x() <= 0) p1.setX(1);
		if (p1.x() >= w - 1) p1.setX(w - 2);
		if (p1.y() <= 0) p1.setY(1);
		if (p1.y() >= h - 1) p1.setY(h - 2);
		if (p2.x() <= 0) p2.setX(1);
		if (p2.x() >= w - 1) p2.setX(w - 2);
		if (p2.y() <= 0) p2.setY(1);
		if (p2.y() >= h - 1) p2.setY(h - 2);

		list.add(p1, p2);
	}
}

void Mask::getMask()
{
	SET list; getET(list);
	ET list_acti;

	auto rect = (rightbottom - lefttop);
	const int w = rect.x(), h = rect.y();
	map = std::vector<std::vector<int>>(w, std::vector<int>(h, 0));

	int id_inner = 0, id_bound = 0; // 为映射图内的边界点和内点编号			
	for (int i = 0; i < h; i++)
	{
		list.read(list_acti, i);
		for (ET* e1 = &list_acti; e1->next; e1 = e1->next)
		{
			// 将活动链表按 x值 大小排列（冒泡法）
			ET* e_min = e1;
			for (ET* e2 = e1; e2->next; e2 = e2->next)
				if (*(e2->next) < *(e_min->next)) e_min = e2;
			ET* e_tmp = e_min->next;
			e_min->next = e_tmp->next;
			e_tmp->next = e1->next;
			e1->next = e_tmp;
		}

		ET* tmp = list_acti.next;
		bool isinpoly = false;
		while (tmp && tmp->next)
		{
			int jnear = tmp->x, jfar = tmp->next->x;
			for (int j = jnear + 1; j < jfar; j++)
			{
				map[j][i] = ++id_inner;
			}

			tmp = tmp->next;
			isinpoly = !isinpoly;
		}

		for (tmp = &list_acti; tmp->next;)
		{
			tmp->next->uplevel(); //扫描交点上移
			if (tmp->next->ymax <= i || tmp->next->x >= w || tmp->next->x < 0)
			{	// 移除多余点
				ET* tmp_del = tmp->next;
				tmp->next = tmp_del->next;
			}
			else
			{
				tmp = tmp->next;
			}
		}
	}

	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (map[i][j] <= 0)continue;
			assert(i != 0 && i != w && j != 0 && j != h);
			map[i][j - 1] = (map[i][j - 1]) ? map[i][j - 1] : --id_bound;
			map[i - 1][j] = (map[i - 1][j]) ? map[i - 1][j] : --id_bound;
			map[i + 1][j] = (map[i + 1][j]) ? map[i + 1][j] : --id_bound;
			map[i][j + 1] = (map[i][j + 1]) ? map[i][j + 1] : --id_bound;
		}
	}
	n_bound = -id_bound;
	n_inner = id_inner;
}

void Mask::get_laplace_precomposed()
{
	if (n_inner == 0) { printf("no inner points here!\n"); return; }

	const int w = rightbottom.x() - lefttop.x(), h = rightbottom.y() - lefttop.y();
	std::vector<Eigen::Triplet<double>> ijv;
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			int id = map[i][j];
			if (id == 0) continue;
			if (id < 0)
			{
				id = -id - 1 + n_inner;
				ijv.push_back(Eigen::Triplet<double>(id, id, 1));
			}
			else
			{
				id -= 1;
				int id_[4] = { map[i][j - 1], map[i][j + 1], map[i - 1][j], map[i + 1][j] };
				for (int neibor = 0; neibor < 4; neibor++)
				{
					assert(id_[4] != 0);
					id_[neibor] = (id_[neibor] > 0) ? (id_[neibor] - 1) : (-id_[neibor] + n_inner - 1);
					ijv.push_back(Eigen::Triplet<double>(id, id_[neibor], -1));
				}
				ijv.push_back(Eigen::Triplet<double>(id, id, 4));
			}
		}
	}
	A.resize(n_inner + n_bound, n_inner + n_bound);
	A.setFromTriplets(ijv.begin(), ijv.end());

	qr.compute(A);
}

void Mask::get_boundary(Eigen::MatrixXd& boun, const QImage& src, QRect zone) const
{
	assert(zone.bottomRight() - zone.topLeft() == rightbottom - lefttop);
	const int w = rightbottom.x() - lefttop.x(), h = rightbottom.y() - lefttop.y();
	const int w_ = zone.left(), h_ = zone.top();

	boun.resize(n_bound, 3);
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			if (map[i][j] >= 0) continue;
			const QRgb& color = src.pixel(i + w_, j + h_);
			boun.row(-map[i][j] - 1) << qRed(color), qGreen(color), qBlue(color);
		}
	}
}

void Mask::get_boundaryAndInner(Eigen::MatrixXd& pix, const QImage& src, QRect zone) const
{// return (n_inner+n_bound)*3 matrix
	if (n_inner == 0) return;
	assert(zone.bottomRight() - zone.topLeft() == rightbottom - lefttop);
	const int w = rightbottom.x() - lefttop.x(), h = rightbottom.y() - lefttop.y();
	const int w_ = zone.left(), h_ = zone.top();

	pix.resize(n_inner + n_bound, 3);
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			if (map[i][j] == 0) continue;
			if (map[i][j] > 0)
			{
				const QRgb& color = src.pixel(i + w_, j + h_);
				pix.row(map[i][j] - 1) << qRed(color), qGreen(color), qBlue(color);
			}
			else
			{
				const QRgb& color = src.pixel(i + w_, j + h_);
				pix.row(-map[i][j] + n_inner - 1) << qRed(color), qGreen(color), qBlue(color);
			}
		}
}

void Mask::fill(QImage& dst, const Eigen::MatrixXd& color, QRect zone) const
{
	assert(zone.bottomRight() - zone.topLeft() == rightbottom - lefttop);
	const int w = rightbottom.x() - lefttop.x(), h = rightbottom.y() - lefttop.y();
	const int w_ = zone.left(), h_ = zone.top();

	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			if (map[i][j] <= 0) continue;
			const int id = map[i][j] - 1;
			int r(color(id, 0)), g(color(id, 1)), b(color(id, 2));
			r = r > 0 ? r : 0; r = r < 255 ? r : 255;
			g = g > 0 ? g : 0; g = g < 255 ? g : 255;
			b = b > 0 ? b : 0; b = b < 255 ? b : 255;
			dst.setPixel(i + w_, j + h_, qRgb(r, g, b));
		}
}

Mask& Mask::operator=(const Mask& rhs)
{ // 仅用于初始化，qr预分解模板会清空
	this->lefttop = rhs.lefttop;
	this->rightbottom = rhs.rightbottom;
	this->kpoints = rhs.kpoints;
	this->map = rhs.map;
	this->n_bound = rhs.n_bound;
	this->n_inner = rhs.n_inner;
	return *this;
}

void Mask::restore()
{
	n_bound = n_inner = 0;
	lefttop = QPoint(0, 0);
	rightbottom = QPoint(100000, 100000);
	kpoints.clear();
	map.clear();
}

void Mask::poission_editing(QImage& dst, const QImage& src, const QRect& zone_dst, const QRect& zone_src)
{
	if (n_inner == 0) return;
	if (qr.info() != Eigen::Success) { get_laplace_precomposed(); }

	// color
	Eigen::MatrixXd color, color_boundary;
	get_boundaryAndInner(color, src, zone_src);
	get_boundaryAndInner(color_boundary, dst, zone_dst);

	color = A * color;
	for (int i = 0; i < n_bound; i++)
		for (int c = 0; c < 3; c++)
			color(i + n_inner, c) = color_boundary(i + n_inner, c);

	color = qr.solve(color);
	fill(dst, color, zone_dst);
}

int Mask::getWidth()
{
	return rightbottom.x() - lefttop.x() + 1;
}

int Mask::getHeight()
{
	return rightbottom.y() - lefttop.y() + 1;
}

QPoint Mask::getLefttop() const
{
	return lefttop;
}

void Mask::Draw(QPainter& painter)
{
	const auto& points = kpoints.data();
	painter.drawPolygon(points, kpoints.size());
}
