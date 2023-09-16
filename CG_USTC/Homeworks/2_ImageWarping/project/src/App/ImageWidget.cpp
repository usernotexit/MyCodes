#include "ImageWidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>
#include <assert.h>

#include <Eigen/Dense>
#include <ANN/ANN.h>

using std::cout;
using std::endl;

ImageWidget::ImageWidget(void)
{
	ptr_image_ = new QImage();
	ptr_image_backup_ = new QImage();
	warpingActivate = false;
	angle = 0.;
	mode = kRBF_ann; // default warping mode
}


ImageWidget::~ImageWidget(void)
{
}

// 一些小功能
bool InRect(const QPoint& pos, const QRect& rec)
{
	return true;
}

void ImageWidget::paintEvent(QPaintEvent *paintevent)
{
	QPainter painter;
	painter.begin(this);

	// Draw background
	painter.setBrush(Qt::lightGray);
	QRect back_rect(0, 0, width(), height());
	painter.drawRect(back_rect);

	// Draw image
	QRect rect = QRect( (width()-ptr_image_->width())/2, (height()-ptr_image_->height())/2, ptr_image_->width(), ptr_image_->height());
	painter.drawImage(rect, *ptr_image_); 

	// Draw control points
	const int n_points = psrc.size();
	for (int i = 0; i < n_points; i++) 
	{
		QPoint topleft((width() - ptr_image_->width()) / 2, (height() - ptr_image_->height()) / 2);
		painter.setPen(QPen(Qt::blue, 2));
		painter.drawLine(psrc[i] + topleft, pdst[i] + topleft);
		painter.setPen(QPen(Qt::red, 4));
		painter.drawPoint(psrc[i] + topleft);
		painter.setPen(QPen());
	}
	painter.end();
}

void ImageWidget::addWarpingPoint(QPoint psrc_new, QPoint pdst_new)
{
	psrc.push_back(psrc_new);
	pdst.push_back(pdst_new);
	assert(psrc.size() == pdst.size());
}

void ImageWidget::delWarpingPoint()
{
	if (!psrc.empty())
	{
		psrc.pop_back();
		pdst.pop_back();
	}
}

void ImageWidget::mousePressEvent(QMouseEvent* event)
{
	if (Qt::LeftButton == event->button())
	{
		if (warpingActivate)
		{
			auto pos_img = event->pos() - QPoint(width()/2, height()/2) + QPoint(ptr_image_->width() / 2, ptr_image_->height() / 2);
			psrc.push_back(pos_img);
			pdst.push_back(pos_img);
		}
	}
}

void ImageWidget::mouseMoveEvent(QMouseEvent* event)
{
	if (warpingActivate)
	{
		assert(pdst.size() > 0);
		pdst.back() = event->pos() - QPoint(width() / 2, height() / 2) + QPoint(ptr_image_->width()/2, ptr_image_->height()/2);
		update();
	}
}

void ImageWidget::mouseReleaseEvent(QMouseEvent* event)
{
	if (warpingActivate)
	{
		pdst.back() = event->pos() - QPoint(width() / 2, height() / 2) + QPoint(ptr_image_->width() / 2, ptr_image_->height() / 2);
		assert(psrc.size() == pdst.size());
		Warp();
		update();
	}
}

void ImageWidget::Open()
{
	// Open file
	QString fileName = QFileDialog::getOpenFileName(this, tr("Read Image"), ".", tr("Images(*.bmp *.png *.jpg)"));

	// Load file
	if (!fileName.isEmpty())
	{
		ptr_image_->load(fileName);
		*(ptr_image_backup_) = *(ptr_image_);
	}

	//ptr_image_->invertPixels(QImage::InvertRgb);
	//*(ptr_image_) = ptr_image_->mirrored(true, true);
	//*(ptr_image_) = ptr_image_->rgbSwapped();
	cout<<"image size: "<<ptr_image_->width()<<' '<<ptr_image_->height()<<endl;
	update();
}

void ImageWidget::Save()
{
	SaveAs();
}

void ImageWidget::SaveAs()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), ".", tr("Images(*.bmp *.png *.jpg)"));
	if (filename.isNull())
	{
		return;
	}	

	ptr_image_->save(filename);
}

void ImageWidget::Invert()
{
	for (int i=0; i<ptr_image_->width(); i++)
	{
		for (int j=0; j<ptr_image_->height(); j++)
		{
			QRgb color = ptr_image_->pixel(i, j);
			ptr_image_->setPixel(i, j, qRgb(255-qRed(color), 255-qGreen(color), 255-qBlue(color)) );
		}
	}

	// equivalent member function of class QImage
	// ptr_image_->invertPixels(QImage::InvertRgb);
	update();
}

void ImageWidget::Mirror(bool ishorizontal, bool isvertical)
{
	QImage image_tmp(*(ptr_image_));
	int width = ptr_image_->width();
	int height = ptr_image_->height();

	if (ishorizontal)
	{
		if (isvertical)
		{
			for (int i=0; i<width; i++)
			{
				for (int j=0; j<height; j++)
				{
					ptr_image_->setPixel(i, j, image_tmp.pixel(width-1-i, height-1-j));
				}
			}
		} 
		else			//仅水平翻转			
		{
			for (int i=0; i<width; i++)
			{
				for (int j=0; j<height; j++)
				{
					ptr_image_->setPixel(i, j, image_tmp.pixel(width-1-i, j));
				}
			}
		}
		
	}
	else
	{
		if (isvertical)		//仅垂直翻转
		{
			for (int i=0; i<width; i++)
			{
				for (int j=0; j<height; j++)
				{
					ptr_image_->setPixel(i, j, image_tmp.pixel(i, height-1-j));
				}
			}
		}
	}
	// equivalent member function of class QImage
	//*(ptr_image_) = ptr_image_->mirrored(true, true);
	update();
}

void ImageWidget::TurnGray()
{
	for (int i=0; i<ptr_image_->width(); i++)
	{
		for (int j=0; j<ptr_image_->height(); j++)
		{
			QRgb color = ptr_image_->pixel(i, j);
			int gray_value = (qRed(color)+qGreen(color)+qBlue(color))/3;
			ptr_image_->setPixel(i, j, qRgb(gray_value, gray_value, gray_value) );
		}
	}

	update();
}

void ImageWidget::Restore()
{
	*(ptr_image_) = *(ptr_image_backup_);
	warpingActivate = false;
	psrc.clear(); pdst.clear();
	angle = 0.;
	update();
}

void ImageWidget::setWarpingActivate()
{
	warpingActivate = !warpingActivate;
}

/* Image warping algorthms */
auto RBF = [](double px, double py, double qx, double qy) {
	const double e = 2.7182818284, r_2 = 2500;
	double sigma_2 = (px - qx) * (px - qx) + (py - qy) * (py - qy);//pow(abs(px - qx), 2) + pow(abs(py - qy), 2);
	return pow(e, -sigma_2/ r_2);
};

auto SIGMA = [](double px, double py, double qx, double qy, double mu) {
	double sigma_2 = (px - qx) * (px - qx) + (py - qy) * (py - qy);
	return pow(sigma_2, -mu/2);
};

void ImageWidget::Warp()
{
	switch (mode)
	{
	case ImageWidget::kRBF_ann:
		Warp_RBF_ann();
		break;
	case ImageWidget::kRBF_std:
		Warp_RBF_std();
		break;
	case ImageWidget::kIDW_ann:
		Warp_IDW();
		break;
	default:
		cout << "Error!" << endl;
		break;
	}
}

void ImageWidget::setRBFstd()
{
	mode = kRBF_std;
	*(ptr_image_) = *(ptr_image_backup_);
	Warp();
}

void ImageWidget::setRBFann()
{
	mode = kRBF_ann;
	Warp();
}

void ImageWidget::setIDWann()
{
	mode = kIDW_ann;
	Warp();
}

void ImageWidget::Warp_RBF_ann()
{
	/// RBF algorithm: f(x) = Ax+b+\sum_j [a_j R(||x-p_j||)] 
	// 计算RBF权重
	int n_psrc = psrc.size();
	Eigen::MatrixXd R(n_psrc, n_psrc);
	Eigen::MatrixXd A(2,2); Eigen::Vector2d b;
	Eigen::MatrixXd Q(n_psrc,2), P(n_psrc,2), B(n_psrc,2);

	A.setIdentity();
	b.setZero();
	B.setZero();

	for (int i = 0; i < n_psrc; i++)
		for (int j = 0;j < n_psrc;j++)
			R(i,j) = RBF(psrc[i].x(), psrc[i].y(), psrc[j].x(), psrc[j].y());	

	for (int i = 0;i < n_psrc;i++)
	{
		P(i, 0) = psrc[i].x(); P(i, 1) = psrc[i].y();
		Q(i, 0) = pdst[i].x(), Q(i, 1) = pdst[i].y();
	}
	Eigen::MatrixXd alpha = R.colPivHouseholderQr().solve(Q - P);

	// 计算各点新坐标
	const unsigned int w = ptr_image_backup_->width(), h = ptr_image_backup_->height();
	Eigen::MatrixXd coords_old(w * h,2), B_(w * h,2);
	Eigen::MatrixXd R_(w * h, n_psrc);

	for (int i = 0;i < w;i++)
		for (int j = 0;j < h;j++)
		{
			coords_old.row(i * h + j) = Eigen::Vector2d(i, j); // ?
			for (int k = 0; k < n_psrc; k++)
				R_(i * h + j, k) = RBF(psrc[k].x(), psrc[k].y(), i, j);
		}
	B_.setZero();
	Eigen::MatrixXd coords_new = R_ * alpha + coords_old * A; coords_new += B_;

	// 用ann算法填补空洞
	const size_t n_point = w*h;
	constexpr size_t dim = 2;
	constexpr size_t K = 4;//取最近邻四个点求平均，目前暂未启用

	ANNpointArray pts_Arr = annAllocPts(n_point, dim);
	for (size_t i = 0; i < n_point; i++)
		for (size_t j = 0; j < dim; j++)
			pts_Arr[i][j] = coords_new(i, j);
	/* 此例可行，问题竟然在值的分配？
	const size_t ww = n_point;
	ANNpointArray pts_arr2 = annAllocPts(ww, dim);
	for (size_t i = 0;i < ww;i++)
		for (size_t j = 0;j < dim;j++)
			pts_arr2[i][j] = 1.3 + i;
	ANNbd_tree tree2(pts_arr2, ww, dim);
	printf("%d", sizeof(tree2)); */
	
	// 似乎在初始化tree时出现堆栈溢出的问题
	// 2023/8/9 更新：问题已解决，原因应该是点集过密导致沿x轴y轴划分时无限递归
	ANNbd_tree tree(pts_Arr, n_point, dim, w, ANN_KD_SUGGEST, ANN_BD_SUGGEST);//, 4, ANN_KD_STD, ANN_BD_SUGGEST);

	ANNpoint queryPt = annAllocPt(dim);
	ANNidx idxArr[K];
	ANNdist distArr[K];

	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			queryPt[0] = i; queryPt[1] = j;
			tree.annkSearch(queryPt, K, idxArr, distArr);
			
			int inv_i = idxArr[0] / h, inv_j = idxArr[0] % h;//求余、整除

			QColor color;
			color = ptr_image_backup_->pixel(QPoint(inv_i, inv_j));
			ptr_image_->setPixelColor(QPoint(i, j), color);
		}
	}
	update();
}

void ImageWidget::Warp_RBF_std()
{
	/// RBF algorithm: f(x) = Ax+b+\sum_j [a_j R(||x-p_j||)] 
	// 计算RBF权重
	int n_psrc = psrc.size();
	Eigen::MatrixXd R(n_psrc, n_psrc);
	Eigen::MatrixXd A(2, 2); Eigen::Vector2d b;
	Eigen::MatrixXd Q(n_psrc, 2), P(n_psrc, 2), B(n_psrc, 2);

	A.setIdentity();
	b.setZero();
	B.setZero();

	for (int i = 0; i < n_psrc; i++)
		for (int j = 0;j < n_psrc;j++)
			R(i, j) = RBF(psrc[i].x(), psrc[i].y(), psrc[j].x(), psrc[j].y());

	for (int i = 0;i < n_psrc;i++)
	{
		P(i, 0) = psrc[i].x(); P(i, 1) = psrc[i].y();
		Q(i, 0) = pdst[i].x(), Q(i, 1) = pdst[i].y();
	}
	Eigen::MatrixXd alpha = R.colPivHouseholderQr().solve(Q - P);

	// 计算各点新坐标
	const unsigned int w = ptr_image_backup_->width(), h = ptr_image_backup_->height();
	Eigen::MatrixXd coords_old(w * h, 2), B_(w * h, 2);
	Eigen::MatrixXd R_(w * h, n_psrc);

	for (int i = 0;i < w;i++)
		for (int j = 0;j < h;j++)
		{
			coords_old.row(i * h + j) = Eigen::Vector2d(i, j); // ?
			for (int k = 0; k < n_psrc; k++)
				R_(i * h + j, k) = RBF(psrc[k].x(), psrc[k].y(), i, j);
		}
	B_.setZero();
	Eigen::MatrixXd coords_new = R_ * alpha + coords_old * A; coords_new += B_;
	Eigen::MatrixXd tmp = R_ * alpha;

#ifdef _DEBUG
	for (int i = 0; i < 20; i++)
	{
		std::cout << R_.row(i) << std::endl;
	}
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	for (int i = 0; i < 20; i++)
	{
		std::cout << tmp.row(i) << std::endl;
	}
#endif // DEBUG


	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			int i_new = round(coords_new(i * h + j, 0)), j_new = round(coords_new(i * h + j, 1));
			QColor color = ptr_image_backup_->pixel(QPoint(i, j));

			ptr_image_->setPixelColor(QPoint(i_new, j_new), color);
		}
	}
	update();
}


void ImageWidget::Warp_IDW()
{
	/// IDW algorithm: f(x) = \sum_j w_j f_j(x) 
	/// f_j(p) = q_j + D_j(p - p_j)
	int n_psrc = psrc.size();
	Eigen::MatrixXd D(2, 2); Eigen::Vector2d b;
	Eigen::MatrixXd Q(n_psrc, 2), P(n_psrc, 2);
	Eigen::MatrixXd W(1, n_psrc);

	// 计算idw权重
	D.setIdentity();// D_j = D = I

	for (int i = 0;i < n_psrc;i++)
	{
		P(i, 0) = psrc[i].x(); P(i, 1) = psrc[i].y();
		Q(i, 0) = pdst[i].x(), Q(i, 1) = pdst[i].y();
	}


	// 计算各点新坐标
	const unsigned int w = ptr_image_backup_->width(), h = ptr_image_backup_->height();
	Eigen::MatrixXd coords_new(w * h, 2);

	coords_new.setZero();
	for (int i = 0;i < w;i++)
		for (int j = 0;j < h;j++)
		{
			double sigma_sum = 0;

			for (int k = 0; k < n_psrc;k++)
				sigma_sum += SIGMA(i, j, psrc[k].x(), psrc[k].y(), 3);
			for (int k = 0; k < n_psrc;k++)
			{
				Eigen::MatrixXd X(1, 2); 
				X(0, 0) = i; X(0, 1) = j;
				coords_new.row(i * h + j) += SIGMA(i, j, psrc[k].x(), psrc[k].y(), 3) / sigma_sum * (Q.row(k) + X - P.row(k));
			}
		}

	// 用ann算法填补空洞
	const size_t n_point = w * h;
	constexpr size_t dim = 2;
	constexpr size_t K = 4;//取最近邻四个点求平均，目前暂未启用

	ANNcoord* const coords = coords_new.data();
	ANNpointArray pts_Arr = annAllocPts(n_point, dim);
	for (size_t i = 0; i < n_point; i++)
		for (size_t j = 0; j < dim; j++)
			pts_Arr[i][j] = coords_new(i, j);

	ANNbd_tree tree(pts_Arr, n_point, dim, w, ANN_KD_STD, ANN_BD_SUGGEST);//, 4, ANN_KD_STD, ANN_BD_SUGGEST);

	ANNpoint queryPt = annAllocPt(dim);
	ANNidx idxArr[K];
	ANNdist distArr[K];

	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			queryPt[0] = i; queryPt[1] = j;
			tree.annkSearch(queryPt, K, idxArr, distArr);

			if (distArr[1] > 3) 
			{
				ptr_image_->setPixelColor(QPoint(i, j), QColor(255,255,255));
				continue; 
			}
			int inv_i = idxArr[0] / h, inv_j = idxArr[0] % h;//求余、整除
			QColor color = ptr_image_backup_->pixel(QPoint(inv_i, inv_j));
			ptr_image_->setPixelColor(QPoint(i, j), color);
		}
	}
	update();
}