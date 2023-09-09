#pragma once
#include <QWidget>
#include <qevent.h>

QT_BEGIN_NAMESPACE
class QImage;
class QPainter;
QT_END_NAMESPACE

class ImageWidget :
	public QWidget
{
	Q_OBJECT

public:
	ImageWidget(void);
	~ImageWidget(void);

protected:
	void paintEvent(QPaintEvent *paintevent);
	void addWarpingPoint(QPoint psrc, QPoint pdst);
	void delWarpingPoint(); // undo

public:
	void mousePressEvent(QMouseEvent* event);
	void mouseMoveEvent(QMouseEvent* event);
	void mouseReleaseEvent(QMouseEvent* event);

public slots:
	// File IO
	void Open();												// Open an image file, support ".bmp, .png, .jpg" format
	void Save();												// Save image to current file
	void SaveAs();												// Save image to another file

	// Image processing
	void Invert();												// Invert pixel value in image
	void Mirror(bool horizontal=false, bool vertical=true);		// Mirror image vertically or horizontally
	void TurnGray();											// Turn image to gray-scale map
	void Restore();												// Restore image to origin
	void setWarpingActivate();									// Change the edit status
	void Warp();
	void Warp_RBF_std();		// RBF ÷±Ω”ÃÓ≥‰
	void Warp_RBF_ann();	// ”√annÃÓ≥‰∑Ïœ∂								
	void Warp_IDW();	// IDW

private:
	QImage		*ptr_image_;				// image 
	QImage		*ptr_image_backup_;
	std::vector<QPoint> psrc;
	std::vector<QPoint> pdst;
	double angle; // rotation
	bool warpingActivate;
};

