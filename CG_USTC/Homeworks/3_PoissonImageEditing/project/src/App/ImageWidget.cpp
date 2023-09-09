#include "ImageWidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>
#include "ChildWindow.h"

using std::cout;
using std::endl;


ImageWidget::ImageWidget(ChildWindow* relatewindow)
{
	image_ = new QImage();
	image_backup_ = new QImage();

	draw_status_ = kNone;
	is_choosing_ = false;
	is_pasting_ = false;

	point_start_ = QPoint(0, 0);
	point_end_ = QPoint(0, 0);

	source_window_ = NULL;
}

ImageWidget::~ImageWidget(void)
{
}

int ImageWidget::ImageWidth()
{
	return image_->width();
}

int ImageWidget::ImageHeight()
{
	return image_->height();
}

void ImageWidget::set_draw_status_to_choose()
{
	draw_status_ = kChoose;	
}

void ImageWidget::set_draw_status_to_paste()
{
	draw_status_ = kPaste;
}

QImage* ImageWidget::image()
{
	return image_;
}

void ImageWidget::set_source_window(ChildWindow* childwindow)
{
	source_window_ = childwindow;
}

void ImageWidget::paintEvent(QPaintEvent* paintevent)
{
	QPainter painter;
	painter.begin(this);

	// Draw background
	painter.setBrush(Qt::lightGray);
	QRect back_rect(0, 0, width(), height());
	painter.drawRect(back_rect);

	// Draw image
	QRect rect = QRect(0, 0, image_->width(), image_->height());
	painter.drawImage(rect, *image_);

	// Draw choose region
	painter.setBrush(Qt::NoBrush);
	painter.setPen(Qt::red);
	painter.drawRect(point_start_.x(), point_start_.y(),
		point_end_.x() - point_start_.x(), point_end_.y() - point_start_.y());

	painter.end();
}

void ImageWidget::mousePressEvent(QMouseEvent* mouseevent)
{
	if (Qt::LeftButton == mouseevent->button())
	{
		switch (draw_status_)
		{
		case kChoose:
			is_choosing_ = true;
			point_start_ = point_end_ = mouseevent->pos();
			mask_selected = Mask(point_start_);
			break;

		case kPaste:
		{
			is_pasting_ = true;

			// Start point in object image
			int xpos = mouseevent->pos().rx();
			int ypos = mouseevent->pos().ry();

			// Start point in source image
			int xsourcepos = source_window_->imagewidget_->point_start_.rx();
			int ysourcepos = source_window_->imagewidget_->point_start_.ry();

			// Width and Height of rectangle region
			int w = source_window_->imagewidget_->point_end_.rx()
				- source_window_->imagewidget_->point_start_.rx() + 1;
			int h = source_window_->imagewidget_->point_end_.ry()
				- source_window_->imagewidget_->point_start_.ry() + 1;

			// Poission Editing with rectangle region
			// editing zone
			QRect zone_dst(xpos, ypos, w, h), zone_src(xsourcepos, ysourcepos, w, h);
			const QImage* img_src = source_window_->imagewidget_->image();

			Mask& mask_ = source_window_->imagewidget_->mask_selected;//.get_boundaryAndInner()

			*(image_) = *(image_backup_);
			mask_.poission_editing(*image(), *img_src, zone_dst, zone_src);
		}

		update();
		break;

		default:
			break;
		}
	}
}

void ImageWidget::mouseMoveEvent(QMouseEvent* mouseevent)
{
	switch (draw_status_)
	{
	case kChoose:
		// Store point position for rectangle region
		if (is_choosing_)
		{
			point_end_ = mouseevent->pos();
		}
		break;

	case kPaste:
		// Paste rectangle region to object image
		if (is_pasting_)
		{
			// Start point in object image
			int xpos = mouseevent->pos().rx();
			int ypos = mouseevent->pos().ry();

			// Start point in source image
			int xsourcepos = source_window_->imagewidget_->point_start_.rx();
			int ysourcepos = source_window_->imagewidget_->point_start_.ry();

			// Width and Height of rectangle region
			int w = source_window_->imagewidget_->point_end_.rx()
				- source_window_->imagewidget_->point_start_.rx() + 1;
			int h = source_window_->imagewidget_->point_end_.ry()
				- source_window_->imagewidget_->point_start_.ry() + 1;

			// Poission Editing with rectangle region
			// editing zone
			QRect zone_dst(xpos, ypos, w, h), zone_src(xsourcepos, ysourcepos, w, h);
			const QImage* img_src = source_window_->imagewidget_->image();
			
			Mask& mask_ = source_window_->imagewidget_->mask_selected;//.get_boundaryAndInner()

			*(image_) = *(image_backup_);
			mask_.poission_editing(*image(), *img_src, zone_dst, zone_src);
		}

	default:
		break;
	}

	update();
}

void ImageWidget::mouseReleaseEvent(QMouseEvent* mouseevent)
{
	switch (draw_status_)
	{
	case kChoose:
		if (is_choosing_)
		{
			point_end_ = mouseevent->pos();
			is_choosing_ = false;
			draw_status_ = kNone;
			
			// copied zone in source image
			int x1 = point_start_.x(), y1 = point_start_.y(),
				x2 = point_end_.x(), y2 = point_end_.y();
			int xpos = x1 < x2 ? x1 : x2, ypos = y1 < y2 ? y1 : y2,
				w = x1>x2 ? x1 - x2 : x2 - x1, h = y1>y2 ? y1 - y2 : y2 - y1;

			QRect zone(xpos, ypos, w, h);
			mask_selected.add(QPoint(x1, y2));
			mask_selected.add(QPoint(x2, y2));
			mask_selected.add(QPoint(x2, y1));
			mask_selected.getMask();
			mask_selected.get_laplace_precomposed();
		}

	case kPaste:
		if (is_pasting_)
		{
			is_pasting_ = false;
			draw_status_ = kNone;
		}

	default:
		break;
	}

	update();
}

void ImageWidget::Open(QString filename)
{
	// Load file
	if (!filename.isEmpty())
	{
		image_->load(filename);
		*(image_backup_) = *(image_);
	}

	//	setFixedSize(image_->width(), image_->height());
	//	relate_window_->setWindowFlags(Qt::Dialog);
	//	relate_window_->setFixedSize(QSize(image_->width(), image_->height()));
	//	relate_window_->setWindowFlags(Qt::SubWindow);

		//image_->invertPixels(QImage::InvertRgb);
		//*(image_) = image_->mirrored(true, true);
		//*(image_) = image_->rgbSwapped();
	cout << "image size: " << image_->width() << ' ' << image_->height() << endl;
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

	image_->save(filename);
}

void ImageWidget::Invert()
{
	for (int i = 0; i < image_->width(); i++)
	{
		for (int j = 0; j < image_->height(); j++)
		{
			QRgb color = image_->pixel(i, j);
			image_->setPixel(i, j, qRgb(255 - qRed(color), 255 - qGreen(color), 255 - qBlue(color)));
		}
	}

	// equivalent member function of class QImage
	// image_->invertPixels(QImage::InvertRgb);
	update();
}

void ImageWidget::Mirror(bool ishorizontal, bool isvertical)
{
	QImage image_tmp(*(image_));
	int width = image_->width();
	int height = image_->height();

	if (ishorizontal)
	{
		if (isvertical)
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, height - 1 - j));
				}
			}
		}
		else
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					image_->setPixel(i, j, image_tmp.pixel(i, height - 1 - j));
				}
			}
		}

	}
	else
	{
		if (isvertical)
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, j));
				}
			}
		}
	}

	// equivalent member function of class QImage
	//*(image_) = image_->mirrored(true, true);
	update();
}

void ImageWidget::TurnGray()
{
	for (int i = 0; i < image_->width(); i++)
	{
		for (int j = 0; j < image_->height(); j++)
		{
			QRgb color = image_->pixel(i, j);
			int gray_value = (qRed(color) + qGreen(color) + qBlue(color)) / 3;
			image_->setPixel(i, j, qRgb(gray_value, gray_value, gray_value));
		}
	}

	update();
}

void ImageWidget::Restore()
{
	*(image_) = *(image_backup_);
	point_start_ = point_end_ = QPoint(0, 0);
	update();
}
