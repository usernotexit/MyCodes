#include "minidraw.h"

#include <QToolBar>

MiniDraw::MiniDraw(QWidget* parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	view_widget_ = new ViewWidget();
	Creat_Action();
	Creat_ToolBar();
	Creat_Menu();

	setCentralWidget(view_widget_);
}

void MiniDraw::Creat_Action() {
	Action_About = new QAction(tr("&About"), this);
	connect(Action_About, &QAction::triggered, this, &MiniDraw::AboutBox);

	Action_Line = new QAction(tr("&Line"), this);
	connect(Action_Line, SIGNAL(triggered()), view_widget_, SLOT(setLine()));

	Action_Rect = new QAction(tr("&Rect"), this);
	connect(Action_Rect, &QAction::triggered, view_widget_, &ViewWidget::setRect);

	Action_Circ = new QAction(tr("&Circle"), this);
	connect(Action_Circ, &QAction::triggered, view_widget_, &ViewWidget::setCirc);

	Action_Poly = new QAction(tr("&Polygon"), this);
	connect(Action_Poly, &QAction::triggered, view_widget_, &ViewWidget::setPoly);

	Action_Undo = new QAction(tr("&undo"), this);
	connect(Action_Undo, &QAction::triggered, view_widget_, &ViewWidget::undo);

	Action_ClearAll = new QAction(tr("&clear"), this);
	connect(Action_ClearAll, SIGNAL(triggered()), view_widget_, SLOT(clearAll()));
}

void MiniDraw::Creat_ToolBar() {
	pToolBar = addToolBar(tr("&Main"));
	pToolBar->addAction(Action_About);
	pToolBar->addAction(Action_Line);
	pToolBar->addAction(Action_Rect);
	pToolBar->addAction(Action_Circ);
	pToolBar->addAction(Action_Poly);

	pToolBar->addAction(Action_Undo);
	pToolBar->addAction(Action_ClearAll);
}

void MiniDraw::Creat_Menu() {
	pMenu = menuBar()->addMenu(tr("&Figure Tool"));
	pMenu->addAction(Action_About);
	pMenu->addAction(Action_Line);
	pMenu->addAction(Action_Rect);
	pMenu->addAction(Action_Circ);
	pMenu->addAction(Action_Poly);
	
	editMenu = menuBar()->addMenu(tr("Edit"));
	editMenu->addAction(Action_Undo);
	editMenu->addAction(Action_ClearAll);
}

void MiniDraw::AboutBox() {
	QMessageBox::about(this, tr("About"), tr("MiniDraw by Zegao Chen"));
}

MiniDraw::~MiniDraw() {}
