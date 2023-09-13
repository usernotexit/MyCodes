#include "minidraw.h"
#include <QtWidgets/QApplication>

#include <QcoreApplication>
//#include "signalsandslots3.h"
int main(int argc, char* argv[]) {
	_CrtSetDbgFlag(_CrtSetDbgFlag(_CRTDBG_REPORT_FLAG) | _CRTDBG_LEAK_CHECK_DF);

	QApplication a(argc, argv);
	MiniDraw w;
	w.show();
	//QCoreApplication a(argc, argv);
	//SignalsAndSlots3 w;
	return a.exec();
}
