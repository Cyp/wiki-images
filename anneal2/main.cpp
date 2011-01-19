#include <QApplication>
#include <QSlider>
#include "disp.h"
#include "main.h"
#include <QGridLayout>
#include <QLabel>

int main(int aa, char **bb)
{
    QApplication prog(aa, bb);
    Main m;
    m.show();
    return prog.exec();
}

Main::Main()
{
    QWidget *wdgt = new QWidget(this);
    setCentralWidget(wdgt);
    layout = new QGridLayout(wdgt);
    img = new Disp(100, 100);
    layout->addWidget(img, 1, 0);
    QSlider *slid = new QSlider(Qt::Vertical);
    slid->setRange(0, 27000);
    layout->addWidget(slid, 0, 1, -1, 1);
    slid->show();
    connect(slid, SIGNAL(sliderMoved(int)), img, SLOT(setFac(int)));
    img->show();
    QLabel *lb = new QLabel();
    layout->addWidget(lb, 0, 0);
    connect(img, SIGNAL(facChanged(const QString &)), lb, SLOT(setText(const QString &)));
    //lb->show();
}
