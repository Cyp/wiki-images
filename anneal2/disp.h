#ifndef DISP_H
#define DISP_H

#include <QWidget>
#include <vector>
#include <math.h>
#include <QString>

using std::vector;
using std::swap;


class Disp : public QWidget
{
    Q_OBJECT

    public:
        Disp(int x, int y, QWidget *parent = 0);
        ~Disp();
    public slots:
        virtual void tick();
        virtual void setFac(int nfac) { fac = exp(nfac/1000.); emit facChanged(QString::number(fac)); }
    signals:
        virtual void facChanged(const QString &);
    protected:
        virtual void paintEvent(QPaintEvent *ev);

    private:
        int phase;
        double fac;
        int bitx, bity;
        vector<unsigned int> bitmap;
};

#endif
