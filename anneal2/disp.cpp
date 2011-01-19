#include "disp.h"
#include "graph.h"
#include <QPainter>
#include <QTimer>
#include <math.h>
#include <stdlib.h>
#include <QSlider>
#include <QImage>

Disp::Disp(int x, int y, QWidget *parent) : QWidget(parent)
{
    setWindowTitle("Display");
    QTimer *timer = new QTimer(this);
    bitx = x;
    bity = y;
    setFixedSize(bitx, bity);
    bitmap.resize(bitx*bity);
    int i;
    for(i = 0; i < bitx*bity; ++i)
        bitmap[i] = rand()&0xFFFFFF;
    //int i, j;
    //for(j = 0; j < bity; ++j) for(i = 0; i < bitx; ++i) bitmap[i+j*bitx] = (i < bitx/2 ? 0xFF0000 : 0)+(j < bity/2 ? 0x00FF00 : 0)+(rand()&1)*0xFF;
    //QImage loaded = QImage("pict3895.jpg").convertToFormat(QImage::Format_RGB32).scaled(bitx, bity, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    //bitmap.assign((unsigned int *)loaded.bits(), ((unsigned int *)loaded.bits())+(bitx*bity));
    //bitmap.resize(bitx*(bity+1));
    connect(timer, SIGNAL(timeout()), this, SLOT(tick()));
    timer->start(10);
    phase = 0;//0.6;
    fac = 0.00000001;//10000;
}

Disp::~Disp()
{
    QImage saver((uchar *)&bitmap[0], bitx, bity, QImage::Format_RGB32);
    saver.save("out.png", "PNG");
}

/*
class ColCount
{
    public:
        ColCount() : r(0), g(0), b(0) {}
        ColCount &operator += (int c)
        {
            r += (c>>16)&0xFF;
            g += (c>>8 )&0xFF;
            b += (c    )&0xFF;
            return *this;
        }
        ColCount &operator -= (int c)
        {
            r -= (c>>16)&0xFF;
            g -= (c>>8 )&0xFF;
            b -= (c    )&0xFF;
            return *this;
        }
        int evalc() { return r*g*b; }
        
        int r, g, b;
};
*/

double calEnergy(int a, int b, int dx, int dy)
{
    int rad = dx*dx+dy*dy;
    const int maxRad = 6*6;
    if(rad >= maxRad)
        return 0;
    int r1 = (a>>16)&0xFF,
        g1 = (a>>8 )&0xFF,
        b1 = (a    )&0xFF,
        r2 = (b>>16)&0xFF,
        g2 = (b>>8 )&0xFF,
        b2 = (b    )&0xFF;
    double e = r1*r2+g1*g2+b1*b2;
    e *= ((maxRad/2)-rad)/double(maxRad);
    //e *= ((maxRad)-rad)/double(maxRad);
    return e;
}

void Disp::tick()
{
    phase = (phase+1)%bitx;//0.001;
    //fac *= 0.99;
    //if(fac < 0.001)
    //    fac = 1/0.00000000000000000001;//10000;
    int n;
    for(n = 0; n < 1000; ++n)
    {
        int x1, y1, x2, y2, dx, dy;
        x1 = rand()%bitx;
        y1 = rand()%bity;
        x2 = ((dx = rand()%21-10)+x1+bitx)%bitx;
        y2 = ((dy = rand()%21-10)+y1+bity)%bity;
        if( dx*dx+dy*dy > 100 )
        {
            --n;
            continue;
        }
        //double e = (x2-x1)*(((bitmap[x2+y2*bitx]&0x0000FF)   )-((bitmap[x1+y1*bitx]&0x0000FF)   ))
        //          +(y2-y1)*(((bitmap[x2+y2*bitx]&0x00FF00)>>8)-((bitmap[x1+y1*bitx]&0x00FF00)>>8));
        double e = 0;
        //ColCount cc1, cf1, cc2, cf2;
        //int dx, dy;
        for(dy = -5; dy <= 5; ++dy) for(dx = -5; dx <= 5; ++dx)
        {
            e -= calEnergy(bitmap[x1+y1*bitx], bitmap[(x1+dx+bitx)%bitx+(y1+dy+bity)%bity*bitx], dx, dy);
            e -= calEnergy(bitmap[x2+y2*bitx], bitmap[(x2+dx+bitx)%bitx+(y2+dy+bity)%bity*bitx], dx, dy);
            //cc1 += bitmap[(x1+dx)%bitx+(y1+dy)%bity*bitx];
            //cc2 += bitmap[(x2+dx)%bitx+(y2+dy)%bity*bitx];
        }
        swap(bitmap[x1+y1*bitx], bitmap[x2+y2*bitx]);
        for(dy = -9; dy <= 9; ++dy) for(dx = -9; dx <= 9; ++dx)
        {
            e += calEnergy(bitmap[x1+y1*bitx], bitmap[(x1+dx+bitx)%bitx+(y1+dy+bity)%bity*bitx], dx, dy);
            e += calEnergy(bitmap[x2+y2*bitx], bitmap[(x2+dx+bitx)%bitx+(y2+dy+bity)%bity*bitx], dx, dy);
        }
        if(!(e>0 || (double)rand()/RAND_MAX < exp(e/fac)))
            swap(bitmap[x1+y1*bitx], bitmap[x2+y2*bitx]);
    }
//    bitmap[(rand()^(rand()&rand()))%400*400+(rand()^(rand()&rand()))%400] = 0xFFFF0000;
    update();
}

void Disp::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    QImage img((unsigned char *)&bitmap[0/*phase*/], bitx, bity, QImage::Format_RGB32);
    p.setRenderHint(QPainter::Antialiasing);
    p.setRenderHint(QPainter::SmoothPixmapTransform);
    p.drawImage(QRect(0, 0, width(), height()), img);
    /*int x;
    for(x = 0; x<20; ++x)
    {
        p.drawEllipse(QRectF(width()*(1-sin(phase))/2, height()*(1-cos(phase))/2, width()*(sin(phase)), height()*(cos(phase))));
        //p.drawEllipse(int(width()*(1-sin(phase))/2), int(height()*(1-cos(phase))/2), int(width()*(sin(phase))), int(height()*(cos(phase))));
        phase += .05;
    }
    phase -= 1;
    */
    /*p.drawEllipse(int(width()*(1-sin(phase))/2), int(height()*(1-cos(phase))/2), int(width()*(sin(phase))), int(height()*(cos(phase))));*/
}
