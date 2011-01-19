#include <cstdio>
#include <stdint.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include "../headers/image.h"



void render(Image &img, std::vector<double> const &func, double xa, double xb, Pixel const &p)
{
    double dx = (xb-xa)/(func.size()-1);
    double lx = img.getDx()*.75;  // 1/2 brush width.
    double ly = img.getDy()*.75;  // 1/2 brush height.
    double invRectPix = 1./(1.5*1.5);  // Reciprocal area of brush in pixels.
    double px = 1/img.getDx();
    double py = 1/img.getDy();
    double pdx = dx*px;  // dx in pixels.
    double x = xa+dx*.5;
    //printf("test %zu  %lf %lf\n", func.size(), lx, ly);
    for(unsigned i = 0; i != func.size()-1; ++i)
    {
        double dy = func[i+1]-func[i];
        double pdy = dy*py;  // dy in pixels.
        double y = (func[i]+func[i+1])*-.5;
        img.addRect(x-lx, y-ly, x+lx, y+ly, p*(invRectPix*sqrt(pdx*pdx + pdy*pdy)));
        //if(i == func.size()/2)
        //    printf("%lf %lf %lf %lf %lf\n", x-lx, y-ly, x+lx, y+ly, (invRectPix*sqrt(pdx*pdx + pdy*pdy)));
        x += dx;
    }
}


void setTanh(std::vector<double> &func, std::vector<double> &dfuncdt, double xa, double xb)
{
    double dx = (xb-xa)/(func.size()-1);
    double x = xa;
    for(unsigned i = 0; i != func.size(); ++i)
    {
        func[i] = tanh(x);
        x += dx;
    }
    dfuncdt.assign(func.size(), 0);
}
void setDiTanh(std::vector<double> &func, std::vector<double> &dfuncdt, double xa, double xb)
{
    double dx = (xb-xa)/(func.size()-1);
    double x = xa;
    for(unsigned i = 0; i != func.size(); ++i)
    {
        func[i] = tanh(x+4) - tanh(x-4) - 1;
        x += dx;
    }
    dfuncdt.assign(func.size(), 0);
}
void setDiTanhBump(std::vector<double> &func, std::vector<double> &dfuncdt, double xa, double xb)
{
    dfuncdt.assign(func.size(), 0);
    double dx = (xb-xa)/(func.size()-1);
    double x = xa;
    for(unsigned i = 0; i != func.size(); ++i)
    {
        double v1 = 0.05;  // Velocity of leftmost bump. (Exact velocity = sinh(v1).)
        double v2 = -0.05;  // Velocity of rightmost bump.
        double tanh1 = tanh((x+6)*cosh(v1));
        double tanh2 = tanh((x-6)*cosh(v2));
        func[i] = tanh1 - tanh2 - 1;// + 0.2*exp(-((x+1)*(x+1)));
        dfuncdt[i] = (1-tanh1*tanh1)*-sinh(v1) + (1-tanh2*tanh2)*sinh(v2);
        x += dx;
    }
}

// \mathcal{L} = \partial_\mu\phi\partial^\mu\phi - (\phi^2 - 1)^2
// d²f/dt² = d²f/dx² - 2f(f²-1)²
void stepTime(std::vector<double> &func, std::vector<double> &dfuncdt, double xa, double xb, double dt)
{
    assert(func.size() == dfuncdt.size());
    double dx = (xb-xa)/(func.size()-1);
    double dtinvdx2 = dt/(dx*dx);
    dfuncdt[0] = (func[1] - func[0]) / dt;
    for(unsigned i = 1; i != func.size()-1; ++i)
        dfuncdt[i] += dtinvdx2*(func[i+1] + func[i-1] - 2*func[i]) - 2*dt*func[i]*(func[i]*func[i] - 1);
    dfuncdt[func.size()-1] = -(func[func.size()-1] - func[func.size()-2]) / dx;
    for(unsigned i = 0; i != func.size(); ++i)
        func[i] += dt*dfuncdt[i];
}

int main()
{
    double fxa = -208;  // Function min x.
    double fxb = 208;  // Function max x.
    unsigned fxr = 20000;  // Function resolution.

    double imxa = -8;  // Display min x.
    double imxb = 8;  // Display max x.
    double imya = -2;  // Display min f(x).
    double imyb = 2;  // Display max f(x).

    Image img(801, 201, imxa, imya, imxb, imyb);
    img.addRect(-1000, 0-img.getDx(), 1000, 0+img.getDx(), -Pixel(1, 1, 1));
    img.addRect(-1000, 1-img.getDx(), 1000, 1+img.getDx(), -Pixel(.5, 1, 1)*.5);
    img.addRect(-1000, -1-img.getDx(), 1000, -1+img.getDx(), -Pixel(.5, 1, 1)*.5);
    std::vector<double> func(fxr);
    std::vector<double> dfuncdt(fxr);
    setTanh(func, dfuncdt, fxa, fxb);
    render(img, func, fxa, fxb, -Pixel(1, 1, .75));
    img.save("soliton");

    setDiTanhBump(func, dfuncdt, fxa, fxb);
    for(unsigned frame = 0; frame != 800; ++frame)
    {
        Image imgf(801, 201, imxa, imya, imxb, imyb);
        imgf.addRect(-1000, 0-img.getDx(), 1000, 0+img.getDx(), -Pixel(1, 1, 1));
        imgf.addRect(-1000, 1-img.getDx(), 1000, 1+img.getDx(), -Pixel(.5, 1, 1)*.5);
        imgf.addRect(-1000, -1-img.getDx(), 1000, -1+img.getDx(), -Pixel(.5, 1, 1)*.5);
        render(imgf, func, fxa, fxb, -Pixel(1, 1, .75));
        for(unsigned i = 0; i != 1000; ++i)
            stepTime(func, dfuncdt, fxa, fxb, 1./2000);
        //render(imgf, func, -8, 8, -Pixel(1, .75, 1));
        //render(imgf, dfuncdt, -8, 8, -Pixel(.75, 1, .75)*.25);
        char str[100];
        std::sprintf(str, "disoliton%03u", frame);
        imgf.save(str);
    }
}
