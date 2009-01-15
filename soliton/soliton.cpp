#include <cstdio>
#include <stdint.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <cassert>


struct Pixel
{
    Pixel() { col[0] = 1; col[1] = 1; col[2] = 1; }
    Pixel(double r, double g, double b)
    {
        col[0] = r;
        col[1] = g;
        col[2] = b;
    }
    Pixel const &operator +=(Pixel const &p)
    {
        col[0] += p.col[0];
        col[1] += p.col[1];
        col[2] += p.col[2];
        return *this;
    }
    Pixel operator -() const
    {
        return Pixel(-col[0], -col[1], -col[2]);
    }
    Pixel operator *(float v) const
    {
        return Pixel(col[0]*v, col[1]*v, col[2]*v);
    }

    float col[3];
};

inline int quant(int v)
{
    //return (v>>4)*0x11;
    return v;
}

struct BytePixel
{
    BytePixel() { col[0] = 0xFF; col[1] = 0xFF; col[2] = 0xFF; }
    BytePixel(Pixel const &p)
    {
        col[0] = quant(std::max(0, std::min<int>(255, 256*p.col[0])));
        col[1] = quant(std::max(0, std::min<int>(255, 256*p.col[1])));
        col[2] = quant(std::max(0, std::min<int>(255, 256*p.col[2])));
    }

    uint8_t col[3];
};

class Image
{
    public:
        Image(unsigned x, unsigned y, double x1_, double y1_, double x2_, double y2_)
            : sx(x), sy(y), x1(x1_), y1(y1_), x2(x2_), y2(y2_), facx(sx/(x2-x1)), facy(sy/(y2-y1)), pixx((x2-x1)/sx), pixy((y2-y1)/sy)
        {
            data.resize(x*y);
//#pragma omp parallel for
        }
        void save(std::string const &name)
        {
            std::vector<BytePixel> byteData(data.begin(), data.end());
            std::string filename = name+".ppm";
            FILE *f = fopen(filename.c_str(), "wb");
            if(f == NULL)
                throw std::runtime_error("Couldn't open file \""+filename+"\" for writing.");
            std::fprintf(f, "P6\n%u %u\n255\n", sx, sy);
            std::fwrite(&byteData[0], 3*sx*sy, 1, f);
            std::fclose(f);
        }
        void addRect(double xa, double ya, double xb, double yb, Pixel const &p)
        {
            if(xa > x2 || xb < x1 || ya > y2 || yb < y1)
                return;
            double txa = std::max<double>(0, (xa-x1)*facx);
            double tya = std::max<double>(0, (ya-y1)*facy);
            double txb = std::min<double>(sx, (xb-x1)*facx);
            double tyb = std::min<double>(sy, (yb-y1)*facy);
            unsigned uxa = txa;
            unsigned uya = tya;
            unsigned uxb = std::min<unsigned>(sx-1, txb);
            unsigned uyb = std::min<unsigned>(sy-1, tyb);
            double xstart = txa-uxa;  // Left edge missing in pixels.
            double ystart = tya-uya;  // Top edge missing in pixels.
            double xend = uxb+1-txb;  // Right edge missing in pixels.
            double yend = uyb+1-tyb;  // Bottom edge missing in pixels.
            if(uxb < uxa || uyb < uya)
                return;  // Degenerate rectangle.
            for(unsigned y = uya; y != uyb+1; ++y)
            {
                double yamount = 1;
                if(y == uya)
                    yamount -= ystart;
                if(y == uyb)
                    yamount -= yend;
                for(unsigned x = uxa; x != uxb+1; ++x)
                {
                    double xamount = 1;
                    if(x == uxa)
                        xamount -= xstart;
                    if(x == uxb)
                        xamount -= xend;
                    data[x + y*sx] += p * (xamount*yamount);
                }
            }
        }
        double getDx() const { return pixx; }
        double getDy() const { return pixy; }

    private:
        unsigned sx, sy;
        double x1, y1, x2, y2;
        double facx, facy;
        double pixx, pixy;
        std::vector<Pixel> data;
};


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
