#ifndef IMAGE_H
#define IMAGE_H

#include <stdint.h>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include <stdexcept>
#include <cmath>


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
            if (std::fwrite(&byteData[0], 3*sx*sy, 1, f))
            {}
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
        void addEllipse(double xa, double ya, double xb, double yb, Pixel const &p)
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
            if(uxb < uxa || uyb < uya)
                return;  // Degenerate ellipse.
            for(unsigned y = uya; y != uyb+1; ++y)
            {
                for(unsigned x = uxa; x != uxb+1; ++x)
                {
                    if (hypot(2*x/(txa + txb) - 1, 2*y/(tya + tyb) - 1) > 1)
                        continue;
                    data[x + y*sx] += p;
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

#endif //IMAGE_H
