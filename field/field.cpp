#include <cstdio>
#include <vector>
#include <stdexcept>
#include <complex>
#include <stdint.h>
#include <tr1/cmath>


class Function
{
    public:
        virtual ~Function() {}

        //virtual std::string name() const = 0;
        virtual std::complex<double> eval(double x, double y) const = 0;
};

class Identity : public Function
{
    std::complex<double> eval(double x, double y) const
    {
        return std::complex<double>(x, y);
    }
};

class Propagator : public Function
{
    public:
        Propagator(double mass/*, char const *name*/) : m(mass)/*, n(name)*/ {}
        //std::string name() const { return n/*"propagator"*/; }
        std::complex<double> eval(double x, double y) const
        {
            double s = y*y-x*x;  // y is time, x is a space dimension.
            //double const m = 0.0000000001;//.0001;
            if(s >= -1.e-20 && s <= 1.e-20)
                return 1.e20;
            if(s >= 0)
                return m/(8*M_PI*sqrt(s))*std::complex<double>(std::tr1::cyl_bessel_j(1, m*sqrt(s)), std::tr1::cyl_neumann(1, m*sqrt(s)));
            else
                return m/(4*M_PI*M_PI*sqrt(-s))*std::complex<double>(0, -std::tr1::cyl_bessel_k(1, m*sqrt(-s)));
            return std::complex<double>(x, y);
        }

    private:
        double m;
        //char const *n;
};

struct Pixel
{
    Pixel() { col[0] = 0xFF; col[1] = 0x00; col[2] = 0xFF; }
    Pixel(double r, double g, double b)
    {
        col[0] = std::max(0, std::min(255, int(0.5+255*r)));
        col[1] = std::max(0, std::min(255, int(0.5+255*g)));
        col[2] = std::max(0, std::min(255, int(0.5+255*b)));
    }
    //Pixel(std::complex<double> const &val) { col[0] = 0; col[1] = std::max(0, std::min(255, int(127.5+63.75*real(val)))); col[2] = std::max(0, std::min(255, int(127.5+63.75*imag(val)))); }

    uint8_t col[3];
};

Pixel colourComplex(std::complex<double> const &z)
{
    double r = abs(z);
    double theta = arg(z);

    double y = tanh(r/2);
    double m = y*(1-y);
    double s = m*sin(theta);
    double c = m*cos(theta);

    double const q = .5;
    double const w = sqrt(.75);

    return Pixel(y + c, y - q*c + w*s, y - q*c - w*s);
}

class Image
{
    public:
        Image(unsigned x, unsigned y, Function const &function_, double x1_, double y1_, double x2_, double y2_)
            : sx(x), sy(y), function(&function_), /*filename(function_.name() + ".ppm"),*/ x1(x1_), y1(y1_), x2(x2_), y2(y2_)
        {
            data.resize(x*y);
//#pragma omp parallel for
            for(int py = 0; py < int(y); ++py)
                for(unsigned px = 0; px != x; ++px)
                    data[px+py*sx] = colourComplex(function->eval(x1+(x2-x1)*(px/(sx-1.)), y1+(y2-y1)*(py/(sy-1.))));
        }
        void save(std::string const &name)
        {
            std::string filename = name+".ppm";
            FILE *f = fopen(filename.c_str(), "wb");
            if(f == NULL)
                throw std::runtime_error("Couldn't open file \""+filename+"\" for writing.");
            std::fprintf(f, "P6\n%u %u\n255\n", sx, sy);
            std::fwrite(&data[0], 3*sx*sy, 1, f);
            std::fclose(f);
        }

    private:
        unsigned sx, sy;
        Function const *function;
        //std::string filename;
        double x1, y1, x2, y2;
        std::vector<Pixel> data;
};

int main()
{
    Image identity(1001, 1001, Identity(), -2, -2, 2, 2);
    identity.save("PropagatorColours");

    Image propagator1(1001, 1001, Propagator(.2), -2, -2, 2, 2);
    propagator1.save("FeynmanPropagatorWithMass0.2");
    Image propagator2(1001, 1001, Propagator(2), -2, -2, 2, 2);
    propagator2.save("FeynmanPropagatorWithMass2");
    Image propagator3(1001, 1001, Propagator(20), -2, -2, 2, 2);
    propagator3.save("FeynmanPropagatorWithMass20");
    Image propagator4(1001, 1001, Propagator(200), -2, -2, 2, 2);
    propagator4.save("FeynmanPropagatorWithMass200");
}
