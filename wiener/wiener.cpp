#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <stdint.h>


double unif()  // (0; 1]
{
    return (1.+rand())*(1./RAND_MAX);
}

double gaus(double s, double m)
{
    return sqrt(-2*log(unif()))*cos(2*M_PI*unif())*s + m;
}

typedef std::map<double, double> Function;
double lookupWiener(double x, Function &f)
{
    if(f.empty())
        f[0] = 0;

    double s, m;
    Function::const_iterator b = f.lower_bound(x);
    if(b == f.begin())
    {
        // Less than all known function values. (da = -∞)
        double db = b->first - x;
        s = sqrt(db);
        m = b->second;
    }
    else
    {
        Function::const_iterator a = b;
        --a;
        double da = x - a->first;
        if(b == f.end())
        {
            // Greater than all known function values. (db = ∞)
            s = sqrt(da);
            m = a->second;
        }
        else
        {
            // Somewhere in range of known function values.
            double db = b->first - x;
            s = sqrt(da*db/(da + db));
            m = (da*b->second + db*a->second)/(da + db);
        }
    }

    // Save value and return.
    return f[x] = gaus(s, m);
}

int main()
{
    srand(42);

    Function f;

    double startScaleX = 1, startScaleY = 2;
    double zoomFactor = 1e32;
    double finishScaleX = startScaleX/zoomFactor, finishScaleY = startScaleY/sqrt(zoomFactor);
    unsigned finishFrame = 800;
    for(unsigned frame = 0; frame != finishFrame; ++frame)
    {
        double scaleX = startScaleX*pow(finishScaleX/startScaleX, double(frame)/finishFrame);
        double scaleY = startScaleY*pow(finishScaleY/startScaleY, double(frame)/finishFrame);
        unsigned w = 500, h = 100;
        std::vector<uint8_t> img(w*h, 240);
        for(unsigned ix = 0; ix != w; ++ix)
        {
            double x = (ix - (w/2.))*(2./w*scaleX);
            double y = lookupWiener(x, f);
            f[x/zoomFactor] = y/sqrt(zoomFactor);  // Make the animation cyclic.
            int iy = y*(h/2./scaleY) + (h/2.);
            if(iy >= 0 && unsigned(iy) < h)
                img[ix + iy*w] = 0;
        }
        char fname[100];
        std::sprintf(fname, "wiener%04u", frame);
        FILE *file = std::fopen((std::string() + fname + ".pgm").c_str(), "wb");
        std::fprintf(file, "P5\n%u %u\n255\n", w, h);
        std::fwrite(&img[0], img.size(), 1, file);
        std::fclose(file);
        std::system((std::string() + "convert " + fname + ".pgm " + fname + ".png && rm " + fname + ".pgm").c_str());
    }
}
