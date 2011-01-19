#ifndef MAIN_H
#define MAIN_H

#include <QMainWindow>

class Main : public QMainWindow
{
    Q_OBJECT

    public:
        Main();

    private:
        class QGridLayout *layout;
        class Disp *img;
};

#endif
