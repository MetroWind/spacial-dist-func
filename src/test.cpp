#include <iostream>
#include <vector>

#include "trojectory.h"
#include "pbc.h"

int main()
{
    libmd::Trojectory t;
    t.open("../test/test.xtc", "../test/test.gro");
    t.nextFrame();


    t.close();

    return 0;
}
