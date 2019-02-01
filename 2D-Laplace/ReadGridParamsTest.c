#include "Utils.h"

int main(int argc, char** argv)
{
    GridParams params;

    ReadGridParams("grid_params.dat", &params);
    PrintGridParameters(&params);

    return 0;
}
