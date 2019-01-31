#include "Utils.h"

int main(int argc, char** argv)
{
    Params params;

    ReadGridParams("grid_params.dat", &params);
    PrintParameters(&params);

    return 0;
}
