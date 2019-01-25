#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/* Analytical solutions
 * A series of functions that satisfy Laplace's equation
 */
double func0(double x, double y)
{
    return 1.0;
}

double func1(double x, double y)
{
    return 2 * x * y;
}

double func2(double x, double y)
{
    return x * x - y * y;
}

double func3(double x, double y)
{
    return pow(x, 4) - 6 * pow(x, 2) * pow(y, 2) + pow(y, 4);
}
/* Function map that holds the analytical solutions */
static double (*const funcs[])(double, double) =
{
    func0,
    func1,
    func2,
    func3
};
#define X_MIN -1.0
#define X_MAX 1.0
#define Y_MIN -1.0
#define Y_MAX 1.0
#define NX 64
#define NY 64
#define FUNC_NUM 2

int main(int argc, char**argv)
{
    /* Args
     * ./InitCondition nx ny function_selection
     */
    double dx = (X_MAX - X_MIN) / (NX - 1);
    double dy = (Y_MAX - Y_MIN) / (NY - 1);
    double x = X_MIN;
    FILE* fptr;
    if ((fptr = fopen("initial.dat", "w")) == NULL)
    {
        exit(1);
    }
    double (*func)(double, double) = funcs[FUNC_NUM];

    for (size_t i = 0; i < NX; ++i)
    {
        double y = Y_MIN;
        for (size_t j = 0; j < NX; ++j)
        {
            fprintf(fptr, "%f %f %f\n", x, y, func(x, y));
            y += dy;
        }
        fprintf(fptr, "\n");
        x += dx;
    }

    if(fclose(fptr))
    {
        exit(1);
    }
}
