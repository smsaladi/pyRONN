#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

// Parameters for each model
typedef struct RONNModelTag
{
    double mu[2], sigma[2], disorder_weight;
    int nD, nW;
    std::vector< std::vector<short> >	dbAA;
} RONNModel;

RONNModel read_model_file(std::string filename);

void read_pdf_file(std::string filename, RONNModel *model);

void align(int i, int j, int nW, int rho[], RONNModel *model);

void detect(std::vector<double> & estimate, RONNModel *model);

int callBBF_driver(std::string query, std::string mod_fn, std::string pdf_fn1,
                   double d_weight, std::vector<double> & scores);

int main(int argc, char **argv);
