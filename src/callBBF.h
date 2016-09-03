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
    std::vector<double> w;
    std::vector<int> Length;
    std::vector< std::vector<short> > dbAA;
} RONNModel;

void read_model_record(std::string filename, RONNModel *model);

void read_pdf_record(std::string filename, RONNModel *model);

void align(std::vector<short> &seqAA, int i, int j, int rho[],
           RONNModel *model);

void detect(std::vector<short> &seqAA, std::vector<double> &estimate,
            RONNModel *model);

int predict_seq(std::string query, RONNModel *model,
				std::vector<double> &scores);

RONNModel read_model_data(std::string mod_fn, std::string pdf_fn,
						  double d_weight);
