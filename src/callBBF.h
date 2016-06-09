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
#include <omp.h>

using namespace std;

void model (string filename);

void pdf(string filename);

void align(int i,int j);

void detect(vector<double> & estimate);

int callBBF_driver(string query, string mod_fn, string pdf_fn1, double d_weight, vector<double> & scores);

int main(int argc,char **argv);
