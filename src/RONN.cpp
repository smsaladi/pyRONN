/*=========================================================
**
** This program is designed for running RONNv3.2 for detecting
** disordered residues in proteins.
** Date:	March 2012
** Designer:	Zheng Rong Yang (Exeter University)
** Programmer:	Varun Ramraj (University of Oxford)
** Contact:	varun@strubi.ox.ac.uk
**
** Modified by Shyam Saladi (saladi@caltech.edu, California Institute of Technology), 01-09-15
**
==========================================================*/

// Number of models to use
#define MODCNT 10
#include "callBBF.h"

//@@@@@ the path can be changed
std::string myPath="/ul/saladi/ronn/data";

std::vector<RONNModel> models;

// void predict_seq(char const *seq, float ronn[]) {
//
// }

int main(int argc, char **argv)
{
	std::stringstream *convert;
	FILE *fp;

	//$VR$: VERY IMPORTANT. CONTROLS COST FUNCTION
	double disorder_weight = 0.53;

	std::string fModel, fPDF;
	int r, m, nR;

	double mean;
	std::vector< std::vector<double> > X;

	std::string query;

	std::string line;
    // print the FASTA header
	getline(std::cin, line);
	std::cout << line;
	std::cout << "\n";

	// get the sequence
    while (std::cin >> line)
    {
        query = query + line;
    }
	nR = query.size();

    // Read all models
	for(m = 0; m < MODCNT; m++)
	{
		//sprintf(fModel,"%s/c%d/model.rec",myPath,m);
		convert = new std::stringstream();
		(*convert) << m;
		std::string mString = convert->str();

		fModel = myPath + "/c" + mString + "/model.rec";

	   	//TODO: get rid of FILE* stuff
		if(!(fp = fopen(fModel.c_str(),"r")))
		{
			std::cout << "Can't find " << fModel << std::endl;
			return -1;
		}
		fclose(fp);

		//sprintf(fPDF,"%s/c%d/pdfs.rec",myPath,m);
		fPDF = myPath + "/c" + mString + "/pdfs.rec";

		//TODO: get rid of FILE* stuff
		if(!(fp = fopen(fPDF.c_str(),"r")))
		{
			std::cout << "Can't find " << fPDF << std::endl;
			return -1;
		}
		fclose(fp);

        models.push_back(read_model_data(fModel, fPDF, disorder_weight));
    }

    for(m = 0; m < MODCNT; m++)
    {
        std::vector<double> scores;
        predict_seq(query, &models[m], scores);

        for(r = 0; r < nR; r++)
		{
			X.push_back(std::vector<double>());
			X[r].push_back(scores[r]);
		}
	}

	for(r = 0; r < nR; r++)
	{
		mean = 0;
		for(m = 0; m < MODCNT; m++)
			mean += X[r][m];
		mean /= (double)MODCNT;

		printf("%c\t%f\n", query[r], mean);
	}
	delete convert;

	return 0;
}

void read_all_models()
{
    return;
}
