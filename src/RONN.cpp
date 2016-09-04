/*=========================================================
**
** This program is designed for running RONNv3.2 for detecting
** disordered residues in proteins.
** Date:	March 2012
** Designer:	Zheng Rong Yang (Exeter University)
** Programmer:	Varun Ramraj (University of Oxford)
** Contact:	varun@strubi.ox.ac.uk
**
** Substantial modifications to speed up code as well as facilitate wrapping
** with Python (or any other language). Can now also handle multiple
** FASTA-formatted records in an input file. Path to `data` directory must be
** specified as a command-line argument if not located at `../data` with respect
** to the executable file.
** Date:    Sept 2016
** Author:  Shyam Saladi (California Institute of Technology)
** Contact: saladi@caltech.edu
**
==========================================================*/

// Number of models to use
#define MODCNT 10

#include <string>
#include <limits.h>
#include <unistd.h>

#include "callBBF.h"

std::vector<RONNModel> models;


std::string getexepath()
// find current directory of executable on Linux
// http://stackoverflow.com/a/19535628/2320823
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

int read_all_models(std::string prefix, float disorder_weight)
/*****
	disorder_weight: VERY IMPORTANT. CONTROLS COST FUNCTION
*****/
{
    FILE *fp;
    std::string fModel, fPDF;
    std::stringstream *convert;

	for(int m = 0; m < MODCNT; m++)
	{
		convert = new std::stringstream();
		(*convert) << m;
		std::string mString = convert->str();
		fModel = prefix + "/c" + mString + "/model.rec";

		if(!(fp = fopen(fModel.c_str(), "r")))
		{
			std::cerr << "Can't find " << fModel << std::endl;
			return -1;
		}
		fclose(fp);

		fPDF = prefix + "/c" + mString + "/pdfs.rec";
		if(!(fp = fopen(fPDF.c_str(), "r")))
		{
			std::cerr << "Can't find " << fPDF << std::endl;
			return -1;
		}
		fclose(fp);

        models.push_back(read_model_data(fModel, fPDF, disorder_weight));
    }
    delete convert;

    return 0;
}


std::vector<double> predict_seq(std::string query)
{
    std::vector<double> scores(query.length(), 0.0);
    for(int m = 0; m < MODCNT; m++){
        predict_model(query, &models[m], scores);

        // for(int i = 0; i < query.length(); i++)
        //     printf("intermed\t%c\t%f\n", query[i], scores[i]);
    }

    for(int i = 0; i < query.length(); i++) {
        scores[i] = scores[i]/(double)MODCNT;
    }


    return scores;
}


int main(int argc, char *argv[])
{
	std::string query, line, prefix;

    if(argc == 2)
        prefix = argv[1];
    else if(argc > 2)
    {
        std::cerr << "Command line argument not recognized";
        return 0;
    }
    else
        prefix = getexepath() + "/../data/";

    // read the models
    int status = read_all_models(prefix, 0.53);
    if(status != 0)
        return(status);

    // print the header for the FASTA-formatted file
	getline(std::cin, line);
	std::cout << line;
	std::cout << "\n";

	// get the sequence
    while (std::cin >> line)
        query = query + line;

    // calculate scores
    std::vector<double> scores = predict_seq(query);

    //print results
    for(int i = 0; i < scores.size(); i++)
        printf("%c\t%f\n", query[i], scores[i]);

	return 0;
}
