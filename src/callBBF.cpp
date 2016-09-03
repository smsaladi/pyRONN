/*

The critical parameters for this script are
    Maximum bases is 2000 (parameter name is "maxD")
    Maximum residues is 10000 (parameter name is "maxP")

The script returns the prediction for a sequence under testing using one trained cross-validation model
** Modified by Shyam Saladi (saladi@caltech.edu, California Institute of Technology)

*/
#include "callBBF.h"

static
//               A b C D E F G H I j K L M  N  o P  Q  R  S  T  u V  W  x Y
int	INDEX[25] = {0,0,1,2,3,4,5,6,7,0,8,9,10,11,0,12,13,14,15,16,0,17,18,0,19};

static
int	Dayhoff[20][20] = {
{40, 24, 32, 32, 16, 36, 28, 28, 28, 24, 28, 32, 36, 32, 24, 36, 36, 32,  8, 20},
{24, 80, 12, 12, 16, 20, 20, 24, 12,  8, 12, 16, 20, 12, 16, 32, 24, 24,  0, 32},
{32, 12, 48, 44,  8, 36, 36, 24, 32, 16, 20, 40, 28, 40, 28, 32, 32, 24,  4, 16},
{32, 12, 44, 48, 12, 32, 36, 24, 32, 20, 24, 36, 28, 40, 28, 32, 32, 24,  4, 16},
{16, 16,  8, 12, 68, 12, 24, 36, 12, 40, 32, 16, 12, 12, 16, 20, 20, 28, 32, 60},
{36, 20, 36, 32, 12, 52, 24, 20, 24, 16, 20, 32, 28, 28, 20, 36, 32, 28,  4, 12},
{28, 20, 36, 36, 24, 24, 56, 24, 32, 24, 24, 40, 32, 44, 40, 28, 28, 24, 20, 32},
{28, 24, 24, 24, 36, 20, 24, 52, 24, 40, 40, 24, 24, 24, 24, 28, 32, 48, 12, 28},
{28, 12, 32, 32, 12, 24, 32, 24, 52, 20, 32, 36, 28, 36, 44, 32, 32, 24, 20, 16},
{24,  8, 16, 20, 40, 16, 24, 40, 20, 56, 48, 20, 20, 24, 20, 20, 24, 40, 24, 28},
{28, 12, 20, 24, 32, 20, 24, 40, 32, 48, 56, 24, 24, 28, 32, 24, 28, 40, 16, 24},
{32, 16, 40, 36, 16, 32, 40, 24, 36, 20, 24, 40, 28, 36, 32, 36, 32, 24, 16, 24},
{36, 20, 28, 28, 12, 28, 32, 24, 28, 20, 24, 28, 56, 32, 32, 36, 32, 28,  8, 12},
{32, 12, 40, 40, 12, 28, 44, 24, 36, 24, 28, 36, 32, 48, 36, 28, 28, 24, 12, 16},
{24, 16, 28, 28, 16, 20, 40, 24, 44, 20, 32, 32, 32, 36, 56, 32, 28, 24, 40, 16},
{36, 32, 32, 32, 20, 36, 28, 28, 32, 20, 24, 36, 36, 28, 32, 40, 36, 28, 24, 20},
{36, 24, 32, 32, 20, 32, 28, 32, 32, 24, 28, 32, 32, 28, 28, 36, 44, 32, 12, 20},
{32, 24, 24, 24, 28, 28, 24, 48, 24, 40, 40, 24, 28, 24, 24, 28, 32, 48,  8, 24},
{ 8,  0,  4,  4, 32,  4, 20, 12, 20, 24, 16, 16,  8, 12, 40, 24, 12,  8,100, 32},
{20, 32, 16, 16, 60, 12, 32, 28, 16, 28, 24, 24, 12, 16, 16, 20, 20, 24, 32, 72}};

static
int Blosum62[20][20]={
{ 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-2,-3,-2},
{ 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2},
{-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0, 1,-3,-4,-3},
{-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0, 0,-3,-3,-2},
{-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3},
{ 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3,-2,-2,-2,-2, 0, 1, 0,-2,-3},
{-2,-3, 1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1, 0,-2,-2, 2},
{-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-2, 1,-3,-1},
{-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0, 0,-3,-3,-2},
{-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-2, 3,-2,-1},
{-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1,-2,-1,-1},
{-2,-3, 1, 0,-3, 0,-1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2},
{-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-1, 7,-1,-2,-1, 1,-2,-4,-3},
{-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0, 0,-2,-2,-1},
{-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2},
{ 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2},
{-1,-1, 1, 0,-2, 1, 0,-2, 0,-2,-1, 0, 1, 0,-1, 1, 4,-2,-3,-2},
{ 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2,-2, 4,-3,-1},
{-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-3,-3,11, 2},
{-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7}};


void read_model_record(std::string filename, RONNModel *model)
{
	std::string str, tmp, tmp2;
	std::stringstream *convert;
	std::vector<std::string> filedata;

	std::ifstream fp(filename.c_str());
	if (fp.is_open())
	{
		while (fp.good())
		{
			getline (fp, tmp);
			filedata.push_back(tmp);
		}
		fp.close();
	}
	else
	{
		std::cout << "No model found." << std::endl;
	}

	//first deal with headers
	//line 1: database size (nD)
	//line 2: window length (nW)

	convert = new std::stringstream(filedata[0]);
	(*convert) >> model->nD;
	convert = new std::stringstream(filedata[1]);
	(*convert) >> model->nW;

	std::vector<std::string> seqdata;
	//j for database sequence index

	for(int j = 2; j < (model->nD*2)+2; j++)
	{
		if (j % 2 == 0)
		{
			seqdata.push_back(filedata[j]);
		}
		else
		{
			convert = new std::stringstream(filedata[j]);
			double wtmp = 0.0;
			(*convert) >> wtmp;
			model->w.push_back(wtmp);
		}
	}
	filedata.clear();

	for (int j = 0; j < model->nD; j++)
	{
		str = seqdata[j];
		model->Length.push_back(str.length());

		model->dbAA.push_back(std::vector<short>());

		for (int r = 0; r < str.length(); r++)
			model->dbAA[j].push_back(INDEX[(int)(str[r] - 'A')]);
	}

	delete convert;
}


void read_pdf_record(std::string filename, RONNModel *model)
{
	std::string tmp;
	std::vector<std::string> filedata;
	std::stringstream *convert;

	std::ifstream fp(filename.c_str());
	if (fp.is_open())
	{
		while (fp.good())
		{
			getline (fp, tmp);
			filedata.push_back(tmp);
		}
		fp.close();
	}
	else
	{
		std::cout << "No PDF record!" << std::endl;
	}

	//deal with the four lines sequentially
	convert = new std::stringstream(filedata[0]);
	(*convert) >> model->mu[0];

	convert = new std::stringstream(filedata[1]);
	(*convert) >> model->mu[1];

	convert = new std::stringstream(filedata[2]);
	(*convert) >> model->sigma[0];

	convert = new std::stringstream(filedata[3]);
	(*convert) >> model->sigma[1];

	//std::cout << mu[0] << " " << mu[1] << " " << sigma[0] << " " << sigma[1] << endl;
	//int tmp2 = fscanf(fp,"%f %f %f %f",&mu[0],&mu[1],&sigma[0],&sigma[1]);

	delete convert;
}


void align(std::vector<short> &seqAA, int i, int j, int rho[], RONNModel *model)
// i for query sequence index
// j for database sequence index
{
	int R, score;
	rho[0] = 0;
	rho[1] = -1000000;

	//go though the database sequence for maximised alignment
	for(int r = 0; r <= model->Length[j]-model->nW; r++)
	{
		score = 0;

		//go through the query sequence for one alignment
		for(int w = 0; w < model->nW; w++)
			score += Blosum62[seqAA[i+w]][model->dbAA[j][r+w]];

		if(score > rho[1])
		{
			rho[1] = score;
			R = r;
		}
	}

	for(int w = 0; w < model->nW; w++)
		rho[0] += Blosum62[model->dbAA[j][R+w]][model->dbAA[j][R+w]];

	return;
}


void detect(std::vector<short> &seqAA, std::vector<double> &estimate,
			RONNModel *model)
{
	double y, fOrder, fDisor, pOrder, pDisor;
	std::vector<int> predictTimes;
	std::map < std::map <int, int>, double > yBar;

	int rho[2];

	//string estimate_filename = "estimate.rec";

	for(int i = 0; i < seqAA.size(); i++)
		predictTimes.push_back(0);

	for(int i = 0; i <= seqAA.size()-model->nW; i++)
	{
		y = 0.0;
		for(int j = 0; j < model->nD; j++)
		{
			// search for the maximum alignment between ith query peptide
			// jth database sequence
			align(seqAA, i, j, rho, model);
			y += model->w[j] * exp((double)(rho[1]-rho[0])/(double)rho[0]);
		}

		fOrder = exp(-0.5*pow(y-model->mu[0], 2.0)/model->sigma[0])/(sqrt(M_2_PI*model->sigma[0])); //$VR$: bug fixed by Ron in Feb07
		fDisor = exp(-0.5*pow(y-model->mu[1], 2.0)/model->sigma[1])/(sqrt(M_2_PI*model->sigma[1])); //$VR$: bug fixed by Ron in Feb07

		pDisor = model->disorder_weight*fDisor/((1.0-model->disorder_weight)*fOrder+model->disorder_weight*fDisor);

		for(int r = i; r < i+model->nW; r++)
		{
			std::map <int, int> r_predtimes;
			r_predtimes[r] = predictTimes[r];

			yBar[r_predtimes] = pDisor;
			predictTimes[r]++;
		}
	}

	for(int i = 0; i < seqAA.size(); i++)
	{
		y = 0.0;
		for(int r = 0; r < predictTimes[i]; r++)
		{
			std::map<int, int> tmp_vals;
			tmp_vals[i] = r;
			y += yBar[tmp_vals];
		}
		estimate.push_back((double)y/(double)predictTimes[i]);
	}

	predictTimes.clear();
}


RONNModel read_model_data(std::string mod_fn, std::string pdf_fn,
						  double d_weight)
{
	RONNModel model;
	read_model_record(mod_fn, &model);
	model.disorder_weight = d_weight;
	read_pdf_record(pdf_fn, &model);
	return model;
}


int predict_seq(std::string query, RONNModel *model,
				std::vector<double> &scores)
{
	std::vector<short> seqAA;

	for(int i = 0; i < query.size(); i++)
	{
		seqAA.push_back(INDEX[(int)(query[i]-'A')]);
		// check if valid character
		if(seqAA[i]<0 || seqAA[i]>19)
		{
			printf("seqAA[%d]=%d\n", i, seqAA[i]);
			return(1);
		}
	}

	// printf("starting detect\n");
	detect(seqAA, scores, model);
	// printf("ending detect\n");

	return 0;
}
