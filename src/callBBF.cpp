/*

The script returns the prediction for a sequence under testing using one trained cross-validation model
** Modified by Shyam Saladi (saladi@caltech.edu, California Institute of Technology)

*/
#include "callBBF.h"


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
			getline(fp, tmp);
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
			seqdata.push_back(filedata[j]);
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


int predict_model(std::string query, RONNModel *model,
				  std::vector<double> &scores)
{
	int rho[2], seqlen = scores.size();
	double y, fOrder, fDisor, pDisor;

	// initialize arrays
	std::vector<int> predictTimes(seqlen, 0);
	std::vector<short> seqAA(seqlen, 0);

	std::map < std::map <int, int>, double > yBar;

	//translate protein sequence to numbers
	for(int i = 0; i < seqlen; i++)
	{
		seqAA[i] = INDEX[(int)(query[i]-'A')];
		// check if valid character
		if(seqAA[i]<0 || seqAA[i]>19)
		{
			printf("seqAA[%d]=%d\n", i, seqAA[i]);
			return 1;
		}
	}

	for(int i = 0; i <= seqlen-model->nW; i++)
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

	for(int i = 0; i < seqlen; i++)
	{
		y = 0.0;
		for(int r = 0; r < predictTimes[i]; r++)
		{
			std::map<int, int> tmp_vals;
			tmp_vals[i] = r;
			y += yBar[tmp_vals];
		}

		// add score to running total
		scores[i] += (double)y/(double)predictTimes[i];
	}

	return 0;
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
