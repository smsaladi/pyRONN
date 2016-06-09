/*

The critical parameters for this script are
    Maximum bases is 2000 (parameter name is "maxD")
    Maximum residues is 10000 (parameter name is "maxP")

The script returns the prediction for a sequence under testing using one trained cross-validation model
** Modified by Shyam Saladi (saladi@caltech.edu, California Institute of Technology), 01-09-15

*/


#include "callBBF.h"

using namespace std;

string query;
vector<double> w;
double		mu[2];
double		sigma[2];
double		disorder_weight;

map < map <int, int>, double > yBar;

vector< vector<short> >	dbAA;
vector<short> seqAA;
vector<int>	Length;
vector<int> predictTimes;

int		nD;//number of database sequences
int		nW;
int		rho0;
int		rho1;
//FILE		*fp;
//                         A b C D E F G H I j K L M  N  o P  Q  R  S  T  u V  W  x Y
int		INDEX[25]={0,0,1,2,3,4,5,6,7,0,8,9,10,11,0,12,13,14,15,16,0,17,18,0,19};
int		Dayhoff[20][20]={
{ 40, 24, 32, 32, 16, 36, 28, 28, 28, 24, 28, 32, 36, 32, 24, 36, 36, 32,  8, 20},
{ 24, 80, 12, 12, 16, 20, 20, 24, 12,  8, 12, 16, 20, 12, 16, 32, 24, 24,  0, 32},
{ 32, 12, 48, 44,  8, 36, 36, 24, 32, 16, 20, 40, 28, 40, 28, 32, 32, 24,  4, 16},
{ 32, 12, 44, 48, 12, 32, 36, 24, 32, 20, 24, 36, 28, 40, 28, 32, 32, 24,  4, 16},
{ 16, 16,  8, 12, 68, 12, 24, 36, 12, 40, 32, 16, 12, 12, 16, 20, 20, 28, 32, 60},
{ 36, 20, 36, 32, 12, 52, 24, 20, 24, 16, 20, 32, 28, 28, 20, 36, 32, 28,  4, 12},
{ 28, 20, 36, 36, 24, 24, 56, 24, 32, 24, 24, 40, 32, 44, 40, 28, 28, 24, 20, 32},
{ 28, 24, 24, 24, 36, 20, 24, 52, 24, 40, 40, 24, 24, 24, 24, 28, 32, 48, 12, 28},
{ 28, 12, 32, 32, 12, 24, 32, 24, 52, 20, 32, 36, 28, 36, 44, 32, 32, 24, 20, 16},
{ 24,  8, 16, 20, 40, 16, 24, 40, 20, 56, 48, 20, 20, 24, 20, 20, 24, 40, 24, 28},
{ 28, 12, 20, 24, 32, 20, 24, 40, 32, 48, 56, 24, 24, 28, 32, 24, 28, 40, 16, 24},
{ 32, 16, 40, 36, 16, 32, 40, 24, 36, 20, 24, 40, 28, 36, 32, 36, 32, 24, 16, 24},
{ 36, 20, 28, 28, 12, 28, 32, 24, 28, 20, 24, 28, 56, 32, 32, 36, 32, 28,  8, 12},
{ 32, 12, 40, 40, 12, 28, 44, 24, 36, 24, 28, 36, 32, 48, 36, 28, 28, 24, 12, 16},
{ 24, 16, 28, 28, 16, 20, 40, 24, 44, 20, 32, 32, 32, 36, 56, 32, 28, 24, 40, 16},
{ 36, 32, 32, 32, 20, 36, 28, 28, 32, 20, 24, 36, 36, 28, 32, 40, 36, 28, 24, 20},
{ 36, 24, 32, 32, 20, 32, 28, 32, 32, 24, 28, 32, 32, 28, 28, 36, 44, 32, 12, 20},
{ 32, 24, 24, 24, 28, 28, 24, 48, 24, 40, 40, 24, 28, 24, 24, 28, 32, 48,  8, 24},
{  8,  0,  4,  4, 32,  4, 20, 12, 20, 24, 16, 16,  8, 12, 40, 24, 12,  8,100, 32},
{ 20, 32, 16, 16, 60, 12, 32, 28, 16, 28, 24, 24, 12, 16, 16, 20, 20, 24, 32, 72}};
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


void model(string filename)
{
	string str, tmp, tmp2;
	stringstream *convert;
	vector<string> filedata;
	
	ifstream fp(filename.c_str());
	
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
		cout << "No model found." << endl;
	}
	
	//first deal with headers
	//the first two lines in vector give database size (nD) and window length (nW)
	
	convert = new stringstream(filedata[0]);
	(*convert) >> nD;
	
	convert = new stringstream(filedata[1]);
	(*convert) >> nW;
	
	vector<string> seqdata;
	//j for database sequence index

	for(int j=2; j<(nD*2)+2; j++)
	{	
		if (j % 2 == 0)
		{
			seqdata.push_back(filedata[j]);
		}
		else
		{
			convert = new stringstream(filedata[j]);
			double wtmp = 0.0;
			(*convert) >> wtmp;
			w.push_back(wtmp);
		}
	}
	
	filedata.clear();
	
	for (int j = 0; j < nD; j++)
	{
		str = seqdata[j];
		Length.push_back(str.length());
		
		dbAA.push_back(vector<short>());
	
		for (int r = 0; r < str.length(); r++)
		{
			dbAA[j].push_back(INDEX[(int)(str[r] - 'A')]);
			
		}

		
	}
	
	delete convert;
	
	 // Length[j]=str.length();
	  
	  //for(int r=0;r<strlen(str);r++)
	//	dbAA[j][r]=INDEX[(int)(str[r]-'A')];
		

	
}

void pdf(string filename)
{

	string tmp;
	vector<string> filedata;
	stringstream *convert;
	
	ifstream fp(filename.c_str());
	
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
		cout << "No PDF record!" << endl;
	}

	
	
	//deal with the four lines sequentially
	
	convert = new stringstream(filedata[0]);
	(*convert) >> mu[0];
	
	convert = new stringstream(filedata[1]);
	(*convert) >> mu[1];
	
	convert = new stringstream(filedata[2]);
	(*convert) >> sigma[0];
	
	convert = new stringstream(filedata[3]);
	(*convert) >> sigma[1];

	//cout << mu[0] << " " << mu[1] << " " << sigma[0] << " " << sigma[1] << endl;

	//int tmp2 = fscanf(fp,"%f %f %f %f",&mu[0],&mu[1],&sigma[0],&sigma[1]);

	delete convert;
}


void align(int i,int j)//i for query sequence index and j for database sequence index
{
	int r,w,R,score;
	rho1=-1000000;

	
	for(r=0; r<=Length[j]-nW; r++)//go though the database sequence for maximised alignment
	{
	  
	 
	  score=0;

	  for(w=0; w<nW; w++)//go through the query sequence for one alignment
	  	score+=Blosum62[seqAA[i+w]][dbAA[j][r+w]];
	  
	  
	  if(score>rho1)
	  {
	  	rho1=score; 
	  	R=r; 
	  }
	  

	}

	
	
	rho0=0;
	for(w=0;w<nW;w++)
		rho0+=Blosum62[dbAA[j][R+w]][dbAA[j][R+w]];	 
	
}



void detect(vector<double> & estimate)
{
	
	double y, fOrder, fDisor, pOrder, pDisor;

	//string estimate_filename = "estimate.rec";
	
	for(int i = 0; i < query.length(); i++) 
	{
		predictTimes.push_back(0);
	}
	
	
	for(int i = 0; i <= query.length()-nW; i++)
	{
		y=0.0;
		
		for(int j = 0; j<nD; j++)
		{
			
			align(i,j);//search for the maximum alignment between ith peptide from the query and the jth database sequence
	
			y+=w[j]*exp((double)(rho1-rho0)/(double)rho0);
			
		}
		

			
		fOrder=exp(-0.5*pow(y-mu[0],2.0)/sigma[0])/(sqrt(M_2_PI*sigma[0])); //$VR$: bug fixed by Ron in Feb07
		fDisor=exp(-0.5*pow(y-mu[1],2.0)/sigma[1])/(sqrt(M_2_PI*sigma[1])); //$VR$: bug fixed by Ron in Feb07

		//$VR$: old buggy version pre-Feb07
		//fOrder=exp(-0.5*pow(y-mu[0],2.0)/sigma[0])/(sqrt(6.28)*sigma[0]);
		//fDisor=exp(-0.5*pow(y-mu[1],2.0)/sigma[1])/(sqrt(6.28)*sigma[1]);


		//pOrder=(1.0-disorder_weight)*fOrder/((1.0-disorder_weight)*fOrder+disorder_weight*fDisor);
		pDisor=disorder_weight*fDisor/((1.0-disorder_weight)*fOrder+disorder_weight*fDisor);
		//fprintf(fp,"%c\t%f\t%f\t%f\t%f\t%f\n",query[i],y,mu[0],mu[1],pOrder,pDisor);
		
		
		for(int r=i;r<i+nW;r++)
		{
			map <int, int> r_predtimes;
			r_predtimes[r] = predictTimes[r];
			
			yBar[r_predtimes] = pDisor;
			//yBar[r][predictTimes[r]]=pDisor;
			predictTimes[r]++;
		}
	}
	
	/*
	fp=fopen("estimate.rec","w");
	
	for(i=0;i<strlen(query);i++)
	{
		y=0.0;
		for(r=0;r<predictTimes[i];r++)
		y+=yBar[i][r];
		fprintf(fp,"%c\t%f\n",query[i],y/(double)predictTimes[i]);
	}
	
	fclose(fp);
	
	*/
		

	for(int i = 0; i < query.length(); i++)
	{
		y=0.0;
		
		for(int r = 0; r < predictTimes[i]; r++)
		{
			map<int, int> tmp_vals;
			tmp_vals[i] = r;
			
			//y+=yBar[i][r];
			
			y += yBar[tmp_vals];
			
		}
		
		//fprintf(fp,"%c\t%f\n",query[i],y/(double)predictTimes[i]);
		estimate.push_back((double)y/(double)predictTimes[i]);
	}

}



int callBBF_driver(string myquery, string mod_fn, string pdf_fn1, double d_weight, vector<double> & scores)
{

	int i;
	
	disorder_weight = d_weight;
	query = myquery;

	for(i = 0; i < query.size(); i++)
	{
		seqAA.push_back(INDEX[(int)(query[i]-'A')]);
		if(seqAA[i]<0 || seqAA[i]>19) 
		{ 
			printf("seqAA[%d]=%d\n",i,seqAA[i]); 
			return(1); 
		}
	}
	
	model(mod_fn);
	//printf("Model Worked\n");
	//fflush(stdout);
	
	
	pdf(pdf_fn1);
	//printf("PDF Worked\n");
	//fflush(stdout);
	
	detect(scores);
	//printf("Detect Worked\n");
	//fflush(stdout);
		
	dbAA.clear();
	Length.clear();
	predictTimes.clear();
	w.clear();
	
	return 0;

}

/*
int main(int argc,char **argv)
{
int i;
if(argc<4)
  { printf("Command line is <callBBF sequence_fn model_fn pdf_fn disorder_weight>\n"); return(1); }
if(argc==5) sscanf(argv[4],"%f",&disorder_weight);
if(!(fp=fopen(argv[1],"r"))) { printf("Can't open %s\n",argv[1]); return(1); }
int tmp1 = fscanf(fp,"%s",query);
fclose(fp);
//printf("callBBF: %s\n",query);
for(i=0;i<strlen(query);i++)
  {
  seqAA[i]=INDEX[(int)(query[i]-'A')];
  if(seqAA[i]<0 || seqAA[i]>19) { printf("seqAA[%d]=%d\n",i,seqAA[i]); exit(1); }
  }
model(argv[2]);
printf("Model Worked\n");
fflush(stdout);
pdf(argv[3]);
printf("PDF Worked\n");
fflush(stdout);
detect();
printf("Detect Worked\n");
fflush(stdout);
return(1);
}
*/


