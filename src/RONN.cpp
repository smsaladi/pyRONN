/*=========================================================
**
** This program is designed for running RONNv3.2 for detecting
** disordered residues in proteins.
** Date:	March 2012
** Designer:	Zheng Rong Yang (Exeter University)
** Programmer:	Varun Ramraj (University of Oxford)
** Contact:	varun@strubi.ox.ac.uk
**
==========================================================*/


#include "callBBF.h"

using namespace std;

//@@@@@ the path can be changed
string myPath="/usr/local/RONNv3_2";


int main(int argc,char **argv)
{
	//number of models, characters and the format
	int nM=10;
	int nChar,format;
	
	stringstream *convert;
	FILE *fp;
	
	string qf = "query.dat";

	//$VR$: VERY IMPORTANT. CONTROLS COST FUNCTION
	double disorder_weight=0.53;
	
	string fModel, fPDF, str;
	int r, m, nR;
	
	
	vector< vector<double> > X;
	vector<double> mean;

	string query;

	ifstream ifile;
	string testoutputstring;

	//argument check
	if(argc!=3) 
	{ 
		printf("RONN filename fileFormat(0=plain text;1=FASTA)\n"); 
		return(1); 
	
	}
	
	convert = new stringstream(argv[2]);
	(*convert) >> format;
	
	//sscanf(argv[2],"%d",&format);

	//$VR$: Debug output
	cout << "The cost coefficient or disorder weight is: " << disorder_weight << endl;

	ifstream qfp(argv[1]);
	
	if (qfp.is_open())
	{
		if (format == 0) //plain text
		{
			getline(qfp, query);
		}
		else //fasta
		{
			getline(qfp, query);
			query = ""; 
			getline(qfp, query);
		}
		
		qfp.close();
	}
	else
	{
		cout << "Can't open " << argv[1] << endl;
		return -1;
	}


	nR=query.size();
	
	ofstream fpout(qf.c_str());
	
	if (fpout.is_open())
	{
		fpout << query << endl;
		fpout.close();
	}
	else
	{
		cout << "Can't output query" << endl;
	}
	

	for(m=0; m<nM; m++)
	{
		//sprintf(fModel,"%s/c%d/model.rec",myPath,m);
		convert = new stringstream();
		(*convert) << m;
		string mString = convert->str();
	
		fModel = myPath + "/c" + mString + "/model.rec";
		
	   	//TODO: get rid of FILE* stuff
		if(!(fp=fopen(fModel.c_str(),"r"))) 
		{ 
			cout << "Can't find " << fModel << endl;
			return -1;
		}
		
		fclose(fp);

		//sprintf(fPDF,"%s/c%d/pdfs.rec",myPath,m);
		fPDF = myPath + "/c" + mString + "/pdfs.rec";
		
		//TODO: get rid of FILE* stuff
		if(!(fp=fopen(fPDF.c_str(),"r"))) 
		{ 
			cout << "Can't find " << fPDF << endl; 
			return -1;
		}
		fclose(fp);


		int retval = callBBF_driver(qf.c_str(), fModel, fPDF, disorder_weight);

		if (retval == 0)  
		  printf("................Completing testing ............! Model: %d\n", m);


		//fp=fopen("estimate.rec","r");

		ifstream efp("estimate.rec");
		
		if (efp.is_open())
		{
			
			for(r=0;r<nR;r++)
			{
				//int tmp = fscanf(fp,"%s %f",str,&X[r][m]);
				getline(efp, str);
				string tmpxrm = str.substr(2);
				
				convert = new stringstream(tmpxrm);
				double tmpval = 0.0;
				X.push_back(vector<double>()); 
				
				(*convert) >> tmpval;
				
				X[r].push_back(tmpval);
			}
		
			efp.close();
		}
		else
		{
			cout << "Can't open estimate.rec" << endl;
			return -1;
		}
		

		//int tmp2 = sprintf(str,"mv estimate.rec estimate%d.rec",m);	
		string command = "mv estimate.rec estimate" + mString + ".rec";
		
		if (system(command.c_str()) != 0)
		{
			cerr << "Can't move files" << endl;
			return -1;
		}
		

	}
	  
	fp=fopen("disorder.prb","w");
	
	for(r=0;r<nR;r++)
	{
		mean.push_back(0.0);
		
		for(m=0; m<nM; m++)
			mean[r]+=X[r][m];
		
		mean[r]/=(double)nM;
		
		fprintf(fp,"%c\t%f\n",query[r],mean[r]);
		
		
	}
	  
	  
	fclose(fp);
	
	delete convert;
	
	return 0;
}
