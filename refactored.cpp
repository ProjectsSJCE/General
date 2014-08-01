//Version 21: Greedy algorithm variations     last modified: 7-9-2014
//NOTE: THIS VERSION ASSUMES INPUT FILES ARE SORTED IN DESCENDING ORDER

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <iomanip>
#include <string>
#include <cfloat>
#include <omp.h>
using namespace std;

void bubble_sort(int *randarray, int **order, int NUMclient, int n)
{
    int j, k, tempint;
    for (j=0;j<NUMclient;j++)
			for (k=1;k<(NUMclient-j);k++)
				if(randarray[k-1]<randarray[k])
				{
					tempint = randarray[k-1];
					randarray[k-1] = randarray[k];
					randarray[k] =tempint;
					tempint = order[n][k-1];
					order[n][k-1]=order[n][k];
					order[n][k] = tempint;
				}
		
}

void setup_firm_cutoff(int firmcut[], int NUMfirm, double *gamma, double g[3][3], int m, double *gammap1)
{
    int i;
	for(i=0;i<NUMfirm;i++)
	{                          
		if(i<firmcut[1])
			gamma[i]= g[0][m];
		else if(i<firmcut[2])
			gamma[i]= g[1][m];
		else
			gamma[i]= g[2][m];
		gammap1[i] = gamma[i]+1.;
	}
}

void init_accumulators(double xtable[3][10], long int priceset[3][3][10])
{    
    int a, b, c;
    for(a=0;a<3;a++)        
	    	for(b=0;b<10;b++)
	    	{
	    		xtable[a][b] = 0.;
	    		for(c=0;c<3;c++)
	    			priceset[c][a][b]=0;
	    	}
}

void write_econout_1(int NUMfirm, int NUMclient, int NUMreps, int firmcut[], double g[3][3], int m, std::ofstream& econout)
{

	econout  << setprecision(0) << endl ;
	econout  << setw(18) << "Firms: "
				<< setw(18) << "Clients: "
				<< setw(18) << "Replications: " << endl;
	econout  << setw(18) << NUMfirm
				<< setw(18) << NUMclient
				<< setw(18) << NUMreps << endl;
	econout	<< setw(18) << "Number of Big: "
				<< setw(18) << "Medium: "
				<< setw(18) << "Small firms: " << endl;
	econout	<< setw(18) << firmcut[1]
				<< setw(18) << (firmcut[2]-firmcut[1])
				<< setw(18) << (firmcut[3]-firmcut[2]) << endl;
	econout  << setprecision(2);
	econout  << setw(18) << "Gamma: Big: "
				<< setw(18) << "Medium: "
				<< setw(18) << "Small firms: " << endl;
	econout  << setw(18) << g[0][m]
				<< setw(18) << g[1][m]
				<< setw(18) << g[2][m] << endl;

}

void write_econout_2(double totalprice, double totalcost, std::ofstream& econout)
{
	econout << setw(18) << totalprice;
	econout << setw(18) << totalcost;
	econout << setw(18) <<(totalprice - totalcost);
	econout << endl ;
}

void print_out_table(int firmcut[], double g[3][3], int m, double xtable[3][10], int table[3][10], std::ofstream& tablout)
{
    int a, b;
    for(a=0;a<3;a++)      
    	{
    		tablout << (firmcut[a+1] - firmcut[a]) << '\t' ;
    		tablout << setw(5) << g[a][m] << '\t';
    		for(b=0;b<10;b++)
    		{
    			xtable[a][b] += double(table[a][b]);
    			cout << setw(6) << table[a][b] ;
    			tablout << setw(6) << table[a][b] ;
    		}
    		cout << endl;
    		tablout << endl;
    	}
    	tablout << endl;
}

void initialise(int NUMclient, double *clientscost, double *clientsprice, int NUMfirm, int **related, double *firmstatus, double *clientsize, double part1, double *gammap1, double **prices, double *clpr_ovr_clsz, double w, double fstogamma[], int start, int end)//verify
{
    int i, j;
    double tempdouble;
    
	for (j=start;j<end;j++)
	{
		clientscost[j]  = 999999999999.;     // initialize for tracking minima
		clientsprice[j] = 999999999999.;
		for (i=0;i<NUMfirm;i++)
		{
			if (related[i][j])
				tempdouble = firmstatus[i] - clientsize[j];
			else
				tempdouble = firmstatus[i];
			part1 = pow(tempdouble,gammap1[i]);
			prices[i][j] = w * ((
					pow((tempdouble + clientsize[j]),gammap1[i])
					- part1) / fstogamma[i]);
			if (prices[i][j] <= clientscost[j])	// note minimum price so far
			{
				 clientsprice[j] = clientscost[j];
				 clientscost[j] = prices[i][j];
			}
			else if (prices[i][j] < clientsprice[j])
				clientsprice[j] = prices[i][j];
		}
		clpr_ovr_clsz[j] = clientsprice[j] / clientsize[j]; // price over client size
	}
}

void SetToZero(double array[], int size)
{
    int i;
    for(i=0; i < size; i++)
    {
        array[i] = 0;
    }
}

void initializeFirmAndRelated(double tclpr_ovr_clsz[], double firmstatus[], int **related, int NUMfirm, int NUMclient)
{
    int i, j;
    for (i=0;i<NUMfirm;i++)
	{
		tclpr_ovr_clsz[i]=0.;
		firmstatus[i] = 0.;             // initialize firm status
		for (j=0;j<NUMclient;j++)       // no one is assigned at start
			related[i][j] = 0;
	}
}

void performComputation(int **related, int **order, double clientsmin[], double firmstatus[], double gammap1[], double **prices, int clientsfirm[], double clientsize[], double fstogamma[], int dupprice[], int n, int NUMclient, int NUMfirm, int start, int end)
{
    int j, i, a, b, tempint;
    double part1, w, tempdouble;
    for (k=start;k<end;k++)
	{
        j = order[n][k];                // according to randomly sorted indices
	    clientsmin[j] = 99999999999.;	// initialize before tracking min.
	    tempint = 0;
	    for (i=0;i<NUMfirm;i++)
	    {
		    part1 = pow(firmstatus[i],gammap1[i]);
		    prices[i][j] = w * ((pow((firmstatus[i] + clientsize[j]),gammap1[i]) - part1) / fstogamma[i]);
		    if (prices[i][j] < clientsmin[j])	// note minimum price so far
		    {
			     clientsmin[j] = prices[i][j];
			     clientsfirm[j] = i;
			     dupprice[tempint] = i;
		    }
		    else if (prices[i][j] == clientsmin[j])  // if new price ties old min
		    {
			    tempint++;
			    dupprice[tempint] = i;
		    }
		    if (tempint)  // randomly select from firms which tie above
		    {
			    tempdouble = double(tempint+1)/RAND_MAX;
			    a = tempint;
			    while (a >= tempint)
				    a = int(rand()*tempdouble);
			    clientsfirm[j] = dupprice[a];
		    }
	    }
	    related[clientsfirm[j]][j]=1;
	    firmstatus[clientsfirm[j]] += clientsize[j];
    }
}

void CalculateLowest(int NUMfirm, int clientsfirm[], int clientrival[], double **prices, double clientsprice[], int start, int end)
{
    int i;
    for (j=start;j<end;j++)	// at this point client is at lowest price firm
	{
        for (i=0;i<NUMfirm;i++)
			if ((prices[i][j] < clientsprice[j])     // prices are lowest of
				&& (clientsfirm[j] != i))
					clientrival[j] = i;
//		CalculateLowest(NUMfirm, clientsfirm, clientrival, prices, clientsprice, j);
	}
}

void ObtainRivalTable(long int priceset[3][3][10], int clientrival[], int firmcut[], int clientsfirm[], double tempdouble, int start, int end)
{
    int a, b, c, j;
	for (j=start;j<end;j++)
	{
		//PERFORM MULITTHREADING HERE
		tempdouble = double(NUMclient)/10.;
        if(clientrival[j] < firmcut[1])	// get nearest rival type
			    a=0;
	    else if(clientrival[j] < firmcut[2])
		    a=1;
	    else
		    a=2;
	    if(clientsfirm[j] < firmcut[1])	// get audit firm type
		    b=0;
	    else if(clientsfirm[j] < firmcut[2])
		    b=1;
	    else
		    b=2;
	    c = int(double(j)/tempdouble);	// get proper client decile
	    priceset[a][b][c]++;				// augment proper counter    
//		ObtainRivalTable(priceset, clientrival, firmcut, clientsfirm, NUMclient, j);
	}

}

int main()
{
	int start_s=clock();

	int NUMclient = 5460;  			//number of clients to create

	int NUMfirm = 100;  			// number of firms to create

	int NUMreps = 30;				// number of replications (for random only)

	int NUMgamma = 3;               // number of gamma sets

	int greedymethod = 2; 			// 0: greedy, 1: inverse, 2: random

	int firmcut[]={0,6,18,100}; 	// cutoffs between firm categories

	double g[3][3]={{1.,.8,.6},{1.,.9,.8},{1.,1.,1.}}; // gamma value

	const char * firmfile = " firmsize.txt ";
	const char * clientfile = " jobsize.txt ";

	int   a, b, c, i, j, k, m, n;   // loop counters, i for firms, j for clients

    int clientcut[11];
    int table[3][10];               // table of distribution
	double xtable[3][10];
	long int priceset[3][3][10];	// table of nearest rival summaries
	const char * tabl_out;
	const char * xtab_out;
	const char * firm_out;
	const char * rivl_out;
	const char * clnt_out;
	const char * econ_out;
	double * firmsize;              // holds firm sizes
	double * fstogamma;             // holds firmsize^gamma
    double * clientsize;            // holds client sizes
	double ** prices;               // holds current prices
    double * firmstatus;            // holds firm status (total size of jobs held)
    double * clientsmin;            // cost used in assigning client
    double * clientscost;           // lowest price available to client
    double * clientsprice;          // next lowest price available to client
    double * clpr_ovr_clsz;         // clients price over clients size
    double * tclpr_ovr_clsz;        // total client price over client size within rep
    double * aclpr_ovr_clsz;        // average client price of client size (tot. across reps)
    double * tpricepfirm;           // total price per firm across reps
    double * tcostpfirm;            // total cost per firm across reps
    double * xtcostpfirm;           // accumulator for taking average
	double * tfirmstatus;           // total of firm status across reps
    double * profit;                // profit per firm
    double * tcpf_ovr_fst;          // total cost per firm over firm status (across reps)
    double * tprmarg;               // total profit margin per firm (across reps)
    double w;                       // price of labor
    double * gamma;                 // alpha / beta
    double * gammap1;               // gamma + 1
    double	totalcost;              // total cost
	double	totalprice;             // total prices paid
	double	tempdouble;
	double part1;
	int	* clientsfirm;              // firm assigned to client
    int * clperfirm;                // clients per firm
    int * clientrival;              // clients nearest rival for price
    int ** related;                 // flag array, 1 = related, 0 = not
    int * randarray;                // array of rand()
    int ** order;                   // array of assignment order
    int ** auditfirm;               // auditor for each client, for each replication
    int * dupprice;                 // holds duplicate firm with same price to break ties
    int tempint;
	time_t t1,t2;

	t1 = time(NULL);
	w = 1.0;                        // constant does not affect outcome
	for (b=0;b<11;b++)              // set up table cutoffs
		clientcut[b] = int(b*(NUMclient/10));

	// dynamic memory allocation
	firmsize = new double[NUMfirm];
	fstogamma = new double[NUMfirm];
	clientsize = new double[NUMclient];
	prices = new double*[NUMfirm];
	for (i=0;i<NUMfirm;i++)
		prices[i] = new double[NUMclient];
	firmstatus = new double[NUMfirm];
	tfirmstatus = new double[NUMfirm];
	clientsmin = new double[NUMclient];
	clientscost = new double[NUMclient];
	clientsprice = new double[NUMclient];
	clientsfirm = new int[NUMclient];
	clpr_ovr_clsz = new double[NUMclient];
	tclpr_ovr_clsz = new double[NUMfirm];
	aclpr_ovr_clsz = new double[NUMfirm];
	gamma = new double[NUMfirm];
	gammap1 = new double[NUMfirm];
	clperfirm = new int[NUMfirm];
	clientrival = new int[NUMclient];
	tpricepfirm = new double[NUMfirm];
	tcostpfirm = new double[NUMfirm];
	xtcostpfirm = new double[NUMfirm];
	profit = new double[NUMfirm];
	tcpf_ovr_fst = new double[NUMfirm];
	tprmarg = new double[NUMfirm];
	randarray = new int[NUMclient];
	dupprice = new int[NUMfirm];
	auditfirm = new int*[NUMreps];
	for (n=0;n<NUMreps;n++)
		auditfirm[n]=new int[NUMclient];
	order = new int*[NUMclient];
	for (n=0;n<NUMreps;n++)
		order[n] = new int[NUMclient];
	related = new int*[NUMfirm];
	for (i=0;i<NUMfirm;i++)
		related[i] = new int[NUMclient];

	// set up files and other initializations
	if (greedymethod==2)
	{
        tabl_out = "rga_tabl.out";
		xtab_out = "rga_xtab.out";
		firm_out = "rga_firm.out";
        rivl_out = "rga_rivl.out";
		clnt_out = "rga_clnt.out";
		econ_out = "rga_econ.out";

		// set up arrays for randomizing order of bin assignment
		for (n=0;n<NUMreps;n++)
			for (j=0;j<NUMclient;j++)
				order[n][j]=j;

		for (n=0;n<NUMreps;n++)
		{
			for (j=0;j<NUMclient;j++)
				randarray[j] = rand();

		// bubble sort array of sorted client indices for replicating runs
        bubble_sort(randarray, order, NUMclient, n);
		
        }
    } // end greedymethod 2 initialization.

	// open,read,close files

	ifstream infil;
	infil.open("firmsize.txt");
	for (i=0;i<NUMfirm;i++)
		infil >> firmsize[i];
	infil.close();

	ifstream infile;
	infile.open("jobsize.txt");
	for (j=0;j<NUMclient;j++)
		infile >> clientsize[j];
	infile.close();
    
	ofstream tablout(tabl_out);      // output stream method calls
	ofstream firmout(firm_out);
	ofstream rivlout(rivl_out);
	ofstream clntout(clnt_out);
	ofstream xtabout(xtab_out);
	ofstream econout(econ_out);
	firmout.setf(ios::showpoint);
	firmout.setf(ios::fixed, ios::floatfield);
	firmout.setf(ios::right, ios::adjustfield);
	xtabout.setf(ios::showpoint);
	xtabout.setf(ios::fixed, ios::floatfield);
	xtabout.setf(ios::right, ios::adjustfield);
	clntout.setf(ios::showpoint);
	clntout.setf(ios::fixed, ios::floatfield);
	clntout.setf(ios::right, ios::adjustfield);
	econout.setf(ios::showpoint);
	econout.setf(ios::fixed, ios::floatfield);
	econout.setf(ios::right, ios::adjustfield);

for (m=0;m<NUMgamma;m++)		// start main loop, m is index of gamma set
{
    setup_firm_cutoff(firmcut, NUMfirm, gamma, g, m, gammap1); // set up gamma by firm cutoff

    init_accumulators(xtable, priceset); //initialize accumulators

    SetToZero(tpricepfirm, NUMfirm);
    SetToZero(xtcostpfirm, NUMfirm);
	SetToZero(tfirmstatus, NUMfirm);
	SetToZero(aclpr_ovr_clsz, NUMfirm);
	SetToZero(tprmarg, NUMfirm);

    write_econout_1(NUMfirm, NUMclient, NUMreps, firmcut, g, m, econout); 

	for (i=0;i<NUMfirm;i++)
		fstogamma[i] = pow(firmsize[i],gamma[i]);  // for computation later

 for (n=0;n<NUMreps;n++)
 {
	
    initializeFirmAndRelated(tclpr_ovr_clsz, firmstatus, related, NUMfirm, NUMclient);

    //PERFORM MULITTHREADING HERE ----------------------------------------------------------------------------------------------------
    performComputation(related, order, clientsmin, firmstatus, gammap1, prices, clientsfirm, clientsize, fstogamma, dupprice, n, NUMclient, NUMfirm, k);        
//-----------------------------------------------------------------------------------------------------------------    
    
	// tally distribution of clients
	for(a=0;a<3;a++)
		for(b=0;b<10;b++)
			table[a][b] = 0;

	for(a=0;a<3;a++)
		for(b=0;b<10;b++)
			for(i=firmcut[a];i<firmcut[a+1];i++)
				for(j=clientcut[b];j<clientcut[b+1];j++)
					table[a][b] += related[i][j];


    print_out_table(firmcut, g, m, xtable, table, tablout);

	totalcost = 0.;
	for (j=0;j<NUMclient;j++) // these are prices developed en route
		totalcost += clientsmin[j];
	cout << "accumulated cost during assignment: " << totalcost << endl;
    
//    PERFORM MULITTHREADING ------------------------------------------------------------------------------------------------
    initialise(NUMclient, clientscost, clientsprice, NUMfirm, related, firmstatus, clientsize, part1, gammap1, prices, clpr_ovr_clsz, w, fstogamma, 0, NUMclient);
//----------------------------------------------------------------------------------------------------------------------------

    // PERFORM MULITTHREADING HERE	-------------------------------------------------------------------------------------------------
	CalculateLowest(NUMfirm, clientsfirm, clientrival, prices, clientsprice, 0, NUMclient);
//---------------------------------------------------------------------------------------------------------

    /// YOU WILL HAVE TO MULTITHREAD HERE ---------------------------------------------------------------------------------------------------------
	tempdouble = double(NUMclient)/10.;  // for rivlout output, nearest rival table
	ObtainRivalTable(priceset, clientrival, firmcut, clientsfirm, tempdouble, 0, NUMclient);
//---------------------------------------------------------------------------------------------------------

	for (i=0;i<NUMfirm;i++)
	{
		clperfirm[i] = 0;
		for (j=0;j<NUMclient;j++)       // number of clients per firm
			clperfirm[i] += related[i][j];
		tfirmstatus[i] += firmstatus[i];	// accumulate for average firm status
		tempdouble = firmstatus[i];
		tcostpfirm[i]=                  // get each firm's cost
			(tempdouble<=0.)? 0.: (pow(tempdouble,gammap1[i]) / fstogamma[i]);
		xtcostpfirm[i] += tcostpfirm[i];  // accumulate over n reps
	}

	totalprice = 0.;


	for (j=0;j<NUMclient;j++)
	{             // accumulate totals and firm totals
		totalprice += clientsprice[j];
		tpricepfirm[clientsfirm[j]] += clientsprice[j];
		tclpr_ovr_clsz[clientsfirm[j]] += clpr_ovr_clsz[j];
		auditfirm[n][j] = clientsfirm[j];     // accumulate for clntout
	}

	for (i=0;i<NUMfirm;i++)
	{
				// accumulate average per replication
		profit[i] =  tpricepfirm[i] - xtcostpfirm[i]  ;
		aclpr_ovr_clsz[i] += clperfirm[i] ? (tclpr_ovr_clsz[i]/clperfirm[i]) : 0.;
		tcpf_ovr_fst[i] += (firmstatus[i]>0.) ? (tcostpfirm[i]/firmstatus[i]) : 0.;
		tprmarg[i] += (tpricepfirm[i]>0.) ? (profit[i] / tpricepfirm[i]) : 0.;
	}

    
	if (n==0)
	{
		econout << setw(18) << "total price: " ;
		econout << setw(18) << "total cost: " ;
		econout << setw(18) << "total profit: " ;
		econout << endl ;
	}
    write_econout_2(totalprice, totalcost, econout);

	if (greedymethod==2)
		cout<<"completed replication "<<(n+1)<<" of "<<NUMreps
			 <<" for gamma set " << (m+1) <<endl;

 } // end replication loop (n)

	// output on firms:
	firmout  << endl;
	firmout  << setw(6) << "size"
				<< setw(4) << "idx"
				<< setw(7) << "A(p/s)"
				<< setw(6) << "c/js"
				<< setw(6) << "p/ptn"
				<< setw(6) << "load"
				<< setw(4) << "gam"
				<< setw(10) << "price"
				<< setw(10) << "cost"
				<< setw(10) << "profit"
				<< setw(6) << "p_mr"
				<< endl;
	tempdouble = double(NUMreps);
	for (i=0;i<NUMfirm;i++)
	{
		firmout << setw(6) << setprecision(0) << firmsize[i]
				  << setw(4) << i
				  << setprecision(2)
				  << setw(7) << (aclpr_ovr_clsz[i] / tempdouble)
				  << setw(6)
				  << (tcpf_ovr_fst[i]/tempdouble)
				  << setw(6) << (profit[i] / (firmsize[i]*tempdouble))
				  << setw(6) << (tfirmstatus[i] / (firmsize[i]*tempdouble))
				  << setw(6) << gamma[i]
				  << setw(10) << (tpricepfirm[i]/tempdouble)
				  << setw(10) << (xtcostpfirm[i]/tempdouble)
				  << setw(10) << (profit[i]/tempdouble)
				  << setw(6) << (tprmarg[i]/tempdouble)
				  << endl;
	}
	xtabout << endl << "Per type basis: "<< endl;
	tempdouble = double(NUMreps);
	for(a=0;a<3;a++)      // print out table
	{
		xtabout << setprecision(0);
		xtabout << (firmcut[a+1] - firmcut[a]) << '\t' ;
		xtabout << setprecision(2);
		xtabout << setw(5) << g[a][m] << '\t';
		for(b=0;b<10;b++)
		{
			xtable[a][b] /= tempdouble;
			xtabout << setw(8) << xtable[a][b] ;
		}
		xtabout << endl;
	}
	xtabout << endl << "Per firm basis: " << endl;
	for(a=0;a<3;a++)
	{
		xtabout << setprecision(0);
		xtabout << (firmcut[a+1] - firmcut[a]) << '\t' ;
		xtabout << setprecision(2);
		xtabout << setw(5) << g[a][m] << '\t';
		tempdouble = double(firmcut[a+1]-firmcut[a]);
		for(b=0;b<10;b++)
		{
			xtabout << setw(8) << (xtable[a][b] / tempdouble);
		}
		xtabout << endl;
	}

	for(a=0;a<3;a++)
	{
		rivlout << endl;
		for(b=0;b<3;b++)
			for(c=0;c<10;c++)
				rivlout << setw(7) << priceset[a][b][c];
	}
	rivlout << endl;

	for(j=0;j<NUMclient;j++)
	{
		clntout << setprecision(2)
			<< setw(7) << clientsize[j]
			<< setprecision(0);
		for(n=0;n<NUMreps;n++)
			clntout << setw(4) << auditfirm[n][j];
		clntout << endl;
	}
	clntout << endl;
int stop_s=clock();
cout<<"time: "<<(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000<<endl;
return 0;

} // end main loop (m)()()


	t2 = time(NULL);
	econout	<< setw(18) << "run time (m:s)" << endl
				<< int((t2-t1)/60) << ":" << ((t2-t1)%60) << endl;
	econout.close();
}

