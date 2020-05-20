/********************************/
// Function "get_crs_array" is used to import the necessary database into the object of "crs_coeff" class
/********************************/
#include "main.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "constants.h"

using namespace std;

void get_crs_array(crs_coeff &crs)
{
	
    cout << "initializing s2s arrays...\n";

	FILE *fp;
	char line[500];
	double buf1,buf2,buf3,buf4;
	
    ifstream inputFile;
    string inputLine;
    size_t MAX_LEN = 200;
    char filename[MAX_LEN];
    char dir[MAX_LEN];
    
    double sig_factor;
    double sigD_factor;
    double sigTOT_factor;
    
	//assigning the location of the database based on the path in the file "qct_datbase.input"  
    snprintf(filename, sizeof(filename), "./qct_database.input");
    inputFile.open(filename); 
	if(!inputFile)
	  cout << "there's something wrong with qct_database.input\n"; //Throw exception 
    
    while (!inputFile.eof())
    {
		getline(inputFile, inputLine);
		char cstr[inputLine.size()+1];
		strcpy(cstr,inputLine.c_str()); 
		if (strcmp(cstr,"#directory of qct database") == 0)
		{	
			getline(inputFile, inputLine);
			strcpy(dir,inputLine.c_str()); 
		}
		 
     }
    inputFile.close(); 
    
    //read in database information from the file "info.dat"
    char dirInfo[MAX_LEN];
    strcpy(dirInfo,dir);
	strcat(dirInfo,"info.dat");
    
    inputFile.open(dirInfo);
	if(!inputFile)
	  cout << "there's something wrong with info.input\n";
  
     while (!inputFile.eof())
    {
		getline(inputFile, inputLine);    
		char cstr[inputLine.size()+1];
		strcpy(cstr,inputLine.c_str()); 
		if (strcmp(cstr,"#number of rovibrational energy levels") == 0)
		  {
			getline(inputFile, inputLine);
		    crs.nint = stod(inputLine);
		  }
		if (strcmp(cstr,"#dissociation energy") == 0)
		  {
			getline(inputFile, inputLine);
		    crs.Ediss = stof(inputLine)*CM2J;
		  }
		  if (strcmp(cstr,"#number of coefficients of curve-fitting for total collisions") == 0)
		  {
			getline(inputFile, inputLine);
		    crs.ncoeff_TOT = stod(inputLine);
		  }
		  if (strcmp(cstr,"#number of coefficients of curve-fitting for elastic/dissociation collisions") == 0)
		  {
			getline(inputFile, inputLine);
		    crs.ncoeff = stod(inputLine);
		  }
		  if (strcmp(cstr,"#number of coefficients of curve-fitting for RVT/exchange collisions") == 0)
		  {
			getline(inputFile, inputLine);
		    crs.ncoeff_RVT = stod(inputLine);
		  }
		  
		  if (strcmp(cstr,"#type of database") == 0)
		  {
			  
			getline(inputFile, inputLine);
     		strcpy(crs.database_type,inputLine.c_str()); 
		  }
		  
		  if (strcmp(cstr,"#low/high threshold energy") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.Et_th = stod(inputLine);
		  }
		  if (strcmp(cstr,"#minimum energy") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.Et_min = stod(inputLine);
		  }
		  if (strcmp(cstr,"#maximum energy") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.Et_max = stod(inputLine);
		  }
		  if (strcmp(cstr,"#number of temperatures") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.nTemp = stod(inputLine);
		  }
		  if (strcmp(cstr,"#species of first particle") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.spec1 = stod(inputLine);
		  }
		  if (strcmp(cstr,"#species of second particle") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.spec2 = stod(inputLine);
		  }
		  
		   if (strcmp(cstr,"#coefficient format") == 0)
		  { 
			getline(inputFile, inputLine);
		    strcpy(crs.coeff_format,inputLine.c_str());
		  }
		   if (strcmp(cstr,"#coefficient format of total cross-sections") == 0)
		  { 
			getline(inputFile, inputLine);
		    strcpy(crs.totcrs_coeff_format,inputLine.c_str());
		  }
		  if (strcmp(cstr,"#number of other states") == 0)
		  { 
			getline(inputFile, inputLine);
		    crs.num_os = stod(inputLine);
		  }
    
     }
 
    inputFile.close(); 
    
    // set up the correction factor of RVT transition and dissociation

    if (strcmp(crs.coeff_format,"8 coefficient poly") == 0)
    {
		sig_factor = 1.0e-20/27.0;
		sigD_factor = 1.0e-20/27.0*16.0/3.0;
		}
	else if (strcmp(crs.coeff_format,"new 8 coefficient poly") == 0)
    {
		sig_factor = 1.0e-20;
		sigD_factor = 1.0e-20;
		}
		else if (strcmp(crs.coeff_format,"Kim and Boyd") == 0)
    {
		sig_factor = 1.0;
		sigD_factor =1.0;
		}
		else if (strcmp(crs.coeff_format,"original Kim and Boyd") == 0)
    {
		sig_factor = 1.0;
		sigD_factor =1.0;
		}
		else if (strcmp(crs.coeff_format,"TJPAN") == 0)
    {
		sig_factor = 1.0e-20;
		sigD_factor =1.0e-20;
		}
		else
		{
			mcexit("No such coefficient format");
			}
			
			
	if (strcmp(crs.totcrs_coeff_format,"8 coefficient poly") == 0)
    {
		sigTOT_factor = 1.0e-20/27.0;
		
		}
	else if (strcmp(crs.totcrs_coeff_format,"new 8 coefficient poly") == 0)
    {
		sigTOT_factor = 1.0e-20;
		
		}
		else if (strcmp(crs.totcrs_coeff_format,"Kim and Boyd") == 0)
    {
		sigTOT_factor = 1.0;
		
		}
		else if (strcmp(crs.totcrs_coeff_format,"original Kim and Boyd") == 0)
    {
		sigTOT_factor = 1.0;
		
		}
		else if (strcmp(crs.totcrs_coeff_format,"TJPAN") == 0)
    {
		sigTOT_factor = 1.0e-20;
		
		}
		else
		{
			mcexit("No such coefficient format");
			}		
			
    
    //Allocate the pointers of the object
   
    int nint = crs.nint;
    int ncoeff = crs.ncoeff;
    int ncoeff_RVT = crs.ncoeff_RVT;
    int ncoeff_TOT =crs.ncoeff_TOT;
    int nTemp = crs.nTemp;
    
    crs.deg_array = new double [nint];
    crs.Ei_array = new double [nint];
    crs.Etr_array = new double [ncoeff];
    crs.Ev_array = new double [nint];
    crs.Er_array = new double [nint];
    crs.Temp_array = new double [nTemp];
    
    crs.coeff_RVT_array = new double **[nint];
    crs.coeff_D_array = new double *[nint];    
    crs.coeff_TOT_array = new double *[nint]; 
    
    for (int i=0;i<nint;i++){
		crs.coeff_RVT_array[i]=new double *[nint];
		crs.coeff_D_array[i] = new double [ncoeff];
		crs.coeff_TOT_array[i] = new double [ncoeff_TOT];
		
		for (int j =0;j<nint; j++){
			crs.coeff_RVT_array[i][j]=new double [ncoeff_RVT];
			}
	
		}
	
  
    if (strcmp(crs.database_type,"VT")!= 0)
    {

     // read in elastic cross-sections and degeneracy
	
	char dirErvList[MAX_LEN];
	char dirEtrList[MAX_LEN];
	char dirCrs[MAX_LEN];
	char dirCrsTot[MAX_LEN];
	
	strcpy(dirErvList,dir);
	strcat(dirErvList,"ErvList.dat");
	
	strcpy(dirEtrList,dir);
	strcat(dirEtrList,"EtrList.dat");
	
	strcpy(dirCrs,dir);
	strcat(dirCrs,"crs%s.dat");
	
	strcpy(dirCrsTot,dir);
	strcat(dirCrsTot,"total_crs_fit.dat");

	/*read in energy levels and degeneracy*/
	
	fp = fopen(dirErvList,"r");
    if(fp==NULL)
		mcexit("File [ErvList.dat] not found");
  
    for (int iint=0;iint<nint;iint++)
	{
	  
	  fgets(line,500,fp); 
	  sscanf(line,"%le %le %le %le",&buf1,&buf2,&buf3,&buf4);
	  
	  crs.deg_array[iint] =buf1;
	  crs.Ei_array[iint]=buf2 * CM2J;;
	  crs.Ev_array[iint]=buf3 * CM2J;;
	  crs.Er_array[iint]=buf4 * CM2J;;
	  
	  
	}
    fclose(fp);

  /*read in energy levels and degeneracy*/
	
	fp = fopen(dirEtrList,"r");

    if(fp==NULL)
		mcexit("File [EtrList.dat] not found");
  
    for (int iint=0;iint<ncoeff;iint++)
	{
	  
	  fgets(line,500,fp); 
	  sscanf(line,"%le",&buf1);
	  
	  crs.Etr_array[iint] =buf1;
	  	  
	}
    fclose(fp);

  	// read in other cross-sections 
	int icoeff,iPost,Block;
    
#pragma omp parallel for private(filename,inputFile,inputLine,icoeff,iPost,Block)
    
    for (int iint =0;iint<nint;iint++)
      {
	
	int iPre = iint;
	char groupid[10];
	
	sprintf(groupid,"%04d",iint);
	
	snprintf(filename, sizeof(filename), dirCrs, groupid);
	inputFile.open(filename); 

	if(!inputFile)
	  cout << "there's something wrong with the crs file!\n"; 

	Block =0;
	while (!inputFile.eof())
	{

	    getline(inputFile, inputLine);
	    
	    stringstream ss(inputLine);	    
				
		if (inputLine[0]=='#')
		  {
		    Block ++;
			
		    continue;
		  }
		
		icoeff =0;
		if (Block ==2)
		  {
			getline(ss,inputLine,'\t');
			while(!ss.eof())
			{
				crs.coeff_D_array[iint][icoeff] = stod(inputLine)*sigD_factor;		
				icoeff ++;
				getline(ss,inputLine,'\t');
			}
		 
		  }
		    
		  else 
		  	{  
				getline(ss,inputLine,'\t');
		    	string ipost_str = inputLine;
		    	getline(ss,inputLine,'\t');
				while(!ss.eof())
				{
					iPost = int(stod(ipost_str));	
					crs.coeff_RVT_array[iPre][iPost][icoeff] = stod(inputLine)* sig_factor;
					icoeff ++;
					getline(ss,inputLine,'\t');
				}
			
			}	
		    
		    
		    
		    
		    
		}	    
		    
		    
		  inputFile.close(); 
	}		
	
	if (crs.nTemp > 0)
		{
			
			char dirTempList[MAX_LEN];
	
			strcpy(dirTempList,dir);
			strcat(dirTempList,"tempList.dat");
			
			fp = fopen(dirTempList,"r");
	
    if(fp==NULL)
		mcexit("File [TempList.dat] not found");
  
    for (int iTemp=0;iTemp<nTemp;iTemp++)
	{
	  
	  fgets(line,500,fp); 
	  sscanf(line,"%le",&buf1);
	  
	  crs.Temp_array[iTemp] =buf1;
	  	  
	}
    fclose(fp);
			
			
			
			}
		
		
	}

	if (LOG)
	{
		
		printf("directory of qct database : %s\n",dir);
		printf("type of database : %s\n",crs.database_type);
		printf("number of rovibrational energy levels : %d\n",crs.nint);		
		printf("number of temperatures : %d\n",crs.nTemp);
	
	}
	
}