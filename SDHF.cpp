#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <fftw3.h>
#include <complex>

using namespace std;

/* Description: Compute spectral decomposition of heat flux (SDHF) for group of atoms
   Prerequisites: Rerun dumpfile using LAMMPS to obtain interatomic forces, outputting data as vx vy vz fx fy fz
				  FFTW3 library
*/


int main()
{
    string line,fName;
    int id,i,j,k,n,m,t,c,tTime,nAtoms,nAtomsInitial,yStep,currentTimeStep,typ;
    double x,y,z,vx,vy,vz,tau1,tau2,tau3,tau4,tau5,tau6,PE,KE,r,rMag,fx,fy,fz;    
    double xLo, yLo, zLo, xHi, yHi, zHi, Lx, Ly, Lz;
    
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-12;
    double refVelocity = refLength/refTime; // 1 Angstrom / 1 fs
    double refVolume = refLength*refLength*refLength;
    double refPressure = 100000;
    double refEnergy = 1.60218e-19;
    double refForce = refEnergy / refLength;
    
    //******* CHANGE **************************
    
    int tSkip = 10; // how often is atom data dumped? 

    double deltaT = 0.001; // MD timestep\

    double deltaTSI = tSkip*deltaT*refTime;

    #define nTimeSteps 20000 // no of time steps in dump file: use command line: grep -o 'TIMESTEP' <dumpfilename> | wc -l

    #define Nfft 39999 // must be 2*FVCFSIZE - 1
    string dumpIn;

    string dumpOut;

    // takes dumpfile and output file names as string inputs
    cout << "Type input file: ";
    cin >> dumpIn;

    cout << "Type out file: ";
    cin >> dumpOut;

    ifstream data(dumpIn, ios::in); // your dump file

    //*****************************************
           
    ofstream fileOut1(dumpOut,ios::out);

    cout << "Reading in dump file" << endl; 
       
    int nEvery = floor(nTimeSteps/100);
    int counter = 0;

    cout << dumpIn << endl;
    cout << dumpOut << endl;


    vector< vector<double> > velXAll;
    vector< vector<double> > velYAll;
    vector< vector<double> > velZAll;
    vector< vector<double> > forXAll;
    vector< vector<double> > forYAll;
    vector< vector<double> > forZAll;
    vector<double> filter;
    vector<double> atomIDs;
    
    vector<double> v_intersection;
    int f = 0;

	// go through dump file

    for(t=0;t<nTimeSteps;t++)
    {
        counter++;
        
        for(n=1;n<10;n++)
        {
            if(n == 4)
            {
                data >> nAtoms;
            }
            
            if(n == 2)
            {
                data >> currentTimeStep;
 
                if(counter >= nEvery) // progress update
                {
                    counter = 0;
                    
                    cout << "currentTimeStep = " << currentTimeStep
                        << "; t = " << t << " [ " 
                        << 100*float(t+1)/float(nTimeSteps) 
                        << "% ]" << endl;       
                }
            }
            
			if(n == 6)
            {
                data >> xLo >> xHi;
            }

            if(n == 7)
            {
                data >> yLo >> yHi;
            }

            if(n == 8)
            {
                data >> zLo >> zHi;
            }  
                        
            getline(data,line);
        }

        
        // initialise
        if(t == 0)
        {
            Lx = xHi - xLo;
            Ly = yHi - yLo;
            Lz = zHi - zLo;
            atomIDs.resize(nAtoms);	    
            velXAll.resize(nAtoms);
            velYAll.resize(nAtoms);
            velZAll.resize(nAtoms);
            forXAll.resize(nAtoms);
            forYAll.resize(nAtoms);
            forZAll.resize(nAtoms);
        }
        
        f = 0;
		
		// extract data line by line
        for(n=0;n<nAtoms;n++)
        {
            data>>id>>vx>>vy>>vz>>fx>>fy>>fz;

                velXAll[f].push_back(vx);
                velYAll[f].push_back(vy);
                velZAll[f].push_back(vz);
                forXAll[f].push_back(fx);
                forYAll[f].push_back(fy);
                forZAll[f].push_back(fz);
        f = f + 1;
	}
        getline(data,line);
    }


    cout << "Start FVCF calculation" << endl;
	
    nEvery = floor(nAtoms/100);

    counter = 0;    

      
    vector<double> XforceofOneAtomInTime, YforceofOneAtomInTime, ZforceofOneAtomInTime;
    vector<double> XvelocityOneAtomInTime, YvelocityOneAtomInTime, ZvelocityOneAtomInTime;
    XforceofOneAtomInTime.resize(nTimeSteps, 0.0);
    YforceofOneAtomInTime.resize(nTimeSteps, 0.0);
    ZforceofOneAtomInTime.resize(nTimeSteps, 0.0);
    XvelocityOneAtomInTime.resize(nTimeSteps, 0.0);
    YvelocityOneAtomInTime.resize(nTimeSteps, 0.0);
    ZvelocityOneAtomInTime.resize(nTimeSteps, 0.0);

    vector<double> XforceofOneAtomInTimePadded;
    vector<double> YforceofOneAtomInTimePadded;
    vector<double> ZforceofOneAtomInTimePadded;
    vector<double> XvelocityOneAtomInTimePadded;
    vector<double> YvelocityOneAtomInTimePadded;
    vector<double> ZvelocityOneAtomInTimePadded;

   
    static double in_forx[nTimeSteps][2], in_fory[nTimeSteps][2], in_forz[nTimeSteps][2];
    static double in_velx[nTimeSteps][2], in_vely[nTimeSteps][2], in_velz[nTimeSteps][2];
    static double For_wx[nTimeSteps][2], For_wy[nTimeSteps][2], For_wz[nTimeSteps][2];
    static double Vel_wx[nTimeSteps][2], Vel_wy[nTimeSteps][2], Vel_wz[nTimeSteps][2];
    static double CFor_wx[nTimeSteps][2], CFor_wy[nTimeSteps][2], CFor_wz[nTimeSteps][2];
    static double fftfvx[nTimeSteps][2], fftfvy[nTimeSteps][2], fftfvz[nTimeSteps][2];
    static double fvcfx[Nfft][2], fvcfy[Nfft][2], fvcfz[Nfft][2];

    fftw_plan fftforx, fftfory, fftforz;
    fftw_plan fftvelx, fftvely, fftvelz;
    fftw_plan ifftx, iffty, ifftz;


    vector<double> fvcfXAll;
    vector<double> fvcfYAll;
    vector<double> fvcfZAll;
    fvcfXAll.resize(Nfft,0.0);
    fvcfYAll.resize(Nfft,0.0);
    fvcfZAll.resize(Nfft,0.0);

    vector<double> spectrumX;
    vector<double> spectrumY;
    vector<double> spectrumZ;
    spectrumX.resize(nTimeSteps, 0.0);
    spectrumY.resize(nTimeSteps, 0.0);
    spectrumZ.resize(nTimeSteps, 0.0);
    

     for(n=0;n< nAtoms;n++) // pick one atom
    {

        counter++; 
        if(counter >= nEvery) // progress update
        {
            counter = 0;
            
            cout << "atom = " << n << " [" 
                << 100*float(n+1)/float(nAtoms)
                << "% ]" << endl;             
        }



        XforceofOneAtomInTime = forXAll[n];
        YforceofOneAtomInTime = forYAll[n];
        ZforceofOneAtomInTime = forZAll[n];
        XvelocityOneAtomInTime = velXAll[n];
        YvelocityOneAtomInTime = velYAll[n];
        ZvelocityOneAtomInTime = velZAll[n];
		
		// define complex numbers

        for (i = 0; i < nTimeSteps; i++) {
            in_forx[i][0] = XforceofOneAtomInTime[i];
            in_forx[i][1] = 0.000000;
            in_fory[i][0] = YforceofOneAtomInTime[i];
            in_fory[i][1] = 0.000000;
            in_forz[i][0] = ZforceofOneAtomInTime[i];
            in_forz[i][1] = 0.000000;
            in_velx[i][0] = XvelocityOneAtomInTime[i];
            in_velx[i][1] = 0.000000;
            in_vely[i][0] = YvelocityOneAtomInTime[i];
            in_vely[i][1] = 0.000000;
            in_velz[i][0] = ZvelocityOneAtomInTime[i];
            in_velz[i][1] = 0.000000;
        }

		// compute FFTs

        fftforx = fftw_plan_dft_1d(nTimeSteps, in_forx, For_wx, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftforx);
        fftw_destroy_plan(fftforx);
        fftfory = fftw_plan_dft_1d(nTimeSteps, in_fory, For_wy, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftfory);
        fftw_destroy_plan(fftfory);
        fftforz = fftw_plan_dft_1d(nTimeSteps, in_forz, For_wz, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftforz);
        fftw_destroy_plan(fftforz);

        fftvelx = fftw_plan_dft_1d(nTimeSteps, in_velx, Vel_wx, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftvelx);
        fftw_destroy_plan(fftvelx);
        fftvely = fftw_plan_dft_1d(nTimeSteps, in_vely, Vel_wy, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftvely);
        fftw_destroy_plan(fftvely);
        fftvelz = fftw_plan_dft_1d(nTimeSteps, in_velz, Vel_wz, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftvelz);
        fftw_destroy_plan(fftvelz);

		// compute complex conjugates
		
        for (i = 0; i < nTimeSteps; i++) {
            CFor_wx[i][0] = For_wx[i][0];
            CFor_wx[i][1] = -1*For_wx[i][1];
            CFor_wy[i][0] = For_wy[i][0];
            CFor_wy[i][1] = -1*For_wy[i][1];
            CFor_wz[i][0] = For_wz[i][0];
            CFor_wz[i][1] = -1*For_wz[i][1];
        }
		
		
		// multiply complex numbers
        for (i = 0; i < nTimeSteps; i++) {
            fftfvx[i][0] = deltaTSI*(Vel_wx[i][0] * CFor_wx[i][0] - Vel_wx[i][1] * CFor_wx[i][1]);
            fftfvx[i][1] = deltaTSI*(Vel_wx[i][0] * CFor_wx[i][1] + Vel_wx[i][1] * CFor_wx[i][0]);

            fftfvy[i][0] = deltaTSI*(Vel_wy[i][0] * CFor_wy[i][0] - Vel_wy[i][1] * CFor_wy[i][1]);
            fftfvy[i][1] = deltaTSI*(Vel_wy[i][0] * CFor_wy[i][1] + Vel_wy[i][1] * CFor_wy[i][0]);

            fftfvz[i][0] = deltaTSI*(Vel_wz[i][0] * CFor_wz[i][0] - Vel_wz[i][1] * CFor_wz[i][1]);
            fftfvz[i][1] = deltaTSI*(Vel_wz[i][0] * CFor_wz[i][1] + Vel_wz[i][1] * CFor_wz[i][0]);
        }

        
		// extract real part
            for (i = 0; i < nTimeSteps/2; i++) {
                spectrumX[i] = spectrumX[i] + fftfvx[i][0];
                spectrumY[i] = spectrumY[i] + fftfvy[i][0];
                spectrumZ[i] = spectrumZ[i] + fftfvz[i][0];
            }



}

for (k = 0; k < nTimeSteps/2; k++)
{
    fileOut1 << spectrumX[k] << '\t'
        << spectrumY[k] << '\t'
        << spectrumZ[k] << '\t'
        << endl;
}

    cout << "end" << endl;

    fileOut1.close();
    data.close();

    return 0;
}
