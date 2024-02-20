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

/* Description: Compute vibrational density of states (VDOS) for group of atoms
   Prerequisites: LAMMPS dump file out-putting x y z vx vy vz
				  FFTW3 library
*/

int main()
{
    string line,fName;
    int id,i,j,k,n,m,t,c,tTime,nAtoms,yStep,currentTimeStep,typ;
    double x,y,z,vx,vy,vz,tau1,tau2,tau3,tau4,tau5,tau6,PE,KE,r,rMag,fx,fy,fz;    
    double xLo, yLo, zLo, xHi, yHi, zHi, Lx, Ly, Lz;
    
	// SI units to LAMMPS Metal units conversions
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-12;
    double refVelocity = refLength/refTime; // 1 Angstrom / 1 ps
    double refVolume = refLength*refLength*refLength;
    double refPressure = 100000;
    double refEnergy = 1.60218e-19;
    double refForce = refEnergy / refLength;
    
    //******* CHANGE THESE PARAMETERS **************************
    
    int tSkip = 25; // how often is atom data dumped? 

    double deltaT = 0.001; // MD timestep

    int nTimeSteps = 4000; // no of time steps in dump file: use command line: grep -o 'TIMESTEP' <dumpfilename> | wc -l

    int FVCFSIZE = 4000; // the size of the CCF in terms of no. of time steps
	
    #define Nfft 7999 // must be 2*FVCFSIZE - 1
	
	//*****************************************
	
    string dumpIn;
    string dumpOut;

    // takes dumpfile and output file names as string inputs
    cout << "Type input file: ";
    cin >> dumpIn;

    cout << "Type out file: ";
    cin >> dumpOut;

    ifstream data(dumpIn, ios::in);

    //*****************************************
        
    int T = FVCFSIZE;

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

    cout << "start filtering" << endl;

    // loop below goes through dumpfile, and identifies which atoms remain in the desired domain throughout the specified time period
    for (t = 0; t < nTimeSteps; t++)
    {
        counter++;

        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                if (counter >= nEvery) // progress update
                {
                    counter = 0;

                    cout << "currentTimeStep = " << currentTimeStep
                        << "; t = " << t << " [ "
                        << 100 * float(t + 1) / float(nTimeSteps)
                        << "% ]" << endl;
                }
            }

            if (n == 6)
            {
                data >> xLo >> xHi;
            }

            if (n == 7)
            {
                data >> yLo >> yHi;
            }

            if (n == 8)
            {
                data >> zLo >> zHi;
            }

            getline(data, line);
        }

        if (t == 0)
        {
            Lx = xHi - xLo;
            Ly = yHi - yLo;
            Lz = zHi - zLo;
	    atomIDs.resize(nAtoms);
            filter.resize(nAtoms);
        }
               
        atomIDs = {};
        for (n = 0; n < nAtoms; n++)
        {
            
            data>>id>>typ>>x>>y>>z>>vx>>vy>>vz;
            
            if (t == 0)
            {
                
                filter[n] = id;
                
            }
            
            atomIDs.push_back(id);
        }
        
		// store all atom IDs, identify which are present across all timesteps
        v_intersection = {};
        std::set_intersection(filter.begin(), filter.end(),
            atomIDs.begin(), atomIDs.end(),
            std::back_inserter(v_intersection));

        filter = v_intersection;

        getline(data, line);
    }
    cout << "end filtering" << endl;
    // go back to top of dumpfile, remembering which atom IDs are present across all timesteps
    data.clear();
    data.seekg(0, ios::beg);
    int nAtomsFilter = filter.size();
    cout << nAtomsFilter << endl;
    
    // for the atoms that remain in the domain, extract their velocities in x y z 
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

        
        
        if(t == 0)
        {
            Lx = xHi - xLo;
            Ly = yHi - yLo;
            Lz = zHi - zLo;      
            velXAll.resize(nAtomsFilter);
            velYAll.resize(nAtomsFilter);
            velZAll.resize(nAtomsFilter);
            forXAll.resize(nAtomsFilter);
            forYAll.resize(nAtomsFilter);
            forZAll.resize(nAtomsFilter);
        }
        
        f = 0;
        for(n=0;n<nAtoms;n++)
        {
          data>>id>>typ>>x>>y>>z>>vx>>vy>>vz;

            if (std::count(filter.begin(), filter.end(), id))
            {
                velXAll[f].push_back(vx* refVelocity);
                velYAll[f].push_back(vy* refVelocity);
                velZAll[f].push_back(vz* refVelocity);
                forXAll[f].push_back(vx* refVelocity);
                forYAll[f].push_back(vy* refVelocity);
                forZAll[f].push_back(vz* refVelocity);
                f = f + 1;
            }
        }
        getline(data,line);
    }



    
    
    vector< vector<double> > molfvcfx(nAtoms);
    vector< vector<double> > molfvcfy(nAtoms);
    vector< vector<double> > molfvcfz(nAtoms);
    vector< vector<double> > nCount(nAtoms);

    
    for(n=0;n<nAtoms;n++)
    {
        molfvcfx[n].resize(T, 0.0);
        molfvcfy[n].resize(T, 0.0);
        molfvcfz[n].resize(T, 0.0);
        nCount[n].resize(T, 0.0);
    }    

    nEvery = floor(nAtoms/100);

    counter = 0;    

    double sum, sumX, sumY, sumZ;
    int N = floor(velYAll[0].size()/T);


    cout << "N = " << N << endl;
    
    if(N < 1)
    {
        cout << "ERROR -> Choose nBlock at least equal to nTimeSteps" << endl;
    }

    double vx0, vy0, vz0;
    double fx0, fy0, fz0;

    cout << N << endl;


    
    vector<double> XforceofOneAtomInTime, YforceofOneAtomInTime, ZforceofOneAtomInTime;
    vector<double> XvelocityOneAtomInTime, YvelocityOneAtomInTime, ZvelocityOneAtomInTime;
    XforceofOneAtomInTime.resize(nTimeSteps, 0.0);
    YforceofOneAtomInTime.resize(nTimeSteps, 0.0);
    ZforceofOneAtomInTime.resize(nTimeSteps, 0.0);
    XvelocityOneAtomInTime.resize(nTimeSteps, 0.0);
    YvelocityOneAtomInTime.resize(nTimeSteps, 0.0);
    ZvelocityOneAtomInTime.resize(nTimeSteps, 0.0);

    vector<double> XforceofOneAtomInTimeWindowed, YforceofOneAtomInTimeWindowed, ZforceofOneAtomInTimeWindowed;
    vector<double> XvelocityOneAtomInTimeWindowed, YvelocityOneAtomInTimeWindowed, ZvelocityOneAtomInTimeWindowed;
    XforceofOneAtomInTimeWindowed.resize(nTimeSteps, 0.0);
    YforceofOneAtomInTimeWindowed.resize(nTimeSteps, 0.0);
    ZforceofOneAtomInTimeWindowed.resize(nTimeSteps, 0.0);
    XvelocityOneAtomInTimeWindowed.resize(nTimeSteps, 0.0);
    YvelocityOneAtomInTimeWindowed.resize(nTimeSteps, 0.0);
    ZvelocityOneAtomInTimeWindowed.resize(nTimeSteps, 0.0);

    vector<double> XforceofOneAtomInTimePadded;
    vector<double> YforceofOneAtomInTimePadded;
    vector<double> ZforceofOneAtomInTimePadded;
    vector<double> XvelocityOneAtomInTimePadded;
    vector<double> YvelocityOneAtomInTimePadded;
    vector<double> ZvelocityOneAtomInTimePadded;

    static double in_forx[Nfft][2], in_fory[Nfft][2], in_forz[Nfft][2];
    static double in_velx[Nfft][2], in_vely[Nfft][2], in_velz[Nfft][2];
    static double For_wx[Nfft][2], For_wy[Nfft][2], For_wz[Nfft][2];
    static double Vel_wx[Nfft][2], Vel_wy[Nfft][2], Vel_wz[Nfft][2];
    static double CVel_wx[Nfft][2], CVel_wy[Nfft][2], CVel_wz[Nfft][2];
    static double fftfvx[Nfft][2], fftfvy[Nfft][2], fftfvz[Nfft][2];
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

     for(n=0;n< nAtomsFilter;n++) // pick one atom
    {

        counter++; 
        if(counter >= nEvery) // progress update
        {
            counter = 0;
            
            cout << "atom = " << n << " [" 
                << 100*float(n+1)/float(nAtomsFilter)
                << "% ]" << endl;             
        }

        XforceofOneAtomInTime = forXAll[n];
        YforceofOneAtomInTime = forYAll[n];
        ZforceofOneAtomInTime = forZAll[n];
        XvelocityOneAtomInTime = velXAll[n];
        YvelocityOneAtomInTime = velYAll[n];
        ZvelocityOneAtomInTime = velZAll[n];

        XforceofOneAtomInTimePadded = XforceofOneAtomInTime;
        YforceofOneAtomInTimePadded = YforceofOneAtomInTime;
        ZforceofOneAtomInTimePadded = ZforceofOneAtomInTime;
        XvelocityOneAtomInTimePadded = XvelocityOneAtomInTime;
        YvelocityOneAtomInTimePadded = YvelocityOneAtomInTime;
        ZvelocityOneAtomInTimePadded = ZvelocityOneAtomInTime;

        // 0 pad to make velocity vectors 2n-1 long

        for (i = 1; i < nTimeSteps; i++)
        {
            XforceofOneAtomInTimePadded.push_back(0);
            YforceofOneAtomInTimePadded.push_back(0);
            ZforceofOneAtomInTimePadded.push_back(0);
            XvelocityOneAtomInTimePadded.push_back(0);
            YvelocityOneAtomInTimePadded.push_back(0);
            ZvelocityOneAtomInTimePadded.push_back(0);
        }

        for (i = 0; i < Nfft; i++) {
            in_forx[i][0] = XforceofOneAtomInTimePadded[i];
            in_forx[i][1] = 0.000000;
            in_fory[i][0] = YforceofOneAtomInTimePadded[i];
            in_fory[i][1] = 0.000000;
            in_forz[i][0] = ZforceofOneAtomInTimePadded[i];
            in_forz[i][1] = 0.000000;
            in_velx[i][0] = XvelocityOneAtomInTimePadded[i];
            in_velx[i][1] = 0.000000;
            in_vely[i][0] = YvelocityOneAtomInTimePadded[i];
            in_vely[i][1] = 0.000000;
            in_velz[i][0] = ZvelocityOneAtomInTimePadded[i];
            in_velz[i][1] = 0.000000;
        }

		// compute FFT

        fftforx = fftw_plan_dft_1d(Nfft, in_forx, For_wx, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftforx);
        fftw_destroy_plan(fftforx);
        fftfory = fftw_plan_dft_1d(Nfft, in_fory, For_wy, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftfory);
        fftw_destroy_plan(fftfory);
        fftforz = fftw_plan_dft_1d(Nfft, in_forz, For_wz, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftforz);
        fftw_destroy_plan(fftforz);

        fftvelx = fftw_plan_dft_1d(Nfft, in_velx, Vel_wx, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftvelx);
        fftw_destroy_plan(fftvelx);
        fftvely = fftw_plan_dft_1d(Nfft, in_vely, Vel_wy, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftvely);
        fftw_destroy_plan(fftvely);
        fftvelz = fftw_plan_dft_1d(Nfft, in_velz, Vel_wz, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftvelz);
        fftw_destroy_plan(fftvelz);

		// multiply complex nrs
		
        for (i = 0; i < Nfft/2; i++) {
            fftfvx[i][0] = For_wx[i][0] * Vel_wx[i][0] - For_wx[i][1] * Vel_wx[i][1];
            fftfvx[i][1] = For_wx[i][0] * Vel_wx[i][1] + For_wx[i][1] * Vel_wx[i][0];

            fftfvy[i][0] = For_wy[i][0] * Vel_wy[i][0] - For_wy[i][1] * Vel_wy[i][1];
            fftfvy[i][1] = For_wy[i][0] * Vel_wy[i][1] + For_wy[i][1] * Vel_wy[i][0];

            fftfvz[i][0] = For_wz[i][0] * Vel_wz[i][0] - For_wz[i][1] * Vel_wz[i][1];
            fftfvz[i][1] = For_wz[i][0] * Vel_wz[i][1] + For_wz[i][1] * Vel_wz[i][0];
        }


		// compute absolute of VDOS
        if (std::count(XforceofOneAtomInTime.begin(), XforceofOneAtomInTime.end(), 0)) {
        }
        else {
            for (i = 0; i < Nfft/2; i++) {
                spectrumX[i] = spectrumX[i] + sqrt(fftfvx[i][0]*fftfvx[i][0]+fftfvx[i][1]*fftfvx[i][1]);
                spectrumY[i] = spectrumY[i] + sqrt(fftfvy[i][0]*fftfvy[i][0]+fftfvy[i][1]*fftfvy[i][1]);
                spectrumZ[i] = spectrumZ[i] + sqrt(fftfvz[i][0]*fftfvz[i][0]+fftfvz[i][1]*fftfvz[i][1]);
            }
        }

        XforceofOneAtomInTime.clear();
        YforceofOneAtomInTime.clear();
        ZforceofOneAtomInTime.clear();
        XvelocityOneAtomInTime.clear();
        YvelocityOneAtomInTime.clear();
        ZvelocityOneAtomInTime.clear();
        XforceofOneAtomInTimePadded.clear();
        YforceofOneAtomInTimePadded.clear();
        ZforceofOneAtomInTimePadded.clear();
        XvelocityOneAtomInTimePadded.clear();
        YvelocityOneAtomInTimePadded.clear();
        ZvelocityOneAtomInTimePadded.clear();
}

for (k = 0; k < Nfft/2; k++)
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
