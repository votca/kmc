/* 
 * File:   kmcstandalone.cc
 * Author: kordt
 *
 * Development test version of the KMC VSSM algorithm
 */

#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

#include "/people/thnfs/homes/kordt/votca/include/votca/tools/database.h"
#include "/people/thnfs/homes/kordt/votca/include/votca/tools/statement.h"
#include "/people/thnfs/homes/kordt/votca/include/votca/tools/tokenizer.h"
#include "/people/thnfs/homes/kordt/votca/include/votca/tools/random.h"

using namespace std;
static int verbose = 0; // 0=minimal output, 1=verbose output

int OMPinfo() 
{
    int nthreads, tid, procs, maxt, inpar, dynamic, nested;
    printf("\n||\n|| openMP PARALLEL COMPUTATION STATUS:\n");
    /* Start parallel region */
    #pragma omp parallel private(tid)
    {
        /* Obtain thread number */
        tid = omp_get_thread_num();
        /* Only master thread does this */
        if (tid == 0) 
        {
            // printf("Thread %d getting environment info...\n", tid);
            /* Get environment information */
            procs = omp_get_num_procs();
            nthreads = omp_get_num_threads();
            maxt = omp_get_max_threads();
            inpar = omp_in_parallel();
            dynamic = omp_get_dynamic();
            nested = omp_get_nested();
            /* Print environment information */
            printf("|| Number of processors = %d\n", procs);
            printf("|| Number of threads = %d\n", nthreads);
            // printf("|| Max threads = %d\n", maxt);
            printf("|| In parallel? = %d\n||\n", inpar);
            // printf("|| Dynamic threads enabled? = %d\n", dynamic);
            // printf("|| Nested parallelism supported? = %d\n", nested);
        }
    }  /* Done */
    return nthreads;
}

void progressbar(double fraction)
{
    int totalbars = 50;
    std::cout << "\r";
    for(double bars=0; bars<double(totalbars); bars++)
    {
        if(bars<=fraction*double(totalbars))
        {
            std::cout << "|";
        }
        else
        {
            std::cout << "-";
        }
    }
    std::cout << "  " << int(fraction*1000)/10. <<" %   ";
    std::cout << std::flush;
    if(fraction*100 == 100)
    {
        std::cout << std::endl;
    }
}


void LoadGraph(vector<int> &segment, vector<int> &injsegment, vector<int> &pair_seg1, vector<int> &pair_seg2, vector<double> &pair_rate)
{
    string filename = "/people/thnfs/homes/kordt/votca/votca_ctp_testsuite/norma/work/state.db";
    string injectionpattern = "*";
    
    // Load segments
    votca::tools::Database db;
    db.Open( filename );
    cout << "LOADING GRAPH" << endl << "database file: " << filename << endl;
    votca::tools::Statement *stmt = db.Prepare("SELECT id, name FROM segments;");

    while (stmt->Step() != SQLITE_DONE)
    {
        int id = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        segment.push_back(id);
        if (votca::tools::wildcmp(injectionpattern.c_str(), name.c_str()))
        {
            injsegment.push_back(id);
        }
    }
    delete stmt;
    cout << "segments: " << segment.size() << endl;
    cout << "injectable segments: " << injsegment.size() << endl;
    
    // Load pairs and rates
    stmt = db.Prepare("SELECT seg1 AS 'segment1', seg2 AS 'segment2', rate12e AS 'rate' FROM pairs UNION SELECT seg2 AS 'segment1', seg1 AS 'segment2', rate21e AS 'rate' FROM pairs ORDER BY segment1;");
    while (stmt->Step() != SQLITE_DONE)
    {
        int seg1 = stmt->Column<int>(0);
        int seg2 = stmt->Column<int>(1);
        double rate12 = stmt->Column<double>(2);
        double rate21 = stmt->Column<double>(3);
        pair_seg1.push_back(seg1);
        pair_seg2.push_back(seg2);
        pair_rate.push_back(rate12);
        //pair_seg1.push_back(seg2);
        //pair_seg2.push_back(seg1);
        //pair_rate.push_back(rate21);
    }    
    delete stmt;

    cout << "pairs: " << pair_seg1.size()/2 << endl;
}

void WriteOcc(vector<double> occP, vector<int> segment)
{
    string filename = "/people/thnfs/homes/kordt/votca/votca_ctp_testsuite/norma/work/state.db";
    votca::tools::Database db;
    cout << "Opening for writing " << filename << endl;
	db.Open(filename);
	db.Exec("BEGIN;");
	votca::tools::Statement *stmt = db.Prepare("UPDATE segments SET occPe = ? WHERE id = ?;");  // electron occ. prob., check (think about) this
	for(unsigned int i=0; i<segment.size(); ++i)
        {
	    stmt->Reset();
	    stmt->Bind(1, occP[i]);
	    stmt->Bind(2, segment[i]);
	    stmt->Step();
	}
	db.Exec("END;");
	delete stmt;
}


vector<double> VSSMRunMultiple(vector<int> &segment, vector<int> &injsegment, vector<int> &pair_seg1, vector<int> &pair_seg2, vector<double> &pair_rate, double runtime, unsigned int numberofcharges)
{
    cout << "Algorithm: VSSM for Multiple Charges" << endl;
    double outputfrequency = runtime/10;
    cout << "number of charges: " << numberofcharges << endl;
    cout << "number of nodes: " << segment.size() << endl;
    vector<int> occupied(segment.size(),0);
    vector<double> occP(segment.size(),0.);
    if(numberofcharges>segment.size())
    {
        throw runtime_error("Error in kmcstandalone: Your number of charges is larger than the number of nodes. This conflicts with single occupation.");
    }
    // Injection
    vector<int> position(numberofcharges,0);
    for(unsigned int charge=0; charge<numberofcharges; charge++)
    {
        #pragma omp critical
        {position[charge]= votca::tools::Random::rand_uniform_int(injsegment.size());}
        while(occupied[position[charge]] == 1)
        {   // maybe already occupied?
            #pragma omp critical
            {position[charge]= votca::tools::Random::rand_uniform_int(injsegment.size());}
        }
        occupied[position[charge]] = 1;
        cout << "starting position for charge " << charge+1 << ": segment " << segment[position[charge]] << " (internal position " << position[charge] << ")" << endl;
    }        

    // determine first and last positions in the pairs database
    vector<unsigned int> first_j(segment.size(),0);
    vector<unsigned int> last_j(segment.size(),0);
    first_j[0] = 0;
    last_j[segment.size()-1] = pair_seg1.size()-1;
    for (unsigned int j=1; j<pair_seg1.size(); j++)
    {
        if(pair_seg1[j-1] != pair_seg1[j])
        {
            last_j[pair_seg1[j-1]-1] = j-1;
            first_j[pair_seg1[j]-1] = j;
        }
    }

    double time = 0.;
    int step = 0;
    double normalize = 0.;
    int do_newconf;
    int do_newaffectedcharge;
    double maxprob = 0.;
    double newprob = 0.;
    double nextoutput = outputfrequency;
    
    // unsigned int maxnumberofjumps = MaxNumberOfJumps(segment, pair_seg1);

    
    while(time < runtime)
    {
        // make a list of all possible reactions and rates
        // get reaction with the highest probability and sum up probabilites
        maxprob = 0.;
        normalize = 0;
        // newprob = 0.;
        for(unsigned int charge=0; charge<numberofcharges; charge++)
        {
            # pragma omp master
            { if(verbose >= 1) {cout << "charge " << charge+1 << " at node " << segment[position[charge]] << " - possible jumps: " ;} }
            for (unsigned int j=first_j[position[charge]]; j<=last_j[position[charge]]; j++)
            {
                #pragma omp critical
                {newprob = pair_rate[j] * votca::tools::Random::rand_uniform();}
                
                if(newprob > maxprob)
                {
                    maxprob = newprob;
                    do_newconf = pair_seg2[j];
                    do_newaffectedcharge = charge;
                }
                normalize += newprob;
                  
                if(verbose >= 1) {cout << pair_seg2[j] << "  ";}
            }
            if(verbose >= 1) {cout << endl;}
        }
        if(verbose >= 1) {cout << endl;}
        
        // go forward in time
        double dt = 0;
        if(normalize == 0)
        {   // this should not happen: no possible jumps defined for a node
            throw runtime_error("Error in kmcsingle: Incorrect rates in the database file. All the outgoing rates for the current nodes are 0.");
        }
        else
        {
            double u;
            #pragma omp critical
            {u = votca::tools::Random::rand_uniform();}
            while(u == 0)
            {
                cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                #pragma omp critical
                {u = votca::tools::Random::rand_uniform();}
            }
                
            dt = -1 / normalize * log(u);
        }
        time += dt;
        step += 1;
        // jump!
        # pragma omp master
        {
            if(verbose == 1)
            {
                cout << "EVENT for charge " << do_newaffectedcharge+1 << endl;
                cout << "  old segment " << segment[position[do_newaffectedcharge]] << " (occupation: " << occupied[position[do_newaffectedcharge]] << ")" << endl;
                cout << "  new segment " << do_newconf << " (occupation: " << occupied[do_newconf-1] << ")" << endl;
            }
        }
        if(occupied[do_newconf-1]==0) // unsafe, make check at the beginning
        {
            occupied[position[do_newaffectedcharge]] = 0;
            position[do_newaffectedcharge] = do_newconf-1; // unsafe
            occupied[do_newconf-1] = 1; // unsafe, make check at the beginning
            # pragma omp master
            { if(verbose == 1) {cout << "  charge " << do_newaffectedcharge+1 << " jumps to node " << do_newconf << "." << endl;} }
        }
        else
        {
            # pragma omp master
            { if(verbose == 1) {cout << "no jump for charge " << do_newaffectedcharge+1 << ", node " << do_newconf << " is occupied." << endl;} }
        }

        # pragma omp master
        { 
            if(verbose >= 1) 
            {
                cout << "t = " << time << "   dt = " << dt << "    (step " << step << ")" << "[thread " << omp_get_thread_num() << "]" << endl;
                cout << "normalize=" << normalize << endl;
            }
        }
        # pragma omp master
        {
            if(time > nextoutput)
            {
                nextoutput = time + outputfrequency;
                progressbar(time/runtime);
            }
    	}
            
            if(verbose == 1) {cout << endl << endl;}
        
        // update occupation probability
        occP[position[do_newaffectedcharge]] += dt;
        
    }
    
    cout << endl << "finished KMC simulation after " << step << " steps." << endl << endl;

    // divide by time to get occupation probabilites instead of occupation times
    for(unsigned int j=0; j<occP.size();j++) 
    {
        occP[j] /= time;
    }
    return occP;
}


int main(int argc, char** argv)
{
    double runtime = 1E6;
    int seed  = 23;
    unsigned int numberofcharges = 1;
    
    std::cout << "-----------------------------------" << std::endl;      
    std::cout << "KMC Standalone for testing purposes" << std::endl;
    std::cout << "-----------------------------------" << std::endl << std::endl;      

    // Load Graph
    vector<int> segment(0);
    vector<int> injsegment(0);
    vector<int> pair_seg1(0);
    vector<int> pair_seg2(0);
    vector<double> pair_rate(0);
    LoadGraph(segment, injsegment, pair_seg1, pair_seg2, pair_rate);

    // Initialise random number generator
    cout << endl << "INITIALISING RANDOM NUMBER GENERATOR" << endl;
    srand(seed); // srand expects any integer in order to initialise the random number generator
    votca::tools::Random::init(rand(), rand(), rand(), rand());
    
   
    // VSSM KMC algorithm
    cout << endl << "KMC SIMULATION" << endl;
    // unsigned int numberofthreads = OMPinfo();
    unsigned int numberofthreads = 1;
    cout << "NUMBER OF THREADS " << numberofthreads << endl;
    vector<double> occP(segment.size(),0.);
    vector< vector< double > > occPOneRun ( numberofthreads, vector<double> ( segment.size(), 0. ) );

    (void) omp_set_num_threads(numberofthreads);
    #pragma omp parallel
    {
       // occPOneRun[omp_get_thread_num()] = FRMRunMultiple(segment, injsegment, pair_seg1, pair_seg2, pair_rate, runtime/double(numberofthreads), numberofcharges);
       occPOneRun[omp_get_thread_num()] = VSSMRunMultiple(segment, injsegment, pair_seg1, pair_seg2, pair_rate, runtime/double(numberofthreads),numberofcharges);
       // occPOneRun[omp_get_thread_num()] = VSSMRunSingle(segment, injsegment, pair_seg1, pair_seg2, pair_rate, runtime/double(numberofthreads));
    }

    
    // get mean of multiple runs
    for(unsigned int j=0; j<occP.size();j++) 
    {
        for(unsigned int thread=0; thread<numberofthreads; thread++)
        {
            occP[j] += occPOneRun[thread][j];
        }
    }
    for(unsigned int j=0; j<occP.size();j++) 
    {
         occP[j] /= numberofthreads;
    }
    
    // output occupation probabilites
    for(unsigned int thread=0; thread<numberofthreads; thread++)
    {
        for(unsigned int j=0; j<occPOneRun[thread].size();j++) 
        {
            if(occPOneRun[thread][j] > 0)
            {
                cout << "[thread " << thread+1 << "] "<<"occupation probability " << segment[j] << ": " << occPOneRun[thread][j] << " (internal position: "<< j << ")" << endl;
            }
        }
    }


    // output occupation probabilites
    for(unsigned int j=0; j<occP.size();j++) 
    {
        if(occP[j] > 0)
        {
            cout << "occupation probability " << segment[j] << ": " << occP[j] << " (internal position: "<< j << ")" << endl;
        }
    }
    
    // write occupation probabilites
    WriteOcc(occP, segment);

    return (EXIT_SUCCESS);
}

