/* 
 * File:   kmcstandalone.cc
 * Author: kordt
 *
 * Development test version of the KMC VSSM algorithm
 */

#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include "/people/thnfs/homes/kordt/votca/include/votca/tools/database.h"
#include "/people/thnfs/homes/kordt/votca/include/votca/tools/statement.h"
#include "/people/thnfs/homes/kordt/votca/include/votca/tools/tokenizer.h"
#include "/people/thnfs/homes/kordt/votca/include/votca/tools/random.h"

using namespace std;
int verbose = 0; // 0=minimal output, 1=verbose output


string singlequery(votca::tools::Database db, string statement)
{
    string result;
    votca::tools::Statement *stmt = db.Prepare(statement);
    while (stmt->Step() != SQLITE_DONE) 
    {
        result += stmt->Column<string>(0);
    }
    return result;
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
    votca::tools::Statement *stmt = db.Prepare("SELECT id, name FROM segments");

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
    stmt = db.Prepare("SELECT seg1, seg2, rate12e, rate21e FROM pairs;");
    while (stmt->Step() != SQLITE_DONE)
    {
        int seg1 = stmt->Column<int>(0);
        int seg2 = stmt->Column<int>(1);
        double rate12 = stmt->Column<double>(2);
        double rate21 = stmt->Column<double>(3);
        pair_seg1.push_back(seg1);
        pair_seg2.push_back(seg2);
        pair_rate.push_back(rate12);
        pair_seg1.push_back(seg2);
        pair_seg2.push_back(seg1);
        pair_rate.push_back(rate21);
    }    
    delete stmt;

    cout << "pairs: " << pair_seg1.size()/2 << endl;
}

void VSSMRun(vector<int> segment, vector<int> injsegment, vector<int> pair_seg1, vector<int> pair_seg2, vector<double> pair_rate, double runtime, vector<double> &occP)
{
    // Injection
    int position = votca::tools::Random::rand_uniform_int(injsegment.size());
    cout << "starting simulation at segment " << segment[position] << " (internal position " << position << ")" << endl;

    double time = 0;
    int step = 0;
    double normalize = 0;
    int do_newconf;
    double maxprob = 0;
    while(time < runtime)
    {
        // cout << int(time/runtime *100) << "% done" << endl;
        if(verbose >= 1) {cout << "current position: segment " << segment[position] << endl;}
        // make a list of all possible reactions and rates
        vector<int> newconf(0);
        vector<double> newprob(0);
        if(verbose >= 1) {cout << "possible jumps: " ;}
        for (unsigned int j=0; j<pair_seg1.size(); j++)
        {
            if(pair_seg1[j] == segment[position])
            {
                newconf.push_back(pair_seg2[j]);
                newprob.push_back(pair_rate[j] * votca::tools::Random::rand_uniform());
                if(verbose >= 1) {cout << pair_seg2[j] << " ";}
            }
        }
        if(verbose >= 1) {cout << endl;}
        // this should not happen: no possible jumps defined for a node
        if (newconf.size() == 0)
        {
            cout << "WARNING: no possible jumps found from seqment " << segment[position] << ". The charge is trapped here. Press Enter to continue anyway." << endl;
            cin.get();
        }
        
        // get reaction with the highest probability and sum up probabilites
        maxprob = 0;
        normalize = 0;
        for(unsigned int j=0; j<newprob.size(); j++)
        {
            if(newprob[j] > maxprob)
            {
                maxprob = newprob[j];
                do_newconf = newconf[j];
                normalize += newprob[j];
            }
        }
        double dt = 0;
        // go forward in time
        if(normalize == 0)
        {
            throw runtime_error("Error in kmcsingle: Incorrect rates in the database file. All the outgoing rates for the current segment are 0.");
        }
        else
        {
            double u = votca::tools::Random::rand_uniform();
            while(u == 0)
            {
                cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                u = votca::tools::Random::rand_uniform();
            }
                
            dt = -1 / normalize * log(u);
            if(verbose >= 1) {cout << "DT=" << dt << endl;}
        }
        time += dt;
        step += 1;
        progressbar(time/runtime);
        if(verbose >= 1) 
        {
            cout << "t = " << time << "   dt = " << dt << "    (step " << step << ")" << endl << endl;
            cout << "normalize=" << normalize << endl;
        }
        // update occupation probability
        occP[position] += dt;
        // jump!
        for(unsigned int j=0; j<=segment.size(); j++)
        {
            if(segment[j] == do_newconf)
            {
                position = j;
                break;
            }
        }
    }
    
    cout << endl << "finished KMC simulation after " << step << " steps." << endl << endl;

    // divide by time to get occupation probabilites instead of occupation times
    for(unsigned int j=0; j<occP.size();j++) 
    {
        occP[j] /= time;
    }
}

int main(int argc, char** argv)
{
    double runtime = 1E5;
    int seed  = 23;
    
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
    vector<double> occP(segment.size(),0.);
    VSSMRun(segment, injsegment, pair_seg1, pair_seg2, pair_rate, runtime, occP);
    
    // output occupation probabilites
    for(unsigned int j=0; j<occP.size();j++) 
    {
        if(occP[j] > 0)
        {
            cout << "occupation probability " << segment[j] << ": " << occP[j] << " (internal position: "<< j << ")" << endl;
        }
    }

    return (EXIT_SUCCESS);
}

