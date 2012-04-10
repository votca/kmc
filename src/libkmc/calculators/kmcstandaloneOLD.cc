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

using namespace std;


string singlequery(votca::tools::Database db, string statement)
{
    string result;
    votca::tools::Statement *stmt = db.Prepare(statement);
    while (stmt->Step() != SQLITE_DONE) {
        result += stmt->Column<string>(0);
    }
    return result;
}

int main(int argc, char** argv)
{
    //float runtime = 60*60;
    //float seedNo  = 23;
    string filename = "/people/thnfs/homes/kordt/votca/votca_ctp_testsuite/norma/work/state.db";
    string injectionpattern = "%";
    
    std::cout << "KMC Standalone for testing purposes" << std::endl;
    
    // Determine number of segments and injectable segments
    votca::tools::Database db;
    db.Open( filename );
    cout << "Loading graph from " << filename << endl;
    string queryresult = singlequery(db, "SELECT COUNT(id) FROM segments");
    int no_segment = atoi(queryresult.c_str());
    cout << "Found " << no_segment << " segments.";
    
    string statement = "SELECT COUNT(id) FROM segments WHERE name LIKE '";
    statement += injectionpattern;
    statement += "'";
    queryresult = singlequery(db, statement);
    int no_injsegment = atoi(queryresult.c_str());
    cout << " Out of those, " << no_injsegment << " segments match the injection pattern." << endl;
    
    // Load segments
    int segment[no_segment];
    votca::tools::Statement *stmt = db.Prepare("SELECT id, name FROM segments");

    int i= 0;
    while (stmt->Step() != SQLITE_DONE) {
        int id = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        segment[i] = id;
        i++;
    }
    delete stmt;
    
    // Load injectable segments
    int injsegment[no_injsegment];
    statement = "SELECT id, name FROM segments WHERE name LIKE '";
    statement += injectionpattern;
    statement += "'";
    stmt = db.Prepare(statement);

    i= 0;
    while (stmt->Step() != SQLITE_DONE) {
        int id = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        injsegment[i] = id;
        i++;
    }
    delete stmt;
    
    // Determine number of pairs
    queryresult = singlequery(db, "SELECT COUNT(id) FROM pairs");
    int no_pair = 2 * atoi(queryresult.c_str());
    
    // Load pairs and rates
    int pair_seg1[no_pair];
    int pair_seg2[no_pair];
    double pair_rate[no_pair];
    stmt = db.Prepare("SELECT seg1, seg2, rate12e FROM pairs;");
    i = 0;
    while (stmt->Step() != SQLITE_DONE)
    {
        pair_seg1[i] = stmt->Column<int>(0);
        pair_seg2[i] = stmt->Column<int>(1);
        pair_rate[i] = stmt->Column<double>(2);
        i++;
    }    
    delete stmt;
    stmt = db.Prepare("SELECT seg2, seg1, rate21e FROM pairs;");
    while (stmt->Step() != SQLITE_DONE)
    {
        pair_seg1[i] = stmt->Column<int>(0);
        pair_seg2[i] = stmt->Column<int>(1);
        pair_rate[i] = stmt->Column<double>(2);
        i++;
    }    
    delete stmt;

    return (EXIT_SUCCESS);
}

