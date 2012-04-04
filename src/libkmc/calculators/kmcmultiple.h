/*
 * Copyright 2009-2012 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_KMC_MULTIPLE_H
#define	__VOTCA_KMC_MULTIPLE_H

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include "node.h"

using namespace std;
using namespace votca::kmc;



class KMCMultiple : public KMCCalculator 
{
public:
    KMCMultiple() {};
   ~KMCMultiple() {};

    void Initialize(const char *filename, Property *options );
    bool EvaluateFrame();

protected:
	    void LoadGraph();
            void RunKMC(void);
            void WriteOcc(void);

            map<int , node_t *> _nodes_lookup;
            vector<node_t *> _nodes;
            map<int , node_t *> _injection_lookup;
            vector<node_t *> _injection;
            string _injection_name;
            double _runtime;
            double _dt;
            int _seed;
            string _filename; // HACK
            
};

void KMCMultiple::Initialize(const char *filename, Property *options )
{
    cout << "---- KMCMultiple::Initialise ---- " << endl;
        // check that input in the file specified via -o optionsfile exists
        if (options->exists("options.kmcmultiple.runtime")) {
	    _runtime = options->get("options.kmcmultiple.runtime").as<double>();
	}
	else {
	    throw invalid_argument("Error in kmcsingle: total run time is not provided");
        }

    	if (options->exists("options.kmcmultiple.outputtime")) {
	    _dt = options->get("options.kmcmultiple.outputtime").as<double>();
	}
	else {
	    throw invalid_argument("Error in kmcmultiple: output frequency is not provided");
        }

    	if (options->exists("options.kmcmultiple.seed")) {
	    _seed = options->get("options.kmcmultiple.seed").as<int>();
	}
	else {
	    throw invalid_argument("Error in kmcmultiple: seed is not provided");
        }
   	if (options->exists("options.kmcmultiple.injection")) {
	    _injection_name = options->get("options.kmcmultiple.injection").as<string>();
	}
        else {
	    _injection_name = "*";
           throw invalid_argument("Error in kmcmultiple: injection pattern is not provided. The simplest option is setting it to *.");
        }
        // check that the input has no empty fields and is valid
        if (_runtime <= 0 || _dt <= 0) {
           throw invalid_argument("Error in kmcmultiple: run time and frequency ('outputtime') must be greater than 0.");
        }
        if (_injection_name.length() < 1) {
           throw invalid_argument("Error in kmcmultiple: injection pattern is empty. (The simplest option is setting it to '*'.)");
        }
        
        _filename = filename;

       srand(_seed);
       votca::tools::Random::init(rand(), rand(), rand(), rand());

}

bool KMCMultiple::EvaluateFrame()
{
    cout << "---- KMCMultiple::EvaluateFrame ----" << endl;

    LoadGraph();
    RunKMC();
    return true;
}

void KMCMultiple::LoadGraph()
{
    cout << "-------- KMCMultiple::LoadGraph ----" << endl;

    Database db;
    db.Open( _filename );
    cout << " Loading graph from " << _filename << endl;
    Statement *stmt = db.Prepare("SELECT id, name FROM segments;");

    while (stmt->Step() != SQLITE_DONE) {
        int id = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        node_t *n =new node_t(id);
        _nodes.push_back(n);
        _nodes_lookup[id] = _nodes.back();
        if (wildcmp(_injection_name.c_str(), name.c_str())) {
            _injection.push_back(n);
            _injection_lookup[id] = _injection.back();
        }
    }
    //delete stmt;
    cout << "  -Nodes: " << _nodes.size() << endl;
    cout << "  -Seed:" << _seed << endl;
    if(_seed > _nodes.size()){ throw invalid_argument ("Error in kmcmultiple: Seed is outside the range of nodes. Please specify an existing seed in your input file."); }
    cout << "  -Injection Points: " << _injection.size() << endl;

    delete stmt;

    int links = 0;
    stmt = db.Prepare("SELECT seg1, seg2, rate12e, rate21e, drX, drY, drZ FROM pairs;"); // electron rates, check (think about) this
    while (stmt->Step() != SQLITE_DONE) {
        node_t *n1 = _nodes_lookup[stmt->Column<int>(0)];
        node_t *n2 = _nodes_lookup[stmt->Column<int>(1)];
        double rate12 = stmt->Column<double>(2);
        double rate21 = stmt->Column<double>(3);
        vec r = vec(stmt->Column<double>(4), stmt->Column<double>(5), stmt->Column<double>(6));
        // cout << "adding event " << n2 << " | " << rate12 << " | " << r << endl;
        // cout << "adding event " << n1 << " | " << rate21 << " | " << r << endl;
        n1->AddEvent(new link_t(n2, rate12, r));
        n2->AddEvent(new link_t(n1, rate21, -r));
        links += 2;
        
        if ( votca::tools::globals::verbose ) {
            cout << "rate12=" << rate12 << endl;
            cout << "rate21=" << rate21 << endl;
            cout << "r=" << r << endl;
        }
    }
    delete stmt;
    cout << "  -Links: " << links << endl;

}

void KMCMultiple::RunKMC(void)
{
    cout << "-------- KMCMultiple::RunKMC ---- " << endl;
    double t = 0;

    srand(_seed);
    votca::tools::Random::init(rand(), rand(), rand(), rand());

    // cout << " seed:size:site " << _seed << ":" << _injection.size() << ":" << Random::rand_uniform_int(_injection.size()) << endl;
    current=_injection[Random::rand_uniform_int(_injection.size())];
    cout <<"       Starting simulation at node: "<<current->_id-1<<endl;

}

void KMCMultiple::WriteOcc()
{
    cout << "---- KMCMultiple::WriteOcc ---- " << endl;
    Database db;
    cout << "     Opening for writing " << _filename << endl;
	db.Open(_filename);
	db.Exec("BEGIN;");
	Statement *stmt = db.Prepare("UPDATE segments SET occPe = ? WHERE id = ?;");  // electron occ. prob., check (think about) this
	for(int i=0; i<_nodes.size(); ++i) {
		stmt->Reset();
		stmt->Bind(1, _nodes[i]->_occ/_runtime);
		stmt->Bind(2, _nodes[i]->_id);
		stmt->Step();
	}
	db.Exec("END;");
	delete stmt;
}

#endif	/* __VOTCA_KMC_MULTIPLE_H */
