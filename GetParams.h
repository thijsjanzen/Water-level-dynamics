#ifndef GETPARAMS_H_INCLUDED
#define GETPARAMS_H_INCLUDED

#include <sstream>
#include <string>
#include <vector>
#include "Gillespie.h"
#include "Beta.h"

struct p;

class GetParams {
public:
	GetParams();
	void readFromIni( const char * filename );

	//GLOBAL VARIABLES
	int global_ID;
	int species_ID;
	double maxT;
	double period_waterlevel;

	double extinctionRate;
	
	double SympatricRate_low;
	double SympatricRate_high;
	double AllopatricRate;
	
	double Jiggle;


	int replicates;

    double pi;
	
	double millionTime;

//	int waterModel;

	template <typename T>
	void readNameValuePair( std::stringstream& ss, std::string iniName, T& value );

	std::vector< double > vec_epsilon;
	std::vector< double > emp_values;
    std::vector< std::vector< double > > sd_values;


	int numberParticles;

	std::vector<spec_point> single_outcome;

	int generateData;
	

	int t;
	int numberExtinctions;

	int numAllo;
	int numExtinct;

    int waterModel;


	double sigma;
	bool store_for_recheck;
};

#endif //GETPARAMS_H_INCLUDED