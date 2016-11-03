
#include "GetParams.h"
#include <fstream>
#include <iostream>

//This cpp file was adapted from
//H.Hildenbrandt 2007




GetParams::GetParams() {
	global_ID = 0;
	species_ID = 1;

	period_waterlevel = 500;

	maxT = 5.45;

    pi = 3.14159265;

	numberParticles = 100;

	sigma = 0.05;

	store_for_recheck = false;

	period_waterlevel = -100;
	SympatricRate_high = -100;
	SympatricRate_low = -100;
	AllopatricRate = -100;
	extinctionRate = -100;
	
    std::vector<double> sd_values_real = {0.156467, 2.332945, 1.067379, 68.19059 };
    std::vector<double> sd_values_norm = {0.2287997, 1.801797, 1.427837, 122.6273 };
    std::vector<double> sd_values_10 = {0.2091974,  2.698178, 1.167943, 77.22117 };

    sd_values.push_back(sd_values_norm);
    sd_values.push_back(sd_values_real);
    sd_values.push_back(sd_values_10);

    waterModel = 0;

}

void GetParams::readFromIni( const char * filename ) {
	
	//locate file and tranfer text to stringstream
	std::ifstream ifs( filename ); 
	std::stringstream ss;

	if( ifs ) { //only for succesfully created ifstream: otherwise null-pointer?
		ss << ifs.rdbuf(); //config.ini content is transferred to stringstream: easier to search?
	}
	else {
		throw "Can't locate file";
		//std::cout << "Can't locate file" << std::endl;
	}
	while( ss.good() ) {	
				readNameValuePair( ss,  "numberParticles",numberParticles);
				readNameValuePair( ss,  "generateData",generateData);
                readNameValuePair( ss,  "waterModel", waterModel);
	}
}

template <typename T>
void GetParams::readNameValuePair( std::stringstream& ss, std::string iniName, T& value ) {
	std::string name;
	char sign;
	ss >> name; //>> copies ss content to string until white space is encountered
	if( name != iniName )
		throw "expect parameter";
	ss >> sign; //copies ss content to character 
	if( sign != '=' )
		throw "text format of ini file is not compatible";
    ss >> value;
	std::cout << iniName << ": " << value << std::endl;
//	std::cerr << iniName << ": " << value << std::endl;
}


