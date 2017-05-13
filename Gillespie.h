#ifndef GILLESPIE_H
#define GILLESPIE_H


#include <utility>    // std::move, std::forward
#include <vector>
#include <istream>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include "randomc.h"

struct spec_point
{
	spec_point(double id, double t): ID(id), time(t)
	{
	}

	double ID;
	double time;

	bool operator<(const spec_point& other) const { return time < other.time; }

	bool operator==(const spec_point& other) const
	{ 
		if(time == other.time) return true; 
		else return false;
	}
	bool operator!=(const spec_point& other) const {return !((*this) == other);}
};




struct species
{
	species()
	{
		ID = -1;
		death_time = -1;
	}
	species(int& id);
	species(const species& parent_species, int& id_count, double b_time);

	
	species(const species& other): ID(other.ID),
		parent(other.parent),
		birth_time(other.birth_time),
		death_time(other.death_time),
		extant_offspring(other.extant_offspring),
		checked(other.checked),
		alloSpeciated(other.alloSpeciated)
	{
	}



	int ID;
	int parent;
	double birth_time; //time of birth;
	double death_time; //time of death;

	bool extant_offspring;
	bool checked;

	bool alloSpeciated;


	void updateHistory(double t);
	bool check_has_viable_offspring(std::vector<species>& v);

	bool operator< (const species& other) const {return       ID < other.ID; }
	bool operator> (const species& other) const {return       ID > other.ID;}
    bool operator<=(const species& other) const {return !operator> (other);} 
	bool operator>=(const species& other) const {return !operator< (other);}
	
	bool operator==(const species& other) const { return ID == other.ID; }
	bool operator!=(const species& other) const {return !((*this) == other);}

	species& operator=(const species& other);

	
};

struct sort_by_birthtime
{
	inline bool operator() (const species& struct1, const species& struct2)
	{
		return (struct1.birth_time < struct2.birth_time);
	}
};



struct point
{
	point(int l, int t);
	int L;
	int T;
};

struct theta
{
	theta();
	theta(double E, double Ss_l, double Ss_h, double As, double J, int m);
	theta(const theta& other);
	theta& operator=(const theta& other);



	double extinct, sym_spec_low, sym_spec_high, allo_spec, jiggle;
    int model;
	
	void changeParams(int time);
    void changeParams();
	void assignParams();

    bool withinPrior();
    void getRandomCombo();

};

struct particle
{
	theta T;
	double weight;
  	std::vector<double> F;
	int lins;
	double specT;
	int numberAlloEvents;
	int numberExtinctions;

    void calculateWeight(const std::vector<particle>& v);
    void calculateWeight(const std::vector< double >& M,
                         const std::vector< std::vector< particle > >& v);


	particle(theta params,std::vector<double>fit, int L, double sT, int nA, int nE);
	particle();

	particle(const particle& other);
	particle& operator=(const particle& other);
};

struct allo_pair
{
	int ID;
	int index_a;
	int index_b;

	allo_pair(int id, int index_A) : ID(id), index_a(index_A)
	{
		index_b = -1;
	}

	allo_pair(int id, int index_A, int index_B) : ID(id),
		index_a(index_A),
		index_b(index_B)
	{
	}

};

struct pair
{
	int ID;
	double dist;
};

struct newick_node
{
	newick_node(const species& S);
	newick_node();

	bool extant;
	int ID;
	int parent;
	double branch_length;

	bool checked;
	int left;
	int right;

	std::string composeString(const std::vector<newick_node>& v) const;
	int getNumberKids(std::vector<newick_node>& v);
};

std::istream& operator >> (std::istream& is, theta& t);
std::ostream& operator << (std::ostream& os, const theta& t);

std::ostream& operator << (std::ostream& os, const particle& p);
std::istream& operator >> (std::istream& is, particle& p);

int drawEvent(double E, double W, double S, double A);

void extinction(std::vector<species>& v, std::vector<species>& extinct_species, double time, int wLevel);
void extinction2(std::vector<species>& v, std::vector<species>& extinct_species, double time, int waterLevel);

void waterLevelChange(std::vector<species>& v, int& wLevel);
void symp(std::vector<species>& v, int i, std::vector<species>& e, int& id_count, double time);

void Symp_speciation(std::vector<species>& v, int& id_count, std::vector<species>& extinct_species, 
					 double time, double waterTime, std::vector<double>& specTimes, int wLevel);
void Allo_speciation(std::vector<species>& v, int& id_count, double time, 
					  double water_time, const std::vector<allo_pair>& p, std::vector<double>& specTimes);

void updateEpsilonVectors();
std::string generateFileName(int i);
void readParticles(int time, std::vector<particle>& particles);
void progressBar(double percent);
void normalizeWeights(std::vector< particle >& p, std::vector<double>& weights);
double calculateWeight(int i, const std::vector<particle>& v, const theta& T);
long double calcK(double theta_i, double theta_j, double sigma);
void convertPointToInt(const std::vector<point>& p, std::vector<double>& R);
void readRealDataPoint(std::vector<spec_point>& R);
void setupVectors(std::vector<spec_point>& real_data);

theta getFromPrevious(const std::vector<double>& weights, const std::vector<particle>& particles);
theta getFromPrevious2(const std::vector<particle>& particles);
theta getFromPrevious3(const std::vector<particle>& particles, double maxWeight);



bool withinLimits(std::vector<double>& F, int t, const std::vector<spec_point>& L, const std::vector<spec_point>& real, int model);

template<typename T> 
void removeDuplicates(std::vector<T>& vec);
void normalizeWeights(std::vector< particle >& p, std::vector<double>& weights);
void divideWeightsbyMax(std::vector< particle >& p);

double get_min_time();


std::vector<double> writeMean(const std::vector< std::vector< int > >& outcomes);
void writeAll(const std::vector< std::vector< int > >& outcomes);
void doSingle(const std::vector<spec_point>& real_data);

void writeWater(const std::vector<std::vector< int> >& W);

template<typename Iterator>  
void bubbleSort(Iterator first, Iterator last);

template<typename T>  
double calcMean(const std::vector<T> v);

std::vector<spec_point> calculateLineages(const std::vector<species>& pop);
void add_histories(std::vector<spec_point>& O, const std::vector<spec_point>& add);
void insert_to_O(std::vector<spec_point>& O, const spec_point& add);
bool present(const std::vector<spec_point>& O, const spec_point& match);

double calcDiff(const std::vector<spec_point>& R, const std::vector< spec_point >& sim, const int& time);
double calculateFit(const std::vector<spec_point>& R, const std::vector<spec_point>& L);
int getLin(const std::vector<spec_point>& v, double t);
int countLineages(const std::vector<spec_point>& v, double time);

std::vector<spec_point> calculateLineages_new(const std::vector<species>& pop, const std::vector<species>& extinct_species);
std::vector<spec_point>  calculateLineages_noextinct(const std::vector<species>& allSpecies);
std::vector<spec_point>  calculateLineages_withextinct(std::vector<species>& allSpecies);

bool sortOnTime(const species& left, const species& right);

std::vector<spec_point> run(const std::vector<spec_point>& real, const std::vector<double>& W, int& id_count, std::vector<species>& allSpecies, int& lins, double& meanSpecTime);

std::vector<double> doRun(const std::vector<spec_point>& real, int& L, int iter, int part_num, std::vector<double>& W, double& beta, double& mST,std::vector<spec_point>& Lins_for_NLTT, int waterModel);

std::vector<spec_point> sumBranches(const std::vector<spec_point>& b1, const std::vector<spec_point>& b2);
void writeOutcomes(const std::vector< std::vector< spec_point > >& v);

void writeLineage(const std::vector<spec_point>& L, int iter, int part_number, double fit);
//void write_newick_file(const std::vector<species>& s1, const std::vector<species>& s2, const std::vector<spec_point>& b1, const std::vector<spec_point>& b2, int iter, int part_number);
void write_newick_file(const std::vector<species>& s1, const std::vector<species>& s2, const std::vector<spec_point>& b1, const std::vector<spec_point>& b2, int iter, int part_number, int jiggle);

void append_newick_file(const std::vector<species>& s1, const std::vector<species>& s2, const std::vector<spec_point>& b1, const std::vector<spec_point>& b2, int iter);

std::vector<double> generateWaterLevelChanges(int model);
int find_parent( int ID, const std::vector<species>& v);

int find_Allo(int i, const std::vector<species>& v);
std::vector<int> findOffspring(int ID, const std::vector<species>& v);
std::vector<int> findOffspring(int ID, const std::vector<newick_node>& v);

template <typename T> 
void sortit(std::vector<T>& v);


void makeFiles(const std::vector<spec_point>& L);

std::vector<newick_node> generateNodeList(const std::vector<species>& v);
double calcNLTT(const std::vector<spec_point>& R, const std::vector<spec_point>& L);
double calcBr(const std::vector<species>& S1, const std::vector<species>& S2);
double calcGamma(const std::vector<spec_point>& L);
bool file_exists(const std::string& name);

void jiggle(std::vector<spec_point>& V, double sigma);
void macstart(const char * argv[]) ;

#endif
