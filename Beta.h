#ifndef BETA_H
#define BETA_H

#include <vector>
#include "Gillespie.h"


struct p  //this struct holds the number of branches and tips
{
	p(): i(-1),n(-1)
	{
	}


public:
	int i;
	int n;
};

struct val
{
	val();
	val(double valX, const std::vector<p>& C);
	val(const val& other);
	double x;
	double f;
	bool operator<(const val& other) const;
	bool operator>(const val& other) const;
	bool operator<=(const val& other) const;
	bool operator>=(const val& other) const;
	val operator=(const val& other);

	double calcLL(double beta, const std::vector<p>& C);
};



double calculateBeta(const std::vector<species>& v1, const std::vector<species>& v2);

double gammaln(double d);
double findMinBeta(const std::vector<p>& C);
void updateSimplex(std::vector<val>& simplex, const std::vector<p>& C);
bool sortByF(const val& left, const val& right);







#endif