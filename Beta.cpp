#include "Gillespie.h"
#include "Beta.h"
#include "GetParams.h"

#include <map>
#include <boost/math/special_functions/gamma.hpp>

extern GetParams P;

///////////////////////////////////////////////////////////////////////////////////////////////////

double calculateBeta(const std::vector<species>& v1, const std::vector<species>& v2)
{
    if(v1.size() == 1 && v2.size() == 1) return 100.0;
	
	
	std::vector< newick_node > node_list1 = generateNodeList(v1);
	std::vector< newick_node > node_list2 = generateNodeList(v2);

	for(std::vector<newick_node>::iterator it = node_list1.begin(); it != node_list1.end(); ++it)
	{
		(*it).getNumberKids(node_list1);
	}

	for(std::vector<newick_node>::iterator it = node_list2.begin(); it != node_list2.end(); ++it)
	{
		(*it).getNumberKids(node_list2);
	}

	std::vector< p > chances;

	int extant1 = 0;

	BOOST_FOREACH(newick_node i, node_list1)
	{
		if(i.extant == false)
		{
			p newP;
			newP.i = i.left;
			newP.n = i.left + i.right;
			if(newP.i > (0.5 * newP.n)) newP.i = newP.n - newP.i;
			chances.push_back(newP);
		}
		else extant1++;
	}

	int extant2 = 0;

	BOOST_FOREACH(newick_node i, node_list2)
	{
		if(i.extant == false)
		{
			p newP;
			newP.i = i.left;
			newP.n = i.left + i.right;
			if(newP.i > (0.5 * newP.n)) newP.i = newP.n - newP.i;
			chances.push_back(newP);
		}
		else extant2++;
	}

	p newP;
	newP.i = extant1;
	newP.n = extant1 + extant2;
	if(newP.i > (0.5 * newP.n)) newP.i = newP.n - newP.i;
	chances.push_back(newP);
	
	//P.C.clear();
	//P.C = chances;

	/*for(int i = 0; i < node_list1; ++i)
	{
		node_list1[i].
	}*/


	double b = findMinBeta(chances);

	return b;
}

int newick_node::getNumberKids(std::vector<newick_node>& v)
{	
	if(extant == true) return 1;

	if(checked == true) return left + right;

	std::vector<int> index = findOffspring(ID,v);

	left = v[index[0]].getNumberKids(v);
	right = v[index[1]].getNumberKids(v);

	checked = true;
	return left + right;
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////  HERE THE SIMPLEX FMINSEARCH ALGORITHM STARTS ///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// A SEPARATE WORKINGEXAMPLE CAN BE FOUND ///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////// IN THE FMINSEARCH PROJECT /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double findMinBeta(const std::vector<p>& C)
{
	if(C.size() == 1) return 100.0; //only two branches, perfectly balanced
	
	
	double initBeta = 1.0;

	double tolx = 1e-4;
	double tolf = 1e-5;
	int maxIter = (int)(1.25 * 10000);

	std::vector< val > simplex(2);

	simplex[0] = val(initBeta,C);

	simplex[1] = val(initBeta * 1.05,C);

	for(int i = 0; i < maxIter; ++i)
	{
		updateSimplex(simplex,C);

		//check tolerances and quit??
		double diffX = simplex[1].x - simplex[0].x; 
		if(diffX < 0) diffX *= -1.0; //abs difference

		double diffF = simplex[1].f - simplex[0].f; 
		if(diffF < 0) diffF *= -1.0; //abs difference

		if(diffX < tolx && diffF < tolf) break;

		//if(simplex[0].x > 40 && simplex[1].x > 40) break;
	}

	//return simplex[0].x;	
	//double beta = simplex[0].x;
	//return (1.0 * (beta + x) / (beta + 3)); //this returns "B", as defined in Jones 2011, which has properties B = (0,1);
	return simplex[0].x;
}

double val::calcLL(double beta, const std::vector<p>& C)
{
	if(beta <= -2) return 1e6+std::rand(); //  a large number, but never exactly the same number, to avoid getting stuck <-2
//	if(C.size() < 1) return 1e6 + std::rand();

	double loglik = 0.0;
	int maxval = 1 + C.back().n;

	//if(P.C.size() == 1) return 0.0; //P.C[0].n == 2 --> loglik = log(1.0) = 0.0;

	std::vector<double> xn(maxval,0);
	xn[3] = 0.5;

	std::vector<double> sn(maxval,0);
	sn[3] = 2 * exp( gammaln(1 + 1 + beta) + gammaln(1 + 2 + beta) - gammaln(1 + 1) - gammaln(2 + 1));

	for(int i = 4; i < maxval; ++i)
	{
		sn[i] = 1.0 / i * ( i + 1 + 2 * beta + 2.0 * (i - 1  + beta) / (i-1) * xn[i-1]) * sn[i-1];
		xn[i] = (i - 1 + beta)* 1.0 / (i - 1) * xn[i-1] * 1.0 * sn[i-1] / sn[i];
	}

	for(std::size_t i = 0; i < C.size(); ++i)
	{
		int val = C[i].n;
		if(val > (int)sn.size())
		{
			
			std::cout << "\n" << val << "\t" << sn.size() << "\n";
		}

		if(val == 2) loglik += log(1.0);
		if(val == 3) loglik += log(0.5);
		if(val > 3)
		{
			double x = C[i].i;
			double y = val - x;
			loglik += gammaln(1 + x + beta) + gammaln(1 + y + beta) - gammaln(x + 1) - gammaln(y + 1) - log(sn[val]);
		}
	}
	xn.clear();
	sn.clear();
	return -1.0 * loglik;
}


val::val() //default constructor
{
	x = 0.0;
	f = 42;
}

val::val(double valX, const std::vector<p>& C)
{
	x = valX;
	f = calcLL(x,C);
}

val::val(const val& other) //copy constructor
{
	x = other.x;
	f = other.f;
}

val val::operator=(const val& other)
{
	if(this == &other)
		return *this;

	x = other.x;
	f = other.f;
	return *this;
}

bool val::operator<(const val& other) const
{
	if(f < other.f) return true;
	else return false;
}

bool val::operator<=(const val& other) const
{
	if(f <= other.f) return true;
	else return false;
}

bool val::operator>=(const val& other) const
{
	if(f >= other.f) return true;
	else return false;
}

bool val::operator>(const val& other) const
{
	if(f > other.f) return true;
	else return false;
}

bool sortByF(const val& left, const val& right)
{
	if(left < right) return true;
	else return false;
}

void reverse(std::vector<val>& simplex)
{
	val temp = simplex[0];
	simplex[0] = simplex[1];
	simplex[1] = temp;
	return;
}

void updateSimplex(std::vector<val>& simplex, const std::vector<p>& C)
{
	static double rho = 1.0;
	static double chi = 2.0;
	static double gam = 0.5;
	static double sig = 0.5;

	//1. Order
	if(simplex[1].f < simplex[0].f) reverse(simplex); //we do not need to explicitly sort, because the simplex only consists of 2 points!

	//2. reflect, compute the reflection point:
	double mean_x = simplex[0].x; //mean_x is simply first x, because we only have 2 points
	double x_n_and_one = simplex[1].x; //the other x

	double xr = (1 + rho) * mean_x - rho * x_n_and_one;

	val fr(xr,C); //reflection can never happen because the following can never be true: if(simplex[0] <= fr < simplex[0], which is a constraint of optimizing a single variable

	if(fr < simplex[0]) //expand!
	{
		double xe = (1 + rho * chi) * mean_x - rho * chi * x_n_and_one;
		val fe(xe,C); 

		if(fe < fr)
		{
			simplex[1] = fe; //accept xe and terminate iteration //expand
			return;
		}
		else {//accept xr and terminate iteration       //reflect
			simplex[1] = fr;
			return;
		}
	}

	if(fr >= simplex[0]) //contract
	{
		bool shrink = false;
		//a. Outside:
		if(fr >= simplex[0] && fr < simplex[1])
		{
			double xc = (1+ rho * gam) * mean_x - rho * gam * x_n_and_one;
			val fc(xc,C);
			if(fc <= fr) //accept xc
			{
				simplex[1] = fc;
				return;
			} 
			else{shrink = true;}  //move to step 5, a shrink
		}
		if(fr >= simplex[1]) //b. inside contraction
		{
			double xcc = (1 - gam) * mean_x + gam * x_n_and_one;
			val fcc(xcc,C);
			if(fcc < simplex[1])
			{
				simplex[1] = fcc;
				return;
			}
			else shrink = true;
		}

		if(shrink == true) //shrinking
		{
			double new_X = simplex[0].x + sig * (simplex[1].x - simplex[0].x);
			val fv(new_X,C);
			simplex[1] = fv;
			return;
		}
	}

	return;
}

double gammaln(double d)
{
	static std::map<double, double> cache;

	// look up function result
	double ret;
	std::map<double, double>::const_iterator it = cache.find(d);
	if (it != cache.end())
	{
		ret = (*it).second;
	}
	else
	{
		ret = boost::math::lgamma(d);
		if(cache.size() < 50000) cache[d] = ret;
	}
	return ret;
}

