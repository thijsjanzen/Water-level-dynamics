
#include "Gillespie.h"
#include "GetParams.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// CLASS PARAMS MEMBER FUNCTIONS /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

extern GetParams P;

theta::theta()
{
	extinct = -100.0;
	sym_spec_low = -100.0;
	sym_spec_high = -100.0;
	allo_spec = -100.0;
	jiggle = -100.0;
    model = -1;
}


theta::theta(double E, double Ss_l, double Ss_h, double As, double J, int m)
{
	extinct = E;
	sym_spec_low = Ss_l;
    sym_spec_high = Ss_h;
	allo_spec = As;
	jiggle = J;
    model = m;
}

void update_val(double& val, double max)
{
	val += normal(0.0, P.sigma);
	if(val > max) val = max;
}

void update_val(double& val, double max, double min)
{
    val += normal(0.0,P.sigma);
    if(val > max) val = max;
    if(val < min) val = min;
}

void theta::changeParams() {

    double min = -3;
    double max = 2;

    for(int whichParam = 0; whichParam < 5; ++whichParam)
    {
        switch(whichParam)
        {
            case 0: update_val(extinct,    max, min);
                break;
            case 1: update_val(sym_spec_low,   max, min);
                break;
            case 2: update_val(sym_spec_high,   max, min);
                break;
            case 3: update_val(allo_spec,  max, min);
                break;
            case 4: update_val(jiggle, max, min);
                break;
            default: break;
        }
    }

    int r = random_number(4);
    if(r == 0) {
        model++;
    }
    if(r == 1) {
        model--;
    }
    if(model > 2) model -= 3;
    if(model < 0) model += 3;


   /* double r = uniform();
    if(r < 0.05) {
        model++;
    }
    if(r > 0.95) {
        model--;
    }*/

    if(model > 2) model -= 3;
    if(model < 0) model += 3;
}


void theta::assignParams()
{
	P.extinctionRate = pow(10.0,extinct);	
	P.SympatricRate_low = pow(10.0,sym_spec_low);
	P.SympatricRate_high = pow(10.0,sym_spec_high);
	P.AllopatricRate = pow(10.0,allo_spec);
    P.Jiggle = pow(10.0,jiggle);

	return;
}


void theta::getRandomCombo()
{
	extinct = -3 + 5 * uniform();
    sym_spec_low = -3 + 5 * uniform();
	sym_spec_high = -3 + 5 * uniform();
	allo_spec = -3 + 5 * uniform();
	jiggle = -3 + 3 * uniform();
    model = random_number(3);

	return;
}

bool theta::withinPrior()
{
    if(extinct < -3) return false;
    if(extinct >  2) return false;

    if(sym_spec_low < -3) return false;
    if(sym_spec_low >  2) return false;

    if(sym_spec_high < -3) return false;
    if(sym_spec_high >  2) return false;

    if(allo_spec < -3) return false;
    if(allo_spec >  2) return false;

    if(jiggle < -3) return false;
    if(jiggle >  0) return false;

    if(model < 0) return false;
    if(model > 2) return false;

    return true;
}




theta getFromPrevious(const std::vector<double>& weights, const std::vector<particle>& particles)
{
	const double r = uniform();

	int min = 0;
	int max = (int)weights.size()- 1;
	int med = (int)((max+min)*0.5);

	while((max-min) > 1)
	{
		if(weights[med] >= r) max = med;
		else min = med;

		med = (int)((max+min)*0.5);
	}

	return particles[med].T;
}

theta getFromPrevious2(const std::vector<particle>& particles)
{
	for(int i = 0; i < 1e6; ++i)    {
		int index = random_number((int)particles.size());
		if(uniform() < particles[index].weight)
		{
			return particles[index].T;
		}
	}
	
	//this code should never be reached, but just in case:
	int index = random_number((int)particles.size());
	return particles[index].T;
}

theta getFromPrevious3(const std::vector<particle>& particles, double maxWeight)
{
    for(int i = 0; i < 1e6; ++i)    {
        int index = random_number((int)particles.size());
        if(uniform() < particles[index].weight / maxWeight)
        {
            return particles[index].T;
        }
    }

    //this code should never be reached, but just in case:
    int index = random_number((int)particles.size());
    return particles[index].T;
}

theta::theta(const theta& other) //copy constructor
{
	extinct = other.extinct;
	sym_spec_low = other.sym_spec_low;
	sym_spec_high = other.sym_spec_high;
	allo_spec = other.allo_spec;
	jiggle = other.jiggle;
    model = other.model;
}

theta& theta::operator=(const theta& other)
{
	if(this == &other) return *this;

	extinct = other.extinct;
	sym_spec_low = other.sym_spec_low;
	sym_spec_high = other.sym_spec_high;
	allo_spec = other.allo_spec;
	jiggle = other.jiggle;
    model = other.model;

	return *this;
}


std::istream& operator >> (std::istream& is, theta& t)
{
	is >> t.extinct;
	is >> t.sym_spec_low;
	is >> t.sym_spec_high;
	is >> t.allo_spec;
	is >> t.jiggle;
    is >> t.model;
	return is;
}

std::ostream& operator << (std::ostream& os, const theta& t)
{
	os << t.extinct << '\t';
	os << t.sym_spec_low << '\t';
	os << t.sym_spec_high << '\t';
	os << t.allo_spec << '\t';
	os << t.jiggle << '\t';
    os << t.model << '\t';
	return os;
}

particle::particle(theta params, std::vector<double> fit, int L, double sT, int nA, int nE)
{
	T = params;
	F = fit;
	lins = L;
	specT = sT;
	numberAlloEvents = nA;
	numberExtinctions = nE;
	weight = -1;
}

particle::particle()
{
	std::vector<double> tempF(4,1e6);
	F = tempF;
	weight = -1;
	T = theta();
	lins = 0;
	specT = 0;
    numberAlloEvents = 0;
    numberExtinctions = 0;
}

particle::particle(const particle& other)
{
	T = other.T;
	F = other.F;
	weight = other.weight;
	lins = other.lins;
	specT = other.specT;
	numberAlloEvents = other.numberAlloEvents;
	numberExtinctions = other.numberExtinctions;
}

particle& particle::operator=(const particle& other)
{
	if(this == &other) return *this;

	T = other.T;
	F = other.F;
	weight = other.weight;
	lins = other.lins;
	specT = other.specT;
	numberAlloEvents = other.numberAlloEvents;
	numberExtinctions = other.numberExtinctions;

	return *this;
}

int isnan(double x) { return x != x; }

std::istream& operator >> (std::istream& is, particle& p)
{

	is >> p.T.extinct;
	is >> p.T.sym_spec_low;
	is >> p.T.sym_spec_high;
	is >> p.T.allo_spec;
	is >> p.T.jiggle;
    is >> p.T.model;
	
	is >> p.F[0];
	is >> p.F[1];
	is >> p.F[2];
	is >> p.lins;
	is >> p.specT;
	is >> p.numberAlloEvents;
	is >> p.numberExtinctions;
	is >> p.weight;

	return is;
}

std::ostream& operator << (std::ostream& os, const particle& p)
{
	os << p.T.extinct << '\t';
	os << p.T.sym_spec_low << '\t';
	os << p.T.sym_spec_high << '\t';
	os << p.T.allo_spec << '\t';
	os << p.T.jiggle << '\t';
    os << p.T.model << '\t';
	
	os << p.F[0] << '\t';
	os << p.F[1] << '\t';
	os << p.F[2] << '\t';
	os << p.lins << '\t';
	os << p.specT << '\t';
	os << p.numberAlloEvents << '\t';
	os << p.numberExtinctions << '\t';
	os << p.weight;
	

	return os;
}

