// Gillespie.cpp : Defines the entry point for the console application.
//

#include "GetParams.h"
#include "Gillespie.h"
#include "Beta.h"
#include <math.h> 
#include <random>

/////// PARAMETER CONTAINER //////////////////////////////////////////////
GetParams P;
//////////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[])
{
    macstart(argv);
	
	std::cout << "Starting reading config file\n";
	P.readFromIni("config.ini"); //read the config file
	
	
	std::random_device rdev;
	unsigned int chosen_seed = rdev();

	set_seed(chosen_seed);

	std::ofstream seedFile("seed.txt",std::ios::app);
	seedFile << chosen_seed << "\n";
	seedFile.close();
	

	std::vector< spec_point> real_data; //vector to contain the empirical data
	std::vector<double> W; //vector containing the water level
	std::vector<double> real; //log converted empirical data

	updateEpsilonVectors();

	std::vector<particle> particles;
	
	boost::timer t1;
	std::cout << "starting to do SMC!\n";
	std::size_t t = 0;
    int numberAccepted = 0;
    int numberProposed = 0;

	
	P.emp_values.resize(4,0);

    // read old acceptance rate
    std::ifstream acceptFile("acceptFile.txt");
    if(file_exists("acceptFile.txt")) {
        std::cout << "found old acceptance file, reading acceptance file: \n";
        while(!acceptFile.eof()) {
            std::string s_iter;
            acceptFile >> s_iter;
            t = atoi(s_iter.c_str());

            std::string s_numberAccepted;
            acceptFile >> s_numberAccepted;
            numberAccepted = atoi(s_numberAccepted.c_str());

            std::string s_numberProposed;
            acceptFile >> s_numberProposed;
            numberProposed = atoi(s_numberProposed.c_str());

            std::string s_temp;
            acceptFile >> s_temp;
            std::cout << t << "\t" << numberAccepted << "\t" << numberProposed << "\t" << s_temp << "\n";
        }
    }


	std::string f_name = "summaryStatistics.txt";
	if(file_exists(f_name)) {
		
		std::ifstream inputFile("summaryStatistics.txt");
		for(std::size_t i = 0; i < 4; ++i)  {
			double temp;
			inputFile >> temp;
			P.emp_values[i] = temp;
		}
			
		double temp2;
		inputFile >> temp2;
		if(temp2 > 1) P.maxT = temp2;
		
		for(int k = 500; k > -1; --k) {
			std::string f_name = generateFileName(k);
			if(file_exists(f_name)) {
				t = k;
				std::cout << "reading previously finished run of time \t" << t << "\n";
				particles.clear();
				readParticles((int)t,particles);
				std::cout << "done reading\n";
				break;
			}
		}

       		if((int)particles.size() == P.numberParticles) {
			++t; //the previously finished run was fully done!
		}
	}

    
	if(P.generateData == 1) {
		std::string f_name_alt = "lampro_alt.txt";
		if (!file_exists(f_name_alt)) //if lampro_alt already exists, a previous run has generated data already!
        {
            std::cout << "Generating data to fit to\n";
			int tries = 1;
			bool data_found = false;
			while(tries < 1e6 && !data_found)
			{
				theta Params = getRandomCombo();
				Params.assignParams();
				tries++;
			
				for(int j = 0; j < 100;++j)
				{
					int lins;
					double beta = 2;
					int t = 0;

					double meanSpecTime = 0;
					std::vector<spec_point> lineages;

					std::vector<double> F = doRun(real_data, lins,t,j,W,beta,meanSpecTime,lineages, P.waterModel);

				  if(F[1] < 1e6 && F[2] < 1e6 && lins > 10) 
				  {
					std::cout << "we have found a match! starting to generate required files\n";
					makeFiles(lineages);
					F[0] = 0.0; //nLTT always check against 0
					P.emp_values = F;
					

					std::cout << "Fitting to a tree with\n";
					std::cout << "nLTT of  " << F[0] << "\n";
					std::cout << "gamma of " << F[1] << "\n";
					std::cout << "mBr of   " << F[2] << "\n"; 
					std::cout << "#lins   " <<  F[3] << "\n";

					std::ofstream statisticsFile("summaryStatistics2.txt");
					for(std::size_t index = 0; index < F.size(); ++index)
					{
						statisticsFile << F[index] << "\n";
						P.emp_values[index] = F[index];
					}
					statisticsFile << P.maxT << "\n";
					statisticsFile.close();
					data_found = true;
					break;
				  }

				}
			}
        } else {
            std::ifstream alternativeSummary("summaryStatistics2.txt");
            for(std::size_t index = 0; index < P.emp_values.size(); ++index)
            {
                alternativeSummary >> P.emp_values[index];
            }
            alternativeSummary >> P.maxT;
            alternativeSummary.close();
        }
	}

	setupVectors(real_data,real); //this function reads either lampro_alt if recheck == 1, or lampro4.txt if else.
	
	for(; t < 10; ++t) //smaller than 10, to not let travis run until eternity
	{	
		// std::vector<particle> previousParticles = particles;
		P.t = (int)t;
	//	if(P.t > 9) break;

		if((int)particles.size() != 0 && (int)particles.size() != P.numberParticles) numberAccepted = (int)particles.size();
		particles.clear();

        std::vector< std::vector< particle > > previousParticles(3); // per model
        std::vector< double > M(3,0); //holds the frequencies of each model.
        std::vector< double > maxWeights(3,-1);
		
		if(t!=0)
		{
			readParticles((int)(t-1),particles);
			
			if(particles.empty())
			{
				std::ofstream partfile("PARTICLES EMPTY.txt");
				partfile << "could not read particles file...";
				partfile.close();
				break;
			}

            for(int i = 0; i < particles.size(); ++i) {
                int modelType = particles[i].T.model;
                if(modelType >= 0 && modelType < (int)previousParticles.size()) {
                    previousParticles[modelType].push_back(particles[i]);
                }
            }


            double sum = 0;
            for(int i = 0; i < 3; ++i) {
                M[i] = (int)previousParticles[i].size();
                sum += M[i];
            }
            for(int i = 0; i < 3; ++i) {
                M[i] = M[i] / sum;
            }


            // normalize the weights per model
            if(!previousParticles.empty()) {
                for(int i = 0; i < 3; ++i) {
                    std::vector < particle > temp = previousParticles[i];
                    double sum = 0.0;
                    for(auto it = temp.begin(); it != temp.end(); ++it) {
                        sum += (*it).weight;
                    }
                    for(auto it = temp.begin(); it != temp.end(); ++it) {
                        (*it).weight = (*it).weight / sum;
                        if(maxWeights[i] < (*it).weight) {
                            maxWeights[i] = (*it).weight;
                        }
                    }
                    previousParticles[i] = temp;
                }
            }

		}




		std::string f_name = generateFileName((int)t);
		std::ofstream out_particle(f_name.c_str(), std::ios::app);
		out_particle.precision(30);

		boost::timer time_of_iteration;
		for(; numberAccepted < P.numberParticles; numberProposed++)
		{
			theta Params;
			if(t == 0) {
				Params = getRandomCombo(); 
			}   else    {

                int chosen = 0;
                double r = uniform();
                for(int i = 0; i < 3; ++i) {
                    r -= M[i];
                    if(r <= 0) {
                        chosen = i;
                        break;
                    }
                }

				Params = getFromPrevious(previousParticles[chosen], maxWeights[chosen]);
				Params.changeParams();
			}
			
			Params.assignParams();
			
			int lins = -1;
			double beta = 100;
			double meanSpecTime = 0;

			std::vector<spec_point> lineages;

			//////////////////////////////////////////////////////////////////////
			std::vector<double> F(4,1e6);
			if(t != 0) F = doRun(real_data, lins,(int)t,numberProposed,W,beta,meanSpecTime,lineages, Params.model);
			///////////////////////////////////////////////////////////////////////

			if(withinLimits(F,(int)t,lineages,real_data, Params.model) || t == 0)
			{
				particle add(Params,F,lins,meanSpecTime,P.numAllo,P.numExtinct);

                if(t == 0) {
                    add.weight = 1.0;
                }   else {
					add.calculateWeight(M,previousParticles);
				}
				
				out_particle << add << "\n"; out_particle.flush();

				numberAccepted++;

				std::cout << "iteration " << t << "\t" << numberAccepted << " of " << P.numberParticles << " accepted\n";
			}

            if(numberProposed % 1000 == 0) {
                // update acceptance rate and write to file in case another simulation will have to restart after this one
                double acceptRate = 1.0 * numberAccepted/ numberProposed;
                std::ofstream acceptFile("acceptFile.txt",std::ios::app);
                acceptFile << t << "\t" << numberAccepted << "\t" << numberProposed << "\t" << acceptRate << "\n";
                acceptFile.close();
                if(1.0 * numberProposed / numberAccepted > 1e6 && numberProposed > 1e6) {
                    // accept rate drops below 1 in a million
                    std::ofstream stopFile("stopFile.txt");
                    stopFile << "Acceptance rate below: " << numberAccepted  << " in " << numberProposed << "particles\n";
                    stopFile.close();
                    return 0;
                }
            }
		}

		std::cout << "\niteration " << t << "  " << "iterations: " << numberProposed << "  " << "acceptance rate  " << 100.0 * P.numberParticles / numberProposed << "%  " << time_of_iteration.elapsed() << " " << t1.elapsed() << "\n";
        std::ofstream acceptRate("acceptRate.txt",std::ios::app);
        acceptRate << t << "\t" << numberProposed << "\t" << P.numberParticles << "\t" << 100.0 * P.numberParticles / numberProposed << "\t" << time_of_iteration.elapsed() << "\n";
        acceptRate.close();
        numberAccepted = 0;
        numberProposed = 0;

		out_particle.flush();
		out_particle.close();
	}

	return 0;
}


template <typename T>
double calcMean(const std::vector<T>& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = 1.0 * sum / v.size();
    return(mean);
}



void jiggle(std::vector<species>& S)    {
	std::vector<int> jiggled((int)S.size(),0);
	for(int i = 0; i < (int)S.size(); ++i)
	{
		double lowerLim = 0.0;
		double upperLim = P.maxT;
		double bTime = S[i].birth_time;
		if(bTime > 0 && jiggled[i] == 0)
		{
			for(int j = 0; j < (int)S.size(); ++j)
			{
				if(S[j].ID == S[i].parent) {
					lowerLim = S[j].birth_time;
				}
				
				if(S[j].parent == S[i].ID) {
					if(S[j].birth_time < upperLim) {
						upperLim = S[j].birth_time;
					}
				}
			}
			
			double newBtime = bTime + normal(0,P.Jiggle);
			while(newBtime < lowerLim || newBtime > upperLim)
			{
				double stdev = P.Jiggle;
				if( (bTime - lowerLim) < stdev) stdev = bTime - lowerLim;
				if( (upperLim - bTime) < stdev) stdev = upperLim - bTime;
				newBtime = bTime + normal(0,stdev);
			}
			
			//find other species with same birth time and parent
			int comrade = -1;
			for(int j = 0; j < (int)S.size(); ++j)
			{
				if(j != i) {
					if(S[j].parent == S[i].parent)
					{
						if(S[j].birth_time == S[i].birth_time)
						{
							comrade = j;
							break;
						}
					}
				}
			}
			
			if(comrade > 0) {
				for(int j = 0; j < (int)S.size(); ++j)
				{
					if(S[j].parent == S[comrade].ID) {
						if(S[j].birth_time < upperLim) {
							upperLim = S[j].birth_time;
						}
					}
				}
				newBtime = bTime + normal(0,P.Jiggle);
				while(newBtime < lowerLim || newBtime > upperLim)
				{
					double stdev = P.Jiggle;
					if( (bTime - lowerLim) < stdev) stdev = bTime - lowerLim;
					if( (upperLim - bTime) < stdev) stdev = upperLim - bTime;
					newBtime = bTime + normal(0,stdev);
				}
				S[comrade].birth_time = newBtime;
				jiggled[comrade] = 1;
			}
			
			S[i].birth_time = newBtime;
		}
		jiggled[i] = 1;
	}

	return;
}




std::vector<double> doRun(const std::vector<spec_point>& real, int& L, int iter, int part_num, std::vector<double>& W, double& beta, double& mST,std::vector<spec_point>& Lins_for_NLTT, int waterModel)
{
    std::vector<double> summaryStats = {1e6,1e6,1e6,1e6}; //default values very, very high
	int idCount = 0;
	W.clear();
	W = generateWaterLevelChanges(waterModel);
	L = 0;
	int lins = 0;
	std::vector<species> s1;
	int numAlloToReturn = 0;
	int numExtinctToReturn = 0;

	static int counter = 0;
	counter++;

	double meanSpecTime = 0;;
	mST = 0;

	std::vector<spec_point> branch1 = run(real,W,idCount,s1,lins,meanSpecTime);
	for(int i = 0; branch1.empty() && i < 100; ++i)
	{ 
		lins = 0;
		branch1 = run(real,W,idCount,s1,lins,meanSpecTime);
	}
	L = lins; 
	if(lins > 1) 
	{
			mST = meanSpecTime;
	} else {
		mST = 0;
	}
	numAlloToReturn += P.numAllo;
	numExtinctToReturn += P.numExtinct;

	int b1Extinct = P.numExtinct;
	
	if(branch1.empty()) { 
		summaryStats[3] = L;
		return summaryStats;
	}

	std::vector<species> s2;

	std::vector<spec_point> branch2 = run(real,W,idCount,s2,lins,meanSpecTime);
	for(int i = 0; branch2.empty() && i < 100; ++i)
	{
		    lins = L;	    
			branch2 = run(real,W,idCount,s2,lins,meanSpecTime);
	}
	int b2Extinct = P.numExtinct - b1Extinct;
	if((lins-L) > 1) { 
		mST += meanSpecTime;
	}
	
	L = lins;
	
	mST *= 0.5;
	if(L == 2) mST = -1;
	numAlloToReturn += P.numAllo;
	numExtinctToReturn += P.numExtinct;
	P.numAllo = numAlloToReturn;
	P.numExtinct = numExtinctToReturn;
	
	if(branch2.empty() || branch1.empty()) { summaryStats[3] = L; return summaryStats;}
	
	
	jiggle(s1);
	jiggle(s2);
	
	if(b1Extinct  == 0) branch1 = calculateLineages_noextinct(s1);
	else
	{
		branch1 = calculateLineages_withextinct(s1);
	}
	
	if(b2Extinct  == 0) branch2 = calculateLineages_noextinct(s2);
	else
	{
		branch2 = calculateLineages_withextinct(s2);
	}
	
	
	std::vector<spec_point> lineages = sumBranches(branch1,branch2);

	L = (int)lineages.back().ID;

	if(L > 2 || iter == 0) {
	
		//if(P.runSingle == true) write_newick_file(s1,s2,branch1,branch2,iter,part_num);

		//fit[0] = calcNLTT(real,lineages);
		Lins_for_NLTT = lineages;
	    summaryStats[1] = calcGamma(lineages);
		summaryStats[2] = calcBr(s1,s2);
		summaryStats[3] = L;


		if( (summaryStats[1] < 1e6 && summaryStats[2] < 1e6) || iter == 0)
		{
			if(iter > 10) {
                    if(withinLimits(summaryStats, iter, lineages, real, waterModel))
					{
						std::string iteration = boost::lexical_cast<std::string>(iter);
						std::string f_name = "lineages_t=" + iteration + ".txt";
						std::ofstream lineagesFile(f_name.c_str(), std::ios::app);	
						for(std::size_t k = 0; k < lineages.size(); ++k) 
						{
							lineagesFile << lineages[k].time << " ";
						}
						lineagesFile << "\n";
						lineagesFile.close();

						append_newick_file(s1,s2,branch1,branch2,iter);
						
						std::string f_name2 = "water_t=" + iteration + ".txt";
						std::ofstream waterFile(f_name2.c_str(), std::ios::app);
						for(std::vector<double>::iterator w = W.begin(); w != W.end(); ++w)
						{
							waterFile << (*w) << "\t";
						}
						waterFile << "\n";
						waterFile.close();
                    }
            }
        }
	}
	return summaryStats;
}

int drawEvent(double E, double S, double A)
{
	double sum = E + S + A;
	double events[3] = {E/sum, S/sum, A/sum};
	double r = uniform();
	for(int i = 0; i < 3; ++i)
	{
		r -= events[i];
		if(r <= 0) return i; 
	}
	return 0;
}

bool onlyInstance(const std::vector<species>& v, int i)
{
	if(v.size() < 2) return true;

	int local_ID = v[i].ID;
	int count = 0;
	
	for(std::vector<species>::const_iterator s = v.begin(); s != v.end(); ++s)
	{
		if((*s).ID == local_ID) count++;

		if(count > 1) return false;
	} 
	return true;
}

void extinction(std::vector<species>& v, std::vector<species>& extinct_species, double time, int wLevel)
{
	
	int i = random_number((int)v.size());
	
	v[i].death_time = time;

	if(wLevel == 0) //low water level, there might be two instances of the same species
	{
		if(onlyInstance(v,i)) extinct_species.push_back(v[i]);
	}
	if(wLevel == 1) //high water level, there is only one instance of this species
	{
		extinct_species.push_back(v[i]);
	}

	v[i] = v.back();
	v.pop_back();	
	
	return;
}


void Symp_speciation(std::vector<species>& v, int& id_count, std::vector<species>& extinct_species, double time, double waterTime, std::vector<double>& specTimes, int wLevel)
{
	int i = random_number((int)v.size());
	if(wLevel == 0)
	{
		if(!onlyInstance(v,i))
		{
			
			v[i].parent = v[i].ID;
			v[i].ID = id_count; id_count++;
			v[i].birth_time = waterTime;
		}
	}

	v[i].death_time = time;

	double initiationTime = waterTime;	if(v[i].birth_time > initiationTime) initiationTime = v[i].birth_time; //this is the case if we have a sympatric speciation event UPON another speciation event
	specTimes.push_back(time - initiationTime);
	

	//generate offspring
	species offspring1 = species(v[i],id_count,time);
	species offspring2 = species(v[i],id_count,time);

	extinct_species.push_back(v[i]);
	v[i] = offspring1;
	v.push_back(offspring2);

	return;
}

void waterLevelChange(std::vector<species>& v, int& wLevel)
{
	if(wLevel == 0) //waterlevel is low, and rises
	{
		removeDuplicates(v);
	}

	if(wLevel == 1) //waterlevel is high and drops
	{
		std::vector<species> copy = v;
		v.insert(v.end(),copy.begin(),copy.end());  //all species are distributed over the two pockets, e.g. doubled
		//std::move(copy.begin(),copy.end(),std::back_inserter(v));
	}

	wLevel = 1 - wLevel;
	return;
}

void Allo_speciation(std::vector<species>& v, int& id_count, double time, double water_time, const std::vector<allo_pair>& p, std::vector<double>& specTimes)
{
	int i = random_number((int)p.size());

	int parent = p[i].index_a;

	species offspring = species(v[parent],id_count,water_time);
	
	v[parent] = offspring;

	specTimes.push_back(time-water_time);
	
	return;
}

int findinP(const std::vector<allo_pair>& p, int ID)
{
	
	int count = 0;
	for(std::vector<allo_pair>::const_iterator it = p.begin(); it != p.end(); ++it)
	{
		if((*it).ID == ID) return count;
		count++;
	}
	
	return -1;
}


void updatePairs2(std::vector<species>& v, std::vector<allo_pair>& p)
{
	p.clear();
	if(v.size() == 2)
	{
		if(v[0].ID == v[1].ID)
		{
			p.push_back(allo_pair(v[0].ID, 0));
			p[0].index_b = 1;
		}
		return;
	}
	
    std::sort(v.begin(),v.end());

	for(std::size_t i = 1; i < v.size(); ++i)
	{
		if(v[i].ID == v[i-1].ID)
		{
			p.push_back(allo_pair(v[i].ID,(int)i-1,(int)i));
			++i;
		}
	}
	return;
}

std::vector<spec_point> run(const std::vector<spec_point>& real, const std::vector<double>& W, int& id_count, 
							 std::vector<species>& allSpecies, int& lins, double& meanSpecTime)
{
	P.numAllo = 0;
	P.numExtinct = 0;
	
	std::vector<double> speciationCompletionTimes;
	id_count+=1e6; //just to be sure

	allSpecies.clear();

	std::vector<spec_point> lineages;

	std::vector<species> pop;
	pop.push_back(species(id_count));
	
	int numberExtinctions = 0;
	
	double time = 0;
	int waterLevel = 1;

	int iter = 0;

	int numberWlevelChanges = 0;
	
	
	std::vector<species> extinct_species;
	std::vector<species> allo_species;

	//int alloSpeciations = 0;
	std::vector<allo_pair> pairs;

	while(time < P.maxT)
	{
		
		iter ++;

		double Pe = P.extinctionRate * pop.size();

		double Ps = 0.0;
		
		if(waterLevel == 0) Ps = P.SympatricRate_low * pop.size();
		if(waterLevel == 1) Ps = P.SympatricRate_high * pop.size();
		
		if(waterLevel == 0) updatePairs2(pop,pairs);

		double Pa = (1 - waterLevel) * P.AllopatricRate * (pairs.size() * 2);

		double rate = Pe + Ps + Pa;

		double timestep = Expon(rate);
		
		time += timestep;
	
		if(time >= W[numberWlevelChanges] && (time-timestep) < W[numberWlevelChanges])
		{
			time = W[numberWlevelChanges];
			numberWlevelChanges++;
			waterLevelChange(pop,waterLevel);
		}
		else
		{
			if(time > P.maxT) break;
			int event_chosen = drawEvent(Pe,Ps,Pa);
		
			double time_of_previous_waterlevelchange = 0.0;
			if(numberWlevelChanges != 0) time_of_previous_waterlevelchange = W[numberWlevelChanges-1];

			switch(event_chosen)
			{
				case 0:
						extinction(pop, extinct_species, time, waterLevel);
						numberExtinctions++;
						break;
				case 1:
						Symp_speciation(pop, id_count, extinct_species, time, time_of_previous_waterlevelchange,speciationCompletionTimes, waterLevel);
						break;
				case 2: 
						Allo_speciation(pop, id_count,time, time_of_previous_waterlevelchange,pairs,speciationCompletionTimes);
						P.numAllo++;
						break;
			}
		}

		if(pop.size() < 1) //everything is extinct
		{
			lins += 0;
			P.numExtinct += extinct_species.size();
			return lineages;
		}

		if(pop.size() > 300) //more than 300 species, unlikely to provide a good fit, but slows down the program considerably
		{
			lins += pop.size();
			P.numExtinct += extinct_species.size();
			return lineages;
		}
		
	}
	


	removeDuplicates(pop); //we remove the duplicates because there might be duplicates due to a low water level stand, these are not interesting (so in effect, we force the simulation to end with high water)

	allSpecies.clear();
	allSpecies = pop;
	//allSpecies.insert(allSpecies.end(), extinct_species.begin(), extinct_species.end());
	P.numExtinct += extinct_species.size();
	std::move(extinct_species.begin(),extinct_species.end(),std::back_inserter(allSpecies));


	if(numberExtinctions == 0) lineages = calculateLineages_noextinct(allSpecies);
	else 
	{
		lineages = calculateLineages_withextinct(allSpecies);
	}
	
	lins += lineages.back().ID;
	meanSpecTime = 0;
	for(std::vector<double>::iterator i = speciationCompletionTimes.begin(); i != speciationCompletionTimes.end(); ++i) meanSpecTime += (*i);

	meanSpecTime *= 1.0 / speciationCompletionTimes.size();

	return lineages;
}


///////////////////////////////////////////////////////////////////////////
////////////////////////// particle stuff /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

long double calcK(double theta_i, double theta_j, double sigma)
{
    static double sqrt_two_pi = sqrt(2 * M_PI);
    double prefactor = 1.0 / (sigma * sqrt_two_pi);

    double numerator = (theta_i - theta_j) * (theta_i - theta_j);
    double denumerator = 2 * sigma * sigma;

    long double output = 1.0 * prefactor * exp(-1.0 * numerator / denumerator);

    return (output);
}

void particle::calculateWeight(const std::vector< double >& M,
                               const std::vector< std::vector< particle > >& v) {

    double sum = 0.0;

    for(int i = 0; i < 3; ++i)  {

        double localSum = 0.0;
        for(int j = 0; j < (int)v[i].size(); ++j)
        {
            long double localProduct = 1.0;
            for(int whichParam = 0; whichParam < 5; ++whichParam) {
                switch(whichParam)
                {
                    case 0: localProduct *= calcK(v[i][j].T.extinct, T.extinct, P.sigma);
                        break;
                    case 1: localProduct *= calcK(v[i][j].T.sym_spec_low, T.sym_spec_low, P.sigma);
                        break;
                    case 2: localProduct *= calcK(v[i][j].T.sym_spec_high, T.sym_spec_high, P.sigma);
                        break;
                    case 3: localProduct *= calcK(v[i][j].T.allo_spec, T.allo_spec, P.sigma);
                        break;
                    case 4: localProduct *= calcK(v[i][j].T.jiggle, T.jiggle, P.sigma);
                        break;
                    default: break;
                }
            }

            localProduct *= v[i][j].weight; // this is already normalized
            localSum += localProduct;
        }

        double prodM = 0.0;
        for(int i = 0; i < 3; ++i) {
            double factor = 0.05;
            if(i == T.model) {
                factor = 0.9;
            }
            prodM += factor * M[i];
        }

        sum += prodM * localSum;
    }

    weight =  (1.0 / sum);
}





void normalizeWeights(std::vector< particle >& p, std::vector<double>& weights)
{
	weights.clear();
	double sum = 0.0;
	//for(std::vector<particle>::iterator i = p.begin(); i != p.end(); ++i) sum += (*i).weight;
	BOOST_FOREACH(particle i, p) sum += i.weight;

	double cumul_weight = 0.0;

	for(std::vector<particle>::iterator i = p.begin(); i != p.end(); ++i)
	{
		cumul_weight += (*i).weight * 1.0 / sum;
		weights.push_back(cumul_weight);
	}
	return;
}

void divideWeightsbyMax(std::vector< particle >& p)
{
	double maxW = -1.0;
	for(std::vector<particle>::iterator i = p.begin(); i != p.end(); ++i)
	{
		if((*i).weight > maxW) maxW = (*i).weight;
	}
	
	for(std::vector<particle>::iterator i = p.begin(); i != p.end(); ++i)
	{
		(*i).weight = 1.0 * (*i).weight / maxW;
	}
	return;
}



void progressBar(double percent)
{
	if(percent == 0) std::cout << "\n";

	int number = (int)(1.0 * percent / 100 * 20);
		
	std::cout << "\r"; //clear line
	for(int i = 0; i <= number; ++i)
	{
		std::cout << char(219);
	}
	std::cout << " " << percent << "%";
	return;
}

void readParticles(int time, std::vector<particle>& particles)
{
	particles.clear();

	std::string f_name = generateFileName(time);
	std::ifstream read_part(f_name.c_str());
	
	std::cout << "attempting to read \t" << f_name << "\n";
	
	if(!read_part.is_open()) 
	{
		std::cout << "No input file of particles found!!!!\n";
		return;
	}
	int counter = 0;
	while(!read_part.eof())
	{
		particle temp;
		read_part >> temp;
		particles.push_back((temp));
		counter++;
		//std::cout << temp << "\n";
		
		if(counter > P.numberParticles + 1) {
			//std::cout << "STUCK IN READING LOOP!!!";
			//exit(0);
			break;
		}
	}
	std::cout << "Done reading particles\n";
	while(particles.size() > P.numberParticles) {
		particles.pop_back();
	}
	return;
}

std::string generateFileName(int i)
{
	std::string time = boost::lexical_cast<std::string>(i);
	std::string f_name = "particles_t=" + time + ".txt";
	return f_name;	
}

void updateEpsilonVectors()
{
	/*double initLTT = 0.2;
	double initGamma = 1;
	double initBr = 1;
	double initNum = 50;
	
	std::vector<double> epsLTT;
	std::vector<double> epsG;
	std::vector<double> epsBr;
	std::vector<double> epsN;

	for(int i = 0; i < 100; ++i)
	{
		epsLTT.push_back( initLTT * exp(-0.25 * i));
		epsG.push_back(   initGamma * exp(-0.25 * i));
		epsBr.push_back( initBr * exp(-0.25 * i));
		epsN.push_back( initNum * exp(-0.25 * i));
	}

	P.vec_epsilon.push_back(epsLTT);
	P.vec_epsilon.push_back(epsG);
	P.vec_epsilon.push_back(epsBr);
	P.vec_epsilon.push_back(epsN);*/


    for(int i = 0; i < 20; ++i)
    {
        P.vec_epsilon.push_back( 3 * exp(-0.25 * i));
    }
	return;
}

///////////////////////////////////////////////////////////////////////////////////
//////////////////////// MEMBER FUNCTIONS /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

species::species(int& id) //constructor
{
	ID = id;  //every species has a unique ID to keep track of it, also over time
	id++;
	birth_time = 0;
	death_time = -1;
	parent = -1;
	alloSpeciated = false;
	checked = false;
}

species& species::operator=(const species& other)
{
	if(this == &other) return *this;

	ID = other.ID;
	birth_time = other.birth_time;
	death_time = other.death_time;
	parent = other.parent;
	alloSpeciated = other.alloSpeciated;
	extant_offspring = other.extant_offspring;
	checked = other.checked;

	return *this;
}

point::point(int l, int t) //container to store LTT data, where lineage and time is linked.
{
	L = l; T =t;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void setupVectors(std::vector<spec_point>& real_data, std::vector<double>& real)
{
	readRealDataPoint(real_data); //read the empirical data from a txt file
	P.millionTime = get_min_time(); //this is the total time in millions of years of the empirical data
}

void readRealDataPoint(std::vector<spec_point>& R)
{
	R.clear();
	std::vector<int> lin;
	std::vector<double> time;
	std::ifstream in;

	if(file_exists("lampro_alt.txt"))
	{
		in.open("lampro_alt.txt");
	} else {
		in.open("lampro3.txt");
	}
	
	
	if(!in.is_open())
	{
		std::cout << "No input file of LTT data found!!!!\n";
		return;
	}

	while(!in.eof())
	{
		int l;
		double t;
		in >> l;
		in >> t;	
		time.push_back(t);
		lin.push_back(l);
	}
	//recalculate to maxT timesteps:
	for(std::size_t i = 0; i < lin.size(); ++i)
	{
		double tim = (P.maxT - P.maxT * 1.0 * time[i] / time[0]);
		//point temp(lin[i],tim);
		spec_point temp(lin[i],tim);
		R.push_back(temp);
	}
	//add final "closing" point
	R.push_back(spec_point(lin.back(),P.maxT));
	return;
}


bool withinLimits(std::vector<double>& F, int t, const std::vector<spec_point>& L, const std::vector<spec_point>& real, int model)
{
    if(t == 0) {
        return true;
    }

    if(L.empty()) {
        return false;
    }

    if(F[0] == 1e6) F[0] = calcNLTT(real,L);

    double euclidDist = 0.0;
    for(int i = 0; i < 4; ++i)  {
        double diff = P.emp_values[i] - F[i];
        diff *= diff;
        diff *= 1.0 / P.sd_values[model][i];
        euclidDist += diff;
    }
    euclidDist = sqrt(euclidDist);

    if(euclidDist < P.vec_epsilon[t]) {
        return true;
    }
    return false;
}





template<typename T> 
void removeDuplicates(std::vector<T>& vec)
{
	std::sort(vec.begin(), vec.end()); //first we have to sort, and we sort on the ID of species
	vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); //species with an identical ID are removed (as we assume that they have had the same history)
}

double get_min_time()
{
	std::ifstream input;
	if(P.generateData == 0) input.open("lampro3.txt");
	else input.open("lampro_alt.txt");
	double lin, T;
	input >> lin;
	input >> T;
	input.close();
	return -T;
}


int getLin2(const std::vector<spec_point>& v, double t)
{
   int min = 0;
   int max = (int)v.size() - 1;
   int med = (int)((max+min)*0.5); 
   while((max-min) > 1)
   {
	   if(v[med].time > t) max = med;
        else min = med;

        med = (int)((max+min)*0.5);
   }
   return (int)v[med].ID;
}



int getLin(const std::vector<spec_point>& v, double t)
{
	for(std::vector<spec_point>::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		if( (*it).time == t) return (int)(*it).ID;

		if(it != v.begin())
		{
			if((*it).time > t) return  (int)(*(it-1)).ID;
		}
	}

	return 0;
}

double calcDiff(const std::vector<spec_point>& R, const std::vector< spec_point >& sim, const double& time)
{
	//lookup timepoint in both vectors, and check the number of lineages at that time
	double lin_R = getLin2(R,time);
	double lin_S = getLin2(sim,time);

	double outcome = (lin_R - lin_S);
	outcome *= outcome; //squared distance, faster and easier than absolute distance

	return  outcome;
}

species::species(const species& parent_species, int& id_count, double b_time)
{
	ID = id_count;
	id_count++;
	death_time = -1;
	birth_time = b_time;
	parent = parent_species.ID;
	alloSpeciated = false;
	checked = false;
}

int findParent(const std::vector<species>& v, int parent_id)
{
	int count = 0;
	BOOST_FOREACH(species i, v) 
	{
		if(i.ID == parent_id) return count;

		count++;
	}
	return 0;
}

int countLin(const std::vector<species>& v, double time)
{
	int count = 0;
	BOOST_FOREACH(species i, v)
	{
		if( i.birth_time <= time)
		{
			if(i.death_time == -1 || i.death_time > time)
			{
				count++;
			}
		}
	}

	return count;
}

bool sortOnTime(const species& left, const species& right)
{
	return left.birth_time < right.birth_time;
}

void purgeOutput(std::vector<spec_point>& v)
{
	std::vector<spec_point> new_vec;

	int previous_ID = -1;

	BOOST_FOREACH(spec_point i, v)
	{
		if(i.ID != previous_ID)
		{
			new_vec.push_back(i);
			previous_ID = (int)i.ID;
		}
	}

	v = new_vec;
	return;
}

newick_node::newick_node()
{
	ID = -1;
	branch_length = -1;
	extant = false;
	parent = -1;
}

std::vector<int> findChildren(int ID, const std::vector<newick_node>& v)
{
	std::vector<int> output;
	int count = 0;
	BOOST_FOREACH(newick_node i, v)
	{
		if(i.parent == ID) output.push_back(count);
		count++;
	}
	return  output;
}

std::string newick_node::composeString(const std::vector<newick_node>& v) const
{
	std::string output;
	static std::string haakjeOpen("(");
	static std::string haakjeSluit(")");
	static std::string komma(",");
	static std::string dubbelp(":");

	std::string s_ID = boost::lexical_cast<std::string>(ID);
	std::string s_BL = boost::lexical_cast<std::string>(branch_length);

	if(extant)
	{
		output += s_ID;
		output += dubbelp;
		output += s_BL;
	}
	else
	{
		std::vector<int> children = findChildren(ID,v);

		if((int)children.size() == 0) 
		{
			std::cout << "error!\t";
			
		}
		output += haakjeOpen;
		for(std::size_t i = 0; i < children.size(); ++i)
		{
			output += v[children[i]].composeString(v);
			if(i != children.size()-1) output += komma;
		}
		output += haakjeSluit;
		
		if(parent != -1)
		{		
			output += dubbelp;
			output += s_BL;
		}
	}
	return output;
}

int getIndex(int ID, const std::vector<species>& v)
{
	for(std::size_t i = 0; i < v.size(); ++i)
	{
		if(v[i].ID == ID) return (int)i;
	}

	return -1;
}

int findYoungest(const std::vector<species>& v)
{
	double min = -1;
	int index = -1;
	int count = 0;

	BOOST_FOREACH(species i,v)
	{
		if(i.birth_time > min)
		{

			min = i.birth_time;
			index = count;
		}
		count++;
	}
	return index;
}

int find_parent( int ID, const std::vector<species>& v)
{
	for(std::size_t i = 0; i < v.size(); ++i)
	{
		if(v[i].ID == ID) return (int)i;
	}
	return -1;
}

int findOther( int youngest, const std::vector<species>& v)
{
	int parent = v[youngest].parent;

	for(int i = 0; i < (int)v.size(); ++i)
	{
		if(youngest != i)
		{
			if(v[i].parent == parent)
			{
				if(v[i].birth_time == v[youngest].birth_time)
				return i;
			}
		}
	}
	return -1;
}

newick_node::newick_node(const species& S)
{
	if(S.death_time == -1) 
	{
		extant = true;
		branch_length = P.maxT - S.birth_time;
	}
	else
	{
		extant = false;
		branch_length = S.death_time - S.birth_time;
	}

	parent = S.parent;
	ID = S.ID;
}

void updateReferences(int oldID, int newID, std::vector<newick_node>& v)
{
	for(std::size_t i = 0; i < v.size(); ++i)
	{
		if( v[i].parent == oldID) v[i].parent = newID;
	}
}

void updateReferences(int oldID, int newID, std::vector<species>& v, double time)
{
	BOOST_FOREACH(species i, v)
	{
		if( i.parent == oldID)
		{
			if(i.birth_time >= time) i.parent = newID;
		}
	}
}

std::vector<newick_node> generateNodeList(const std::vector<species>& v)
{
	std::vector<newick_node> node_list;

	std::vector<species> extant;
	std::vector<species> extinct; 

	int maxID = 0;

	BOOST_FOREACH(species i, v)
	{
		if(i.ID > maxID) maxID = i.ID;

		if(i.death_time == -1) extant.push_back(i);
		else extinct.push_back(i);
	}
	maxID++;

	//we have all extant species, that are connected through the extinct species to the ancestor (parent ID = -1);
	while(!extant.empty())
	{
		//we find the instance with the most recent branching time, e.g. the most recent branching moment
		int youngest = findYoungest(extant);
		int other = findOther(youngest,extant);
		
		if(extant[youngest].parent == -1) //we have reached the root!
		{
			newick_node root(extant[0]);
			node_list.push_back(root);
			break;
		}

		if(other != -1) //two species branched off from one other species
		{
			newick_node temp1(extant[youngest]);
			newick_node temp2(extant[other]); 

			node_list.push_back(temp1);
			node_list.push_back(temp2);

			//remove both offspring, and add the parent
			int parent = find_parent(extant[youngest].parent,extinct);

			if(parent != -1)
			{
				extant[youngest] = extinct[parent];
				extant[other] = extant.back();
				extant.pop_back();

				extinct[parent] = extinct.back();
				extinct.pop_back();
			}
			else
			{
				//the parent is already in the extant community
				std::vector<species> temp;
				for(int i = 0; i < (int)extant.size(); ++i)
				{
					if(i != youngest && i != other) temp.push_back(extant[i]);
				}
				extant = temp;
			}
		}
		else
		{
			//either the species branched off from an extant species, or not:
			int parent = find_parent(extant[youngest].parent,extant);

			if(parent == -1)
			{
				///////////////////////////////////
				/////////   1   /          1 /
				////////       /    -->     /
				////////   2  /          1 /
				///////////////////////////////////////
				parent = find_parent(extant[youngest].parent,extinct);
				int oldID = extant[youngest].ID;
				int newID = extinct[parent].ID;

				double time = extant[youngest].birth_time;

				extinct[parent].death_time = extant[youngest].death_time;
				extant[youngest] = extinct[parent];
				extinct[parent] = extinct.back();
				extinct.pop_back();

				updateReferences(oldID,newID,extinct,time);
				updateReferences(oldID,newID,extant,time);
				updateReferences(oldID,newID,node_list);
			}
			else //I don't think this can happen, but just in case
			{
				if(extant[parent].death_time == extant[youngest].birth_time)
				{
					///////////////////////////////////
					/////////   1   /           
					////////       /
					////////   2  /        
					///////////////////////////////////////
					newick_node temp(extant[youngest]);
					node_list.push_back(temp);

					extant[youngest] = extant.back();
					extant.pop_back();
				}
				else
				{
					///////////////////////////////////
					/////////   1   /           1  /
					////////       /\    -->      /\
					////////   1  /  \ 2       3 /  \ 2
					///////////////////////////////////////
					newick_node temp1(extant[youngest]); //number 2
					newick_node temp2(extant[parent]); //future number 3

					temp1.parent = extant[parent].ID;
					temp2.parent = extant[parent].ID; //they both have the same parents

					maxID++;
					temp2.ID = maxID;  //make it number 3
					temp2.branch_length = P.maxT - extant[youngest].birth_time; //adjust branch length (shorten it)
					if(extant[parent].death_time != -1) temp2.branch_length = extant[parent].death_time - extant[youngest].birth_time;


					extant[parent].death_time = extant[youngest].birth_time;

					updateReferences(extant[parent].ID, maxID, extant, extant[youngest].birth_time); //adjust all downstream references
					updateReferences(extant[parent].ID, maxID, extinct, extant[youngest].birth_time);
					updateReferences(extant[parent].ID, maxID, node_list);

					node_list.push_back(temp1); //add
					node_list.push_back(temp2);

					extant[youngest] = extant.back();
					extant.pop_back();
				}
			}
		}
	}


	return node_list;
}


std::string writeTREE2(const std::vector<species>& v)
{
	std::vector<species> extant;
	std::vector<species> extinct; 

	int maxID = 0;

	BOOST_FOREACH(species i, v)
	{
		if(i.ID > maxID) maxID = i.ID;

		if(i.death_time == -1) extant.push_back(i);
		else extinct.push_back(i);
	}
	maxID++;


	if(extant.size() ==1) //there is only one, or two species
	{
		std::string s_ID = boost::lexical_cast<std::string>(extant[0].ID);
		std::string s_BL = boost::lexical_cast<std::string>(P.maxT - extant[0].birth_time);
		std::string output = s_ID + ":" + s_BL;
		return output;
	}
	if(extant.size()==2)
	{
		std::string s_ID1 = boost::lexical_cast<std::string>(extant[0].ID);
		std::string s_BL1 = boost::lexical_cast<std::string>(P.maxT - extant[0].birth_time);
		std::string s_ID2= boost::lexical_cast<std::string>(extant[1].ID);
		std::string s_BL2 = boost::lexical_cast<std::string>(P.maxT - extant[1].birth_time);

		std::string output = "(" + s_ID1 + ":" + s_BL1 + "," + s_ID2 + ":" + s_BL2 + ")";
		return output;
	}

	std::vector< newick_node > node_list = generateNodeList(v);


	std::string core = node_list.back().composeString(node_list);
	return core;
}

std::vector<spec_point>  calculateLineages_noextinct(const std::vector<species>& allSp)
{
	std::vector<spec_point> output;
	std::vector<species> allSpecies = allSp;

	std::sort(allSpecies.begin(), allSpecies.end(), sortOnTime);

	BOOST_FOREACH(species i, allSpecies)
	{
		double time = i.birth_time;
		double previously_checked_time = -10;
		if(!output.empty()) previously_checked_time = output.back().time;

		if(time != previously_checked_time)	
		{
			int lins = countLin(allSpecies,time);
			output.push_back(spec_point(lins,time));
		}
	}

	purgeOutput(output);

	return output;
}

std::vector<int> find_indices(const std::vector<species>& v, int ID)
{
	std::vector<int> indices;
	int count = 0;

	BOOST_FOREACH(species i, v)
	{
		if(i.parent == ID) indices.push_back(count);
		
		count++;
	}

	return indices;
}

bool species::check_has_viable_offspring(std::vector<species>& v)
{
	if(checked == true) return extant_offspring;

	std::vector<int> offspring = find_indices(v,ID); //find the positions of the offspring;
	extant_offspring = false;

	for(std::size_t i = 0; i < offspring.size(); ++i) //ofspring is of size 2 (or 1), so using iterators is useless here
	{
		int index = offspring[i];
		
		if(v[index].death_time == -1) extant_offspring  = true;  //the offspring is extant
		else   //the offspring has died, but might have given birth to other species that have survived
		{
			if(v[index].checked == true) extant_offspring  = v[index].extant_offspring;
			else extant_offspring = v[index].check_has_viable_offspring(v);
		}

		if(extant_offspring == true) break;
	}

	checked = true;

	return extant_offspring;
}

std::vector<spec_point>  calculateLineages_withextinct(std::vector<species>& allSpecies)
{
	std::vector<spec_point> output;
	std::vector<species> lineages;

	BOOST_FOREACH(species i, allSpecies)
	{
	    if(i.death_time == -1) lineages.push_back(i); //the species is an extant species and is added to the list to create the ltt plot
		else
		{
			i.check_has_viable_offspring(allSpecies);
			if(i.extant_offspring == true) lineages.push_back(i);			
		}
	}
	allSpecies = lineages;
	std::vector<species> filler;
	output = calculateLineages_noextinct(allSpecies);
	output.push_back(spec_point(output.back().ID,P.maxT));

	return output;
}

int countLineages(const std::vector<spec_point>& v, double time)
{
	for(std::size_t i = 0; i < v.size(); ++i)
	{
		if(v[i].time == time) return (int)v[i].ID;

		if(i > 0)
		{
			if(v[i].time > time && v[i-1].time < time) return (int)v[i-1].ID;
		}
	}
	return (int)v.back().ID;
}



std::vector<spec_point> sumBranches(const std::vector<spec_point>& b1, const std::vector<spec_point>& b2)
{
	std::vector<spec_point> output;
	std::vector<double> times;
	BOOST_FOREACH(spec_point i, b1) times.push_back(i.time);
	BOOST_FOREACH(spec_point i, b2) times.push_back(i.time);

	removeDuplicates(times); //to remove duplicates, the list also gets sorted, how convenient!

	int L1, L2;

	BOOST_FOREACH(double i, times)
	{
		L1 = countLineages(b1,i);
		L2 = countLineages(b2,i);

		spec_point add(L1+L2,i);
		output.push_back(add);
	}


	output.push_back(spec_point(output.back().ID,P.maxT));
	removeDuplicates(output);
	return output;
}

void writeLineage(const std::vector<spec_point>& L, int iter, int part_number, double fit)
{
	std::string s_iter = boost::lexical_cast<std::string>(iter);
	std::string s_part_number = boost::lexical_cast<std::string>(part_number);

	std::string f_name = "LTT_" + s_iter + "_" + s_part_number + ".txt";	
	std::ofstream out_file(f_name.c_str());

	BOOST_FOREACH(spec_point i, L)
	{
		out_file << i.time << "\t" << i.ID << "\n";
	}

	out_file.close();


	std::string f_name2 = "values_" + s_iter + "_" + s_part_number + ".txt";
	std::ofstream values_file(f_name2.c_str());

	values_file << fit << "\t"; 
	values_file << P.extinctionRate << "\t"; 
	values_file << P.SympatricRate_low << "\t";
	values_file << P.SympatricRate_high << "\t";
	values_file << P.AllopatricRate << "\t"; 
    values_file << P.period_waterlevel << "\t";

	values_file << L.back().ID;
	values_file.close();


	return;
}

void write_newick_file(const std::vector<species>& s1, const std::vector<species>& s2, const std::vector<spec_point>& b1, const std::vector<spec_point>& b2, int iter, int part_number, int jiggle)
{
	double b_1,b_2;
	if(b1.size() == 1)	b_1 = P.maxT - b1[0].time;
	else b_1 = b1[1].time;

	if(b2.size() ==1) b_2 = P.maxT - b2[0].time;
	else b_2 = b2[1].time;


	std::string left = writeTREE2(s1);
	std::string right = writeTREE2(s2);
	std::string BL_left = boost::lexical_cast<std::string>(b_1);
	std::string BL_right = boost::lexical_cast<std::string>(b_2);

	std::string s_iter = boost::lexical_cast<std::string>(iter);
	std::string s_part_number = boost::lexical_cast<std::string>(part_number);
	std::string s_jiggle = boost::lexical_cast<std::string>(jiggle);
	
	std::string f_name = "NEWICK_" + s_iter + "_" + s_part_number + "_" + s_jiggle +  ".tre";
	std::ofstream out_file(f_name.c_str());

	std::string output = "(";
    output += left;
	output += ":";
	output += BL_left;
	output += ",";
	output += right;
	output += ":";
	output += BL_right;
	output += ");";
	
	out_file << output;
	out_file.close();

	return;
}

void append_newick_file(const std::vector<species>& s1, const std::vector<species>& s2, const std::vector<spec_point>& b1, const std::vector<spec_point>& b2, int iter)
{
	double b_1,b_2;
	if(b1.size() == 1)	b_1 = P.maxT - b1[0].time;
	else b_1 = b1[1].time;

	if(b2.size() ==1) b_2 = P.maxT - b2[0].time;
	else b_2 = b2[1].time;


	std::string left = writeTREE2(s1);
	std::string right = writeTREE2(s2);
	std::string BL_left = boost::lexical_cast<std::string>(b_1);
	std::string BL_right = boost::lexical_cast<std::string>(b_2);

	std::string s_iter = boost::lexical_cast<std::string>(iter);

	std::string f_name = "NEWICK_" + s_iter + ".tre";
	std::ofstream out_file(f_name.c_str(), std::ios::app);

	std::string output = "(";
    output += left;
	output += ":";
	output += BL_left;
	output += ",";
	output += right;
	output += ":";
	output += BL_right;
	output += ");";
	
	out_file << output << "\n";
	out_file.close();

	return;
}

std::vector<double> generateWaterLevelChanges(int model)
{
	std::vector<double> output;
	
	if(model == 0)
	{
		//constant high water
		output.push_back(P.maxT * 2);
	}
	
	if(model == 1)
	{
		//0,35,40,169,193,262,295,363,293,550,1100
		output.push_back(P.maxT - 1.1);			//0
		output.push_back(P.maxT - 0.55);		//1
		output.push_back(P.maxT - 0.393);		//0
		output.push_back(P.maxT - 0.363);		//1
		output.push_back(P.maxT - 0.295);		//0
		output.push_back(P.maxT - 0.262);		//1
		output.push_back(P.maxT - 0.193);		//0
		output.push_back(P.maxT - 0.169);		//1
		output.push_back(P.maxT - 0.04);		//0
		output.push_back(P.maxT - 0.035);		//1
		output.push_back(P.maxT);				//1
	}
	
	if(model == 2)
	{ //literature values for the past 1 million years, non-literature values before that
		
		int waterLevel = 1;
		double time = 0;
		double tempTime = time; // + Expon(10); //HARDCODED RATE
		while(tempTime < (P.maxT - 1.1))
		{
			tempTime += Expon(10);
			if(tempTime > (P.maxT - 1.1)) break;
			time = tempTime;
			output.push_back(time);
			waterLevel = 1 - waterLevel;
		}
		if(waterLevel == 0) output.pop_back(); //the waterlevel has to be high.
		
		output.push_back(P.maxT - 1.1);			//0
		output.push_back(P.maxT - 0.55);		//1
		output.push_back(P.maxT - 0.393);		//0
		output.push_back(P.maxT - 0.363);		//1
		output.push_back(P.maxT - 0.295);		//0
		output.push_back(P.maxT - 0.262);		//1
		output.push_back(P.maxT - 0.193);		//0
		output.push_back(P.maxT - 0.169);		//1
		output.push_back(P.maxT - 0.04);		//0
		output.push_back(P.maxT - 0.035);		//1
	}

	return output;
}

int find_Allo(int i, const std::vector<species>& v)
{
	int ID = v[i].ID;
	for(int j = 0; j < (int)v.size(); ++j)
	{
		if(v[j].ID == ID && i != j) return j;
	}
	return -1;
}

std::vector<int> findOffspring(int ID, const std::vector<species>& v)
{
	std::vector<int> output;
	int counter = 0;
	BOOST_FOREACH(species i, v)
	{
		if(i.parent == ID) output.push_back(counter);
		counter++;
	}

	return output;
}






void makeFiles(const std::vector<spec_point>& L)
{

	std::ofstream out("config_generated.ini");
	out << "numberParticles = " << P.numberParticles << "\n";
	out << "generateData " << P.generateData << "\n";
	out << "maxT = " << P.maxT << "\n";
	out << "waterlevel_period = " << log10(P.period_waterlevel) << "\n";
	out << "sympatricRate_low = " << log10(P.SympatricRate_low) << "\n";
	out << "sympatricRate_high = " << log10(P.SympatricRate_high) << "\n";
	out << "allopatricRate = " << log10(P.AllopatricRate) << "\n";
	out << "extinctionRate = " << log10(P.extinctionRate) << "\n";
	out << "Jiggle = " << log10(P.Jiggle);

	out.close();

	std::ofstream ltt_out("lampro_alt.txt");
	for(std::size_t i = 0; i < L.size(); ++i)
	{
		ltt_out << L[i].ID << "\t" << P.maxT - L[i].time << "\n";
	}
	ltt_out.close();
}


std::vector<int> findOffspring(int ID, const std::vector<newick_node>& v)
{
	std::vector<int> output;

	int counter = 0;
	BOOST_FOREACH(newick_node i, v)
	{
		if(i.parent == ID) output.push_back(counter);
		counter++;
	}

	return output;
}



double calcNLTT(const std::vector<spec_point>& R, const std::vector<spec_point>& L)
{
	int maxLinR = (int)R.back().ID;
	int maxLinL = (int)L.back().ID;

	
	std::vector<double> times;
	for(int i = 0; i < (int)R.size(); ++i)
	{
		times.push_back(R[i].time);
	}
	for(int j = 0; j < (int)L.size(); ++j)
	{
		times.push_back(L[j].time);
	}
	removeDuplicates(times);
	
	double sum = 0;
	
	for(int t = 1; t < (int)times.size(); ++t)
	{
		double dt = 1.0 * times[t] / P.maxT - 1.0 * times[t-1] / P.maxT;
		int numLinR = countLineages(R, times[t-1]);
		int numLinL = countLineages(L, times[t-1]);
		double absDiff = 1.0 * numLinR / maxLinR - 1.0 * numLinL / maxLinL;
		if(absDiff < 0) absDiff *= -1.0;
		sum += absDiff * dt;
	}
	return sum;
}

double calcGamma(const std::vector<spec_point>& L)
{
	if(L.size() == 2) {
		return 1e6;
	}
	
	std::vector<double> b_times;
	b_times.push_back(L[0].time);
	for(int i = 1; i < (int)L.size(); ++i)
	{
		int reps = L[i].ID - L[i-1].ID;
		for(int r = 0; r < reps; ++r)
		{
			b_times.push_back(L[i].time);
		}
	}
	int N = 1+(int)b_times.size();
	std::vector<double> bt;
	for(int i = 0; i < (int)b_times.size(); ++i)
	{
		double diff = P.maxT - b_times[i];
		if(diff != 0) {
			bt.push_back(P.maxT - b_times[i]);
		}
	}
	std::sort(bt.begin(),bt.end());
	
	
	std::vector<double> g;
	g.push_back(bt[0]);
	for(int i = 1; i < (int)bt.size(); ++i)
	{
		g.push_back(bt[i]-bt[i-1]);
	}
	
	std::reverse(g.begin(),g.end());
	
	//ST <- sum((2:N) * g)
	std::vector<int> values;
	for(int i = 2; i <= N; ++i) values.push_back(i);
	double ST = 0;
	for(int i = 0; i < (int)g.size(); ++i)
	{
		double add = g[i] * values[i];
		ST += add;
	}
	//stat <- sum(cumsum((2:(N - 1)) * g[-(N - 1)]))/(N - 2)
	
	values.clear();
	
	for(int i = 2; i <= (N-1); ++i) values.push_back(i);
	std::vector<double> toCumSum;
	
	for(int i = 0; i < (int)g.size(); ++i)
	{
		if(i != (N))
		{
			toCumSum.push_back(g[i] * values[i]);
		}
	}
	std::vector<double> cumSum(1,0);
	for(int i = 0; i < (int)toCumSum.size(); ++i)
	{
		double add = toCumSum[i] + cumSum[i];
		cumSum.push_back(add);
	}
	cumSum.pop_back();
	double stat = 0.0;
	for(int i = 0; i < (int)cumSum.size(); ++i)
	{
		stat += cumSum[i];
	}
	stat *= 1.0 / (N-2);
	
	double s = ST * sqrt(1.0 / (12 * (N-2)));
	double m = 1.0 * ST / 2;
	double gamma = (stat - m) / s;
	return gamma;
}

double calcBr(const std::vector<species>& S1, const std::vector<species>& S2)
{
	std::vector< newick_node > node_list1 = generateNodeList(S1);
	std::vector< newick_node > node_list2 = generateNodeList(S2);

	double sum = 0.0;
	for(std::size_t i = 0; i < node_list1.size(); ++i) sum += node_list1[i].branch_length;
	for(std::size_t i = 0; i < node_list2.size(); ++i) sum += node_list2[i].branch_length;

	sum *= 1.0 / (node_list1.size() + node_list2.size());
	return(sum);
}


bool file_exists(const std::string& name) 
{
    std::ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}


void macstart(const char * argv[])  {
    std::cout << "\n\n\n";
#ifdef __APPLE__
    {
        char *dirsep = strrchr(argv[0], '/');
        if (dirsep != NULL) *dirsep = 0;
        int changeDir = chdir(argv[0]);
        std::cout << "Changing Dir: " << changeDir << "\n";
        std::string cwd = getcwd(NULL, 0);
        std::cout << cwd << "\n";
        std::cout << "Starting simulation\n";
    }
#endif
}
