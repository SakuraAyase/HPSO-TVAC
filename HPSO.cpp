#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi 3.1415926535
#define E  2.71828182845904523536

using namespace std;

double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}
class Particle
{
public:
	double fitness;
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;
	double pBestFitness;

	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;
		this->pBestFitness = best_fitness;
	}
};

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c1s, double c2s, 
		double c1f, double c2f,double dt, double percent)
	{
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c1i = c1s;
		this->c2i = c2s;
		this->c1f = c1f;
		this->c2f = c2f;

		this->dt = dt;
		this->percent = percent;
		particles.resize(m);
		gBest.resize(dim);
	}


	double fitnessFunction(Particle particle)
	{
		double result = 0.0;
		for (int i = 0; i < dim; i++)
		{
			double x = particle.position[i];
			result += pow(x, 2);
		}

		//result = -20 * exp(-0.2*sqrt(temp1 / dim)) + 20 + E - exp(temp2 / dim);

		return result;
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(this->min, this->max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].position[j] = particles[i].position[j];

		}
		particles[i].fitness = fitnessFunction(particles[i]);
		particles[i].pBestFitness = fitnessFunction(particles[i]);
	}

	void initialAllParticles()
	{
		initialParticles(0);
		gBest_fitness = particles[0].pBestFitness;
		for (int i = 0; i < dim; i++)
			gBest[i] = particles[0].pBest[i];

		for (int i = 1; i < m; i++)
		{
			initialParticles(i);
			if (particles[i].pBestFitness < gBest_fitness)
			{
				gBest_fitness = particles[i].pBestFitness;
				for (int j = 0; j < dim; j++)
				{
					gBest[j] = particles[i].pBest[j];
				}
			}
		}
	}

	void inertiaWeight()
	{
		//w = randDouble(0.4, 0.6);
		double t = T / ((double)Tmax);
		w = 0.9;
	}

	void coefficient()
	{
		c1 = (c1f - c1i)*T / Tmax + c1i;
		c2 = (c2f - c2i)*T / Tmax + c2i;
	}
	void updateParticle(int i)
	{
		for (int j = 0; j < dim; j++)
		{
			double last_position = particles[i].position[j];
			double range = percent * (max - min);

			particles[i].velocity[j] = c1 * randDouble(0, 1) * (particles[i].pBest[j] - particles[i].position[j])
				+ c2 * randDouble(0, 1) * (gBest[j] - particles[i].position[j]);
			if (particles[i].velocity[j] == 0)
			{
				if (randDouble(0, 1) < 0.5)
					particles[i].velocity[j] = randDouble(0, 1)*randDouble(0, range);
				else
					particles[i].velocity[j] = -randDouble(0, 1)*randDouble(0, range);
			}

			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = range;

			if (particles[i].velocity[j] < -range)
				particles[i].velocity[j] = -range;


			particles[i].position[j] += dt * particles[i].velocity[j];

			

			if (particles[i].position[j] > max)
				particles[i].position[j] = max;
			if (particles[i].position[j] < min)
				particles[i].position[j] = min;




		}
		particles[i].fitness = fitnessFunction(particles[i]);
		if (particles[i].fitness < particles[i].pBestFitness)
		{
			particles[i].pBestFitness = particles[i].fitness;
			for (int j = 0; j < dim; j++)
			{
				particles[i].pBest[j] = particles[i].position[j];
			}
		}

	}


	void updateAllParticles()
	{
		coefficient();
		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
			if (particles[i].pBestFitness < gBest_fitness)
			{
				gBest_fitness = particles[i].pBestFitness;
				for (int j = 0; j < dim; j++)
				{
					gBest[j] = particles[i].pBest[j];
				}
			}
		}
		T++;
	}

	double getFitness()
	{
		return gBest_fitness;
	}
private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	double w;
	double max;
	double min;
	double c1;
	double c2;
	double c1i;
	double c2i;
	double c1f;
	double c2f;


	double dt;//时间步长
	double percent;

	vector<double>gBest;
	double gBest_fitness;

	vector<Particle> particles;


};

void run(vector<double>& result1)
{
	int dim = 30;
	int m = 20;
	int Tmax = 2000;
	double max = 100;
	double min = -100;
	double c1s = 2.5;
	double c2s = 0.5;
	double c1f = 0.5;
	double c2f = 2.5;

	double dt = 1.0;
	double percent = 0.2;

	PSO pso = PSO(dim, m, Tmax, max, min, c1s, c2s, c1f, c2f, dt, percent);
	pso.initialAllParticles();

	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		
		cout << ":";
		//fitness.push_back(pso.getFitness());
		fitness.push_back(pso.getFitness());
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;
		pso.updateAllParticles();
	}



	result1 = fitness;
}


int main()
{


	int times = 5;
	int interval = 10;
	vector<double> result1;

	run(result1);

	for (int i = 1; i < times; i++)
	{
		vector<double> result1_temp;
		run(result1_temp);
		for (int j = 0; j < result1_temp.size(); j++)
		{
			result1[j] += result1_temp[j];
		}
	}
	for (int j = 0; j < result1.size(); j++)
	{
		result1[j] /= times;
	}

	for (int j = 0; j < result1.size(); j++)
	{
		if (j%interval == 0)
			cout << result1[j] << " ";
	}

	system("pause");
}