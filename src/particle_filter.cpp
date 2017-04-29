/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

#ifndef NUM_PARTICLES
  #define NUM_PARTICLES (500)
#endif /* NUM_PARTICLES */



/**
 * Efficient way to compare distances: https://en.wikibooks.org/wiki/Algorithms/Distance_approximations
 * */


void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = NUM_PARTICLES;
  std::default_random_engine engine;
  std::normal_distribution<double_t> x_rand(x,std[0]);
  std::normal_distribution<double_t> y_rand(y,std[1]);
  std::normal_distribution<double_t> th_rand(theta,std[2]);
  Particle p;
  for (int i = 0; i < num_particles; ++i)
  {
    particles.push_back(Particle{0, x_rand(engine), y_rand(engine), th_rand(engine), 1.0f});
  }

	is_initialized = 1;
}


/* Define what is zero in this implementation  */
#ifndef EPSILON
#define EPSILON ((double_t) 1e-3)
#endif /* EPSILON */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  /* Add noise to velocity and yaw rate */
  std::default_random_engine eng;


  /* Using the bycicle motion model update the particle position or each particle */
  /* Bycicle motion model for Yaw rate != 0:
   * x_f = x_0 + (v/thetha_dot)*(sin(theta_0 + theta_dot *dt) - sin(theta_0))
   * y_f = y_0 + (v/thetha_dot)*(cos(theta_0) - sin(theta_0 + theta_dot *dt))
   * theta_f = theta_0 + theta_dot *dt
   *
   * For yaw rate == 0
   * x_f = x_0 + v*dt * cos(theta_o)
   * x_f = y_0 + v*dt * sin(theta_o)
   * theta_f = 0;
   *
   * */

  for(std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); ++p)
  {
    if(fabs(yaw_rate) < EPSILON)
    {
      p->x = p->x + velocity * delta_t * std::cos(p->theta);
      p->y = p->y + velocity * delta_t * std::sin(p->theta);
      p->theta = 0;
    }
    else
    { /* Yaw rate is different from zero */
      p->x = p->x + ((velocity/yaw_rate) * (std::sin(p->theta + yaw_rate * delta_t) - std::sin(p->theta)));
      p->y = p->y + ((velocity/yaw_rate) * (std::cos(p->theta) - std::cos(p->theta + yaw_rate * delta_t)));
      p->theta += yaw_rate * delta_t;
    }
    /*Adding noise*/
    std::normal_distribution<double_t> noise_x(p->x, std_pos[0]);
    std::normal_distribution<double_t> noise_y(p->y, std_pos[1]);
    std::normal_distribution<double_t> noise_th(p->theta, std_pos[2]);
    /* Additive Gaussian noise */
    p->x += noise_x(eng);
    p->y += noise_y(eng);
    p->theta += noise_th(eng);
  }
}
/***
 * Compare each observation with each predicted(trasnformed) landmark and assign the id
 * of the landmark to the predicted observation
 * @param predicted
 * @param observations
 * reference:  http://www.geeksforgeeks.org/closest-pair-of-points/
 * http://www.geeksforgeeks.org/closest-pair-of-points-onlogn-implementation/
 * This algorithm can be improved, but it implies sorting the points, not needed now
 * we'll keep this O(n^2) implementaion
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations)
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  for(std::vector<LandmarkObs>::iterator obs_it = observations.begin(); obs_it != observations.end(); ++obs_it)
  {
    /* Initialize with a large number */
    double_t min_dist = std::numeric_limits<double_t>::max();
    for(std::vector<LandmarkObs>::iterator l_it = predicted.begin(); l_it != predicted.end(); ++l_it)
    {
      /*get the distance between points */
      // is it a better way instead of using sqrt()?
      // https://en.wikibooks.org/wiki/Algorithms/Distance_approximations

      if(dist(obs_it->x,obs_it->y,l_it->x,l_it->y) < min_dist)
      {
        obs_it->id = l_it->id;
      }
    }
  }
}
/***
 * paso 1: trasnformar las mediciones del auto a coordenadas globales tomando como x0,y0 las coordenadas de la n part
 * paso 2: Por cada medicion de cada particula ya transformada asociar un landmark a ella si es que esta en el rango
 *         del sensor ie < sensor_range
 * paso 3: Encontrar el peso de la particula usando la distribucion Gaussiana multivariada
 * @param sensor_range
 * @param std_landmark
 * @param observations
 * @param map_landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks)
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  /* Transform the  observations to global coordinates */
  /* $\displaystyle \begin{pmatrix}x \cos\theta - y \sin\theta + x_t x \sin\theta + y \cos\theta + y_t \end{pmatrix} .$ */
/*
  std::vector<LandmarkObs>
  for (int i = 0; i < observations.size(); ++i) {

  }*/
}

/***
 *
 */

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
