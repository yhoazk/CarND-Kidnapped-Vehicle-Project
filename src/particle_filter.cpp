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
    p.weight = 1.0;
    p.x = x_rand(engine);
    p.y = y_rand(engine);
    p.theta = th_rand(engine);
    p.id = i;
    particles.push_back(p);
  }

	is_initialized = 1;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  /* Add noise to velocity and yaw rate */
  std::default_random_engine eng;
  std::normal_distribution<double_t> yaw_noise();
  std::normal_distribution<double_t> velocity_noise();

  /* Using the bycicle motion model update the particle position or each particle */
  /* Bycicle motion model:
   * x_f = x_0 + (v/thetha_dot)*(sin(theta_0 + theta_dot *dt) - sin(theta_0))
   * y_f = y_0 + (v/thetha_dot)*(cos(theta_0) - sin(theta_0 + theta_dot *dt))
   * theta_f = theta_0 + theta_dot *dt
   * */
  for(std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); ++p)
  {
    p->x = p->x + ((velocity/yaw_rate) * (std::sin(p->theta + yaw_rate * delta_t) - std::sin(p->theta)));
    p->y = p->y + ((velocity/yaw_rate) * (std::cos(p->theta) - std::cos(p->theta + yaw_rate * delta_t)));
    p->theta += yaw_rate * delta_t;
  }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
}

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
