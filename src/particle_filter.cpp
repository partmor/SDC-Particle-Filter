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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std_pos[]) {

  num_particles = 100;

  //random number engine
  default_random_engine gen;

  // standard deviations for x, y, and theta
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  // create normal distributions for x, y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);


  for (int i = 0; i < num_particles; ++i) {
    Particle particle;
    particle.id = i + 1;

    // take sample from distributions
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    // add created particle to the list of particles
    particles.push_back(particle);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

  //random number engine
  default_random_engine gen;

  // iterate over the vector of particles
  for(vector<Particle>::size_type i = 0; i != particles.size(); i++) {

    // retrieve particle current state
    Particle particle = particles[i];
    double x0 = particle.x;
    double y0 = particle.y;
    double theta0 = particle.theta;

    // predict new state of the particle: apply motion model
    double xp, yp, thetap;
    // avoid division by zero
    if (fabs(yaw_rate) > 0.001) {
      xp = x0 + velocity/yaw_rate * ( sin (theta0 + yaw_rate*delta_t) - sin(theta0));
      yp = y0 + velocity/yaw_rate * ( cos(theta0) - cos(theta0 + yaw_rate*delta_t) );
      thetap = theta0 + yaw_rate * delta_t;
    }
    else {
      xp = x0 + velocity * delta_t * cos(theta0);
      yp = y0 + velocity * delta_t * sin(theta0);
      thetap = theta0;
    }

    // NOTE: the method is receiving std of x, y, theta given for GPS.
    // std should be related to velocity and yaw_rate, and be applied on these,
    // not directly on position: issue on Udacity's side.

    // create normal distributions for x, y and theta
    normal_distribution<double> dist_x(xp, std_pos[0]);
    normal_distribution<double> dist_y(yp, std_pos[1]);
    normal_distribution<double> dist_theta(thetap, std_pos[2]);

    // add random noise and update particle
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particles[i] = particle;
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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
