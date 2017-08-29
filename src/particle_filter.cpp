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

  // create normal distributions for x, y and theta
  normal_distribution<double> dist_x(x, std_pos[0]);
  normal_distribution<double> dist_y(y, std_pos[1]);
  normal_distribution<double> dist_theta(theta, std_pos[2]);


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

  double min_dist_meas;
  double dist_meas;
  int index_min_dist_meas;

  // iterate over actual predicted landmarks
  for(vector<LandmarkObs>::size_type i = 0; i != predicted.size(); i++){
    // initialize minimum distance to landmark
    min_dist_meas = dist(predicted[i].x, predicted[i].y, observations[0].x, observations[0].y);
    index_min_dist_meas = 0;
    // iterate over lidar measurements looking for nearest measurement to
    // predicted landmark i
    for(vector<LandmarkObs>::size_type j = 1; j != observations.size(); j++){
      dist_meas = dist(predicted[i]. x,predicted[i].y, observations[j].x, observations[j].y);
      if(dist_meas < min_dist_meas){
        min_dist_meas = dist_meas;
        index_min_dist_meas = j;
      }
    }
    // association via id
    observations[index_min_dist_meas].id = predicted[i].id;
  }
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

  // store particle weights into the `weight` vector
  for(vector<Particle>::size_type i = 0; i != particles.size(); i++){
    weights.push_back(particles[i].weight);
  }

  // discrete distribution based on particle weights
  discrete_distribution<> dd(weights.begin(), weights.end());

  //random number engine
  default_random_engine gen;

  // resample vector of particles
  vector<Particle> temp(particles); // copy
  for(vector<Particle>::size_type i = 0; i != particles.size(); i++){
    int sample_index = dd(gen);
    particles[i] = temp[sample_index];
  }

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
