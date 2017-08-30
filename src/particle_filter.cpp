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

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {

  weights.clear();

  for(vector<Particle>::size_type i = 0; i != particles.size(); i++){

    // set weight to 1:
    double weight = 1.0;

    // extract particle position
    double xp = particles[i].x;
    double yp = particles[i].y;
    double thetap = particles[i].theta;

    for(vector<LandmarkObs>::size_type j = 0; j != observations.size(); j++){

      // extract observation components (in the particle's FoR)
      double xc = observations[j].x;
      double yc = observations[j].y;

      // transform observation to map global FoR
      double xm = xp + xc * cos(thetap) - yc * sin(thetap);
      double ym = yp + xc * sin(thetap) + yc * cos(thetap);

      double min_dist_ij;
      double mu [2];
      for(vector<Map::single_landmark_s>::size_type k = 0; k != map_landmarks.landmark_list.size(); k++){

        // extract landmark coordinates (in map FoR)
        double xl = map_landmarks.landmark_list[k].x_f;
        double yl = map_landmarks.landmark_list[k].y_f;

        // calculate between transformed observation and landmark (in map FoR)
        double dist_ijk = dist(xm, ym, xl, yl);

        // identify nearest landmark to the observation "j"
        if(k == 0 || dist_ijk < min_dist_ij){
          min_dist_ij = dist_ijk;
          mu[0] = xl;
          mu[1] = yl;
        }
      }
      weight *= bivariateNormalDist(xm, ym, mu, std_landmark);
    }
    particles[i].weight = weight;
    weights.push_back(weight);
  }
}

double ParticleFilter::bivariateNormalDist(double x, double y, double mu[], double sigma[]){
  double ex = pow(x - mu[0], 2) / (2 * pow(sigma[0], 2));
  double ey = pow(y - mu[1], 2) / (2 * pow(sigma[1], 2));
  double c = 2 * M_PI * sigma[0] * sigma[1];
  return 1 / c * exp(-(ex + ey));
}

void ParticleFilter::resample() {

  // discrete distribution based on particle weights
  discrete_distribution<> dd(weights.begin(), weights.end());

  //random number engine
  default_random_engine gen;

  // resample vector of particles
  vector<Particle> sampled_particles;
  for(vector<Particle>::size_type i = 0; i != particles.size(); i++){
    int sample_index = dd(gen);
    sampled_particles.push_back(particles[sample_index]);
  }
  // update particles
  particles = sampled_particles;
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
