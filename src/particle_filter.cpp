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
    // initialize weight list
    weights.push_back(1.0);
  }

  // mark particle filter as initialized
  is_initialized = true;
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

    // extract particle position
    double xp = particles[i].x;
    double yp = particles[i].y;
    double thetap = particles[i].theta;

    // set weight to 1:
    double weight = 1.0;

    // vectors to store the association details for each particle
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;

    // landmarks in sensor range
    vector<Map::single_landmark_s> landmarks_ir;

    for(vector<Map::single_landmark_s>::size_type k = 0; k != map_landmarks.landmark_list.size(); k++){
      Map::single_landmark_s landmark = map_landmarks.landmark_list[k];
      if(dist(xp, yp, landmark.x_f, landmark.y_f) < sensor_range){
        landmarks_ir.push_back(landmark);
      }
    }

    for(vector<LandmarkObs>::size_type j = 0; j != observations.size(); j++){

      // extract observation components (in the particle's FoR)
      double xc = observations[j].x;
      double yc = observations[j].y;

      // transform observation to map global FoR
      double xm = xp + xc * cos(thetap) - yc * sin(thetap);
      double ym = yp + xc * sin(thetap) + yc * cos(thetap);

      // initialize minimum distance to infinity
      double min_dist = numeric_limits<double>::infinity();
      int nearest_landmark_id;
      double nearest_landmark_x, nearest_landmark_y;
      for(vector<Map::single_landmark_s>::size_type k = 0; k != landmarks_ir.size(); k++){

        // extract landmark coordinates (in map FoR)
        double xl = landmarks_ir[k].x_f;
        double yl = landmarks_ir[k].y_f;

        // calculate between transformed observation and landmark (in map FoR)
        double dist_obs_landmark = dist(xm, ym, xl, yl);

        // identify nearest landmark to the observation "j"
        if(dist_obs_landmark < min_dist){
          min_dist = dist_obs_landmark;
          nearest_landmark_id = landmarks_ir[k].id_i;
          nearest_landmark_x = xl;
          nearest_landmark_y = yl;
        }

      }

      // store associations: nearest landmark id, and respective (x,y) coordinates in global FoR
      associations.push_back(nearest_landmark_id);
      sense_x.push_back(nearest_landmark_x);
      sense_y.push_back(nearest_landmark_y);

      // mean of bivariate normal distribution: nearest landmark position
      // evaluation of distribution @ (x,y): particle's observation in global FoR
      double mu [2] = { nearest_landmark_x, nearest_landmark_y };
      weight *= bivariateNormalDist(xm, ym, mu, std_landmark);
    }

    // set the associations for particle "i"
    particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
    // set weight for particle "i"
    particles[i].weight = weight;
    // add particle weight to weight list
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
