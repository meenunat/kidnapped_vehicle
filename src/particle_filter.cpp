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
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles= 10;
	std::default_random_engine gen;
	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);
	for(int i=0; i<num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen) ;
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(particle.weight);
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	default_random_engine gen;
	for (int i = 0; i < particles.size(); i++)
	{
		double new_x;
		double new_y;
		double new_theta;
		if (fabs(yaw_rate)<1e-5){
			double v_dt = velocity*delta_t;
			new_x = particles[i].x + (v_dt*cos(particles[i].theta));
			new_x = particles[i].y + (v_dt*sin(particles[i].theta)); 
			new_theta = particles[i].theta;
		}
		else{
			double v_dt = velocity/yaw_rate;
			new_x = particles[i].x + (v_dt)*(sin(particles[i].theta+(yaw_rate*delta_t)) - sin(particles[i].theta));
			new_y = particles[i].y + (v_dt)*(cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate*delta_t)));
			new_theta = particles[i].theta + (yaw_rate*delta_t);
		}
		normal_distribution<double> N_x(new_x, std_pos[0]);
		normal_distribution<double> N_y(new_y, std_pos[1]);
		normal_distribution<double> N_theta(new_theta, std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);

	} 
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (auto& obs :observations) {
	    double min_dist = numeric_limits<float>::max();
    
	    for (const auto& pred : predicted) {
	      	
      		//get distance between current/predicted landmarks
      		 double curr_dist = dist(obs.x, obs.y, pred.x, pred.y);

	        // find the predicted landmark nearest the current observed landmark
		if (curr_dist < min_dist) {
		        min_dist = curr_dist;
       			obs.id = pred.id;
      		}
    	    }

        }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> observations, const Map map_landmarks) {
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

	
	const double landmark_x = std_landmark[0];
	const double landmark_y = std_landmark[1];

	weights.clear();

	for (auto& p:particles){

		p.weight = 1.0;
		map<int, LandmarkObs> landmarksInRange;
		
		for (const auto& landmark : map_landmarks.landmark_list){
			if ((landmark.x_f>=p.x-sensor_range && landmark.x_f<=p.x+sensor_range)&&(landmark.y_f>=p.y-sensor_range && landmark.y_f<=p.y+sensor_range))		
			{
				
				LandmarkObs pred;
				pred.id = landmark.id_i;
				pred.x = landmark.x_f;
				pred.y = landmark.y_f;
				landmarksInRange[pred.id] = pred;

			}
		}

	vector<LandmarkObs> trans_observations;
	for (const auto& obs: observations)
	{
		LandmarkObs trans_obs;
		trans_obs.x = p.x +(obs.x*cos(p.theta)) - (obs.y*sin(p.theta));
		trans_obs.y = p.y +(obs.x*sin(p.theta)) + (obs.y*cos(p.theta));
		trans_obs.id = -1;
		

	        double min_dist = numeric_limits<float>::max();
		for (const auto &L: landmarksInRange){
			auto &pred = L.second;
			double curr_dist = dist(trans_obs.x, trans_obs.y, pred.x, pred.y);
			if(curr_dist < min_dist)
			{
				min_dist = curr_dist;
				trans_obs.id = pred.id;
			}
		}

		trans_observations.push_back(trans_obs);
	}

		
	vector<int> associations;
	vector<double> sense_x;
	vector<double> sense_y;

	for (const auto& obs: trans_observations)
	{
		LandmarkObs& association = landmarksInRange.find(obs.id)->second;
		associations.push_back(association.id);
		sense_x.push_back(association.x);
		sense_y.push_back(association.y);
		
		double diff_x = association.x - obs.x;
		double diff_y = association.y - obs.y;
		double p1 = (diff_x*diff_x)/(2*landmark_x*landmark_x);
		double p2 = (diff_y*diff_y)/(2*landmark_y*landmark_y);
		
		double weight = exp(-(p1 + p2)) / ( 2 * M_PI * landmark_x * landmark_y);
		p.weight *= weight;

	}

	SetAssociations(p, associations, sense_x, sense_y);
	weights.push_back(p.weight);
}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	vector<Particle> resample_particles;
	
	for (int i=0; i<num_particles; i++)
	{
		resample_particles.push_back(particles[distribution(gen)]);	
	}

	particles = resample_particles;
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
