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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if (is_initialized) {
	  return;
	} 
	else {
		num_particles= 100;

		default_random_engine gen;

		normal_distribution<double> N_x(x, std[0]);
		normal_distribution<double> N_y(y, std[1]);
		normal_distribution<double> N_theta(theta, std[2]);
		for(int i=0; i<num_particles; i++)
		{
			Particle p;
			p.id = i;
			p.x = N_x(gen) ;
			p.y = N_y(gen);
			p.theta = N_theta(gen);
			p.weight = 1.0;
			weights.push_back(p.weight);
			particles.push_back(p);
		}
	
		is_initialized = true;
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	default_random_engine gen;

	normal_distribution<double> N_x(0.0, std_pos[0]);
	normal_distribution<double> N_y(0.0, std_pos[1]);
	normal_distribution<double> N_theta(0.0, std_pos[2]);

	for (unsigned int i = 0; i < num_particles; i++)
	{
		//Particle p = particles[i]; 
		double noise_x = N_x(gen);
		double noise_y = N_y(gen);
		double noise_theta = N_theta(gen);

		if (abs(yaw_rate)<0.00001){
			double v_dt = velocity*delta_t;
			particles[i].x += (v_dt*cos(particles[i].theta)) + noise_x;
			particles[i].y += (v_dt*sin(particles[i].theta)) + noise_y; 
			particles[i].theta += noise_theta;
		}
		else{
			double v_dt = velocity/yaw_rate;
			particles[i].x += (v_dt)*(sin(particles[i].theta+(yaw_rate*delta_t)) - sin(particles[i].theta)) + noise_x;
			particles[i].y += (v_dt)*(cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate*delta_t))) + noise_y;
			particles[i].theta += (yaw_rate*delta_t) + noise_theta;
		}

	} 
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	for (int i = 0; i < observations.size(); i++) {
	    int near_id = -1;
	    double min_dist = numeric_limits<double>::max();

	    for (int j = 0; j < predicted.size(); j++) {
	      double delta_x = predicted[j].x - observations[i].x;
	      double delta_y = predicted[j].y - observations[i].y;
	      double curr_dist = delta_x * delta_x + delta_y * delta_y;

	      if (curr_dist < min_dist) {
		near_id = j;
	        min_dist = curr_dist;
	       }
    	    }
	    observations[i].id = near_id;
  	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	double d_0 = sqrt(2.0 * M_PI * std_x * std_y);
	double d_1 = 2 * std_x * std_x;
	double d_2 = 2 * std_y * std_y;

	for(int i=0;i<num_particles;i++){

	    double p_x = particles[i].x;
    	    double p_y = particles[i].y;
	    double p_theta = particles[i].theta;

	    vector<LandmarkObs> landmarks_in_range;
	    vector<LandmarkObs> transformed_landmarks;
	    
	    //Transformation
    
	    for (int j = 0; j < observations.size(); j++) {

		LandmarkObs transformed_landmark;		
		transformed_landmark.id = observations[j].id;

	        transformed_landmark.x = p_x + (observations[j].x * cos(p_theta)) - (observations[j].y * sin(p_theta));
		transformed_landmark.y = p_y + (observations[j].y * cos(p_theta)) + (observations[j].x * sin(p_theta));

	        transformed_landmarks.push_back(transformed_landmark);
    	     }
			
	     //Landmarks within sensor Range	
	     for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
	
		      int landmark_id = map_landmarks.landmark_list[j].id_i;
      		      double landmark_x = map_landmarks.landmark_list[j].x_f;
		      double landmark_y = map_landmarks.landmark_list[j].y_f;
	
		      double dist_x = landmark_x - p_x;
	      	      double dist_y = landmark_y - p_y;
		      double dist = sqrt(dist_x * dist_x + dist_y * dist_y);

		      if (dist < sensor_range) {
        			LandmarkObs landmarksInRange;
        			landmarksInRange.id = landmark_id;
			        landmarksInRange.x = landmark_x;
			        landmarksInRange.y = landmark_y;

			        landmarks_in_range.push_back(landmarksInRange);
      			}
    		}
             
	     dataAssociation(landmarks_in_range, transformed_landmarks);
	     double weight = 1.0;
 

	     for (int j = 0; j < transformed_landmarks.size(); j++) {
		
		      int obs_id = transformed_landmarks[j].id;
		      double obs_x = transformed_landmarks[j].x;
		      double obs_y = transformed_landmarks[j].y;

		      double pred_x = landmarks_in_range[obs_id].x;
		      double pred_y = landmarks_in_range[obs_id].y;

      		      double diff_x = obs_x - pred_x;
		      double diff_y = obs_y - pred_y;

		      double p1 = diff_x * diff_x / d_1;
		      double p2 = diff_y * diff_y / d_2;

		      weight *= exp(-(p1 + p2)) / d_0;
    		}

	    if (weight == 0) {
		      particles[i].weight = 0.00001;
		      weights[i] = 0.00001;
	    } 
	    else {
		      particles[i].weight = weight;
		      weights[i] = weight;
    	    }
  	}

}





//	weights.clear();


/*	for (auto& p:particles){

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
	cout << "Inside update" <<p.weight<< endl;

}*/


//}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> index(weights.begin(), weights.end());
	vector<Particle> resample_particles;
	
	for (int i=0; i<num_particles; i++)
	{
		
		const int idx = index(gen);
		Particle p;
		p.id = idx;
		p.x = particles[idx].x;
		p.y = particles[idx].y;
		p.theta = particles[idx].theta;
		p.weight = 1.0;
		resample_particles.push_back(p);	
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
