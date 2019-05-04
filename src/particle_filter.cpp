/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);
    
  num_particles = 20;  // TODO: Set the number of particles

  //Sample from the x,y and z normal distrubtions and concatenate to create particles
  for (int i=0; i<num_particles; i++) {
    Particle P = {0};
    P.id = i;
    P.x = dist_x(gen);
    P.y = dist_y(gen);
    P.theta = dist_theta(gen);
    P.weight = 1;

    //push particle to vector
    weights.push_back(1); 
    particles.push_back(P); 
  }

  //set initialization flag
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  double pred_x, pred_y, pred_theta;
  std::default_random_engine gen;  
    
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);
    
  for (int i=0; i<num_particles; i++)
  {
      //if yaw_rate is too small, it may cause error
      if( fabs(yaw_rate) > 0.0001 ){
        pred_theta = particles[i].theta + delta_t*yaw_rate;
        pred_x = particles[i].x + velocity/yaw_rate*(sin(pred_theta) -sin(particles[i].theta));
        pred_y = particles[i].y + velocity/yaw_rate*(-cos(pred_theta) +cos(particles[i].theta));
      }else{
        pred_x = particles[i].x + (velocity * delta_t * cos(particles[i].theta));
        pred_y = particles[i].y + (velocity * delta_t * sin(particles[i].theta));
        pred_theta = particles[i].theta;
      }
          
    particles[i].x = dist_x(gen) + pred_x;
    particles[i].y = dist_y(gen) + pred_y;
    particles[i].theta = dist_theta(gen) + pred_theta;
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    
  for(LandmarkObs &obsv: observations){      
    double min = std::numeric_limits<double>::max();
      
      for(LandmarkObs &pred: predicted){
        double d = dist(obsv.x, obsv.y, pred.x, pred.y);
        if(d < min){
          min = d;
          // changing the id of the observation to match the id of the landmark predicted
          obsv.id = pred.id;
        }
     }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // get all observations in the range of the particle
  for (int i=0; i<num_particles; i++){  
    //check range, get close landmarks
    vector<LandmarkObs> predictions;  
    
    for(unsigned int landmark = 0 ;landmark<map_landmarks.landmark_list.size();landmark++){
      double d = dist(map_landmarks.landmark_list[landmark].x_f, 
                      map_landmarks.landmark_list[landmark].y_f, particles[i].x, particles[i].y); 
      if(d < sensor_range){
        LandmarkObs new_obs = {};
        new_obs.id = map_landmarks.landmark_list[landmark].id_i;
        new_obs.x  = map_landmarks.landmark_list[landmark].x_f;
        new_obs.y  = map_landmarks.landmark_list[landmark].y_f;
        predictions.push_back(new_obs);
      }
    }
      
    vector<LandmarkObs> observations_trans;
    // transform each observation to the map coordinate system
    for (unsigned int j=0;j<observations.size();j++){
      LandmarkObs obs_trans;
      obs_trans.id = observations[j].id;
      obs_trans.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) 
                    - (sin(particles[i].theta) * observations[j].y);
      obs_trans.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) 
                    + (cos(particles[i].theta) * observations[j].y);
      observations_trans.push_back(obs_trans);
    }
    
    dataAssociation(predictions,observations_trans);
    // after this step, each observation should have x,y in world position and should also have a pair in the observations
    // Now the important thing is to update the weights
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
    //compute weights
    double weight=1;
    double std_2_pi = 2.0*M_PI*std_landmark[0]*std_landmark[1];
    double std_2_x  = 2.0*std_landmark[0]*std_landmark[0];
    double std_2_y = 2.0*std_landmark[1]*std_landmark[1];
    
     // iteration the vector obs_trans where each is an struct of the landmark type 
    for(LandmarkObs &obsv : observations_trans)
    {
      // remeber id = 1 is the landmarklist 0 for the map
      Map::single_landmark_s mark = map_landmarks.landmark_list[obsv.id -1];
      double e_x = std::pow(obsv.x - mark.x_f, 2);
      double e_y = std::pow(obsv.y - mark.y_f, 2);     
      double e = (e_x/std_2_x + e_y/std_2_y);

      double ee = exp(-e);
      double w = ee/std_2_pi;
      
      //prod of all weights
      weight *= w;
      //record association
      associations.push_back(obsv.id);
      sense_x.push_back(obsv.x);
      sense_y.push_back(obsv.y); 
    }
    particles[i].weight = weight;
      
    //insert into weight vector
    weights[i]= weight;
    //update particle's associations
    SetAssociations(particles[i], associations, sense_x, sense_y); 
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::discrete_distribution<> distribution(weights.begin(), weights.end());
  vector<Particle> n_particles;
 
  //resample  
  for(int i=0; i < num_particles; i++){      
    const int index = distribution(gen);  
    n_particles.push_back(particles[index]);
  }

  //new set of particles
  particles = n_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}