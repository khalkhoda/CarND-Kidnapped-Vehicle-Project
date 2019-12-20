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
#include <cstring>
#include <stdio.h>
#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::cout;
using std::endl;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles

  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // Create a normal (Gaussian) distribution centered at x
  normal_distribution<double> dist_x(x, std_x);
  // Create a normal (Gaussian) distribution centered at y
  normal_distribution<double> dist_y(y, std_y);
  // Create a normal (Gaussian) distribution centered at theta
  normal_distribution<double> dist_theta(theta, std_theta);

  for (int i = 0; i < num_particles; ++i)
  {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p);
    weights.push_back(p.weight);
    // std::cout << "p.x, p.y, p.theta: " << p.x << "," << p.y << "," << p.theta << std::endl;
  }

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

   double std_x = std_pos[0];
   double std_y = std_pos[1];
   double std_theta = std_pos[2];
   double r;
   double delta_theta;
   double new_theta;
   double deltaDist = velocity * delta_t;
   std::default_random_engine gen;

   if (yaw_rate != 0.0){
     r = velocity / yaw_rate;
     delta_theta = yaw_rate*delta_t;
   }

   for (int i = 0; i < num_particles; ++i)
   {
     Particle &p = particles[i];
     double theta = p.theta;
     if (abs(yaw_rate) > 0.001)
     {
       new_theta = theta + delta_theta;
       p.x = p.x + r * ( sin(new_theta) - sin(theta));
       p.y = p.y + r * (cos(theta) - cos(new_theta));
       p.theta = new_theta;
     }
     else
     {
       p.x = p.x + deltaDist * cos(theta);
       p.y = p.y + deltaDist * sin(theta);
     }
     //
     /*Add noise*/
     // Create a normal (Gaussian) distribution for x
     normal_distribution<double> dist_x(p.x, std_x);
     // Create a normal (Gaussian) distribution for y
     normal_distribution<double> dist_y(p.y, std_y);
     // Create a normal (Gaussian) distribution for theta
     normal_distribution<double> dist_theta(p.theta, std_theta);

     /* Draw a random sample around the predicted particle to simulate
     noise in the control signals*/
     p.x = dist_x(gen);
     p.y = dist_y(gen);
     p.theta = dist_theta(gen);

   }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predObss,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predObss measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */


   LandmarkObs predLandmark;
   double x_meas, x_pred;
   double y_meas, y_pred;
   double measToPredDistTmp, measToPredDistNearest;
   // double x_associated, y_associated;
   for (unsigned int i_obs = 0; i_obs < observations.size(); i_obs++)
   {
     measToPredDistNearest = 1000000.0; // Any relatively big number

     // Get a reference to next observation to work with
     LandmarkObs &landmarkMeas = observations[i_obs];
     landmarkMeas.id = -1;
     x_meas = landmarkMeas.x;
     y_meas = landmarkMeas.y;
     // x_associated = -1.0;
     // y_associated = -1.0;
     /*Loop through all predicted measurements to choose the nearest one to
     the actual measurements as a solution to the association problem*/
     for (unsigned int j_pred = 0; j_pred < predObss.size(); j_pred++)
     {
       predLandmark = predObss[j_pred];
       x_pred = predLandmark.x;
       y_pred = predLandmark.y;
       // printf ("Predicted(%d): %4.2f, %4.2f ------------ \n",
       // predLandmark.id, x_pred, y_pred);
       // printf ("prediction(%d): %4.2f, %4.2f\n", predLandmark.id, x_pred, x_pred);
       measToPredDistTmp = dist(x_meas, y_meas, x_pred, y_pred);
       if (measToPredDistTmp < measToPredDistNearest)
       {
         // printf ("measToPredDistNearest: %4.2f\n", measToPredDistTmp);
         measToPredDistNearest = measToPredDistTmp;
         // x_associated = x_pred;
         // y_associated = y_pred;
         landmarkMeas.id = predLandmark.id;
       }
     }
     // printf ("measurement[%d] (%4.2f, %4.2f) is associated with prediction (%4.2f, %4.2f)\n",
     // observations[i_obs], x_meas, y_meas, x_associated, y_associated);
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   vector<LandmarkObs> &observations,
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

   for (int i = 0; i < num_particles; ++i)
   {
     Particle &p = particles[i];
     p.associations.clear();
     p.sense_x.clear();
     p.sense_y.clear();
     // printf ("p[%d]: %4.2f, %4.2f, %4.2f\n", i, p.x, p.y, p.theta);

     /* Predict the sensor measurements to all the map landmarks within the
     sensor range for each particle.*/

     /* Update the weights and resample */

     auto start = high_resolution_clock::now();
     vector<LandmarkObs> predObss;
     predictObss(map_landmarks, sensor_range, p, predObss);
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
     std::cout << "UpdateWeights: Predict measurements step time: "<<duration.count() << " microseconds"<<std::endl;

     // printf("predObss landmarks vector size:  %d\n", predObss.size());
     // for (int i_pred = 0; i_pred < predObss.size(); i_pred++){
     //   printf("predObss landmark (x, y) = (%4.2f, %4.2f)\n", predObss[i_pred].x, predObss[i_pred].y);
     // }
     // printf("measured landmarks vector size:  %d\n", observations.size());
     // for (int i = 0; i < observations.size(); i++){
     //   printf("predObss landmark (x, y) = (%4.2f, %4.2f)\n", observations[i].x, observations[i].y);
     // }
     // Assign observations to map landmarks
     // printf ("Calling Data association\n");
     // p.x = 4.0;
     // p.y = 6.0;
     // p.theta = M_PI/2;
     // printf ("p[%d]: %4.2f, %4.2f, %4.2f\n", i, p.x, p.y, p.theta);

     /* Transform observations from vehicle to map coordinates*/
     std::vector<LandmarkObs> transObss;
     transformObss(observations, p, transObss);

     /* Associate measurements with map landmarks using
     the predicted measurements*/
     dataAssociation(predObss, transObss);

     /*Copy the association result to the current particle structure
     This will be used in the simulator to visualize the association result
     of the best particle*/
     for (unsigned int i_obs = 0; i_obs < transObss.size(); i_obs++){
       LandmarkObs &transObs = transObss[i_obs];
       p.associations.push_back(transObs.id);
       p.sense_x.push_back(transObs.x);
       p.sense_y.push_back(transObs.y);
     }

     /* Calculate the new weight of the current particle
     using a mult-variate Gaussian
     printf ("Calling calculate_weight\n");*/
     p.weight = calculate_weight(predObss, transObss, std_landmark);
     // printf ("Weight of p[%d] is %4.2f\n", i, p.weight);
     weights[i] = p.weight;

   }

   // Normalize particles weights
   // Find sum of weights
   double sum_weights = 0.0;
   for (int i = 0; i < num_particles; ++i)
   {
     sum_weights += weights[i];
   }
   // printf ("Sum of weights is %4.4f\n", sum_weights);
   for (int i = 0; i < num_particles; ++i)
   {
     particles[i].weight = particles[i].weight / sum_weights;
     weights[i] = particles[i].weight;
   }

}

double ParticleFilter::calculate_weight(const vector<LandmarkObs>& predObss,
                        const vector<LandmarkObs>& observations,
                        double std_landmark[])
{
  double par_weight = 1.0;
  double std_x = std_landmark[0];
  double std_y = std_landmark[1];
  for (unsigned int i = 0; i < observations.size(); i++)
  {
    LandmarkObs landmarkMeas = observations[i];
    // landmarkMeas.id, landmarkMeas.x, landmarkMeas.y);
    for (unsigned int j = 0; j < predObss.size(); j++)
    {
      LandmarkObs predLandmark = predObss[j];
      // printf ("predLandmark(%d) (%4.2f, %4.2f)\n",
      // predLandmark.id, predLandmark.x, predLandmark.y);
      if (predLandmark.id == landmarkMeas.id)
      {
        /* The association of the predicted measurement and
        the actual measurement  is found.*/
        double tmp = multiv_prob(std_x, std_y,
                                  landmarkMeas.x, landmarkMeas.y,
                                  predLandmark.x, predLandmark.y);
        // printf ("multiv_prob of (%4.2f, %4.2f) & (%4.2f, %4.2f) is %4.9f\n",
        //        landmarkMeas.x, landmarkMeas.y,
        //        predLandmark.x, predLandmark.y,
        //        tmp);
        par_weight *= tmp;
      }
    }
  }
  return par_weight;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   std::default_random_engine generator;
   // std::cout << "discrete_distribution initializing " <<std::endl;
   std::discrete_distribution<int> d (weights.begin(), weights.end());
   // Set of new particles
   std::vector<Particle> particles_resampled;
   for (int i = 0; i < num_particles; i++)
   {
     // std::cout << "discrete_distribution done " << d(generator) <<std::endl;
     // std::cout << "Particles[i](x, y): " << particles[i].x<< ","<<particles[i].y <<std::endl;
     particles_resampled.push_back(particles[d(generator)]);
   }
   //Copy back the new resampled particles to the particles member
   for (int i = 0; i < num_particles; i++)
   {
     particles[i] = particles_resampled[i];
   }

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

void ParticleFilter::transformObss(vector<LandmarkObs> const &observations,
                                   Particle const &p,
                                   vector<LandmarkObs> &transObss){
  for (unsigned int i_obs = 0; i_obs < observations.size(); i_obs++){
    LandmarkObs const &obs = observations[i_obs];
    LandmarkObs transObs;
    // printf ("Obs (In vehicle)(%d): %4.2f, %4.2f\n", obs.id, obs.x, obs.y);
    transObs.x = p.x + cos(p.theta)*obs.x - sin(p.theta)*obs.y;
    transObs.y = p.y + sin(p.theta)*obs.x + cos(p.theta)*obs.y;
    transObss.push_back(transObs);
    // printf ("Obs (In world)(%d): %4.2f, %4.2f\n", obs.id, transObs.x, transObs.y);
  }
}

void ParticleFilter::predictObss(const Map &map_landmarks,
                                 double sensor_range,
                                 Particle const &p,
                                 vector<LandmarkObs> &predObss){
  for (unsigned int i_landmark = 0; i_landmark < map_landmarks.landmark_list.size(); i_landmark++)
  {
    double x_f = (double)map_landmarks.landmark_list[i_landmark].x_f;
    double y_f = (double)map_landmarks.landmark_list[i_landmark].y_f;
    int lm_id  = (double)map_landmarks.landmark_list[i_landmark].id_i;
    // printf ("map_landmark_GT(%3d): %4.2f, %4.2f\n", lm_id, x_f, y_f);
    if (dist(x_f, y_f, p.x, p.y) < sensor_range){
      LandmarkObs predLandmark;
      predLandmark.x = x_f; // In world coordinates
      predLandmark.y = y_f; // In world coordinates
      predLandmark.id = lm_id;
      predObss.push_back(predLandmark);
    }
  }
}
