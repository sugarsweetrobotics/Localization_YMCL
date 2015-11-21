// -*- C++ -*-
/*!
 * @file  Localization_YMCL.h
 * @brief Localization YMCL Component
 * @date  $Date$
 *
 * GPLv3
 *
 * $Id$
 */

#ifndef LOCALIZATION_YMCL_H
#define LOCALIZATION_YMCL_H
#include "MobileRobotStub.h"
#include <rtm/idl/BasicDataTypeSkel.h>
#include <rtm/idl/ExtendedDataTypesSkel.h>
#include <rtm/idl/InterfaceDataTypesSkel.h>

#include <rtm/Manager.h>
#include <rtm/DataFlowComponentBase.h>
#include <rtm/CorbaPort.h>
#include <rtm/DataInPort.h>
#include <rtm/DataOutPort.h>
// Service implementation headers
// <rtc-template block="service_impl_h">

// </rtc-template>

// Service Consumer stub headers
// <rtc-template block="consumer_stub_h">


// </rtc-template>

using namespace RTC;
#define _YMCL_USE_STATIC
#include "ymcl.h"

#include <opencv2/opencv.hpp>
/*!
 * @class Localization_YMCL
 * @brief Localization YMCL Component
 *
 * Localization RT-component using Monte Carlo Localization.
 */
class Localization_YMCL
  : public RTC::DataFlowComponentBase
{
 public:
  /*!
   * @brief constructor
   * @param manager Maneger Object
   */
  Localization_YMCL(RTC::Manager* manager);

  /*!
   * @brief destructor
   */
  ~Localization_YMCL();

  // <rtc-template block="public_attribute">
  
  // </rtc-template>

  // <rtc-template block="public_operation">
  
  // </rtc-template>

  /***
   *
   * The initialize action (on CREATED->ALIVE transition)
   * formaer rtc_init_entry() 
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
   virtual RTC::ReturnCode_t onInitialize();

  /***
   *
   * The finalize action (on ALIVE->END transition)
   * formaer rtc_exiting_entry()
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onFinalize();

  /***
   *
   * The startup action when ExecutionContext startup
   * former rtc_starting_entry()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onStartup(RTC::UniqueId ec_id);

  /***
   *
   * The shutdown action when ExecutionContext stop
   * former rtc_stopping_entry()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onShutdown(RTC::UniqueId ec_id);

  /***
   *
   * The activated action (Active state entry action)
   * former rtc_active_entry()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
   virtual RTC::ReturnCode_t onActivated(RTC::UniqueId ec_id);

  /***
   *
   * The deactivated action (Active state exit action)
   * former rtc_active_exit()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
   virtual RTC::ReturnCode_t onDeactivated(RTC::UniqueId ec_id);

  /***
   *
   * The execution action that is invoked periodically
   * former rtc_active_do()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
   virtual RTC::ReturnCode_t onExecute(RTC::UniqueId ec_id);

  /***
   *
   * The aborting action when main logic error occurred.
   * former rtc_aborting_entry()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onAborting(RTC::UniqueId ec_id);

  /***
   *
   * The error action in ERROR state
   * former rtc_error_do()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onError(RTC::UniqueId ec_id);

  /***
   *
   * The reset action that is invoked resetting
   * This is same but different the former rtc_init_entry()
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onReset(RTC::UniqueId ec_id);
  
  /***
   *
   * The state update action that is invoked after onExecute() action
   * no corresponding operation exists in OpenRTm-aist-0.2.0
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onStateUpdate(RTC::UniqueId ec_id);

  /***
   *
   * The action that is invoked when execution context's rate is changed
   * no corresponding operation exists in OpenRTm-aist-0.2.0
   *
   * @param ec_id target ExecutionContext Id
   *
   * @return RTC::ReturnCode_t
   * 
   * 
   */
  // virtual RTC::ReturnCode_t onRateChanged(RTC::UniqueId ec_id);


 protected:
  // <rtc-template block="protected_attribute">
  
  // </rtc-template>

  // <rtc-template block="protected_operation">
  
  // </rtc-template>

  // Configuration variable declaration
  // <rtc-template block="config_declare">
  /*!
   * Initial pose of robot
   * - Name:  initial_pose_x
   * - DefaultValue: 0
   */
  float m_initial_pose_x;
  /*!
   * Initial pose of robot
   * - Name:  initial_pose_y
   * - DefaultValue: 0
   */
  float m_initial_pose_y;
  /*!
   * Initial pose of robot
   * - Name:  initial_pose_phi
   * - DefaultValue: 0
   */
  float m_initial_pose_phi;
  /*!
   * minimum value range of initial particle position (This is
   * added to initial position)
   * - Name:  initial_particle_min_x
   * - DefaultValue: -0.3
   */
  float m_initial_particle_min_x;
  /*!
   * maximum value range of initial particle position (This is
   * added to initial position)
   * - Name:  initial_particle_max_x
   * - DefaultValue: 0.3
   */
  float m_initial_particle_max_x;
  /*!
   * minimum value range of initial particle position (This is
   * added to initial position)
   * - Name:  initial_particle_min_y
   * - DefaultValue: -0.3
   */
  float m_initial_particle_min_y;
  /*!
   * maximum value range of initial particle position (This is
   * added to initial position)
   * - Name:  initial_particle_max_y
   * - DefaultValue: 0.3
   */
  float m_initial_particle_max_y;
  /*!
   * minimum value range of initial particle position (This is
   * added to initial position)
   * - Name:  initial_particle_min_phi
   * - DefaultValue: -0.3
   */
  float m_initial_particle_min_phi;
  /*!
   * maximum value range of initial particle position (This is
   * added to initial position)
   * - Name:  initial_particle_max_phi
   * - DefaultValue: 0.3
   */
  float m_initial_particle_max_phi;
  /*!
   * Initial particle count
   * - Name:  initial_particle_count
   * - DefaultValue: 1000
   */
  long m_initial_particle_count;

  float m_initial_particle_std_xy;
  float m_initial_particle_std_phi;


  /*!
   * Timeout threshold [sec] of odometry port.
   * - Name:  poseTimeOut
   * - DefaultValue: 3.0
   */
  float m_poseTimeOut;
  /*!
   * Type name for Motion model.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_model
   * - DefaultValue: Gausian
   */
  std::string m_motion_model;
  /*!
   * Parameter of motion_model.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_alpha1
   * - DefaultValue: 0.01
   */
  float m_motion_alpha1;
  /*!
   * Parameter of motion_model.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_alpha2
   * - DefaultValue: 0.05729
   */
  float m_motion_alpha2;
  /*!
   * Parameter of motion_model.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_alpha3
   * - DefaultValue: 0.01745
   */
  float m_motion_alpha3;
  /*!
   * Parameter of motion_model.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_alpha4
   * - DefaultValue: 0.05
   */
  float m_motion_alpha4;
  /*!
   * Parameter for Motion model.
   * if motion_model is Gausian, this is Sigma_min
   * if motion_model is Thrun, this is additional noise
   * parameter.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_std_XY
   * - DefaultValue: 0.01
   */
  float m_motion_std_XY;
  /*!
   * Parameter for Motion model.
   * if motion_model is Gausian, this is Sigma_min
   * if motion_model is Thrun, this is additional noise
   * parameter.
   * http://www.mrpt.org/tutorials/programming/odometry-and-motion
   * -models/probabilistic_motion_models/
   * - Name:  motion_std_PHI
   * - DefaultValue: 0.003490
   */
  float m_motion_std_PHI;
  /*!
   * The selected method to compute an observation likelihood.
   * http://www.mrpt.org/tutorials/programming/maps-for-localizati
   * on-slam-map-building/occupancy_grids/#42_Observations_likelih
   * ood
   * - Name:  LM_likelihoodMethod
   * - DefaultValue: lmLikelihoodField_Thrun
   */
  std::string m_LM_likelihoodMethod;
  /*!
   * Enables the usage of a cache of likelihood values (for LF
   * methods), if set to true (default=true).
   * - Name:  LM_enableLikelihoodCache
   * - DefaultValue: true
   */
  std::string m_LM_enableLikelihoodCache;
  /*!
   * [LikelihoodField] The laser range "sigma" used in
   * computations; Default value = 0.35
   * - Name:  LM_LF_decimation
   * - DefaultValue: 5
   */
  int m_LM_LF_decimation;
  /*!
   * The laser range "sigma" used in computations; Default value
   * = 0.35
   * - Name:  LM_LF_stdHit
   * - DefaultValue: 0.35
   */
  float m_LM_LF_stdHit;
  /*!
   * Ratios of the hit/random components of the likelihood;
   * Default values=0.05
   * - Name:  LM_LF_zRandom
   * - DefaultValue: 0.05
   */
  float m_LM_LF_zRandom;
  /*!
   * Set this to "true" ot use an alternative method, where the
   * likelihood of the whole range scan is computed by
   * "averaging" of individual ranges, instead of by the
   * "product". Default = false
   * - Name:  LM_LF_alternateAverageMethod
   * - DefaultValue: false
   */
  std::string m_LM_LF_alternateAverageMethod;
  /*!
   * Ratios of the hit/random components of the likelihood;
   * Default values=0.95
   * - Name:  LM_LF_zHit
   * - DefaultValue: 0.95
   */
  float m_LM_LF_zHit;
  /*!
   * The max. distance for searching correspondences around each
   * sensed point default  0.3
   * - Name:  LM_LF_maxCorrsDistance
   * - DefaultValue: 0.3
   */
  float m_LM_LF_maxCorrsDistance;
  /*!
   * The max. range of the sensor (def=81meters)
   * - Name:  LM_LF_maxRange
   * - DefaultValue: 81
   */
  float m_LM_LF_maxRange;
  /*!
   * The exponent in the MI likelihood computation. Default value
   * = 2.5
   * - Name:  LM_MI_exponent
   * - DefaultValue: 2.5
   */
  float m_LM_MI_exponent;
  /*!
   * [MI] The ratio for the max. distance used in the MI
   * computation and in the insertion of scans, e.g. if set to
   * 2.0 the MI will use twice the distance that the update
   * distance. def=1.5
   * - Name:  LM_MI_ratio_max_distance
   * - DefaultValue: 1.5
   */
  float m_LM_MI_ratio_max_distance;
  /*!
   * The scan rays decimation: at every N rays, one will be used
   * to compute the MI: def=10
   * - Name:  LM_MI_skip_rays
   * - DefaultValue: 10
   */
  float m_LM_MI_skip_rays;
  /*!
   * The power factor for the likelihood (default=5)
   * - Name:  LM_consensus_pow
   * - DefaultValue: 5
   */
  float m_LM_consensus_pow;
  /*!
   * The down-sample ratio of ranges (default=1, consider all the
   * ranges)
   * - Name:  LM_consensus_takeEachRange
   * - DefaultValue: 1
   */
  int m_LM_consensus_takeEachRange;
  /*!
   * [rayTracing] The laser range sigma. def=1.0
   * - Name:  LM_rayTracing_stdHit
   * - DefaultValue: 1.0
   */
  float m_LM_rayTracing_stdHit;
  /*!
   * One out of "rayTracing_decimation" rays will be simulated
   * and compared only: set to 1 to use all the sensed ranges.
   * def=10
   * - Name:  LM_rayTracing_decimation
   * - DefaultValue: 10
   */
  int m_LM_rayTracing_decimation;
  /*!
   * If true (default), the rayTracing method will ignore
   * measured ranges shorter than the simulated ones.
   * - Name:  LM_rayTracing_useDistanceFilter
   * - DefaultValue: true
   */
  std::string m_LM_rayTracing_useDistanceFilter;
  /*!
   * Sequential Importance Resampling – SIR (pfStandardProposal)
   * Standard proposal distribution + weights according to
   * likelihood function.
   * Auxiliary Particle Filter – APF (pfAuxiliaryPFStandard)
   * This method was introduced by Pitt and Shephard in 1999 [1]
   * Optimal Sampling (pfOptimalProposal)
   * Use the exact optimal proposal distribution (where
   * available!, usually this will perform approximations). In
   * the case of the RBPF-SLAM implementation, this method
   * follows [2]
   * Approximate Optimal Sampling (pfAuxiliaryPFOptimal)
   * Use the optimal proposal and a auxiliary particle filter
   * (see paper [3] ).
   * See :
   * http://www.mrpt.org/tutorials/programming/statistics-and-baye
   * s-filtering/particle_filter_algorithms/
   * - Name:  PF_algorithm
   * - DefaultValue: pfStandardProposal
   */
  std::string m_PF_algorithm;
  /*!
   * Setting of resampling method.
   * See:
   * http://www.mrpt.org/tutorials/programming/statistics-and-baye
   * s-filtering/resampling_schemes/
   * - Name:  PF_resamplingMethod
   * - DefaultValue: prMultinomial
   */
  std::string m_PF_resamplingMethod;
  /*!
   * The resampling of particles will be performed when ESS (in
   * range [0,1]) < BETA (default is 0.5)
   * - Name:  PF_BETA
   * - DefaultValue: 0.5
   */
  float m_PF_BETA;
  /*!
   * An optional step to "smooth" dramatic changes in the
   * observation model to affect the variance of the particle
   * weights, eg weight*=likelihood^powFactor (default=1 = no
   * effects).
   * - Name:  PF_powFactor
   * - DefaultValue: 1.0
   */
  float m_PF_powFactor;
  /*!
   * The initial number of particles in the filter (it can change
   * only if adaptiveSampleSize=true) (default=1)
   * - Name:  PF_sampleSize
   * - DefaultValue: 1
   */
  std::string m_PF_sampleSize;
  /*!
   * A flag that indicates whether the CParticleFilterCapable
   * object should perform adative sample size (default=true).
   * - Name:  PF_adaptiveSampleSize
   * - DefaultValue: true
   */
  std::string m_PF_adaptiveSampleSize;
  /*!
   * Only for PF_algorithm=pfAuxiliaryPFOptimal: If a given
   * particle has a max_likelihood (from the a-priori estimate)
   * below the maximum from all the samples - This is done to
   * assure that the rejection sampling doesn't get stuck in an
   * infinite loop trying to get an acceptable sample.  Default =
   * 15 (in logarithmic likelihood)
   * - Name:  PF_max_loglikelihood_dyn_range
   * - DefaultValue: 15
   */
  double m_PF_max_loglikelihood_dyn_range;
  /*!
   * In the algorithm "CParticleFilter::pfAuxiliaryPFOptimal"
   * (and in "CParticleFilter::pfAuxiliaryPFStandard" only if
   * pfAuxFilterStandard_FirstStageWeightsMonteCarlo = true) the
   * number of samples for searching the maximum likelihood value
   * and also to estimate the "first stage weights" (see papers!)
   * (default=100)
   * - Name:  PF_AuxFilterOptimal_MaximumSearchSamples
   * - DefaultValue: 100
   */
  int m_PF_AuxFilterOptimal_MaximumSearchSamples;
  /*!
   * Only for PF_algorithm==pfAuxiliaryPFStandard:  If false, the
   * APF will predict the first stage weights just at the mean of
   * the prior of the next time step.  If true, these weights
   * will be estimated as described in the papers for the
   * "pfAuxiliaryPFOptimal" method, i.e. through a monte carlo
   * simulation.  In that case,
   * "pfAuxFilterOptimal_MaximumSearchSamples" is the number of
   * MC samples used.
   * - Name:  PF_AuxFilterStandard_FirstStageWeightsMonteCarlo
   * - DefaultValue: false
   */
  std::string m_PF_AuxFilterStandard_FirstStageWeightsMonteCarlo;
  /*!
   * (Default=false) In the algorithm
   * "CParticleFilter::pfAuxiliaryPFOptimal", if set to true, do
   * not perform rejection sampling, but just the most-likely
   * (ML) particle found in the preliminary weight-determination
   * stage.
   * - Name:  PF_AuxFilterOptimal_MLE
   * - DefaultValue: false
   */
  std::string m_PF_AuxFilterOptimal_MLE;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_binSize_PHI
   * - DefaultValue: 0.01
   */
  float m_KLD_binSize_PHI;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_binSize_XY
   * - DefaultValue: 0.01
   */
  float m_KLD_binSize_XY;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_delta
   * - DefaultValue: 0.02
   */
  float m_KLD_delta;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_epsilon
   * - DefaultValue: 0.02
   */
  float m_KLD_epsilon;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_maxSampleSize
   * - DefaultValue: 1000
   */
  int m_KLD_maxSampleSize;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_minSampleSize
   * - DefaultValue: 150
   */
  int m_KLD_minSampleSize;
  /*!
   * Parameters for the KLD adaptive sample size algorithm (see
   * Dieter Fox's papers)
   * - Name:  KLD_minSamplesPerBin
   * - DefaultValue: 0
   */
  double m_KLD_minSamplesPerBin;

  // </rtc-template>

  // DataInPort declaration
  // <rtc-template block="inport_declare">
  RTC::RangeData m_range;
  /*!
   * Range Sensor Data (usually, LiDAR)
   */
  InPort<RTC::RangeData> m_rangeIn;
  RTC::TimedPose2D m_odometry;
  /*!
   * Robot position. Usually, estimated with encoders.
   */
  InPort<RTC::TimedPose2D> m_odometryIn;
  
  // </rtc-template>


  // DataOutPort declaration
  // <rtc-template block="outport_declare">
  RTC::TimedPose2D m_estimatedPose;
  /*!
   * Estimated pose output.
   */
  OutPort<RTC::TimedPose2D> m_estimatedPoseOut;
  
  // </rtc-template>

  // CORBA Port declaration
  // <rtc-template block="corbaport_declare">
  /*!
   */
  RTC::CorbaPort m_mapServerPort;
  
  // </rtc-template>

  // Service declaration
  // <rtc-template block="service_declare">
  
  // </rtc-template>

  // Consumer declaration
  // <rtc-template block="consumer_declare">
  /*!
   * Map Server Interface
   */
  RTC::CorbaConsumer<RTC::OGMapServer> m_mapServer;
  
  // </rtc-template>

 private:
  // <rtc-template block="private_attribute">
  
  // </rtc-template>

  // <rtc-template block="private_operation">
  
  // </rtc-template>

	 OGMap* m_pOGMap;

	 ymcl_t m_ymcl;

	 void setRanger(ranger_data_t* ranger, RangeData& range) {
		 ranger_setconfig(ranger,
			 range.config.minAngle,
			 range.config.maxAngle,
			 range.config.angularRes,
			 range.config.minRange,
			 range.config.maxRange);
		 ranger_setoffset(ranger,
			 range.geometry.geometry.pose.position.x,
			 range.geometry.geometry.pose.position.y,
			 range.geometry.geometry.pose.position.z,
			 range.geometry.geometry.pose.orientation.y,
			 range.geometry.geometry.pose.orientation.r,
			 range.geometry.geometry.pose.orientation.p);
		 ranger_setdata(ranger,
			 &(range.ranges[0]),
			 range.ranges.length());
	 }

	 TimedPose2D m_lastResampledPose;

};


extern "C"
{
  DLL_EXPORT void Localization_YMCLInit(RTC::Manager* manager);
};

#endif // LOCALIZATION_YMCL_H
