// -*- C++ -*-
/*!
 * @file  Localization_YMCL.cpp
 * @brief Localization YMCL Component
 * @date $Date$
 *
 * GPLv3
 *
 * $Id$
 */

#include "Localization_YMCL.h"

#define _USE_MATH_DEFINES
#include <math.h>


IplImage *pMapImage;
IplImage *pWorkImage;
IplImage *pPreviewImage;
RTC::OGMap* pOGMap;
std::string instanceName;

// Module specification
// <rtc-template block="module_spec">
static const char* localization_ymcl_spec[] =
  {
    "implementation_id", "Localization_YMCL",
    "type_name",         "Localization_YMCL",
    "description",       "Localization YMCL Component",
    "version",           "2.0.0",
    "vendor",            "Sugar Sweet Robotics",
    "category",          "Navigation",
    "activity_type",     "PERIODIC",
    "kind",              "DataFlowComponent",
    "max_instance",      "1",
    "language",          "C++",
    "lang_type",         "compile",
    // Configuration variables
    "conf.default.initial_pose_x", "0",
    "conf.default.initial_pose_y", "0",
    "conf.default.initial_pose_phi", "0",
    "conf.default.initial_particle_min_x", "-0.3",
    "conf.default.initial_particle_max_x", "0.3",
    "conf.default.initial_particle_min_y", "-0.3",
    "conf.default.initial_particle_max_y", "0.3",
    "conf.default.initial_particle_min_phi", "-0.3",
    "conf.default.initial_particle_max_phi", "0.3",
    "conf.default.initial_particle_count", "1000",
    "conf.default.poseTimeOut", "3.0",
    "conf.default.motion_model", "Gausian",
    "conf.default.motion_alpha1", "0.01",
    "conf.default.motion_alpha2", "0.05729",
    "conf.default.motion_alpha3", "0.01745",
    "conf.default.motion_alpha4", "0.05",
    "conf.default.motion_std_XY", "0.01",
    "conf.default.motion_std_PHI", "0.003490",
    "conf.default.LM_likelihoodMethod", "lmLikelihoodField_Thrun",
    "conf.default.LM_enableLikelihoodCache", "true",
    "conf.default.LM_LF_decimation", "5",
    "conf.default.LM_LF_stdHit", "0.35",
    "conf.default.LM_LF_zRandom", "0.05",
    "conf.default.LM_LF_alternateAverageMethod", "false",
    "conf.default.LM_LF_zHit", "0.95",
    "conf.default.LM_LF_maxCorrsDistance", "0.3",
    "conf.default.LM_LF_maxRange", "81",
    "conf.default.LM_MI_exponent", "2.5",
    "conf.default.LM_MI_ratio_max_distance", "1.5",
    "conf.default.LM_MI_skip_rays", "10",
    "conf.default.LM_consensus_pow", "5",
    "conf.default.LM_consensus_takeEachRange", "1",
    "conf.default.LM_rayTracing_stdHit", "1.0",
    "conf.default.LM_rayTracing_decimation", "10",
    "conf.default.LM_rayTracing_useDistanceFilter", "true",
    "conf.default.PF_algorithm", "pfStandardProposal",
    "conf.default.PF_resamplingMethod", "prMultinomial",
    "conf.default.PF_BETA", "0.5",
    "conf.default.PF_powFactor", "1.0",
    "conf.default.PF_sampleSize", "1",
    "conf.default.PF_adaptiveSampleSize", "true",
    "conf.default.PF_max_loglikelihood_dyn_range", "15",
    "conf.default.PF_AuxFilterOptimal_MaximumSearchSamples", "100",
    "conf.default.PF_AuxFilterStandard_FirstStageWeightsMonteCarlo", "false",
    "conf.default.PF_AuxFilterOptimal_MLE", "false",
    "conf.default.KLD_binSize_PHI", "0.1",
    "conf.default.KLD_binSize_XY", "0.1",
    "conf.default.KLD_delta", "0.01",
    "conf.default.KLD_epsilon", "0.05",
    "conf.default.KLD_maxSampleSize", "1000",
    "conf.default.KLD_minSampleSize", "300",
    "conf.default.KLD_minSamplesPerBin", "0",
    // Widget
    "conf.__widget__.initial_pose_x", "text",
    "conf.__widget__.initial_pose_y", "text",
    "conf.__widget__.initial_pose_phi", "text",
    "conf.__widget__.initial_particle_min_x", "text",
    "conf.__widget__.initial_particle_max_x", "text",
    "conf.__widget__.initial_particle_min_y", "text",
    "conf.__widget__.initial_particle_max_y", "text",
    "conf.__widget__.initial_particle_min_phi", "text",
    "conf.__widget__.initial_particle_max_phi", "text",
    "conf.__widget__.initial_particle_count", "text",
    "conf.__widget__.poseTimeOut", "text",
    "conf.__widget__.motion_model", "radio",
    "conf.__widget__.motion_alpha1", "text",
    "conf.__widget__.motion_alpha2", "text",
    "conf.__widget__.motion_alpha3", "text",
    "conf.__widget__.motion_alpha4", "text",
    "conf.__widget__.motion_std_XY", "text",
    "conf.__widget__.motion_std_PHI", "text",
    "conf.__widget__.LM_likelihoodMethod", "text",
    "conf.__widget__.LM_enableLikelihoodCache", "radio",
    "conf.__widget__.LM_LF_decimation", "text",
    "conf.__widget__.LM_LF_stdHit", "text",
    "conf.__widget__.LM_LF_zRandom", "text",
    "conf.__widget__.LM_LF_alternateAverageMethod", "radio",
    "conf.__widget__.LM_LF_zHit", "text",
    "conf.__widget__.LM_LF_maxCorrsDistance", "text",
    "conf.__widget__.LM_LF_maxRange", "text",
    "conf.__widget__.LM_MI_exponent", "text",
    "conf.__widget__.LM_MI_ratio_max_distance", "text",
    "conf.__widget__.LM_MI_skip_rays", "text",
    "conf.__widget__.LM_consensus_pow", "text",
    "conf.__widget__.LM_consensus_takeEachRange", "text",
    "conf.__widget__.LM_rayTracing_stdHit", "text",
    "conf.__widget__.LM_rayTracing_decimation", "text",
    "conf.__widget__.LM_rayTracing_useDistanceFilter", "radio",
    "conf.__widget__.PF_algorithm", "radio",
    "conf.__widget__.PF_resamplingMethod", "radio",
    "conf.__widget__.PF_BETA", "text",
    "conf.__widget__.PF_powFactor", "text",
    "conf.__widget__.PF_sampleSize", "text",
    "conf.__widget__.PF_adaptiveSampleSize", "text",
    "conf.__widget__.PF_max_loglikelihood_dyn_range", "text",
    "conf.__widget__.PF_AuxFilterOptimal_MaximumSearchSamples", "text",
    "conf.__widget__.PF_AuxFilterStandard_FirstStageWeightsMonteCarlo", "radio",
    "conf.__widget__.PF_AuxFilterOptimal_MLE", "radio",
    "conf.__widget__.KLD_binSize_PHI", "text",
    "conf.__widget__.KLD_binSize_XY", "text",
    "conf.__widget__.KLD_delta", "text",
    "conf.__widget__.KLD_epsilon", "text",
    "conf.__widget__.KLD_maxSampleSize", "text",
    "conf.__widget__.KLD_minSampleSize", "text",
    "conf.__widget__.KLD_minSamplesPerBin", "text",
    // Constraints
    "conf.__constraints__.motion_model", "(Thrun,Gausian)",
    "conf.__constraints__.LM_likelihoodMethod", "(lmLikelihoodField_Thrun,lmLikelihoodField_II,lmRayTracing,lmCellsDifference,lmConsensus,lmConsensusOWA,lmMeanInformation)",
    "conf.__constraints__.LM_enableLikelihoodCache", "(true,false)",
    "conf.__constraints__.LM_LF_alternateAverageMethod", "(true,false)",
    "conf.__constraints__.LM_rayTracing_useDistanceFilter", "(true,false)",
    "conf.__constraints__.PF_algorithm", "(pfStandardProposal,pfAuxiliaryPFStandard,pfOptimalProposal,pfAuxiliaryPFOptimal)",
    "conf.__constraints__.PF_resamplingMethod", "(prMultinomial,prResidual,prStratified,prSystematic)",
    "conf.__constraints__.PF_AuxFilterStandard_FirstStageWeightsMonteCarlo", "(true,false)",
    "conf.__constraints__.PF_AuxFilterOptimal_MLE", "(true,false)",
    ""
  };
// </rtc-template>

/*!
 * @brief constructor
 * @param manager Maneger Object
 */
Localization_YMCL::Localization_YMCL(RTC::Manager* manager)
    // <rtc-template block="initializer">
  : RTC::DataFlowComponentBase(manager),
    m_rangeIn("range", m_range),
    m_odometryIn("odometry", m_odometry),
    m_estimatedPoseOut("estimatedPose", m_estimatedPose),
    m_mapServerPort("mapServer")

    // </rtc-template>
{
}

/*!
 * @brief destructor
 */
Localization_YMCL::~Localization_YMCL()
{
}



RTC::ReturnCode_t Localization_YMCL::onInitialize()
{
  // Registration: InPort/OutPort/Service
  // <rtc-template block="registration">
  // Set InPort buffers
  addInPort("range", m_rangeIn);
  addInPort("odometry", m_odometryIn);
  
  // Set OutPort buffer
  addOutPort("estimatedPose", m_estimatedPoseOut);
  
  // Set service provider to Ports
  
  // Set service consumers to Ports
  m_mapServerPort.registerConsumer("mapServer", "RTC::OGMapServer", m_mapServer);
  
  // Set CORBA Service Ports
  addPort(m_mapServerPort);
  
  // </rtc-template>

  // <rtc-template block="bind_config">
  // Bind variables and configuration variable
  bindParameter("initial_pose_x", m_initial_pose_x, "0");
  bindParameter("initial_pose_y", m_initial_pose_y, "0");
  bindParameter("initial_pose_phi", m_initial_pose_phi, "0");
  bindParameter("initial_particle_min_x", m_initial_particle_min_x, "-0.3");
  bindParameter("initial_particle_max_x", m_initial_particle_max_x, "0.3");
  bindParameter("initial_particle_min_y", m_initial_particle_min_y, "-0.3");
  bindParameter("initial_particle_max_y", m_initial_particle_max_y, "0.3");
  bindParameter("initial_particle_min_phi", m_initial_particle_min_phi, "-0.3");
  bindParameter("initial_particle_max_phi", m_initial_particle_max_phi, "0.3");
  bindParameter("initial_particle_count", m_initial_particle_count, "1000");
  bindParameter("initial_particle_std_xy", m_initial_particle_std_xy, "0.3");
  bindParameter("initial_particle_std_phi", m_initial_particle_std_phi, "0.3");
  bindParameter("poseTimeOut", m_poseTimeOut, "3.0");
  bindParameter("motion_model", m_motion_model, "Gausian");
  bindParameter("motion_alpha1", m_motion_alpha1, "0.01");
  bindParameter("motion_alpha2", m_motion_alpha2, "0.05729");
  bindParameter("motion_alpha3", m_motion_alpha3, "0.01745");
  bindParameter("motion_alpha4", m_motion_alpha4, "0.05");
  bindParameter("motion_std_XY", m_motion_std_XY, "0.01");
  bindParameter("motion_std_PHI", m_motion_std_PHI, "0.003490");
  bindParameter("LM_likelihoodMethod", m_LM_likelihoodMethod, "lmLikelihoodField_Thrun");
  bindParameter("LM_enableLikelihoodCache", m_LM_enableLikelihoodCache, "true");
  bindParameter("LM_LF_decimation", m_LM_LF_decimation, "5");
  bindParameter("LM_LF_stdHit", m_LM_LF_stdHit, "0.35");
  bindParameter("LM_LF_zRandom", m_LM_LF_zRandom, "0.05");
  bindParameter("LM_LF_alternateAverageMethod", m_LM_LF_alternateAverageMethod, "false");
  bindParameter("LM_LF_zHit", m_LM_LF_zHit, "0.95");
  bindParameter("LM_LF_maxCorrsDistance", m_LM_LF_maxCorrsDistance, "0.3");
  bindParameter("LM_LF_maxRange", m_LM_LF_maxRange, "81");
  bindParameter("LM_MI_exponent", m_LM_MI_exponent, "2.5");
  bindParameter("LM_MI_ratio_max_distance", m_LM_MI_ratio_max_distance, "1.5");
  bindParameter("LM_MI_skip_rays", m_LM_MI_skip_rays, "10");
  bindParameter("LM_consensus_pow", m_LM_consensus_pow, "5");
  bindParameter("LM_consensus_takeEachRange", m_LM_consensus_takeEachRange, "1");
  bindParameter("LM_rayTracing_stdHit", m_LM_rayTracing_stdHit, "1.0");
  bindParameter("LM_rayTracing_decimation", m_LM_rayTracing_decimation, "10");
  bindParameter("LM_rayTracing_useDistanceFilter", m_LM_rayTracing_useDistanceFilter, "true");
  bindParameter("PF_algorithm", m_PF_algorithm, "pfStandardProposal");
  bindParameter("PF_resamplingMethod", m_PF_resamplingMethod, "prMultinomial");
  bindParameter("PF_BETA", m_PF_BETA, "0.5");
  bindParameter("PF_powFactor", m_PF_powFactor, "1.0");
  bindParameter("PF_sampleSize", m_PF_sampleSize, "1");
  bindParameter("PF_adaptiveSampleSize", m_PF_adaptiveSampleSize, "true");
  bindParameter("PF_max_loglikelihood_dyn_range", m_PF_max_loglikelihood_dyn_range, "15");
  bindParameter("PF_AuxFilterOptimal_MaximumSearchSamples", m_PF_AuxFilterOptimal_MaximumSearchSamples, "100");
  bindParameter("PF_AuxFilterStandard_FirstStageWeightsMonteCarlo", m_PF_AuxFilterStandard_FirstStageWeightsMonteCarlo, "false");
  bindParameter("PF_AuxFilterOptimal_MLE", m_PF_AuxFilterOptimal_MLE, "false");
  bindParameter("KLD_binSize_PHI", m_KLD_binSize_PHI, "0.1");
  bindParameter("KLD_binSize_XY", m_KLD_binSize_XY, "0.1");
  bindParameter("KLD_delta", m_KLD_delta, "0.02");
  bindParameter("KLD_epsilon", m_KLD_epsilon, "0.02");
  bindParameter("KLD_maxSampleSize", m_KLD_maxSampleSize, "1000");
  bindParameter("KLD_minSampleSize", m_KLD_minSampleSize, "150");
  bindParameter("KLD_minSamplesPerBin", m_KLD_minSamplesPerBin, "0");
  // </rtc-template>
  
  return RTC::RTC_OK;
}

/*
RTC::ReturnCode_t Localization_YMCL::onFinalize()
{
  return RTC::RTC_OK;
}
*/

/*
RTC::ReturnCode_t Localization_YMCL::onStartup(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/

/*
RTC::ReturnCode_t Localization_YMCL::onShutdown(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/


RTC::ReturnCode_t Localization_YMCL::onActivated(RTC::UniqueId ec_id)
{


	// Initialization
	std::cout << "[RTC::Localization_YMCL] onActivated called." << std::endl;
	m_pOGMap = new OGMap();
	while (m_mapServerPort.get_connector_profiles()->length() == 0) {
		coil::sleep(1);
		std::cout << "[RTC::Localization_YMCL] Waiting for Map Server Connection" << std::endl;
	}

	RTC::ConnectorProfileList& pList = *(m_mapServerPort.get_connector_profiles());
	RTC::RTObjectRef rto = (RTObjectRef)pList[0].ports[0]->get_port_profile()->owner;
	if (std::string((const char*)rto->get_component_profile()->instance_name) == this->getInstanceName()) {
		rto = (RTObjectRef)pList[0].ports[1]->get_port_profile()->owner;
	}
	do {
		coil::sleep(1);
		std::cout << "[RTC::Localization_YMCL] Waiting for Map Server Activation" << std::endl;
		if ((*(rto->get_owned_contexts()))[0]->get_component_state(rto) == RTC::ERROR_STATE) {
			std::cout << "[RTC::Localization_YMCL] Map Server RTC is now in ERROR_STATE" << std::endl;
			return RTC::RTC_ERROR;
		}
	} while ((*(rto->get_owned_contexts()))[0]->get_component_state(rto) != RTC::ACTIVE_STATE);

	RTC::RETURN_VALUE ret = m_mapServer->requestCurrentBuiltMap(m_pOGMap);

	if (ret != RTC::RETVAL_OK) {
		std::cout << "[RTC::Localization_YMCL] Acquiring Map from Server Failed." << std::endl;
	}

	std::cout << "[RTC::Localization_YMCL] Acquiring Map from Server Succeeded." << std::endl;


	std::cout << "[RTC::Localization_YMCL] Waiting For Odometry Input." << std::endl;
	while (!m_odometryIn.isNew()) {
	}
	m_odometryIn.read();
	std::cout << "[RTC::Localization_YMCL] Acquiring Odometry Information Succeeded." << std::endl;


	std::cout << "[RTC::Localization_YMCL] Waiting For Ranger Input." << std::endl;
	while (!m_rangeIn.isNew()) {
	}
	m_rangeIn.read();
	std::cout << "[RTC::Localization_YMCL] Acquiring Ranger Information Succeeded." << std::endl;


	std::cout << "[RTC::Localization_YMCL] Initializing Monte Carlo Localization." << std::endl;


	m_ymcl.param.map.pixel_width = m_pOGMap->config.width;
	m_ymcl.param.map.pixel_height = m_pOGMap->config.height;
	m_ymcl.param.map.resolution = m_pOGMap->config.xScale;
	m_ymcl.param.map.topleft_x = -m_pOGMap->config.origin.position.x;
	m_ymcl.param.map.topleft_y = -m_pOGMap->config.origin.position.y;

	m_ymcl.param.initial.sample_mode = INITIAL_SAMPLE_GAUSIAN;
	m_ymcl.param.initial.sample_size = m_initial_particle_count;

	m_ymcl.param.initial.pose_x = m_initial_pose_x;
	m_ymcl.param.initial.pose_y = m_initial_pose_y;
	m_ymcl.param.initial.pose_phi = m_initial_pose_phi;

	if (m_ymcl.param.initial.sample_mode == INITIAL_SAMPLE_GAUSIAN) {
		m_ymcl.param.initial.std_xy = m_initial_particle_max_x;
		m_ymcl.param.initial.std_phi = m_initial_particle_std_phi;
	}
	else {
		m_ymcl.param.initial.sample_x_min = m_initial_particle_min_x;
		m_ymcl.param.initial.sample_x_max = m_initial_particle_max_x;
		m_ymcl.param.initial.sample_y_min = m_initial_particle_min_y;
		m_ymcl.param.initial.sample_y_max = m_initial_particle_max_y;
		m_ymcl.param.initial.sample_phi_min = m_initial_particle_min_phi;
		m_ymcl.param.initial.sample_phi_max = m_initial_particle_std_phi;
	}

	m_ymcl.param.motion.alpha1 = m_motion_alpha1;
	m_ymcl.param.motion.alpha2 = m_motion_alpha2;
	m_ymcl.param.motion.alpha3 = m_motion_alpha3;
	m_ymcl.param.motion.alpha4 = m_motion_alpha4;
	m_ymcl.param.motion.update_distance = -1;
	m_ymcl.param.motion.update_heading = -1;


	m_ymcl.param.laser.model = LASER_LIKELIHOOD;
	if (m_ymcl.param.laser.model == LASER_LIKELIHOOD) {
		m_ymcl.param.laser.likelihood.zeta_hit = m_LM_LF_zHit;
		m_ymcl.param.laser.likelihood.sigma_hit = m_LM_LF_stdHit;
		m_ymcl.param.laser.likelihood.zeta_rand = m_LM_LF_zRandom;
		m_ymcl.param.laser.likelihood.scan_num = 20;
		m_ymcl.param.laser.likelihood.ranger_max_distance = 10;
	}
	else {
		m_ymcl.param.laser.beam.scan_num = 20;
		m_ymcl.param.laser.beam.zeta_hit = 0.95;
		m_ymcl.param.laser.beam.zeta_rand = 0.01;
		m_ymcl.param.laser.beam.zeta_max = 0.05;
		m_ymcl.param.laser.beam.zeta_short = 0.04;

		m_ymcl.param.laser.beam.sigma_max = 0.01;
		m_ymcl.param.laser.beam.sigma_hit = m_LM_rayTracing_stdHit;
		m_ymcl.param.laser.beam.lambda_short = 0.01;
	}

	m_ymcl.param.sampler.kld_sampling = coil::toBool(m_PF_adaptiveSampleSize, "true", "false", true);
	m_ymcl.param.sampler.kld_bin_size_xy = m_KLD_binSize_XY;
	m_ymcl.param.sampler.kld_bin_size_phi = m_KLD_binSize_PHI;
	m_ymcl.param.sampler.kld_epsilon = m_KLD_epsilon;
	m_ymcl.param.sampler.kld_delta = m_KLD_delta;
	m_ymcl.param.sampler.kld_max_particles = m_KLD_maxSampleSize;
	m_ymcl.param.sampler.kld_min_particles = m_KLD_minSampleSize;
	m_ymcl.param.sampler.random_alpha_fast = 0.8;
	m_ymcl.param.sampler.random_alpha_slow = 0.2;

	m_ymcl.param.sampler.random_mode = RANDOM_SAMPLING_NONE;
	if (m_ymcl.param.sampler.random_mode == RANDOM_SAMPLING_GAUSIAN) {
		m_ymcl.param.sampler.random_sample_std_xy = 0.5;
		m_ymcl.param.sampler.random_sample_std_phi = 0.2;
	}
	else if (m_ymcl.param.sampler.random_mode == RANDOM_SAMPLING_UNIFORM){
		m_ymcl.param.sampler.random_sample_x_min = -0.5;
		m_ymcl.param.sampler.random_sample_x_max = 0.5;
		m_ymcl.param.sampler.random_sample_y_min = -0.5;
		m_ymcl.param.sampler.random_sample_y_max = 0.5;
		m_ymcl.param.sampler.random_sample_phi_min = -0.2;
		m_ymcl.param.sampler.random_sample_phi_max = 0.2;
	}
	m_ymcl.param.sampler.sampling_method = SAMPLING_SYSTEMATIC;
	m_ymcl.param.sampler.resample_update_count = -1.0;
	m_ymcl.param.sampler.resample_distance = 0.3;
	m_ymcl.param.sampler.resample_heading = 0.2;

	m_ymcl.param.ranger.min_angle = m_range.config.minAngle;
	m_ymcl.param.ranger.max_angle = m_range.config.maxAngle;
	m_ymcl.param.ranger.min_distance = m_range.config.minRange;
	m_ymcl.param.ranger.max_distance = m_range.config.maxRange;
	m_ymcl.param.ranger.resolution = m_range.config.angularRes;
	m_ymcl.param.ranger.num_range = m_range.ranges.length();

	// Initialize YMCL object
	ymcl_init(&m_ymcl);


	// Set Map Data
	for (int i = 0; i < m_pOGMap->config.height; i++) {
		for (int j = 0; j < m_pOGMap->config.width; j++) {
			uint32_t cell = m_pOGMap->map.cells[i * m_pOGMap->config.width + j];
			map_cell_t value = cell < 100 ? MAP_OCCUPIED : (cell > 200 ? MAP_EMPTY : MAP_UNKNOWN);
			ymcl_set_map_pixel(&m_ymcl, j, i, value);
//			map.cell[i * m_pOGMap->config.width + j] = m_pOGMap->map.cells[i * m_pOGMap->config.width + j];
		}
	}
	// Set Ranger Data
	setRanger(&m_ymcl.ranger, m_range);

	pose_t pose;
	pose.x = m_odometry.data.position.x;
	pose.y = m_odometry.data.position.y;
	pose.th = normalize_angle(m_odometry.data.heading);
	ymcl_set_initial_pose(&m_ymcl, &pose);

	// Reset & initialize Particle Filter
	ymcl_reset(&m_ymcl);

	/// For Debug
	int windowWidth = 800;
	int windowHeight = 600;
	pPreviewImage = cvCreateImage(cvSize(m_pOGMap->config.width, m_pOGMap->config.height), IPL_DEPTH_8U, 3);
	pMapImage = cvCreateImage(cvSize(m_pOGMap->config.width, m_pOGMap->config.height), IPL_DEPTH_8U, 3);
	pWorkImage = cvCreateImage(cvSize(m_pOGMap->config.width, m_pOGMap->config.height), IPL_DEPTH_8U, 3);

	cvNamedWindow(this->getInstanceName(), CV_WINDOW_AUTOSIZE);

	/// Show map imiage 
	for (int h = 0; h < m_pOGMap->config.height; h++) {
		for (int w = 0; w < m_pOGMap->config.width; w++) {
			int index = h * m_pOGMap->config.width + w;
			int f = m_pOGMap->map.cells[index];
			//int f = 255 * m_ymcl.likelihood_field.distance_cell[index];
			pMapImage->imageData[index * 3 + 0] = f;
			pMapImage->imageData[index * 3 + 1] = f;
			pMapImage->imageData[index * 3 + 2] = f;
		}
	}
	pOGMap = m_pOGMap;
	instanceName = getInstanceName();

	m_lastResampledPose = m_odometry;

	return RTC::RTC_OK;
}



RTC::ReturnCode_t Localization_YMCL::onDeactivated(RTC::UniqueId ec_id)
{
	cvDestroyWindow(this->getInstanceName());
	cvReleaseImage(&pPreviewImage);
	cvReleaseImage(&pWorkImage);
	cvReleaseImage(&pMapImage);

	ymcl_cleanup(&m_ymcl);

	return RTC::RTC_OK;
}


CvPoint& positionToPixel(RTC::OGMap& map, double x, double y) {
	double zoomFactor = 1.0;
	return cvPoint((int)(zoomFactor* (x + map.config.origin.position.x) / map.config.xScale),
		-(int)(zoomFactor* (y + map.config.origin.position.y) / map.config.yScale));
}

void showPose(CvArr* img, RTC::OGMap& map, double x0, double y0, double th, CvScalar color) {
	double d = 0.3;

	CvPoint p0 = positionToPixel(map, x0, y0);

	double x1 = x0 + d * cos(th);
	double y1 = y0 + d * sin(th);
	CvPoint p1 = positionToPixel(map, x1, y1);

	double x2 = x1 + d / 2 * cos(th + M_PI * 3 / 4);
	double y2 = y1 + d / 2 * sin(th + M_PI * 3 / 4);
	CvPoint p2 = positionToPixel(map, x2, y2);

	double x3 = x1 + d / 2 * cos(th - M_PI * 3 / 4);
	double y3 = y1 + d / 2 * sin(th - M_PI * 3 / 4);
	CvPoint p3 = positionToPixel(map, x3, y3);

	cvLine(img, p0, p1, color);
	cvLine(img, p1, p2, color);
	cvLine(img, p2, p3, color);
}


void put_range() {


}

void show_range(const pose_t* pose, const map_t* map, const ranger_data_t* ranger) {

	cvCopy(pMapImage, pWorkImage);
	int j = 0;
	int step = ranger->param.num_range / 10;
	for (int i = 0; i < ranger->param.num_range; i++) {
		int x, y;
		//double distance = ray_casting(pose, map, ranger, ranger->min_angle + ranger->resolution * i, &x, &y);
		double th = pose->th + ranger->param.offset.y + ranger->param.min_angle + ranger->param.resolution * i;
		double sinth = sin(th);
		double costh = cos(th);
		//double a = sinth / costh;
		//int robot_cell_x = llroundf(robot_pose->x / map->resolution) + map->origin_index_x;
		//int robot_cell_y = -llroundf(robot_pose->y / map->resolution) + map->origin_index_y;

		if (ranger->ranges[i] == ranger->param.max_distance) {
			continue; // skip
		}

		pose_t z;
		//int x, y;
		z.x = pose->x + ranger->ranges[i] * cos(th);
		z.y = pose->y + ranger->ranges[i] * sin(th);
		pose_to_cell(&z, map, &x, &y);


		//real_t weight = map->distance_cell[y * map->pixel_width + x];
		real_t weight = 1.0;
		int r = 255 * weight;
		int b = 255 * (1 - weight);
		if (i == j) {
			cvRectangle(pWorkImage, cvPoint(x - 2, y - 2), cvPoint(x + 2, y + 2), CV_RGB(r, 0, b));
			j += step;
		}
		else {
			//cvRectangle(pWorkImage, cvPoint(x - 1, y - 1), cvPoint(x + 1, y + 1), CV_RGB(0, 255, 0));
		}
	}

	showPose(pWorkImage, *pOGMap, pose->x, pose->y, pose->th,
		CV_RGB(0, 255, 0));


	cvResize(pWorkImage, pPreviewImage);

	::cvShowImage(instanceName.c_str(), pPreviewImage);
	::cvWaitKey(1);
}


RTC::ReturnCode_t Localization_YMCL::onExecute(RTC::UniqueId ec_id)
{
	if (m_rangeIn.isNew()) {
		m_rangeIn.read();
		setRanger(&m_ymcl.ranger, m_range);
	}

	if (m_odometryIn.isNew()) {
		m_odometryIn.read();
		pose_t pose;
		pose.x = m_odometry.data.position.x;
		pose.y = m_odometry.data.position.y;
		pose.th = normalize_angle(m_odometry.data.heading);
		ymcl_push_odometry(&m_ymcl, &pose);
		ymcl_get_mean_pose(&m_ymcl, &pose);
		m_estimatedPose.data.position.x = pose.x;
		m_estimatedPose.data.position.y = pose.y;
		m_estimatedPose.data.heading = pose.th;
		setTimestamp<RTC::TimedPose2D>(m_estimatedPose);
		m_estimatedPoseOut.write();

		std::cout << "[[" << m_ymcl.particle_pool.min_weight << "/" << std::endl;

		cvCopy(pMapImage, pWorkImage);
		
		for (int i = 0; i < m_ymcl.particle_pool.num_particles; i++) {
			showPose(pWorkImage, *m_pOGMap, m_ymcl.particle_pool.particles[i].pose.x, m_ymcl.particle_pool.particles[i].pose.y, m_ymcl.particle_pool.particles[i].pose.th,
				CV_RGB(255, 0, 0));
			///CV_RGB((int)(255 * m_ymcl.particle_pool.particles[i].weight), 0, (int)(255 * (1 - m_ymcl.particle_pool.particles[i].weight))));
		}


		showPose(pWorkImage, *m_pOGMap, pose.x,
			pose.y,
			pose.th,
			CV_RGB(0, 0, 255));

		cvResize(pWorkImage, pPreviewImage);

		::cvShowImage(getInstanceName(), pPreviewImage);
		::cvWaitKey(1);

	}

	return RTC::RTC_OK;
}

/*
RTC::ReturnCode_t Localization_YMCL::onAborting(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/

/*
RTC::ReturnCode_t Localization_YMCL::onError(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/

/*
RTC::ReturnCode_t Localization_YMCL::onReset(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/

/*
RTC::ReturnCode_t Localization_YMCL::onStateUpdate(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/

/*
RTC::ReturnCode_t Localization_YMCL::onRateChanged(RTC::UniqueId ec_id)
{
  return RTC::RTC_OK;
}
*/



extern "C"
{
 
  void Localization_YMCLInit(RTC::Manager* manager)
  {
    coil::Properties profile(localization_ymcl_spec);
    manager->registerFactory(profile,
                             RTC::Create<Localization_YMCL>,
                             RTC::Delete<Localization_YMCL>);
  }
  
};


