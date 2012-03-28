/*
 * pset3.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <omp.h>
#include "log.h"
#include "Timer.h"
#include "Options.h"
#include "BatchBinSet.h"
#include "Isotope.h"
#include "Material.h"
#include "Neutron.h"
#include "gnuplot.h"
#include "Region1D.h"
#include "Fissioner.h"


//TODO: variance for k-inf
//TODO: output k-inf as computed by two group xs

int main(int argc, const char **argv) {

	Options options(argc, argv);
	Timer timer;

	log_setlevel(options.getVerbosity());

	/* Get the number of neutrons, bins and batches */
	int num_neutrons = options.getNumNeutrons();
	int num_bins = options.getNumBins();
	int num_batches = options.getNumBatches();
    int num_threads = options.getNumThreads();
	int num_gen;
	int num_alive;

	log_printf(NORMAL, "Beginning two region problem with %d neutrons, "
			"%d bins, %d batches, %d threads...", num_neutrons, num_bins,
			num_batches, num_threads);

	/* Create a handle for plotting with gnuplot */
	gnuplot_ctrl* handle;


	/* Create a set of plotting flux bins for each batch */
	BatchBinSet* total_flux = new BatchBinSet();
	BatchBinSet* fuel_flux = new BatchBinSet();
	BatchBinSet* moderator_flux = new BatchBinSet();

	total_flux->createBinners(1E-6, 1E7, num_bins, num_batches,
							LOGARITHMIC, FLUX_ENERGY, (char*)"all");
	fuel_flux->createBinners(1E-6, 1E7, num_bins, num_batches,
							LOGARITHMIC, FLUX_ENERGY, (char*)"all");
	moderator_flux->createBinners(1E-6, 1E7, num_bins, num_batches,
							LOGARITHMIC, FLUX_ENERGY, (char*)"all");


	/* Create bins to compute total fission and absorption rates */
	BatchBinSet* tot_fiss_rate = new BatchBinSet();
	BatchBinSet* tot_abs_rate = new BatchBinSet();

	tot_fiss_rate->createBinners(1E-7, 1E7, 1, num_batches, EQUAL,
									FISSION_RATE_ENERGY, (char*)"all");
	tot_abs_rate->createBinners(1E-7, 1E7, 1, num_batches, EQUAL,
									ABSORPTION_RATE_ENERGY, (char*)"all");
	float nu_bar = 2.455;	/* CASMO edit for average # neutrons per fission */


	/* Create bins to compute two group cross-sections */
	BatchBinSet* H1_capture_rate_2G = new BatchBinSet();
	BatchBinSet* H1_elastic_rate_2G = new BatchBinSet();
	BatchBinSet* O16_elastic_rate_2G = new BatchBinSet();
	BatchBinSet* ZR90_elastic_rate_2G = new BatchBinSet();

	BatchBinSet* U235_capture_rate_2G = new BatchBinSet();
	BatchBinSet* U235_elastic_rate_2G = new BatchBinSet();
	BatchBinSet* U235_fission_rate_2G = new BatchBinSet();
	BatchBinSet* U238_capture_rate_2G = new BatchBinSet();
	BatchBinSet* U238_elastic_rate_2G = new BatchBinSet();
	BatchBinSet* U238_fission_rate_2G = new BatchBinSet();

	BatchBinSet* two_group_flux = new BatchBinSet();

	float two_group_E_ranges[3] = {0.0, 0.625, 1E7};

	H1_capture_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										CAPTURE_RATE_ENERGY, (char*)"H1");
	H1_elastic_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										ELASTIC_RATE_ENERGY, (char*)"H1");
	O16_elastic_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										ELASTIC_RATE_ENERGY, (char*)"O16");
	ZR90_elastic_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										ELASTIC_RATE_ENERGY, (char*)"ZR90");

	U235_capture_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										CAPTURE_RATE_ENERGY, (char*)"U235");
	U235_elastic_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										ELASTIC_RATE_ENERGY, (char*)"U235");
	U235_fission_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										FISSION_RATE_ENERGY, (char*)"U235");
	U238_capture_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										CAPTURE_RATE_ENERGY, (char*)"U238");
	U238_elastic_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										ELASTIC_RATE_ENERGY, (char*)"U238");
	U238_fission_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
										FISSION_RATE_ENERGY, (char*)"U238");

	two_group_flux->createBinners(two_group_E_ranges, 2, num_batches,
									FLUX_ENERGY, (char*)"all");


	/* Create bins to compute moderator to fuel flux ratios */
	int num_ratios = 13;
	BatchBinSet* fuel_flux_ratio = new BatchBinSet();
	BatchBinSet* moderator_flux_ratio = new BatchBinSet();

	float flux_ratio_E_ranges[14] = {0.0, 0.1, 0.5, 1.0, 6.0, 10.0, 25.0,
									50.0, 100.0, 1000.0, 10000.0, 100000.0,
									500000.0, 10000000.0};

	fuel_flux_ratio->createBinners(flux_ratio_E_ranges, num_ratios,
							num_batches, FLUX_ENERGY, (char*)"all");

	moderator_flux_ratio->createBinners(flux_ratio_E_ranges, num_ratios,
							num_batches, FLUX_ENERGY, (char*)"all");


	/* Create bins to compute the diffusion coefficient for three methods */
	BatchBinSet* coll_rate_2G = new BatchBinSet();
	BatchBinSet* transport_rate_2G = new BatchBinSet();
	BatchBinSet* diffusion_rate_2G = new BatchBinSet();

	coll_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
								COLLISION_RATE_ENERGY, (char*)"all");
	transport_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
								TRANSPORT_RATE_ENERGY, (char*)"all");
	diffusion_rate_2G->createBinners(two_group_E_ranges, 2, num_batches,
								DIFFUSION_RATE_ENERGY, (char*)"all");


	/* 2-region pin cell geometric parameters (units in cm) */
	float r_fuel = 0.4096;
	float r_gap = 0.4178;
	float r_cladding = 0.4750;
	float pitch = 1.26;
	float p2 = pitch * pitch;

	/* 2-region homogenized densities (g/cm^3) and enrichment */
	float rho_fuel = 10.2;
	float rho_cladding = 6.549;
	float rho_coolant = 0.9966;
	float enrichment = 0.03035;

	/* Isotope number densities */
	float N_A = 6.023E23;	/* Avogadro's number (at / mol) */
	float N_U238 = rho_fuel*N_A*(1.0 - enrichment) / ((238.0 *
					(1.0 - enrichment)) + (235.0*enrichment) + (16.0*2.0));
	float N_U235 = rho_fuel*N_A*enrichment / ((238.0 *
					(1.0 - enrichment)) + (235.0*enrichment) + (16.0*2.0));
	float N_O16 = rho_fuel*N_A*2.0 / ((238.0 *
					(1.0 - enrichment)) + (235.0*enrichment) + (16.0*2.0));
	float N_ZR90 = rho_cladding*N_A / 90.0;
	float N_H2O = rho_coolant*N_A / 18.0;
	float N_H1 = rho_coolant*N_A*2.0 / 18.0;

	/* 2-region pin cell volumes (cm^3) */
	float v_fuel = M_PI*r_fuel*r_fuel;
	float v_gap = M_PI*(r_gap*r_gap - r_fuel*r_fuel);
	float v_cladding = M_PI*(r_cladding*r_cladding - r_gap*r_gap);
	float v_coolant = p2 - M_PI*r_cladding*r_cladding;
	float v_moderator = v_gap + v_cladding + v_coolant;
	float v_total = v_fuel + v_moderator;

	/* Compute homogenized moderator number densities using volume weighting */
	N_H2O *= (v_coolant / v_moderator);
	N_H1 *= (v_coolant / v_moderator);
	N_ZR90 *= (v_cladding / v_moderator);

	/* Dancoff factor from CASMO-5 */
	float dancoff = 0.277;

	/* Escape cross-section */
	float sigma_e = 1.0 / (2.0*r_fuel);

	/* Carlvik's two-term rational model */
	float A = (1.0 - dancoff) / dancoff;
	float alpha1 = ((5.0*A + 6.0) - sqrt(A*A + 36.0*A + 36.0)) /
														(2.0*(A+1.0));
	float alpha2 = ((5.0*A + 6.0) + sqrt(A*A + 36.0*A + 36.0)) /
														(2.0*(A+1.0));
	float beta = (((4.0*A + 6.0) / (A + 1.0)) - alpha1) / (alpha2 - alpha1);

	/* Print out the geometry parameters */
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "\t\t\t\tGeometry Parameters (cm)");
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "");
	log_printf(NORMAL, "r_fuel = %f", r_fuel);
	log_printf(NORMAL, "r_gap  = %f", r_gap);
	log_printf(NORMAL, "r_cladding = %f", r_cladding);
	log_printf(NORMAL, "pitch = %f", pitch);
	log_printf(NORMAL, "total cell area = %f", p2);
	log_printf(NORMAL, "v_fuel = %f", v_fuel);
	log_printf(NORMAL, "v_gap = %f", v_gap);
	log_printf(NORMAL, "v_cladding = %f", v_cladding);
	log_printf(NORMAL, "v_coolant = %f", v_coolant);
	log_printf(NORMAL, "v_moderator = %f", v_moderator);
	log_printf(NORMAL, "v_total = %f", v_total);
	log_printf(NORMAL, "");

	/* Print to the console the number densities */
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "\t\t\t\tNumber Densities (at/cm^3)");
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "");
	log_printf(NORMAL, "H1:\t%1.5e", N_H1);
	log_printf(NORMAL, "H2O:\t%1.5e", N_H2O);
	log_printf(NORMAL, "ZR90:\t%1.5e", N_ZR90);
	log_printf(NORMAL, "U235:\t%1.5e", N_U235);
	log_printf(NORMAL, "U238:\t%1.5e", N_U238);
	log_printf(NORMAL, "O16:\t%1.5e", N_O16);
	log_printf(NORMAL, "");

	/* Print to the console the collision probability parameters */
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "\t\t\tTwo Region Collision Probability Parameters");
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "");
	log_printf(NORMAL, "dancoff = %f", dancoff);
	log_printf(NORMAL, "sigma_e = %f", sigma_e);
	log_printf(NORMAL, "A = %f", A);
	log_printf(NORMAL, "alpha1 = %f", alpha1);
	log_printf(NORMAL, "alpha2 = %f", alpha2);
	log_printf(NORMAL, "beta = %f", beta);
	log_printf(NORMAL, "");


	/* Create isotopes*/
	char* delim = (char*)"\t";

	Isotope* H1 = new Isotope();
	H1->setA(1);
	H1->setIsotopeType((char*)"H1");
	H1->loadXS((char*)"pendf/h-1_capture.txt", CAPTURE, delim);
	H1->loadXS((char*)"pendf/h-1_elastic.txt", ELASTIC, delim);
	H1->setElasticAngleType(ISOTROPIC_LAB);
	H1->initializeThermalScattering(1E-6, 15, 1000, 15);

	Isotope* O16 = new Isotope();
	O16->setA(16);
	O16->setIsotopeType((char*)"O16");
	O16->loadXS((char*)"pendf/o-16_elastic.txt", ELASTIC, delim);
	O16->setElasticAngleType(ISOTROPIC_LAB);

	Isotope* ZR90 = new Isotope();
	ZR90->setA(90);
	ZR90->setIsotopeType((char*)"ZR90");
	ZR90->loadXS((char*)"pendf/zr-90_elastic.txt", ELASTIC, delim);
	ZR90->setElasticAngleType(ISOTROPIC_LAB);

	Isotope* U235 = new Isotope();
	U235->setA(235);
	U235->setIsotopeType((char*)"U235");
	U235->loadXS((char*)"pendf/u-235_capture.txt", CAPTURE, delim);
	U235->setOneGroupElasticXS(11.4, ISOTROPIC_LAB);
	U235->loadXS((char*)"pendf/u-235_fission.txt", FISSION, delim);

	Isotope* U238 = new Isotope();
	U238->setA(238);
	U238->setIsotopeType((char*)"U238");
	U238->loadXS((char*)"pendf/u-238_capture.txt", CAPTURE, delim);
	U238->setOneGroupElasticXS(11.3, ISOTROPIC_LAB);
	U238->loadXS((char*)"pendf/u-238_fission.txt", FISSION, delim);


	/* Create Materials */
	Material* moderator = new Material[num_threads];
	Material* fuel = new Material[num_threads];


	/* Create Regions for each thread */
	Region1D* pellet = new Region1D[num_threads];
	Region1D* coolant = new Region1D[num_threads];

	/* Create Fissioners for each thread */
	Fissioner* fissioners = new Fissioner[num_threads];

	/* Create Region class objects for each thread */
	for (int i=0; i < num_threads; i++) {

		/* Initialize Materials for each thread with isotope clones */
		moderator[i].setMaterialName((char*)"moderator");
		fuel[i].setMaterialName((char*)"fuel");

		moderator[i].addIsotope(ZR90->clone(), N_ZR90);
		moderator[i].addIsotope(H1->clone(), N_H1);
		moderator[i].addIsotope(O16->clone(), N_H2O);
		moderator[i].rescaleCrossSections(1E-7, 1E7, 50000, LOGARITHMIC);

		fuel[i].addIsotope(U235->clone(), N_U235);
		fuel[i].addIsotope(U238->clone(), N_U238);
		fuel[i].addIsotope(O16->clone(), N_O16);
		fuel[i].rescaleCrossSections(1E-7, 1E7, 50000, LOGARITHMIC);

		/* Set the two region collision probability parameters */
		pellet[i].setRegionName((char*)"pellet");
		pellet[i].setMaterial(&fuel[i]);
		pellet[i].setAsFuel();
		pellet[i].setOtherPinCellRegion(&coolant[i]);
		pellet[i].setVolume(v_fuel);
		pellet[i].setTwoRegionPinCellParams(sigma_e, beta, alpha1, alpha2);

		coolant[i].setRegionName((char*)"coolant");
		coolant[i].setMaterial(&moderator[i]);
		coolant[i].setAsModerator();
		coolant[i].setOtherPinCellRegion(&pellet[i]);
		coolant[i].setVolume(v_moderator);
		coolant[i].setTwoRegionPinCellParams(sigma_e, beta, alpha1, alpha2);

		/* Set the fissioner class for this thread to have 10MeV maximum and
		 * 5000 sample bins */
		fissioners[i].setEMax(10.0);
		fissioners[i].setNumBins(200);
		fissioners[i].buildCDF();
	}


	/* Run the simulation */
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "\t\t\t\tBeginning Simulation...");
	log_printf(NORMAL, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "");

	timer.start();

	omp_set_num_threads(num_threads);
	#pragma omp parallel shared(total_flux, fuel_flux, moderator_flux,\
							fuel_flux_ratio, moderator_flux_ratio,\
							tot_fiss_rate, tot_abs_rate, U235_capture_rate_2G,\
							U235_elastic_rate_2G, U235_fission_rate_2G,\
							U238_capture_rate_2G, U238_elastic_rate_2G,\
							U238_fission_rate_2G, H1_capture_rate_2G,\
							H1_elastic_rate_2G, O16_elastic_rate_2G,\
							ZR90_elastic_rate_2G, fuel, moderator, \
							pellet, coolant, fissioners)
	{
		/* Loop over batches */
		#pragma omp for private(num_gen, num_alive)
		for (int b=0; b < num_batches; b++) {

			int thread_num = omp_get_thread_num();
			log_printf(NORMAL, "Batch: %d\tThread: %d", b, thread_num);

			/* Set the binns for this batch */
			pellet[thread_num].clearBinners();
			pellet[thread_num].addBinner(total_flux->getBinner(b));
			pellet[thread_num].addBinner(fuel_flux->getBinner(b));
			pellet[thread_num].addBinner(fuel_flux_ratio->getBinner(b));
			pellet[thread_num].addBinner(tot_fiss_rate->getBinner(b));
			pellet[thread_num].addBinner(tot_abs_rate->getBinner(b));
			pellet[thread_num].addBinner(U235_capture_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(U235_elastic_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(U235_fission_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(U238_capture_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(U238_elastic_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(U238_fission_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(O16_elastic_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(two_group_flux->getBinner(b));
			pellet[thread_num].addBinner(coll_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(transport_rate_2G->getBinner(b));
			pellet[thread_num].addBinner(diffusion_rate_2G->getBinner(b));

			coolant[thread_num].clearBinners();
			coolant[thread_num].addBinner(total_flux->getBinner(b));
			coolant[thread_num].addBinner(moderator_flux->getBinner(b));
			coolant[thread_num].addBinner(moderator_flux_ratio->getBinner(b));
			coolant[thread_num].addBinner(tot_fiss_rate->getBinner(b));
			coolant[thread_num].addBinner(tot_abs_rate->getBinner(b));
			coolant[thread_num].addBinner(H1_capture_rate_2G->getBinner(b));
			coolant[thread_num].addBinner(H1_elastic_rate_2G->getBinner(b));
			coolant[thread_num].addBinner(O16_elastic_rate_2G->getBinner(b));
			coolant[thread_num].addBinner(ZR90_elastic_rate_2G->getBinner(b));
			coolant[thread_num].addBinner(two_group_flux->getBinner(b));
			coolant[thread_num].addBinner(coll_rate_2G->getBinner(b));
			coolant[thread_num].addBinner(transport_rate_2G->getBinner(b));
			coolant[thread_num].addBinner(diffusion_rate_2G->getBinner(b));

			/* Initialize all neutrons for this batch and add them to slab 1 */
			for (int n=0; n < num_neutrons; n++) {
				neutron* new_neutron = initializeNewNeutron();
				new_neutron->_x = 0.0;
				new_neutron->_mu = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
				new_neutron->_energy = fissioners[thread_num].emitNeutroneV();
				pellet[thread_num].addNeutron(new_neutron);
			}

			/* Loop over all neutrons until they are all dead */
			num_gen = 1;
			num_alive = num_neutrons;

			while (num_alive > 0) {

				log_printf(DEBUG, "batch = %d, thread = %d, gen = %d, "
						"num_alive = %d", b, thread_num, num_gen, num_alive);

				num_gen++;
				num_alive = 0;

				/* Transfer neutrons between regions based on
				 * two region collision probabilities */
				pellet[thread_num].twoRegionNeutronTransferral();
				coolant[thread_num].twoRegionNeutronTransferral();

				/* Update each region's vector of neutrons with those
				 * neutrons which were just transferred */
				pellet[thread_num].initializeTransferredNeutrons();
				coolant[thread_num].initializeTransferredNeutrons();

				/* Move neutrons within each region */
				pellet[thread_num].moveNeutrons();
				coolant[thread_num].moveNeutrons();

				num_alive = pellet[thread_num].getNumNeutrons() +
							coolant[thread_num].getNumNeutrons();
			}
		}
	}

	log_printf(NORMAL, "");

	/* Stop the timer record the timing split for this simulation */
	timer.stop();
	timer.recordSplit("Pset 4 time (sec)");

	/* Compute batch statistics for total flux and flux in fuel, moderator */
	total_flux->computeScaledBatchStatistics(num_neutrons*v_total);
	fuel_flux->computeScaledBatchStatistics(num_neutrons*v_fuel);
	moderator_flux->computeScaledBatchStatistics(num_neutrons*v_moderator);

	/* Compute batch statistics for total fission and absorption rates */
	tot_fiss_rate->computeScaledBatchStatistics(num_neutrons*v_total);
	tot_abs_rate->computeScaledBatchStatistics(num_neutrons*v_total);

	/* Compute batch statistics for one group cross-sections */
	H1_capture_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	H1_elastic_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	O16_elastic_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	ZR90_elastic_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	U235_capture_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	U235_elastic_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	U235_fission_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	U238_capture_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	U238_elastic_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	U238_fission_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	two_group_flux->computeScaledBatchStatistics(num_neutrons*v_total);
	coll_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	transport_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);
	diffusion_rate_2G->computeScaledBatchStatistics(num_neutrons*v_total);

	/* Compute k-infinity */
	float fiss_rate_mu = tot_fiss_rate->getBatchMu()[0];
	float fiss_rate_var = tot_fiss_rate->getBatchVariance()[0];
	float abs_rate_mu = tot_abs_rate->getBatchMu()[0];
	float abs_rate_var = tot_abs_rate->getBatchVariance()[0];

	float k_inf = fiss_rate_mu * nu_bar / abs_rate_mu;

	float k_inf_var = (fiss_rate_mu*fiss_rate_mu)*abs_rate_var +
						(abs_rate_mu*abs_rate_mu)*fiss_rate_var +
									fiss_rate_var*abs_rate_var;

	float k_inf_std_dev = sqrt(k_inf_var);

	/* Compute moderator to fuel flux ratios */
	fuel_flux_ratio->computeScaledBatchStatistics(num_neutrons*v_fuel);
	moderator_flux_ratio->computeScaledBatchStatistics(num_neutrons*v_moderator);
	fuel_flux_ratio->computeScaledBatchStatistics(num_neutrons*v_fuel);
	moderator_flux_ratio->computeScaledBatchStatistics(num_neutrons*v_moderator);

	/* Print to the console the total fission rate */
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "\t\t\tTotal Fission Rate (Batch Statistics)");
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(NORMAL, "");
	log_printf(RESULT, "Tot fission rate = %1.8f\t\tVariance = %1.8f",
									tot_fiss_rate->getBatchMu()[0],
									tot_fiss_rate->getBatchVariance()[0]);
	log_printf(RESULT, "Tot absorption rate = %f\t\tVariance = %f",
									tot_abs_rate->getBatchMu()[0],
									tot_abs_rate->getBatchVariance()[0]);
	log_printf(RESULT, "k_inf = %f\t\tvariance = %1.8f \t\t 2 sigma = %1.8f", k_inf,
													k_inf_var, k_inf_std_dev);
	log_printf(RESULT, "");

	/* Print to the console the moderator/fuel flux ratios */
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "\t\t\t\tModerator/Fuel Flux Ratios");
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "");

	float ratio;
	for (int i=1; i < num_ratios+1; i++) {

		ratio = moderator_flux_ratio->getBatchMu()[i-1] /
						fuel_flux_ratio->getBatchMu()[i-1];

		log_printf(RESULT, "[%2.e eV - %2.e eV]:\t%f",
				flux_ratio_E_ranges[i-1], flux_ratio_E_ranges[i], ratio);
	}

	log_printf(RESULT, "");


	/* Print to the console the cell-averaged fast to thermal flux ratio */
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "\t\t\tCell-Averaged Fast-to-Thermal Flux Ratio");
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "");
	double* two_group_flux_mu = two_group_flux->getBatchMu();
	double flux1 = two_group_flux_mu[0];
	double flux2 = two_group_flux_mu[1];
	log_printf(RESULT, "Ratio = %f", flux2 / flux1);
	log_printf(RESULT, "");


	/* Print to the console the two group macroscopic cross-sections */
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "\t\t\tTwo Group Macroscopic Cross-Sections (cm^-1)");
	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "");

	float xs1, xs2;

	log_printf(RESULT, "\t\t\t[%1.1f eV - %1.3f eV]\t[%1.3f eV - %1.1e eV]",
						two_group_flux->getBinner(0)->getBinEdges()[0],
						two_group_flux->getBinner(0)->getBinEdges()[1],
						two_group_flux->getBinner(0)->getBinEdges()[1],
						two_group_flux->getBinner(0)->getBinEdges()[2]);

	/* H1 capture */
	xs1 = H1_capture_rate_2G->getBatchMu()[0] / flux1;
	xs2 = H1_capture_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "H1 Capture: \t\t%f\t\t%f", xs1, xs2);

	/* H1 elastic */
	xs1 = H1_elastic_rate_2G->getBatchMu()[0] / flux1;
	xs2 = H1_elastic_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "H1 Elastic: \t\t%f\t\t%f", xs1, xs2);

	/* O16 elastic */
	xs1 = O16_elastic_rate_2G->getBatchMu()[0] / flux1;
	xs2 = O16_elastic_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "O16 Elastic: \t\t%f\t\t%f", xs1, xs2);

	/* ZR90 elastic */
	xs1 = ZR90_elastic_rate_2G->getBatchMu()[0] / flux1;
	xs2 = ZR90_elastic_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "ZR90 Elastic: \t\t%f\t\t%f", xs1, xs2);

	/* U235 capture */
	xs1 = U235_capture_rate_2G->getBatchMu()[0] / flux1;
	xs2 = U235_capture_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "U235 Capture: \t\t%f\t\t%f", xs1, xs2);

	/* U235 elastic */
	xs1 = U235_elastic_rate_2G->getBatchMu()[0] / flux1;
	xs2 = U235_elastic_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "U235 Elastic: \t\t%f\t\t%f", xs1, xs2);

	/* U235 fission */
	xs1 = U235_fission_rate_2G->getBatchMu()[0] / flux1;
	xs2 = U235_fission_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "U235 Fission: \t\t%f\t\t%f", xs1, xs2);

	/* U238 capture */
	xs1 = U238_capture_rate_2G->getBatchMu()[0] / flux1;
	xs2 = U238_capture_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "U238 Capture: \t\t%f\t\t%f", xs1, xs2);

	/* U238 elastic */
	xs1 = U238_elastic_rate_2G->getBatchMu()[0] / flux1;
	xs2 = U238_elastic_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "U238 Elastic: \t\t%f\t\t%f", xs1, xs2);

	/* U238 fission */
	xs1 = U238_fission_rate_2G->getBatchMu()[0] / flux1;
	xs2 = U238_fission_rate_2G->getBatchMu()[1] / flux2;
	log_printf(RESULT, "U238 Fission: \t\t%f\t\t%f", xs1, xs2);

	log_printf(RESULT, "");


	/* Print to the console the two group macroscopic cross-sections */
	log_printf(RESULT, "*******************************************************"
												"*************************");
	log_printf(RESULT, "\t\t\tTwo Group Diffusion Coefficients");
	log_printf(RESULT, "*******************************************************"
												"*************************");
	log_printf(RESULT, "");
	log_printf(RESULT, "\t\t\t[%1.1f eV - %1.3f eV]\t[%1.3f eV - %1.1e eV]",
						two_group_flux->getBinner(0)->getBinEdges()[0],
						two_group_flux->getBinner(0)->getBinEdges()[1],
						two_group_flux->getBinner(0)->getBinEdges()[1],
						two_group_flux->getBinner(0)->getBinEdges()[2]);

	float sigma_t1, sigma_t2;
	float sigma_tr1, sigma_tr2;
	float D1, D2;

	sigma_t1 = coll_rate_2G->getBatchMu()[0] / flux1;
	sigma_t2 = coll_rate_2G->getBatchMu()[1] / flux2;
	D1 = 1.0 / (3.0 * sigma_t1);
	D2 = 1.0 / (3.0 * sigma_t2);

	log_printf(RESULT, "1/(3*sigma_t):\t\t\t%f\t\t%f", D1, D2);

	sigma_tr1 = transport_rate_2G->getBatchMu()[0] / flux1;
	sigma_tr2  = transport_rate_2G->getBatchMu()[1] / flux2;
	D1 = 1.0 / (3.0 * sigma_tr1);
	D2 = 1.0 / (3.0 * sigma_tr2);

	log_printf(RESULT, "1/(3*sigma_tr):\t\t\t%f\t\t%f", D1, D2);

	D1 = diffusion_rate_2G->getBatchMu()[0] / flux1;
	D2 = diffusion_rate_2G->getBatchMu()[1] / flux2;

	log_printf(RESULT, "Diff coeff:\t\t\t%f\t\t%f", D1, D2);

	log_printf(RESULT, "");


	/* Plot the total neutron flux */
	handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Energy (eV)");
	gnuplot_set_ylabel(handle, (char*)"flux");
	gnuplot_set_xrange(handle, 0.005, 1E7);
	gnuplot_cmd(handle, (char*)"set logscale xy");
	gnuplot_cmd(handle, (char*)"set title \"Normalized Flux\"");
	gnuplot_setstyle(handle, (char*)"lines");
	gnuplot_plot_xy(handle, total_flux->getBinner(0)->getBinCenters(),
			total_flux->getBatchMu(), num_bins, (char*)"Total Flux");
	gnuplot_plot_xy(handle, fuel_flux->getBinner(0)->getBinCenters(),
			fuel_flux->getBatchMu(), num_bins, (char*)"Fuel Flux");
	gnuplot_saveplot(handle, (char*)"flux");
	gnuplot_plot_xy(handle, moderator_flux->getBinner(0)->getBinCenters(),
			moderator_flux->getBatchMu(), num_bins, (char*)"Moderator Flux");
	gnuplot_close(handle);


	/* Free all allocated memory */
	delete [] pellet;
	delete [] coolant;
	delete [] fissioners;

	delete [] moderator;
	delete [] fuel;

	delete total_flux;
	delete fuel_flux;
	delete moderator_flux;
	delete tot_fiss_rate;
	delete tot_abs_rate;
	delete fuel_flux_ratio;
	delete moderator_flux_ratio;
	delete H1_capture_rate_2G;
	delete H1_elastic_rate_2G;
	delete O16_elastic_rate_2G;
	delete U235_capture_rate_2G;
	delete U235_elastic_rate_2G;
	delete U235_fission_rate_2G;
	delete U238_capture_rate_2G;
	delete U238_elastic_rate_2G;
	delete U238_fission_rate_2G;
	delete ZR90_elastic_rate_2G;
	delete two_group_flux;
	delete coll_rate_2G;
	delete transport_rate_2G;
	delete diffusion_rate_2G;

	delete H1;
	delete O16;
	delete ZR90;
	delete U235;
	delete U238;

	log_printf(RESULT, "*******************************************************"
													"*************************");
	log_printf(RESULT, "\t\t\t\tTiming Results");
	log_printf(RESULT, "*******************************************************"
													"*************************");
	timer.printSplits();
}
