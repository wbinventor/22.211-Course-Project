/*
 * Region1D.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Region1D.h"


/**
 * Region1D constructor sets default values
 */
Region1D::Region1D() {

	/* By default the Region1D is infinite (unbounded) */
	_region_type = INFINITE;
	_region_name = (char*)"";
	_left_surface = NULL;
	_right_surface = NULL;

	/* Default volume */
	_volume = 0.0;

	/* Default two region pin cell parameters */
	_fuel = false;
	_moderator = false;
	_other_region = NULL;
	_sigma_e = 0.0;
	_beta = 0.0;
	_alpha1 = 0.0;
	_alpha2 = 0.0;

	/* By default region does not use variance reduction */
	_use_implicit_capture = false;
	_use_forced_collision = false;
	_weight_low = 0.0;
	_weight_avg = 0.0;
}


/**
 * Region1D destructor
 */
Region1D::~Region1D() { }


/**
 * Returns a pointer to the left surface bounding this region
 * @return pointer to the left surface
 */
Surface* Region1D::getLeftSurface() const {
	return _left_surface;
}


/**
 * Returns a pointer to the right surface bounding this region
 * @return pointer to the right surface
 */
Surface* Region1D::getRightSurface() const {
	return _right_surface;
}


/**
 * Returns the volume of this Region1D (cm^3)
 * @return the Region1D's volume
 */
float Region1D::getVolume() {
	return _volume;
}


/**
 * Returns a pointer to the Material filling this Region1D
 * @return a pointer to a Material
 */
Material* Region1D::getMaterial() {
	return _material;
}


/**
 * Return the type of Region1D (BOUNDED or INFINITE)
 * @return the region type
 */
regionType Region1D::getRegionType() {
	return _region_type;
}


/**
 * Return the name, if any, of this region as specified by
 * the user
 * @return a character array representing this region's name
 */
char* Region1D::getRegionName() {
	return _region_name;
}


/**
 * Returns true if this region has been labeled as fuel, false otherwise
 * @return true if fuel, false otherwise
 */
bool Region1D::isFuel() {
	return _fuel;
}


/**
 * Returns true if this region has been labeled as moderator, false otherwise
 * @return true if moderator, false otherwise
 */
bool Region1D::isModerator() {
	return _moderator;
}


/**
 * Returns the other Region1D in the 2 region pin cell model if this Region1D
 * is being used for a simulation employing that model
 * @return a pointer to the other Region1D
 */
Region1D* Region1D::getOtherPinCellRegion() {

	if (_other_region == NULL)
		log_printf(ERROR, "Unable to return %s's other region since it has "
				"not been set yet", _region_name);

	return _other_region;
}


/**
 * Sets this Region1D's name as specified by a user
 * @param region_name the name of this Region1D
 */
void Region1D::setRegionName(char* region_name) {
	_region_name = region_name;
}


/**
 * Set the Material filling this Region1D
 * @param material a pointer to a Material
 */
void Region1D::setMaterial(Material* material) {
	_material = material;
}


/**
 * Assigns this region's left surface
 * @param a pointer to a Surface
 */
void Region1D::setLeftSurface(Surface* surface) {
	_left_surface = surface;
	_region_type = BOUNDED;
}


/**
 * Assigns this region's right surface
 * @param a pointer to a Surface
 */
void Region1D::setRightSurface(Surface* surface) {
	_right_surface = surface;
	_region_type = BOUNDED;
}


/**
 * Sets this Region1D's volume (cm^3)
 * @param volume the volume of this Region1D
 */
void Region1D::setVolume(float volume) {
	_volume = volume;
}


/**
 * Adds a new Binner to this region for tallying
 * @param bins a pointer to a Binner class object
 */
void Region1D::addBinner(Binner* bins) {
	_bin_sets.push_back(bins);
}


/**
 * Sets this Region1D to be the fuel region in a two region pin
 * cell simulation
 */
void Region1D::setAsFuel() {

	if (_moderator)
		log_printf(WARNING, "Setting region %s to be fuel after it was set"
				"to be moderator", _region_name);

	_fuel = true;
	_moderator = false;
	_region_type = TWO_REGION_PIN_CELL;
}



/**
 * Sets this Region1D to be the moderator region in a two region pin
 * cell simulation
 */
void Region1D::setAsModerator() {

	if (_moderator)
		log_printf(WARNING, "Setting region %s to be moderator after it"
				" was set to be fuel", _region_name);

	_moderator = true;
	_fuel = false;
	_region_type = TWO_REGION_PIN_CELL;
}


/**
 * Sets this Region1D's other region (either moderator or fuel if this region
 * is fuel or moderator) for a two region pin cell simulation
 * @param other_region a pointer to the other Region1D
 */
void Region1D::setOtherPinCellRegion(Region1D* other_region) {
	_other_region = other_region;
}



/**
 * Sets the escape cross-section, beta, alpha1 and alpha2 parameters
 * used for the two region pin cell simulation
 * @param sigma_e the escape cross-section
 * @param beta Carlvik's beta parameter
 * @param alpha1 Carlvik's alpha1 parameter
 * @param alpha2 Carlvik's alpha2 parameter
 */
void Region1D::setTwoRegionPinCellParams(float sigma_e, float beta,
											float alpha1, float alpha2) {
	_sigma_e = sigma_e;
	_beta = beta;
	_alpha1 = alpha1;
	_alpha2 = alpha2;
}


/**
 * Sets this Region1D to use implicit capture variance reduction
 * @param weight_cutoff the weight cutoff below which we use Russian
 * Roulette to decide which neutrons survive
 */
void Region1D::useImplicitCapture(float weight_low, float weight_avg) {
	_use_implicit_capture = true;
	_weight_low = weight_low;
	_weight_avg = weight_avg;
	return;
}


/**
 * Sets this Region1D to use forced collision variance reduction
 * @param weight_cutof the weight cutoff below which we use Russian
 * Roulette to decide which neutrons survive
 */
void Region1D::useForcedCollision(float weight_low, float weight_avg) {
	_use_forced_collision = true;
	_weight_low = weight_low;
	_weight_avg = weight_avg;
	return;
}


/**
 * This method plays Russian Roulette with a neutron by checking if it
 * has a weight below the low weight and if so, it kills the neutron
 * with a probability equal to the neutron's weight divided by weight
 * average
 * @param neutron a pointer to the neutron of interest
 * @return a boolean to kill (true) or not (false)
 */
bool Region1D::playRussianRoulette(neutron* neutron) {

	bool kill_neutron = false;

	log_printf(DEBUG, "Playing russian roulette in region %s with neutron "
			"with weight = %f against cutoff= %f", _region_name,
			neutron->_weight, _weight_low);

	/* Only play Roulette if the weight of the neutron is below the
	 * cutoff */
	if (neutron->_weight < _weight_low) {
		float test = float(rand()) / RAND_MAX;

		if (test > (neutron->_weight / _weight_avg))
			kill_neutron = true;
	}


	return kill_neutron;
}


/**
 * This function computes the two-region fuel-to-fuel collision probability for
 * a two-region pin cell simulation. It uses Carlvik's two-term rational model
 * and assumes that the escape cross-section (_sigma_e), _beta, _alpha1 and
 * _alpha2 have all been set or else it throws an error
 * @param energy_index index into the material's energy grid
 * @return the fuel-to-fuel collision probability at that energy
 */
float Region1D::computeFuelFuelCollisionProb(int energy_index) {

	float p_ff;

	/* If this is the fuel region, we can compute p_ff directly */
	if (_fuel) {

		/* Check that all necessary parameters to compute p_ff have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to compute a fuel-fuel collision "
					"probability for region %s since beta, sigma_e, "
					"alpha1, or alpha2 for this region has not yet been set",
					_region_name);

		float sigma_tot_fuel = _material->getTotalMacroXS(energy_index);

		p_ff = ((_beta*sigma_tot_fuel) / (_alpha1*_sigma_e + sigma_tot_fuel)) +
			((1.0 - _beta)*sigma_tot_fuel / (_alpha2*_sigma_e + sigma_tot_fuel));
	}

	/* If this is the moderator region, we ask fuel region to compute p_ff */
	else {

		if (_other_region == NULL)
			log_printf(ERROR, "Unable to compute fuel-fuel collision "
					"probability for region %s since other region has"
					" not been set", _region_name);

		p_ff = _other_region->computeFuelFuelCollisionProb(energy_index);
	}

	return p_ff;
}


/**
 * This function computes the two-region moderator-to-fuel collision
 * probability for a two-region pin cell simulation. It uses Carlvik's
 * two-term rational model and assumes that the escape cross-section
 * (_sigma_e), _beta, _alpha1 and _alpha2 have all been set or else it
 * throws an error
 * @param energy_index index into the material's energy grid
 * @return the moderator-to-fuel collision probability at that energy
 */
float Region1D::computeModeratorFuelCollisionProb(int energy_index) {

	float p_mf;

	/* If this is the fuel region, we can compute p_mf directly */
	if (_fuel) {

		/* Check that all necessary parameters to compute p_mf have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to compute a moderator-fuel collision "
					"probability for region %s since beta, sigma_e, alpha1, "
					"or alpha2 for this region has not yet been set",
					_region_name);

		float p_ff = computeFuelFuelCollisionProb(energy_index);
		float p_fm = 1.0 - p_ff;

		float tot_sigma_f = _material->getTotalMacroXS(energy_index);
		float tot_sigma_mod = _other_region->getMaterial()->getTotalMacroXS(energy_index);
		float v_mod = _other_region->getVolume();

		p_mf = p_fm*(tot_sigma_f*_volume) / (tot_sigma_mod*v_mod);
	}

	/* If this is the moderator region, we ask fuel region to compute p_mf */
	else {

		if (_other_region == NULL)
			log_printf(ERROR, "Unable to compute moderator-fuel collision "
					"probability for region %s since other region has"
					" not been set", _region_name);

		p_mf = _other_region->computeModeratorFuelCollisionProb(energy_index);
	}

	return p_mf;
}


/**
 * Clear this region's vector of Binner class object pointers
 */
void Region1D::clearBinners() {
	_bin_sets.clear();
}


/**
 * Check if this region contains a certain x value
 * @param x the value to check
 * @return if contained (true), otherwise (false)
 */
bool Region1D::contains(float x) {

	if  (_region_type != BOUNDED)
		return true;
	else if (x >= _left_surface->getX() && x <= _right_surface->getX())
		return true;

	return false;
}


/**
 * Checks if a point is on either the left or right walls
 * @param x the value to check
 * @return if on boundary (true), otherwise false
 */
bool Region1D::onBoundary(float x) {

	if (_region_type == INFINITE)
		return false;
	else
		return (_left_surface->onSurface(x) || _right_surface->onSurface(x));
}


/**
 * Adds a new neutron to this region
 * @param neutron a pointer to a neutron
 */
void Region1D::addNeutron(neutron* neutron) {

	if (!contains(neutron->_x))
		log_printf(ERROR, "Cannot add a neutron to region %s"
				"since it does not contain it (bounded by %f and %f)",
				_region_name, _left_surface->getX(), _right_surface->getX());

	_neutrons.push_back(neutron);
}


/**
 * This method adds a neutron to the Region1D's vector of tranfser
 * neutrons - neutrons which have been transferred to this region during
 * a two-region pin cell simulation
 * @param neutron a pointer to the neutron to transfer
 */
void Region1D::transferNeutron(neutron* neutron) {
	_transfer_neutrons.push_back(neutron);
}


/**
 * This method loops over all neutrons which have been transferred to this
 * Region1D in a two-region pin cell simulation and adds each neutron to
 * the vector of neutrons within this Region1D. This is an intermediate
 * step between moving neutrons between regions depending on the collision
 * probabilities and actually colliding neutrons within each region
 */
void Region1D::initializeTransferredNeutrons() {

	neutron* curr;
	std::vector<neutron*>::iterator iter1;

	if (_region_type != TWO_REGION_PIN_CELL)
		log_printf(ERROR, "Cannot initialize transferred neutrons in region %s"
				"since the parameters for two region pin cell model have not "
				"been initialized yet", _region_name);

	/* Iterate over all neutrons transferred to this Region1D */
	for (iter1 = _transfer_neutrons.begin(); iter1 !=
										_transfer_neutrons.end(); ++iter1) {

		curr = (*iter1);

		/* Adds the transferred neutron to the neutrons to collide */
		addNeutron(curr);
		iter1 = _transfer_neutrons.erase(iter1);
		--iter1;
	}

}


/**
 * This method loops over all this Region1D's neutrons and either computes
 * the two-region collision probabilities and either keeps the neutron
 * in this Region1D or transfers for its next collision in the 2nd Region1D
 * in the pin cell simulation
 */
void Region1D::twoRegionNeutronTransferral() {

	neutron* curr;
	std::vector<neutron*>::iterator iter1;
	int energy_index;

	if (!_material->isRescaled())
		log_printf(ERROR, "Region %s is unable to transfer neutrons since "
				"it's Material %s has not been rescaled",
				_region_name, _material->getMaterialName());

	if (_region_type != TWO_REGION_PIN_CELL)
		log_printf(ERROR, "Cannot initialize transferred neutrons in region %s"
				"since the parameters for two region pin cell model have not "
				"been initialized yet", _region_name);

	/* Iterate over all neutrons that are inside this Region1D */
	for (iter1 = _neutrons.begin(); iter1 != _neutrons.end(); ++iter1) {

		curr = (*iter1);

		/* find index into the material's energy grid for this neutron */
		energy_index = _material->getEnergyGridIndex(curr->_energy);

		float p_ff = computeFuelFuelCollisionProb(energy_index);
		float p_mf = computeModeratorFuelCollisionProb(energy_index);
		float test = float(rand()) / RAND_MAX;

		/* If this Region1D is the fuel */
		if (_fuel) {

			/* If test is larger than p_ff, move to moderator */
			if (test > p_ff) {

				_other_region->transferNeutron(curr);
				iter1 = _neutrons.erase(iter1);
				--iter1;
			}
		}

		/* If this Region1D is the moderator */
		else {

			/* If test is larger than p_mf, move to fuel */
			if (test < p_mf) {

				_other_region->transferNeutron(curr);
				iter1 = _neutrons.erase(iter1);
				--iter1;
			}
		}
	}

	return;
}


/**
 * This method contains all of the neutron physics for the simulation
 * It moves neutrons by sampling their path traveled and either collides
 * them with a isotope in the region, or places them on one of the bounding
 * surfaces
 */
void Region1D::moveNeutrons() {

	log_printf(DEBUG, "Inside moveNeutrons method for region %s with "
			"%d neutrons", _region_name, getNumNeutrons());

	if (!_material->isRescaled())
		log_printf(ERROR, "Region %s is unable to move neutrons since it's "
				"Material %s has not been rescaled",
				_region_name, _material->getMaterialName());

	int energy_index;
	float sigma_t;
	float sigma_a;
	float path_length;
	float new_x;
	neutron* curr;
	Isotope* isotope;
	Binner* curr_bin;
	collisionType collision_type;
	std::vector<Binner*>::iterator iter3;
	std::vector<neutron*>::iterator iter1;

	/**************************************************************************
	 ************************  MONTE CARLO KERNEL  ****************************
	 *************************************************************************/
	/* Iterate over all neutrons that are inside this Region1D */
	for (iter1 = _neutrons.begin(); iter1 != _neutrons.end(); ++iter1) {
		curr = (*iter1);

		log_printf(DEBUG, "Neutron x = %f, mu = %f, energy = %f", curr->_x,
													curr->_mu, curr->_energy);

		/* Compute total sigma_t from all isotopes */
		energy_index = _material->getEnergyGridIndex(curr->_energy);
		sigma_t = _material->getTotalMacroXS(energy_index);

		/* Compute path length and the new x position of this neutron */
		path_length = -log(float(rand()) / RAND_MAX) / sigma_t;
		new_x = curr->_x + curr->_mu * path_length;

		log_printf(DEBUG, "sigma_t = %f, path_length = %f, new_x = %f",
											sigma_t, path_length, new_x);

		/* The neutron collided within this region */
		if (contains(new_x)) {

			/* Figure out which isotope the neutron collided in */
			isotope = _material->sampleIsotope(curr->_energy);

			log_printf(DEBUG, "Neutron collided in isotope: %s",
											isotope->getIsotopeType());

			/* Figure out the collision type */
			collision_type = isotope->getCollisionType(curr->_energy);

			/* Update the neutron position */
			curr->_x = new_x;

			log_printf(DEBUG, "collision type = %d, new_x = %f",
												collision_type, new_x);

			/* Tally neutron for each bin set */
			for (iter3 = _bin_sets.begin(); iter3 != _bin_sets.end(); ++iter3){

				curr_bin = *iter3;
				char* isotopes = curr_bin->getIsotopes();

				/* This tally type is for all isotopes */
				if (strcmp(isotopes, "all") == 0) {

					/* Spatial flux tally for all isotopes */
					if (curr_bin->getTallyType() == FLUX_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight /
																sigma_t);

					/* Energy flux tally */
					else if (curr_bin->getTallyType() == FLUX_ENERGY)
						curr_bin->weightedTally(curr->_energy,
												curr->_weight / sigma_t);

					/* Spatial capture rate tally */
					else if (curr_bin->getTallyType()==CAPTURE_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
						_material->getCaptureMacroXS(energy_index) / sigma_t);

					/* Energy absorption rate tally */
					else if (curr_bin->getTallyType()==CAPTURE_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
						_material->getCaptureMacroXS(energy_index) / sigma_t);

					/* Spatial absorption rate tally */
					else if (curr_bin->getTallyType()==ABSORPTION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
						_material->getAbsorbMacroXS(energy_index) / sigma_t);

					/* Energy absorption rate tally */
					else if (curr_bin->getTallyType()==ABSORPTION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
						_material->getAbsorbMacroXS(energy_index) / sigma_t);

					/* Spatial elastic rate tally */
					else if (curr_bin->getTallyType()==ELASTIC_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
						_material->getElasticMacroXS(energy_index) / sigma_t);

					/* Energy elastic rate tally */
					else if (curr_bin->getTallyType()==ELASTIC_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
						_material->getElasticMacroXS(energy_index) / sigma_t);

					/* Spatial fission rate tally */
					else if (curr_bin->getTallyType()==FISSION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
						_material->getFissionMacroXS(energy_index) / sigma_t);

					/* Energy fission rate tally */
					else if (curr_bin->getTallyType()==FISSION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
						_material->getFissionMacroXS(energy_index) / sigma_t);

					/* Spatial transport rate tally */
					else if (curr_bin->getTallyType()==TRANSPORT_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
						_material->getTransportMacroXS(energy_index) / sigma_t);

					/* Energy transport rate tally */
					else if (curr_bin->getTallyType()==TRANSPORT_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
						_material->getTransportMacroXS(energy_index) / sigma_t);

					/* Spatial collision rate tally */
					else if (curr_bin->getTallyType()==COLLISION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight);

					/* Energy collision rate tally */
					else if (curr_bin->getTallyType()==COLLISION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight);

					/* Spatial diffusion rate tally */
					else if (curr_bin->getTallyType()==DIFFUSION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
							1.0 / (3.0 *_material->getTransportMacroXS(curr->_energy)
																* sigma_t));

					/* Energy diffusion rate tally */
					else if (curr_bin->getTallyType()==DIFFUSION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
							1.0 / (3.0 *_material->getTransportMacroXS(energy_index)
																* sigma_t));
				}

				/* This tally type is for one isotope type */
				else {

					/* Spatial flux tally for all isotopes */
					if (curr_bin->getTallyType() == FLUX_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight
															/ sigma_t);

					/* Energy flux tally */
					else if (curr_bin->getTallyType() == FLUX_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight
																/ sigma_t);

					/* Spatial capture rate tally */
					else if (curr_bin->getTallyType() == CAPTURE_RATE_SPATIAL) {
						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_c =
						isotope->getCaptureXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_c / sigma_t);
					}

					/* Energy capture rate tally */
					else if (curr_bin->getTallyType() == CAPTURE_RATE_ENERGY) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_c =
						isotope->getCaptureXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_c / sigma_t);
					}

					/* Spatial absorption rate tally */
					else if (curr_bin->getTallyType() == ABSORPTION_RATE_SPATIAL) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_a =
						isotope->getAbsorbXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_a / sigma_t);
					}

					/* Energy absorption rate tally */
					else if (curr_bin->getTallyType() == ABSORPTION_RATE_ENERGY) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_a =
						isotope->getAbsorbXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_a / sigma_t);
					}

					/* Spatial elastic rate tally */
					else if (curr_bin->getTallyType() == ELASTIC_RATE_SPATIAL) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_e =
						isotope->getElasticXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_e / sigma_t);
					}

					/* Energy elastic rate tally */
					else if (curr_bin->getTallyType() == ELASTIC_RATE_ENERGY) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_e =
						isotope->getElasticXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_e / sigma_t);
					}

					/* Spatial fission rate tally */
					else if (curr_bin->getTallyType() == FISSION_RATE_SPATIAL) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_f =
						isotope->getFissionXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_f / sigma_t);
					}

					/* Energy fission rate tally */
					else if (curr_bin->getTallyType() == FISSION_RATE_ENERGY) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_f =
						isotope->getFissionXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_f / sigma_t);
					}

					/* Spatial fission rate tally */
					else if (curr_bin->getTallyType() == TRANSPORT_RATE_SPATIAL) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_tr =
						isotope->getTransportXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_tr / sigma_t);
					}

					/* Energy fission rate tally */
					else if (curr_bin->getTallyType() == TRANSPORT_RATE_ENERGY) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_tr =
						isotope->getTransportXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_tr / sigma_t);
					}

					/* Spatial collision rate tally */
					else if (curr_bin->getTallyType() == COLLISION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight);

					/* Energy collision rate tally */
					else if (curr_bin->getTallyType() == COLLISION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight);

					/* Spatial diffusion rate tally */
					else if (curr_bin->getTallyType() == DIFFUSION_RATE_SPATIAL) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_tr =
						isotope->getTransportXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
										1.0 / (3.0 * isotope_sigma_tr * sigma_t));
					}

					/* Energy diffusion rate tally */
					else if (curr_bin->getTallyType() == DIFFUSION_RATE_ENERGY) {

						Isotope* isotope = _material->getIsotope(isotopes);
						float num_density = _material->getIsotopeNumDensity(isotopes);
						float isotope_sigma_tr =
						isotope->getTransportXS(energy_index) * num_density * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
										1.0 / (3.0 * isotope_sigma_tr * sigma_t));
					}
				}
			}

			/******************************************************************
			 **********************  VARIANCE REDUCTION  **********************
			 *****************************************************************/
			/* Implicit capture variance reduction */
			if (_use_implicit_capture) {

				log_printf(DEBUG, "Applying implicit capture to neutron with "
						"weight = %f", curr->_weight);

				/* If we are using implicit capture, force the collision to be
				 * a scattering collision */
				if (_use_implicit_capture) {

					/* Sample a isotope with a scattering cross-section */
					while (isotope->getScatterXS(curr->_energy) == 0)
						isotope = _material->sampleIsotope(curr->_energy);

					while (collision_type == CAPTURE)
						collision_type =
									isotope->getCollisionType(curr->_energy);
				}

				/* Get the total absorption cross-section for this Region1D */
				sigma_a = _material->getAbsorbMacroXS(energy_index);

				/* Reduce weight by sigma_a / sigma_t */
				curr->_weight *= (1.0 - (sigma_a / sigma_t));

				log_printf(DEBUG, "New weight = %f", curr->_weight);
			}

			/* Forced collision variance reduction */
			else if (_use_forced_collision) {

				log_printf(DEBUG, "Applying forced collision to neutron with "
						"weight = %f", curr->_weight);

				/* Find distance to nearest wall along neutron's trajectory */
				float d;
				if (curr->_mu <= 0)
					d = (curr->_x - _left_surface->getX()) / curr->_mu;
				else
					d = (_right_surface->getX() - curr->_x) / curr->_mu;

				/* Compute exponential term in forced collision prob dist */
				float px = exp(-sigma_t * d);

				/* Create a new uncollided neutron that hits next surface
				 * in the region along the collided neutron's trajectory */
				neutron* new_neutron = initializeNewNeutron();
				new_neutron->_energy = curr->_energy;
				new_neutron->_mu = curr->_mu;
				new_neutron->_weight = curr->_weight * px;

				/* Figure out which surface to put uncollided neutron on */
				if (curr->_mu <= 0)
					_left_surface->addNeutron(new_neutron);
				else
					_right_surface->addNeutron(new_neutron);

				/* Update collided neutron's properties */
				curr->_x += -curr->_mu * log(1.0 - float(rand()/RAND_MAX)
												* (1.0 - px)) / sigma_t;
				curr->_weight *= (1.0 - px);
			}

			/******************************************************************
			 ***********************  COLLISIONS  *****************************
			 *****************************************************************/
			/* Check if collision type was absorption and if so kill neutron */
			if (collision_type == CAPTURE && !_use_implicit_capture) {

				log_printf(DEBUG, "Capture type collision");

				iter1 = _neutrons.erase(iter1);
				--iter1;
				delete curr;

				continue;
			}

			/* Check if collision type was fission and if so kill neutron */
			if (collision_type == FISSION && !_use_implicit_capture) {

				log_printf(DEBUG, "Fission type collision");

				iter1 = _neutrons.erase(iter1);
				--iter1;
				delete curr;
				continue;
			}

			/* Check if collision was general scattering */
			else if (collision_type == SCATTER) {

				log_printf(DEBUG, "Normal scatter type collision");

				/* Isotropic in lab */
				if (isotope->getScatterAngleType() == ISOTROPIC_LAB) {
					float A = float(isotope->getA());
					float energy = curr->_energy * (1.0 -
						(1.0 - isotope->getAlpha())*(float(rand())/RAND_MAX));
					float mu_cm = ((2 * energy / curr->_energy) -
									isotope->getAlpha() - 1.0) / (1.0 -
														isotope->getAlpha());
					curr->_mu = (1.0 + A*mu_cm) / sqrt(A*A + 1.0 + 2*mu_cm*A);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy = energy;
				}

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = isotope->getA();
					float mu_l = (1.0 + A*mu_cm) / (sqrt(A*A + 2*A*mu_cm+1.0));
					float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy *= (A*A + 2*A*mu_cm + 1.0) /
															((A+1.0)*(A+1.0));
				}

				log_printf(DEBUG, "Updated mu= %f, energy = %f",
										curr->_mu, curr->_energy);
			}

			/* Check if collision type was inelastic scattering */
			else if (collision_type == INELASTIC) {

				log_printf(DEBUG, "Inelastic scatter type collision");

				/* Isotropic in lab */
				if (isotope->getInelasticAngleType() == ISOTROPIC_LAB)
					curr->_mu = (float(rand()) / RAND_MAX) * 2.0 - 1.0;

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = isotope->getA();
					float mu_l = (1.0 + A*mu_cm) / (sqrt(A*A + 2*A*mu_cm+1.0));
					float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);
				}

				/* Update the neutron energy */
				curr->_energy =
						isotope->getInelasticScatterEnergy(curr->_energy);

				log_printf(DEBUG, "Updated mu= %f, energy = %f",
									curr->_mu, curr->_energy);
			}

			/* Check if collision type was elastic scattering */
			else if (collision_type == ELASTIC) {

				log_printf(DEBUG, "Elastic scatter type collision");

				/* Isotropic in lab */
				if (isotope->getElasticAngleType() == ISOTROPIC_LAB) {
					float A = float(isotope->getA());
					float energy = curr->_energy * (1.0 -
						(1.0 - isotope->getAlpha())*(float(rand())/RAND_MAX));
					float mu_cm = ((2 * energy / curr->_energy) -
									isotope->getAlpha() - 1.0) / (1.0 -
														isotope->getAlpha());
					curr->_mu = (1.0 + A*mu_cm) / sqrt(A*A + 1.0 + 2*mu_cm*A);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy = energy;

					//CORRECTION: xs only extend to 1E-5
					if (curr->_energy < 1E-5)
						curr->_energy = 1.1E-5;

				}

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = isotope->getA();
					float mu_l = (1.0 + A*mu_cm)/(sqrt(A*A + 2*A*mu_cm + 1.0));
					float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy *= (A*A + 2*A*mu_cm + 1.0) /
															((A+1.0)*(A+1.0));
				}

				log_printf(DEBUG, "Updated mu= %f, energy = %f",
													curr->_mu, curr->_energy);
			}

			/* Play Russian Roulette with neutron */
			if (curr->_weight < _weight_low && !_use_forced_collision) {
				if (playRussianRoulette(curr)) {
					log_printf(DEBUG, "Killing neutron");
					iter1 = _neutrons.erase(iter1);
					delete curr;
					--iter1;
				}
				else
					curr->_weight = _weight_avg;
			}
		}

		/******************************************************************
		 ********************  SURFACE INTERSECTIONS  *********************
		 *****************************************************************/
		/* If the new_x crossed the left surface, place it on the surface */
		else if (new_x < _left_surface->getX()) {
			_left_surface->addNeutron(curr);
			iter1 = _neutrons.erase(iter1);
			--iter1;
		}

		/* If the new_x crossed the right surface, place it on the surcace */
		else {
			_right_surface->addNeutron(curr);
			iter1 = _neutrons.erase(iter1);
			--iter1;
		}

		/* If we removed the final neutron from the vector, break loop */
		if (iter1 == _neutrons.end())
			break;
	}
}


/**
 * Returns the number of neutrons inside of this region
 * @return the number of neutrons contained by this region
 */
int Region1D::getNumNeutrons() {
	return _neutrons.size();
}
