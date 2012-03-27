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

	/* Default number density and volume */
	_tot_num_density = 0.0;
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
Region1D::~Region1D() {

	std::map<char*, std::pair<float, Material*> >::iterator iter;
	Material* curr;

	for (iter = _materials.begin(); iter != _materials.end(); ++iter) {
		curr = iter->second.second;
		delete curr;
	}

}


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
 * Returns the total number density for all materials within
 * this region
 * @return the total number density (at/cm^3)
 */
float Region1D::getTotalNumberDensity() {
	return _tot_num_density;
}


/**
 * Returns the volume of this Region1D (cm^3)
 * @return the Region1D's volume
 */
float Region1D::getVolume() {
	return _volume;
}

/**
 * Returns the total macroscopic cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic cross-section (cm^-1)
 */
float Region1D::getTotalMacroXS(float energy) {

	float sigma_t = 0;

	/* Increment sigma_t for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_t += iter->second.second->getTotalXS(energy)
											* iter->second.first * 1E-24;

	return sigma_t;
}


/**
 * Returns the total microscopic cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic cross-section (barns)
 */
float Region1D::getTotalMicroXS(float energy) {

	float sigma_t = 0;

	/* Increment sigma_t for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_t += iter->second.second->getTotalXS(energy);

	return sigma_t;
}


/**
 * Returns the total macroscopic capture cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic capture cross-section (cm^-1)
 */
float Region1D::getCaptureMacroXS(float energy) {

	float sigma_c = 0;

	/* Increment sigma_a for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_c += iter->second.second->getCaptureXS(energy) *
											iter->second.first * 1E-24;

	return sigma_c;
}


/**
 * Returns the total microscopic capture cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic capture cross-section (barns)
 */
float Region1D::getCaptureMicroXS(float energy) {

	float sigma_a = 0;

	/* Increment sigma_a for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_a += iter->second.second->getCaptureXS(energy);

	return sigma_a;
}


/**
 * Returns the total macroscopic scattering cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic scattering cross-section (cm^-1)
 */
float Region1D::getScatterMacroXS(float energy){

	float sigma_s = 0;

	/* Increment sigma_s for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter) {
		sigma_s += iter->second.second->getScatterXS(energy)
												* iter->second.first * 1E-24;
		sigma_s += iter->second.second->getInelasticXS(energy)
												* iter->second.first * 1E-24;
		sigma_s += iter->second.second->getElasticXS(energy)
												* iter->second.first * 1E-24;
	}

	return sigma_s;
}


/**
 * Returns the total microscopic scattering cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic scattering cross-section (barns)
 */
float Region1D::getScatterMicroXS(float energy){

	float sigma_s = 0;

	/* Increment sigma_s for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter) {
		sigma_s += iter->second.second->getScatterXS(energy);
		sigma_s += iter->second.second->getInelasticXS(energy);
		sigma_s += iter->second.second->getElasticXS(energy);
	}

	return sigma_s;
}


/**
 * Returns the total macroscopic elastic scattering cross-section within
 * this Region1D at some energy
 * @param energy energy of interest (eV)
 * @return the total elastic macroscopic scattering cross-section (cm^-1)
 */
float Region1D::getElasticMacroXS(float energy) {

	float sigma_s = 0;

	/* Increment sigma_s for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_s += iter->second.second->getElasticXS(energy)
												* iter->second.first * 1E-24;

	return sigma_s;
}


/**
 * Returns the total macroscopic elastic scattering cross-section within
 * this Region1D at some energy
 * @param energy energy of interest (eV)
 * @return the total elastic microscopic scattering cross-section (barns)
 */
float Region1D::getElasticMicroXS(float energy) {

	float sigma_s = 0;

	/* Increment sigma_s for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_s += iter->second.second->getElasticXS(energy);

	return sigma_s;
}


/**
 * Returns the total macroscopic inelastic scattering cross-section within
 * this Region1D at some energy
 * @param energy energy of interest (eV)
 * @return the total inelastic macroscopic scattering cross-section (cm^-1)
 */
float Region1D::getInelasticMacroXS(float energy) {

	float sigma_s = 0;

	/* Increment sigma_s for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_s += iter->second.second->getInelasticXS(energy)
												* iter->second.first * 1E-24;

	return sigma_s;
}


/**
 * Returns the total macroscopic inelastic scattering cross-section within
 * this Region1D at some energy
 * @param energy energy of interest (eV)
 * @return the total inelastic microscopic scattering cross-section (barns)
 */
float Region1D::getInelasticMicroXS(float energy) {

	float sigma_s = 0;

	/* Increment sigma_s for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_s += iter->second.second->getInelasticXS(energy);

	return sigma_s;
}


/**
 * Returns the total macroscopic fission cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic fission cross-section (cm^-1)
 */
float Region1D::getFissionMacroXS(float energy) {

	float sigma_f = 0;

	/* Increment sigma_f for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_f += iter->second.second->getFissionXS(energy) *
											iter->second.first * 1E-24;

	return sigma_f;
}


/**
 * Returns the total microscopic fission cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic fission cross-section (barns)
 */
float Region1D::getFissionMicroXS(float energy) {

	float sigma_f = 0;

	/* Increment sigma_f for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_f += iter->second.second->getFissionXS(energy);

	return sigma_f;
}


/**
 * Returns the total macroscopic absorption cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic absorption cross-section (cm^-1)
 */
float Region1D::getAbsorbMacroXS(float energy) {
	float sigma_a = 0;

	/* Increment sigma_a for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_a += iter->second.second->getAbsorbXS(energy) *
											iter->second.first * 1E-24;

	return sigma_a;
}


/**
 * Returns the total microscopic absorption cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic absorption cross-section (barns)
 */
float Region1D::getAbsorbMicroXS(float energy) {
	float sigma_a = 0;

	/* Increment sigma_a for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_a += iter->second.second->getAbsorbXS(energy);

	return sigma_a;
}


/**
 * Returns the total microscopic transport cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic transport cross-section (barns)
 */
float Region1D::getTransportMicroXS(float energy) {
	float sigma_tr = 0;

	/* Increment sigma_a for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_tr += iter->second.second->getTransportXS(energy);

	return sigma_tr;
}


/**
 * Returns the total macroscopic transport cross-section within this Region1D
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic transport cross-section (cm^-1)
 */
float Region1D::getTransportMacroXS(float energy) {
	float sigma_tr = 0;

	/* Increment sigma_a for each material */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		sigma_tr += iter->second.second->getTransportXS(energy) *
										iter->second.first * 1E-24;

	return sigma_tr;
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
 * Adds a new material to this region
 * @param material a pointer to a material class object
 * @param num_density the number density (at/cm^3) for this material
 */
void Region1D::addMaterial(Material* material, float num_density) {

	/* Creates a pair between the number density and material pointer */
	std::pair<float, Material*> new_pair = std::pair<float, Material*>
													(num_density, material);

	std::pair<char*, std::pair<float, Material*> > new_material =
							std::pair<char*, std::pair<float, Material*> >
									(material->getIsotopeType(), new_pair);

	/* Inserts the material and increments the total number density */
	_materials.insert(new_material);
	_tot_num_density += num_density;

	return;
}


/**
 * Sets this Region1D's name as specified by a user
 * @param region_name the name of this Region1D
 */
void Region1D::setRegionName(char* region_name) {
	_region_name = region_name;
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
 * Samples a material for a collision with a probability based on the
 * ratios of each material's total cross-section to the total cross-section
 * of all material's in this Region1D
 * @return a pointer to the chosen material
 */
Material* Region1D::sampleMaterial(float energy) {

	float sigma_t = getTotalMacroXS(energy);
	float sigma_t_ratio = 0.0;
	float new_sigma_t_ratio = 0.0;
	float test = float(rand()) / RAND_MAX;

	/* Loop over all materials */
	std::map<char*, std::pair<float, Material*> >::iterator iter;
	Material* material = NULL;
	for (iter =_materials.begin(); iter !=_materials.end(); ++iter){

		new_sigma_t_ratio += (iter->second.second->getTotalXS(energy) *
										iter->second.first * 1E-24) / sigma_t;

		if (test >= sigma_t_ratio && ((test <= new_sigma_t_ratio) ||
							fabs(test - new_sigma_t_ratio) < 1E-5)) {
			material = iter->second.second;
			break;
		}
		sigma_t_ratio = new_sigma_t_ratio;
	}

	if (material == NULL)
		log_printf(ERROR, "Unable to find material type in region %s"
				" moveNeutron method, test = %1.20f, new_num_density_ratio "
				"= %1.20f", _region_name, test, new_sigma_t_ratio);

	return material;
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
 * @param energy a neutrons' energy (eV)
 * @return the fuel-to-fuel collision probability at that energy
 */
float Region1D::computeFuelFuelCollisionProb(float energy) {

	float p_ff;

	/* If this is the fuel region, we can compute p_ff directly */
	if (_fuel) {

		/* Check that all necessary parameters to compute p_ff have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to compute a fuel-fuel collision "
					"probability for region %s since beta, sigma_e, "
					"alpha1, or alpha2 for this region has not yet been set",
					_region_name);

		float sigma_tot_fuel = getTotalMacroXS(energy);

		p_ff = ((_beta*sigma_tot_fuel) / (_alpha1*_sigma_e + sigma_tot_fuel)) +
			((1.0 - _beta)*sigma_tot_fuel / (_alpha2*_sigma_e + sigma_tot_fuel));
	}

	/* If this is the moderator region, we ask fuel region to compute p_ff */
	else {

		if (_other_region == NULL)
			log_printf(ERROR, "Unable to compute fuel-fuel collision "
					"probability for region %s since other region has"
					" not been set", _region_name);

		p_ff = _other_region->computeFuelFuelCollisionProb(energy);
	}

	return p_ff;
}


/**
 * This function computes the two-region moderator-to-fuel collision
 * probability for a two-region pin cell simulation. It uses Carlvik's
 * two-term rational model and assumes that the escape cross-section
 * (_sigma_e), _beta, _alpha1 and _alpha2 have all been set or else it
 * throws an error
 * @param energy a neutrons' energy (eV)
 * @return the moderator-to-fuel collision probability at that energy
 */
float Region1D::computeModeratorFuelCollisionProb(float energy) {

	float p_mf;

	/* If this is the fuel region, we can compute p_mf directly */
	if (_fuel) {

		/* Check that all necessary parameters to compute p_mf have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to compute a moderator-fuel collision "
					"probability for region %s since beta, sigma_e, alpha1, "
					"or alpha2 for this region has not yet been set",
					_region_name);

		float p_ff = computeFuelFuelCollisionProb(energy);
		float p_fm = 1.0 - p_ff;

		float tot_sigma_f = getTotalMacroXS(energy);
		float tot_sigma_mod = _other_region->getTotalMacroXS(energy);
		float v_mod = _other_region->getVolume();

		p_mf = p_fm*(tot_sigma_f*_volume) / (tot_sigma_mod*v_mod);
	}

	/* If this is the moderator region, we ask fuel region to compute p_mf */
	else {

		if (_other_region == NULL)
			log_printf(ERROR, "Unable to compute moderator-fuel collision "
					"probability for region %s since other region has"
					" not been set", _region_name);

		p_mf = _other_region->computeModeratorFuelCollisionProb(energy);
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

//	log_printf(NORMAL, "Transferring neutrons...");

	neutron* curr;
	std::vector<neutron*>::iterator iter1;

	if (_region_type != TWO_REGION_PIN_CELL)
		log_printf(ERROR, "Cannot initialize transferred neutrons in region %s"
				"since the parameters for two region pin cell model have not "
				"been initialized yet", _region_name);

	/* Iterate over all neutrons that are inside this Region1D */
	for (iter1 = _neutrons.begin(); iter1 != _neutrons.end(); ++iter1) {

		curr = (*iter1);

		float p_ff = computeFuelFuelCollisionProb(curr->_energy);
		float p_mf = computeModeratorFuelCollisionProb(curr->_energy);
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
 * them with a material in the region, or places them on one of the bounding
 * surfaces
 */
void Region1D::moveNeutrons() {

	log_printf(DEBUG, "Inside moveNeutrons method for region %s with "
			"%d neutrons", _region_name, getNumNeutrons());

	float sigma_t;
	float sigma_a;
	float path_length;
	float new_x;
	neutron* curr;
	Material* material;
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

		/* Compute total sigma_t from all materials */
		sigma_t = getTotalMacroXS(curr->_energy);

		/* Compute path length and the new x position of this neutron */
		path_length = -log(float(rand()) / RAND_MAX) / sigma_t;
		new_x = curr->_x + curr->_mu * path_length;

		log_printf(DEBUG, "sigma_t = %f, path_length = %f, new_x = %f",
											sigma_t, path_length, new_x);

		/* The neutron collided within this region */
		if (contains(new_x)) {

			/* Figure out which material the neutron collided in */
			material = sampleMaterial(curr->_energy);

			log_printf(DEBUG, "Neutron collided in material: %s",
											material->getIsotopeType());

			/* Figure out the collision type */
			collision_type = material->getCollisionType(curr->_energy);

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
									getCaptureMacroXS(curr->_energy) / sigma_t);

					/* Energy absorption rate tally */
					else if (curr_bin->getTallyType()==CAPTURE_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
									getCaptureMacroXS(curr->_energy) / sigma_t);

					/* Spatial absorption rate tally */
					else if (curr_bin->getTallyType()==ABSORPTION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
									getAbsorbMacroXS(curr->_energy) / sigma_t);

					/* Energy absorption rate tally */
					else if (curr_bin->getTallyType()==ABSORPTION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
									getAbsorbMacroXS(curr->_energy) / sigma_t);

					/* Spatial elastic rate tally */
					else if (curr_bin->getTallyType()==ELASTIC_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
									getElasticMacroXS(curr->_energy) / sigma_t);

					/* Energy elastic rate tally */
					else if (curr_bin->getTallyType()==ELASTIC_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
									getElasticMacroXS(curr->_energy) / sigma_t);

					/* Spatial fission rate tally */
					else if (curr_bin->getTallyType()==FISSION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
									getFissionMacroXS(curr->_energy) / sigma_t);

					/* Energy fission rate tally */
					else if (curr_bin->getTallyType()==FISSION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
									getFissionMacroXS(curr->_energy) / sigma_t);

					/* Spatial transport rate tally */
					else if (curr_bin->getTallyType()==TRANSPORT_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
								getTransportMacroXS(curr->_energy) / sigma_t);

					/* Energy transport rate tally */
					else if (curr_bin->getTallyType()==TRANSPORT_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
								getTransportMacroXS(curr->_energy) / sigma_t);

					/* Spatial collision rate tally */
					else if (curr_bin->getTallyType()==COLLISION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight);

					/* Energy collision rate tally */
					else if (curr_bin->getTallyType()==COLLISION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight);

					/* Spatial diffusion rate tally */
					else if (curr_bin->getTallyType()==DIFFUSION_RATE_SPATIAL)
						curr_bin->weightedTally(curr->_x, curr->_weight *
							1.0 / (3.0 *	getTransportMacroXS(curr->_energy)
																* sigma_t));

					/* Energy diffusion rate tally */
					else if (curr_bin->getTallyType()==DIFFUSION_RATE_ENERGY)
						curr_bin->weightedTally(curr->_energy, curr->_weight *
							1.0 / (3.0 *	getTransportMacroXS(curr->_energy)
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
						float isotope_sigma_c =
									_materials.at(isotopes).second->getCaptureXS(
						curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_c / sigma_t);
					}

					/* Energy capture rate tally */
					else if (curr_bin->getTallyType() == CAPTURE_RATE_ENERGY) {

						float isotope_sigma_c =
							_materials.at(isotopes).second->getCaptureXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_c / sigma_t);
					}

					/* Spatial absorption rate tally */
					else if (curr_bin->getTallyType() == ABSORPTION_RATE_SPATIAL) {

						float isotope_sigma_a =
									_materials.at(isotopes).second->getAbsorbXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_a / sigma_t);
					}

					/* Energy absorption rate tally */
					else if (curr_bin->getTallyType() == ABSORPTION_RATE_ENERGY) {

						float isotope_sigma_a =
							_materials.at(isotopes).second->getAbsorbXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_a / sigma_t);
					}

					/* Spatial elastic rate tally */
					else if (curr_bin->getTallyType() == ELASTIC_RATE_SPATIAL) {

						float isotope_sigma_e =
									_materials.at(isotopes).second->getElasticXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_e / sigma_t);
					}

					/* Energy elastic rate tally */
					else if (curr_bin->getTallyType() == ELASTIC_RATE_ENERGY) {

						float isotope_sigma_e =
							_materials.at(isotopes).second->getElasticXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;


						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_e / sigma_t);
					}

					/* Spatial fission rate tally */
					else if (curr_bin->getTallyType() == FISSION_RATE_SPATIAL) {

						float isotope_sigma_f =
									_materials.at(isotopes).second->getFissionXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_f / sigma_t);
					}

					/* Energy fission rate tally */
					else if (curr_bin->getTallyType() == FISSION_RATE_ENERGY) {

						float isotope_sigma_f =
							_materials.at(isotopes).second->getFissionXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_energy, curr->_weight *
													isotope_sigma_f / sigma_t);
					}

					/* Spatial fission rate tally */
					else if (curr_bin->getTallyType() == TRANSPORT_RATE_SPATIAL) {

						float isotope_sigma_tr =
									_materials.at(isotopes).second->getTransportXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
													isotope_sigma_tr / sigma_t);
					}

					/* Energy fission rate tally */
					else if (curr_bin->getTallyType() == TRANSPORT_RATE_ENERGY) {

						float isotope_sigma_tr =
							_materials.at(isotopes).second->getTransportXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

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

						float isotope_sigma_tr =
									_materials.at(isotopes).second->getTransportXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

						curr_bin->weightedTally(curr->_x, curr->_weight *
										1.0 / (3.0 * isotope_sigma_tr * sigma_t));
					}

					/* Energy diffusion rate tally */
					else if (curr_bin->getTallyType() == DIFFUSION_RATE_ENERGY) {

						float isotope_sigma_tr =
							_materials.at(isotopes).second->getTransportXS(
							curr->_energy) * _materials.at(isotopes).first * 1E-24;

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

					/* Sample a material with a scattering cross-section */
					while (material->getScatterXS(curr->_energy) == 0)
						material = sampleMaterial(curr->_energy);

					while (collision_type == CAPTURE)
						collision_type =
									material->getCollisionType(curr->_energy);
				}

				/* Get the total absorption cross-section for this Region1D */
				sigma_a = getAbsorbMacroXS(curr->_energy);

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
				if (material->getScatterAngleType() == ISOTROPIC_LAB) {
					float A = float(material->getA());
					float energy = curr->_energy * (1.0 -
						(1.0 - material->getAlpha())*(float(rand())/RAND_MAX));
					float mu_cm = ((2 * energy / curr->_energy) -
									material->getAlpha() - 1.0) / (1.0 -
														material->getAlpha());
					curr->_mu = (1.0 + A*mu_cm) / sqrt(A*A + 1.0 + 2*mu_cm*A);

					if (curr->_energy < 4 && material->usesThermalScattering())
						curr->_energy =
						material->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy = energy;
				}

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = material->getA();
					float mu_l = (1.0 + A*mu_cm) / (sqrt(A*A + 2*A*mu_cm+1.0));
					float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);

					if (curr->_energy < 4 && material->usesThermalScattering())
						curr->_energy =
						material->getThermalScatteringEnergy(curr->_energy);
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
				if (material->getInelasticAngleType() == ISOTROPIC_LAB)
					curr->_mu = (float(rand()) / RAND_MAX) * 2.0 - 1.0;

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = material->getA();
					float mu_l = (1.0 + A*mu_cm) / (sqrt(A*A + 2*A*mu_cm+1.0));
					float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);
				}

				/* Update the neutron energy */
				curr->_energy =
						material->getInelasticScatterEnergy(curr->_energy);

				log_printf(DEBUG, "Updated mu= %f, energy = %f",
									curr->_mu, curr->_energy);
			}

			/* Check if collision type was elastic scattering */
			else if (collision_type == ELASTIC) {

				log_printf(DEBUG, "Elastic scatter type collision");

				/* Isotropic in lab */
				if (material->getElasticAngleType() == ISOTROPIC_LAB) {
					float A = float(material->getA());
					float energy = curr->_energy * (1.0 -
						(1.0 - material->getAlpha())*(float(rand())/RAND_MAX));
					float mu_cm = ((2 * energy / curr->_energy) -
									material->getAlpha() - 1.0) / (1.0 -
														material->getAlpha());
					curr->_mu = (1.0 + A*mu_cm) / sqrt(A*A + 1.0 + 2*mu_cm*A);

					if (curr->_energy < 4 && material->usesThermalScattering())
						curr->_energy =
						material->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy = energy;
				}

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = material->getA();
					float mu_l = (1.0 + A*mu_cm)/(sqrt(A*A + 2*A*mu_cm + 1.0));
					float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);

					if (curr->_energy < 4 && material->usesThermalScattering())
						curr->_energy =
						material->getThermalScatteringEnergy(curr->_energy);
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


/**
 * This function plots the total macroscopic cross-sections for
 * the materials specified with a variable argument list of character
 * arrays. The final argument must be NULL so that this function knows
 * when to stop loop over materials.
 * @param start_energy the starting energy for the plot (eV)
 * @param end_energy the ending energy for the plot (eV)
 * @param num_energies the number of energies to plot (eV)
 * @param isotopes a variable length parameter list of character arrays
 * of material types
 */
void Region1D::plotMacroscopicCrossSections(float start_energy,
		float end_energy, int num_energies, char* isotopes, ...) {

	/* Allocate memory for energies and xs values */
	float* energies = logspace(start_energy, end_energy, num_energies);
	float* xs_values = new float[num_energies];

	/* Initialize variable parameters data structures of different materials */
	va_list xs_types;
	va_start(xs_types, isotopes);
	Material* material;
	float num_density;
	char* i;

	/* Create title and filename for plot */
	std::stringstream title;
	std::stringstream filename;
	title << "set title \"" << _region_name;
	title << " Macroscopic Total Cross-sections\"";
	filename << _region_name << "_macro_xs";

	/* Initialize the plot */
	gnuplot_ctrl* handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Energy (eV)");
	gnuplot_set_ylabel(handle, (char*)"Cross-section (cm^-1)");
	gnuplot_cmd(handle, (char*)title.str().c_str());
	gnuplot_cmd(handle, (char*)"set logscale xy");
	gnuplot_setstyle(handle, (char*)"lines");

	/* Loop through each material */
	for (i=isotopes; i != NULL; i=va_arg(xs_types, char*)) {

		material = _materials.at(i).second;
		num_density = _materials.at(i).first;

		/* Load xs_values vector */
		for (int j=0; j < num_energies; j++)
			xs_values[j] = material->getTotalXS(energies[j])
											* num_density * 1E-24;

		/* Plot the cross-section */
		gnuplot_saveplot(handle, (char*)filename.str().c_str());
		gnuplot_plot_xy(handle, energies, xs_values, num_energies, i);
	}

	gnuplot_close(handle);
	va_end(xs_types);

	delete [] energies;
	delete [] xs_values;

	return;
}



/**
 * This function plots the total microscopic cross-sections for
 * the materials specified with a variable argument list of character
 * arrays. The final argument must be NULL so that this function knows
 * when to stop loop over materials.
 * @param start_energy the starting energy for the plot (eV)
 * @param end_energy the ending energy for the plot (eV)
 * @param num_energies the number of energies to plot (eV)
 * @param isotopes a variable length parameter list of character arrays
 * of material types
 */
void Region1D::plotMicroscopicCrossSections(float start_energy,
		float end_energy, int num_energies, char* isotopes, ...) {

	/* Allocate memory for energies and xs values */
	float* energies = logspace(start_energy, end_energy, num_energies);
	float* xs_values = new float[num_energies];

	/* Initialize variable parameters data structures of different materials */
	va_list xs_types;
	va_start(xs_types, isotopes);
	Material* material;
	char* i;

	/* Create title and filename for plot */
	std::stringstream title;
	std::stringstream filename;
	title << "set title \"" << _region_name;
	title << " Microscopic Total Cross-sections\"";
	filename << _region_name << "_micro_xs";

	/* Initialize the plot */
	gnuplot_ctrl* handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Energy (eV)");
	gnuplot_set_ylabel(handle, (char*)"Cross-section (cm^-1)");
	gnuplot_cmd(handle, (char*)title.str().c_str());
	gnuplot_cmd(handle, (char*)"set logscale xy");
	gnuplot_setstyle(handle, (char*)"lines");

	/* Loop through each material */
	for (i=isotopes; i != NULL; i=va_arg(xs_types, char*)) {

		material = _materials.at(i).second;

		/* Load xs_values vector */
		for (int j=0; j < num_energies; j++)
			xs_values[j] = material->getTotalXS(energies[j]);

		/* Plot the cross-section */
		gnuplot_saveplot(handle, (char*)filename.str().c_str());
		gnuplot_plot_xy(handle, energies, xs_values, num_energies, i);
	}

	gnuplot_close(handle);
	va_end(xs_types);

	delete [] energies;
	delete [] xs_values;

	return;
}
