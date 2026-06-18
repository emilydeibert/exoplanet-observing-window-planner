"""
Author: Emily Deibert
Last Modified: 12 December 2024
Description: Calculates timing windows for transits and pre-/post-eclipse, generates optional plots, and outputs the windows in an OT-accessible format.
"""

from astroplan import (PrimaryEclipseConstraint, is_event_observable, is_observable, is_always_observable, AtNightConstraint, AltitudeConstraint, LocalTimeConstraint, EclipsingSystem, FixedTarget, PeriodicEvent, AirmassConstraint, MoonSeparationConstraint)
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.coordinates import SkyCoord
from astroplan.plots import plot_airmass
import matplotlib.pyplot as plt
from astroplan import Observer
from astropy import units as u
from astropy.time import Time
import datetime as dt
import numpy as np
import argparse
import batman

###############################################################################################################
########################################## SET UP COMMAND LIND FLAGS ##########################################
###############################################################################################################
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--planet", help='Planet name, formatted for Exoplanet Archive (e.g., WASP-121 b)', type=str)
parser.add_argument("-t", "--type", help='Type of window to generate (transit or eclipse)', type=str)
parser.add_argument("-plot", "--plot", action=argparse.BooleanOptionalAction)
parser.add_argument("-e", "--exp_time", help='Exposure time for observations, in seconds', type=float)
parser.add_argument("-length", "--length", help='Minimum time for eclipse observation, in hours', type=float, default=5)
parser.add_argument("-a", "--airmass", help='Maximum airmass for observations', type=float, default=2)
###############################################################################################################


#######################################################################################################
########################################## SET UP PARAMETERS ##########################################
#######################################################################################################
### GEMINI SOUTH ###
gemini = Observer.at_site('gemini south', timezone='America/Santiago')

# ### LOWELL ###
# gemini = Observer.at_site('lowell', timezone='America/Phoenix')

### 2025A SEMESTER ###
start_of_semester = gemini.twilight_evening_nautical(Time('2026-08-01T00:00:00', format='isot', scale='utc'), which='next')
end_of_semester = gemini.twilight_morning_nautical(Time('2027-02-01T00:00:00', format='isot', scale='utc'), which='next')

### GHOST OBSERVATIONS ###
overhead = 22. * u.second

# ### LOWELL OBSERVATIONS ###
# overhead = 45. * u.second

transit_baseline = 60. * u.minute # on either side of the transit
window_flexibility = 20. * u.minute # flexibility on either side of the window

min_op = 0.3 # minimum orbital phase for pre-eclipse windows; maximum will be defined by right before secondary eclipse
max_op = 0.7 # similar to above, for post-eclipse windows
#######################################################################################################


######################################################################################################
########################################## DEFINE FUNCTIONS ##########################################
######################################################################################################
def orbital_phase(T0, t, P0, type):
	"""
	Calculates orbital phases given a list of observing times

	Parameters
	----------
	T0:
	t:
	P0:
	type: str
		'transit' or 'eclipse'
	"""
	
	nt = np.ceil(((t[0] - T0) / P0).value).astype(int)
	t_mid = T0 + nt * P0
	op = (t - t_mid) / P0

	if type == 'transit':
		op = op.value
	elif type == 'eclipse':
		op = op.value + 1
	
	return op


def occult_eclipse_model(t0, per, inc, rp, a, t, u_list, type, fp = 0.000204, ecc = 0., w = 90., limb_dark = 'quadratic', ts = 0.5):
	"""
	Calculates transit or eclipse depth, given an input list of orbital phases.
	NOTE: accurate limb-darkening parameters or planet flux ratio are not really necessary for the purposes of this code.
	"""
	
	params = batman.TransitParams()
	params.t0 = t0
	params.per = per
	params.inc = inc
	params.rp = rp
	params.a = a
	params.ecc = ecc
	params.w = w
	params.limb_dark = limb_dark
	params.u = u_list
	params.fp = fp
	params.t_secondary = ts

	if type == 'transit':
		m = batman.TransitModel(params, np.array(t, dtype=float))
		
	elif type == 'eclipse':
		m = batman.TransitModel(params, np.array(t, dtype=float), transittype="secondary")
		
	flux = m.light_curve(params)

	return flux


def setup_astroplan(planet):
	"""
	Sets up astroplan.
	"""

	star = planet.split(' ')[0]

	planet_table = NasaExoplanetArchive.query_object(planet, select='*', table='pscomppars')
	
	# get relevant planet properties
	epoch = Time(planet_table['pl_tranmid'].value, format='jd')
	period = planet_table['pl_orbper'].value * u.day
	transit_duration = planet_table['pl_trandur'].value * u.hour

	astroplanet = EclipsingSystem(primary_eclipse_time = epoch, orbital_period = period, duration = transit_duration)
	target = FixedTarget.from_name(star)

	return astroplanet, target, planet_table

######################################################################################################

##########################################################################################################
########################################## PLAN TRANSIT WINDOWS ##########################################
##########################################################################################################
def get_semester_windows(start_of_semester, end_of_semester, planet, n_test = 100):
	"""
	Gets all windows (not necessarily observable) in a given semester.

	Parameters
	----------
	start_of_semester: astropy Time object
		astropy time object corresponding to start of semester
	end_of_semester: astropy Time object
		astropy time object corresponding to end of semester
	planet: astroplan object
		planed defined in astroplanet

	Returns
	----------
	ing_egr: list 
		list of ingresses and egresses in the semester
	"""
	
	ing_egr = planet.next_primary_ingress_egress_time(start_of_semester, n_eclipses=n_test)
	
	while ing_egr[-1][-1] < end_of_semester:
		n_test += 25
		ing_egr = planet.next_primary_ingress_egress_time(start_of_semester, n_eclipses=n_test)
		
	ing_egr = ing_egr[np.where(ing_egr[:,1] < end_of_semester)[0]]

	return ing_egr


def get_observable_transit_windows(start_of_semester, end_of_semester, planet, star, type, site, baseline, max_airmass):
	"""
	Check which windows are actually observable, given some constraints.

	Parameters
	----------

	Returns
	----------
	observable_windows: list
		list of windows that are observable in the semester
	"""

	all_windows = get_semester_windows(start_of_semester, end_of_semester, planet)
	constraints = [AtNightConstraint.twilight_nautical(), AirmassConstraint(max_airmass), MoonSeparationConstraint(min = 30 * u.degree)]

	all_windows_plus_baseline = np.copy(all_windows)
	all_windows_plus_baseline[:,0] = all_windows_plus_baseline[:,0] - (baseline * (3/3))
	all_windows_plus_baseline[:,1] = all_windows_plus_baseline[:,1] + (baseline * (3/3))
	
	observable = is_event_observable(constraints, site, star, times_ingress_egress = all_windows_plus_baseline)[0]

	observable_windows = all_windows_plus_baseline[observable]

	return observable_windows


def check_window_flexibility(observable_windows, flexibility, site, star, max_airmass):
	"""
	Take in the observable windows, then check if there's space to add flexibility before/after the window.
	"""

	flexible_windows = []

	constraints = [AtNightConstraint.twilight_nautical(), AirmassConstraint(max_airmass), MoonSeparationConstraint(min = 30 * u.degree)]

	for window in observable_windows:

		start_flexible = window[0] - flexibility 
		end_flexible = window[1] + flexibility
		temp = is_event_observable(constraints, site, star, times_ingress_egress = np.asarray([[start_flexible, end_flexible]]))[0]
		if temp == True:
			# 15 mins pre/post are still observable
			flexible_windows.append([start_flexible, end_flexible])
		else:
			temp = is_event_observable(constraints, site, star, times_ingress_egress = np.asarray([[window[0], end_flexible]]))[0]
			if temp == True:
				# 15 mins post are still observable. DON'T ADD THE PRE YET...; the pre can be safely added as well since
				# the 15 minutes are needed for acquisition anyways
				flexible_windows.append([window[0], end_flexible])
			else:
				temp = is_event_observable(constraints, site, star, times_ingress_egress = np.asarray([[start_flexible, window[1]]]))[0]
				if temp:
					# only the 15 minutes pre is possible
					flexible_windows.append([start_flexible, window[1]])
				else:
					# no flexibility is possible
					flexible_windows.append([window[0], window[1]])

	return np.asarray(flexible_windows)


def plot_single_transit_observation(planet_params, star, site, ing_egr_list, flexible_list, baseline, overhead, exp_time):
	"""
	Generate the plots for a single window.
	"""

	start_flexible = Time(flexible_list[0], format='jd', scale='utc')
	end_flexible = Time(flexible_list[1], format='jd', scale='utc')

	starting_time = Time(ing_egr_list[0], format='jd', scale='utc')
	ending_time = Time(ing_egr_list[1], format='jd', scale='utc')

	total_obs = ((ending_time - starting_time).value * u.day).to(u.hour)
	total_flexible = ((end_flexible - start_flexible).value * u.day).to(u.hour)
	acquisition = Time(flexible_list[0] - (8 * u.minute), format='jd', scale='utc')
	
	cadence = exp_time + overhead
	num_obs = int((total_obs / cadence).decompose())

	jds_observation = np.linspace(starting_time, ending_time, num_obs)

	jds_observation = Time(jds_observation, format='jd', scale='utc', precision=9)

	jds_acquisition = np.arange(acquisition.value, start_flexible.value, cadence.to(u.day).value)
	jds_acquisition = Time(jds_acquisition, format='jd', scale='utc', precision=9)

	jds_flexible = Time(np.linspace(start_flexible, end_flexible, int((total_flexible / cadence).decompose())),
		format='jd', scale='utc', precision=9)

	if planet_params['pl_name'] == 'WASP-121 b':
		ratdor = 3.754
	elif planet_params['pl_name'] == 'WASP-103 b':
		ratdor = 3.013
	else:
		ratdor = planet_params['pl_ratdor']

	#print('OKAY TO HERE')

	op = orbital_phase(Time(planet_params['pl_tranmid'].value, format='jd'), jds_observation, 
					   planet_params['pl_orbper'].value * u.day, type='transit')

	#print('OKAY TO HERE')

	# print("a =", ratdor, type(ratdor))
	# print("rp =", planet_params['pl_ratror'], type(planet_params['pl_ratror']))
	# print("inc =", planet_params['pl_orbincl'].value)

	occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), 
		float(np.asarray(planet_params['pl_ratror']).item()), 
		float(np.asarray(ratdor).item()),
		op, [0.25, 0.25], type = 'transit')

	print('OKAY TO HERE')


	# op_flexible = orbital_phase(Time(planet_params['pl_tranmid'].value, format='jd'), jds_flexible, 
	# 				   planet_params['pl_orbper'].value * u.day, type='transit')
	# occult_flexible = occult_eclipse_model(0, 1, planet_params['pl_orbincl'].value, planet_params['pl_ratror'], ratdor,
	# 							 op_flexible, [0.25, 0.25], type = 'transit')
	# op_acquisition = orbital_phase(Time(planet_params['pl_tranmid'].value, format='jd'), jds_acquisition, 
	# 				   planet_params['pl_orbper'].value * u.day, type='transit')

	op_flexible = orbital_phase(Time(planet_params['pl_tranmid'].value, format='jd'), jds_flexible,
		planet_params['pl_orbper'].value * u.day, type='transit')

	occult_flexible = occult_eclipse_model(
		0, 1,
		float(np.asarray(planet_params['pl_orbincl'].value).item()),
		float(np.asarray(planet_params['pl_ratror']).item()),
		float(np.asarray(ratdor).item()),
		op_flexible,
		[0.25, 0.25],
		type='transit')

	op_acquisition = orbital_phase(Time(planet_params['pl_tranmid'].value, format='jd'), jds_acquisition,
		planet_params['pl_orbper'].value * u.day, type='transit')

	if np.isnan(occult).all():
		if planet_params['pl_name'] == 'WASP-100 b':
			ratdor = 4.97
		elif planet_params['pl_name'] == 'WASP-78 b':
			ratdor = 3.52
		elif planet_params['pl_name'] == 'WASP-103 b':
			ratdor = 3.013
		elif planet_params['pl_name'] == 'WASP-121 b':
			ratdor = 3.754
		else:
			print('You have an error with a planet parameter...')
		occult = occult_eclipse_model(0, 1, planet_params['pl_orbincl'].value, planet_params['pl_ratror'], ratdor, op, [0.25, 0.25], type = 'transit')
		occult_flexible = occult_eclipse_model(0, 1, planet_params['pl_orbincl'].value, planet_params['pl_ratror'], ratdor, op_flexible, [0.25, 0.25], type = 'transit')

	airmass = site.altaz(jds_observation, star).secz
	airmass_flexible = site.altaz(jds_flexible, star).secz
	
	fig = plt.figure(figsize=(10, 4))
	ax_transit = fig.add_subplot(121) 
	ax_transit.plot(op, occult, '.', zorder=10)
	ax_transit.plot(op_flexible, occult_flexible, 'o', zorder=1)
	ax_transit.plot(op_acquisition, np.ones(len(op_acquisition)), 'o', zorder=1, color='red')
	ax_transit.axvspan(op_acquisition[0], op_flexible[0], color='lightgrey', alpha=0.5)

	ax_airmass = fig.add_subplot(122)
	ax_airmass.plot(np.asarray(jds_observation.value), np.asarray(airmass), '.', zorder=10)

	ax_airmass.plot(np.asarray(jds_flexible.value), np.asarray(airmass_flexible), 'o', zorder=2)

	ax_airmass.axvline(site.twilight_evening_nautical(starting_time, which='previous').value, color='k', linestyle=':')
	ax_airmass.axvline(site.twilight_morning_nautical(starting_time, which='next').value, color='k', linestyle=':')

	ax_airmass.axvspan(acquisition.value, jds_flexible[0].value, color='lightgrey', alpha=0.5)

	ax_transit.set_xlabel('Orbital Phase')
	ax_airmass.set_xlabel('Julian Date')
	ax_airmass.set_ylabel('Airmass')
	ax_airmass.invert_yaxis()
	ax_airmass.set_title(starting_time.isot.split('T')[0])
	
	plt.show()

	print('UT Start of Night: ', starting_time.isot.split('T')[0])
	print('Number of Observations: ', num_obs)
	print('Number of OOT Observations: ', len(np.where(occult == 1.)[0]))
	print(' ')
		
	return


##########################################################################################################
########################################## PLAN ECLIPSE WINDOWS ##########################################
##########################################################################################################
def find_eclipse_windows(start_time, end_time, planet_params, star, phase_lim, max_airmass, eclipseLength, site = gemini):

	eclipse_windows = []

	# find the orbital phase corresponding to the start and end of the secondary eclipse
	op_test = np.linspace(0.25, 0.75, 1000)

	if planet_params['pl_name'] == 'WASP-103 b':
		ratdor = 3.013
		occult = occult_eclipse_model(0, 1, planet_params['pl_orbincl'].value, planet_params['pl_ratror'].value, ratdor, op_test, [0.25, 0.25], type = 'eclipse')
	elif planet_params['pl_name'] == 'WASP-100 b':
		ratdor = 4.97
		occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), float(np.asarray(planet_params['pl_ratror']).item()), ratdor, op_test, [0.25, 0.25], type = 'eclipse')
	elif planet_params['pl_name'] == 'WASP-78 b':
		ratdor = 3.52
		occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), float(np.asarray(planet_params['pl_ratror'].value).item()), ratdor, op_test, [0.25, 0.25], type = 'eclipse')
	elif planet_params['pl_name'] == 'WASP-121 b':
		ratdor = 3.754
		occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), float(np.asarray(planet_params['pl_ratror']).item()), float(np.asarray(ratdor).item()), op_test, [0.25, 0.25], type = 'eclipse')
	else:
		occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), float(np.asarray(planet_params['pl_ratror']).item()), float(np.asarray(planet_params['pl_ratdor']).item()),
								 op_test, [0.25, 0.25], type = 'eclipse')
	if np.isnan(occult).all():
		if planet_params['pl_name'] == 'WASP-100 b':
			ratdor = 4.97
		elif planet_params['pl_name'] == 'WASP-78 b':
			ratdor = 3.52
		elif planet_params['pl_name'] == 'WASP-103 b':
			ratdor = 3.013
		elif planet_params['pl_name'] == 'WASP-121 b':
			ratdor = 3.754
		else:
			print('You have an error with a planet parameter...')
		occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), float(np.asarray(planet_params['pl_ratror']).item()), float(np.asarray(ratdor).item()), op_test, [0.25, 0.25], type = 'eclipse')
	
	eclipse_phases = np.where(occult != occult[0])[0]
	#print(eclipse_phases)
	if len(eclipse_phases) > 1:
		pre_eclipse = op_test[eclipse_phases[0]]
		post_eclipse = op_test[eclipse_phases[-1]]
		#print(pre_eclipse)
		#print(post_eclipse)
	else:
		for i in range(10):
			print('WARNING WARNING')
		#pre_eclipse = 0.47
		#post_eclipse = 0.53


	if phase_lim > 0.5:
		# you are doing post-eclipse phases
		lims = [post_eclipse, phase_lim]

	elif phase_lim < 0.5:
		# you are doing pre-eclipse phases
		lims = [phase_lim, pre_eclipse]

	start = start_time - 1. * u.day

	while start < end_time:
		start = site.twilight_evening_nautical(start + 12 * u.hour, which='next')
		end = site.twilight_morning_nautical(start, which='next')

		op_day_temp = orbital_phase(Time(planet_params['pl_tranmid'][0].value, format='jd'), Time([start.value, end.value], format='jd', scale='utc'), planet_params['pl_orbper'].value * u.day, type = 'transit')
		if op_day_temp[0] < 0:
			op_day_temp += 1

		if (lims[0] > min(op_day_temp) or lims[0] < max(op_day_temp)) or (lims[1] > min(op_day_temp) or lims[1] < max(op_day_temp)):
			# the phases you're trying to observe are present on this date

			# now need to check if at least 5 hours are observable...
			jds_night = Time(np.arange(start.value, end.value, ((1 * u.minute).to(u.day)).value), format='jd', scale='utc')
			airmass = site.altaz(jds_night, star).secz
			acceptable_airmass = np.where((airmass < max_airmass) & (airmass > 1))[0]
			jds_above_altitude = jds_night[acceptable_airmass]
			if len(jds_above_altitude) > 1:
				time_above_altitude = (jds_above_altitude[-1] - jds_above_altitude[0]).to(u.hour)

				ops_above_altitude = orbital_phase(Time(planet_params['pl_tranmid'][0].value, format='jd'), jds_above_altitude, planet_params['pl_orbper'].value * u.day, type = 'transit')
				if ops_above_altitude[0] < 0:
					ops_above_altitude += 1

				if time_above_altitude > eclipseLength:
					# the target is visible at a good airmass for more than 5 hours

					# check obs
					earliest_obs = max([lims[0], ops_above_altitude[0]])
					latest_obs = min([lims[-1], ops_above_altitude[-1]])

					# find the orbital phases corresponding to this...?
					suitable_obs = np.where((ops_above_altitude > earliest_obs) & (ops_above_altitude < latest_obs))[0]
					# find the jds corresponding to this...?
					suitable_jds = jds_above_altitude[suitable_obs]
					if len(suitable_jds) > 1:
						suitable_observation_length = (suitable_jds[-1] - suitable_jds[0]).to(u.hour)

						if suitable_observation_length > eclipseLength:

							start_of_window = suitable_jds[0]
							window_length = suitable_observation_length

							if window_length < eclipseLength + (8 * u.minute):
								start_of_window = start_of_window - (8 * u.minute)
								window_length = window_length + (8 * u.minute)

							eclipse_windows.append([start_of_window.isot, window_length])

	# now do a double check to make sure those windows are observable...
	observable_eclipse_windows = []
	constraints = [MoonSeparationConstraint(min = 30 * u.degree)]
	for idx, window in enumerate(eclipse_windows):
		start_time = Time(window[0], format='isot', scale='utc')
		end_time = start_time + window[1]
		check = is_always_observable(constraints, site, star, time_range = [start_time, end_time])
		if check == True:
			observable_eclipse_windows.append(window)

	return observable_eclipse_windows


def plot_single_eclipse_observation(planet_params, star, site, eclipse_window, overhead, exp_time, eclipseLength):
	"""
	Generate the plots for a single window.
	"""

	total_obs = eclipse_window[1]

	starting_time = Time(eclipse_window[0], format='isot', scale='utc')
	ending_time = starting_time + total_obs

	starting_time = Time(starting_time.jd, format='jd', scale='utc')
	ending_time = Time(ending_time.jd, format='jd', scale='utc')

	
	cadence = exp_time + overhead
	num_obs = int((total_obs / cadence).decompose())
	real_num_obs = int((eclipseLength / cadence).decompose())

	jds_observation = np.linspace(starting_time, ending_time, num_obs)
	jds_observation = Time(jds_observation, format='jd', scale='utc', precision=9)

	op = orbital_phase(Time(planet_params['pl_tranmid'].value, format='jd'), jds_observation, 
					   planet_params['pl_orbper'].value * u.day, type='eclipse')
	occult = occult_eclipse_model(0, 1, float(np.asarray(planet_params['pl_orbincl'].value).item()), float(np.asarray(planet_params['pl_ratror']).item()), 
		float(np.asarray(planet_params['pl_ratdor']).item()),
		op, [0.25, 0.25], type = 'eclipse')

	airmass = site.altaz(jds_observation, star).secz
	
	fig = plt.figure(figsize=(10, 4))
	ax_transit = fig.add_subplot(121) 
	ax_transit.plot(op, occult, '.', zorder=10)

	ax_airmass = fig.add_subplot(122)
	ax_airmass.plot(np.asarray(jds_observation.value), np.asarray(airmass), '.', zorder=10)

	ax_airmass.axvline(site.twilight_evening_nautical(starting_time, which='previous').value, color='k', linestyle=':')
	ax_airmass.axvline(site.twilight_morning_nautical(starting_time, which='next').value, color='k', linestyle=':')

	ax_transit.set_xlabel('Orbital Phase')
	ax_airmass.set_xlabel('Julian Date')
	ax_airmass.set_ylabel('Airmass')
	ax_airmass.invert_yaxis()
	ax_airmass.set_title(starting_time.isot.split('T')[0])
	
	plt.show()

	print('UT Start of Night: ', starting_time.isot.split('T')[0])
	print('Number of Observations: ', real_num_obs)
	print('Window Length: ', total_obs)
	print(' ')
		
	return


###########################################################################################################
############################################ FORMAT FOR THE OT ############################################
###########################################################################################################
def format4gemini(t1, total_obs):
	""" 
	Reformats into the format needed for the OT.

	Parameters
	----------
	t1: astropy.Time
		astropy time object for the start of the observation
	total_obs: astropy.unit
		total length of observing window, in hours

	Returns
	----------
	gemstr: str
		Gemini-compatible string
	"""
	t1 = t1
	date = t1.split('T')[0]
	start_time = t1.split('T')[-1].split('.')[0]
	hour = str(int(np.floor(total_obs).value)).zfill(2)
	minute = str(int((total_obs - np.floor(total_obs)).to(u.minute).value)).zfill(2)

	gemstr = date+' '+start_time+' '+hour+':'+minute

	return gemstr


def format_all_transits(windows, fname, save_file = True):
	"""
	Format all transit windows and save into an OT-compatible text file.
	"""

	formatted_windows = []

	for idx, window in enumerate(windows):
		total_obs = ((window[1] - window[0]).value * u.day).to(u.hour)
		t1 = (window[0] - (8 * u.minute)).isot#window[0].isot#
		
		gemstr = format4gemini(t1, total_obs)
		formatted_windows.append(gemstr)

	if save_file:
		with open(fname+'.tw', 'w') as f:
			for line in formatted_windows:
				f.write(f"{line}\n")
	
	return

def format_all_eclipses(windows, fname, save_file = True):
	"""
	Format all transit windows and save into an OT-compatible text file.
	"""

	formatted_windows = []

	for idx, window in enumerate(windows):
		
		gemstr = format4gemini(window[0], window[1])
		formatted_windows.append(gemstr)

	if save_file:
		with open(fname+'.tw', 'w') as f:
			for line in formatted_windows:
				f.write(f"{line}\n")
	
	return
######################################################################################################


##############################################################################################
########################################## RUN CODE ##########################################
##############################################################################################
def main():

	args = parser.parse_args()

	planet, star, planet_params = setup_astroplan(args.planet)

	if args.type == 'transit':
		planet_observable_windows = get_observable_transit_windows(start_of_semester, end_of_semester, planet, star, type=args.type, site=gemini, baseline = transit_baseline, max_airmass = args.airmass)
		planet_flexible_windows = check_window_flexibility(planet_observable_windows, window_flexibility, gemini, star, max_airmass = args.airmass)

	elif args.type == 'eclipse':
		pre_eclipse_windows = find_eclipse_windows(start_of_semester, end_of_semester, planet_params, star, min_op, args.airmass, args.length * u.hour, site = gemini)
		post_eclipse_windows = find_eclipse_windows(start_of_semester, end_of_semester, planet_params, star, max_op, args.airmass, args.length * u.hour, site = gemini)


	savename = args.planet.replace(' ', '')+'_'+args.type

	if args.type == 'transit':
		format_all_transits(planet_flexible_windows, savename)
	elif args.type == 'eclipse':
		format_all_eclipses(pre_eclipse_windows, savename+'_pre')
		format_all_eclipses(post_eclipse_windows, savename+'_post')

	if args.plot == True:
		if args.type == 'transit':
			for idx, window in enumerate(planet_observable_windows):
				plot_single_transit_observation(planet_params, star, gemini, window, planet_flexible_windows[idx], transit_baseline, overhead, args.exp_time * u.second)
		elif args.type == 'eclipse':
			for idx, window in enumerate(pre_eclipse_windows):
				plot_single_eclipse_observation(planet_params, star, gemini, window, overhead, args.exp_time * u.second, args.length * u.hour)
			for idx, window in enumerate(post_eclipse_windows):
				plot_single_eclipse_observation(planet_params, star, gemini, window, overhead, args.exp_time * u.second, args.length * u.hour)

	return


if __name__ == '__main__':
	main()
##############################################################################################

















