import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from atmosphere import time_conversion as TC
from atmosphere.fmcd import call_mcd
import time as T

LAST_RESULTS = None # list that will contain the results from the last MCD call

ALL_VALUES = []     # list that will contain all density, temperature, and pressure
TIMES = []          # list that will contain the times at which the MCD was called  

class mcd_interface:
    """
    This class is used as an interface to the Mars Climate Database.
    """
    def __init__(self, default_inputs=True, save_all_vals=False):
        # The following two arrays contain the times in solar longitude (Ls) at which different files will be loaded
        self.limiting_Ls = np.arange(0, 330.01, 30)
        self.n_species = (56, 65)   # list of the index at which the atmospheric atomic volumetric fractions are to be loaded
        self.n_wind = 25            # index of the vertical wind
        self.Mars_R = 3389.5e3      # Mars radius in [m]
        self.save_all_vals = save_all_vals
        self.species_name = ["CO2", "N2", "Ar", "CO", "O", "O2", "O3", "H", "H2"] # note: He is also accessible if required, at index 77
        self.species_weight = [44.01, 28.0134, 39.948, 28.01, 15.999, 31.999, 48, 1.00784, 1.00784] # [g/mol], molecular weights
        # Load a set of default inputs to the MCD
        if default_inputs:
            self.default_inputs()

    def default_inputs(self):
        # Load a set of default inputs for the MCD interface module
        self.zkey = 2                    # xz is thus the the altitude above the Martian zero datum (Mars geoid)
        self.xz = 125e3                  # [m], distance above the Martian geoid (125km altitude as default)
        self.xlon = 137.4                # [deg], East longitude
        self.xlat = -4.6                 # [deg], Latitude
        self.hireskey = 0                # use the lower resolution grid of 5.625x3.75 deg
        self.localtime = 12              # time of the day in hours
        self.dset = "/cala/jeremie/MCD/MCD5.3/data/"   # path to the MCD data
        self.scena = 1                   # Climatology Scenario, solar EUV average conditions; use 2 for minimum solar conditions, and 3 for maximum
        self.perturkey = 0               # do not additional perturbations
        self.seedin = 1                  # (not) used for the random number generator
        self.gwlength = 0                # (not) used for gravity waves length
        # select which external variables to return
        self.extvarkeys = np.zeros(100)
        self.extvarkeys[self.n_species[0]:self.n_species[1]+1] = 1
        self.extvarkeys[self.n_wind] = 1
        self.datekey = 1                 # 0 = Julian, 1 = Ls
        self.xdate = 0

    def call(self, Ls=None, localtime=None, lat=None, lon=None, h=None, print_results=False):
        global LAST_RESULTS, ALL_VALUES
        # Call the MCD
        if Ls is not None: self.xdate = Ls                      # solar longitude [deg], between 0 and 360
        if localtime is not None: self.localtime = localtime    # local time [hours], between 0 and 24
        if lat is not None: self.xlat = lat                     # latitude [deg], between -90 and 90
        if lon is not None: self.xlon = lon                     # longitude [deg], between -180 and 180
        if h is not None: self.xz = h                           # distance above the Martian geoid [m]
        t0 = T.time()
        # Make the actual call to the MCD (with the dataset already loaded, this takes in the order of 0.05 ms)
        self.pres,self.dens,self.temp,self.zonwind,self.merwind,self.meanvar,self.extvars,self.seedout,self.ier = \
            call_mcd(self.zkey,self.xz,self.xlon,self.xlat,self.hireskey,self.datekey,self.xdate,\
            self.localtime,self.dset,self.scena,self.perturkey,self.seedin,self.gwlength,self.extvarkeys)
        dt = T.time() - t0
        # Print the call time if the file was loaded (call time above 2 sec)
        if dt > 2:
            print("Loading the MCD took %.3f seconds." % dt)
        # Extract the vertical wind
        self.vertwind = self.extvars[self.n_wind]
        # Extract the volumetric ratio of each species in the atmosphere
        self.species_frac = []
        for n in range(*self.n_species):
            self.species_frac.append(self.extvars[n])
        # Save the atomic composition in a dictionnary
        self.species_dict = dict(zip(self.species_name, self.species_frac))
        # Save the atomic composition as a function of density
        mass = np.array(self.species_frac) * np.array(self.species_weight)
        self.species_dens = mass / mass.sum() * self.dens
        self.species_dict_dens = dict(zip(self.species_name, self.species_dens))

        LAST_RESULTS = self
        if self.save_all_vals:
            ALL_VALUES.append([self.dens, self.temp, self.pres, *self.species_frac])

        # Print the results
        if print_results:
            print("Call MCD at Ls=%.2f deg and %.2f hrs, Lat=%.2f deg, Lon=%.2f deg, and h=%.2f km" % \
                (self.xdate, self.localtime, self.xlat, self.xlon, (self.xz-self.Mars_R)/1e3))
            print("Density = %.5e [kg/m3]" % self.dens)
            print("Wind = %.5f E / %.5f N / %.5f vert [m/s]" % (self.zonwind, self.merwind, self.vertwind))
            print("Species [mol/mol] = %s" % self.species_dict)
            # input()

    def density(self, h, lon, lat, time, time_is_JD=True, JD_Tudat=True):
        """
        Return the density at a specific position and time. Used mainly interfaced to Tudat.
        Inputs:
         * h: altitude, in [m]
         * lon: longitude, in [rad]
         * lat: latitude, in [rad]
         * time: Julian date in seconds since J2000 by default. Can be changed (see optional inputs)
        Optional inputs:
         * time_is_JD: boolean specifying whether the input time is a Julian date (true), or a tuple containing the solar longitude and time of day (false)
         * JD_Tudat: boolean specifying whether the input time is a Julian date from Tudat (true, in seconds from J2000), or from the MCD (false, in days since J2023)
        Output:
         * density: float, in [kg/m3]
        """
        self.xz = h
        # Return fixed value if negative absolute altitude
        abs_h = h+2500-373
        if abs_h < 0:
            return 0.02
        self.xlat = np.rad2deg(lat)
        self.xlon = np.rad2deg(lon)
        # If the time is a Julian date, convert it to solar longitude and day of the year
        if time_is_JD:
            Ls, Ds = TC.JD_to_Ls(time, JD_Tudat=JD_Tudat)
        else:
            Ls, Ds = time
        # Convert the day of the year to hour of the day
        self.localtime = Ds % 24
        self.xdate = Ls
        # Call the MCD
        self.call()
        # Return the density
        if self.save_all_vals:
            TIMES.append(time)
        return self.dens

    def wind(self, h, lon, lat, time, time_is_JD=True, JD_Tudat=True):
        """
        Return the density at a specific position and time. Used mainly interfaced to Tudat.
        Inputs:
         * h: altitude, in [m]
         * lon: longitude, in [deg]
         * lat: latitude, in [deg]
         * time: Julian date in seconds since J2000 by default. Can be changed (see optional inputs)
        Optional inputs:
         * time_is_JD: boolean specifying whether the input time is a Julian date (true), or a tuple containing the solar longitude and time of day (false)
         * JD_Tudat: boolean specifying whether the input time is a Julian date from Tudat (true, in seconds from J2000), or from the MCD (false, in days since J2023)
        Output:
         * wind: [float, float, float], in [m/s]. It is a right-handed orthogonal vector in the vertical frame as follows:
           * x component: points to the North, positive towards North
           * y component: points to the East, positive towards East
           * z component: points to the centre of Mars, positive down
        """
        self.xz = self.Mars_R + h   # convert altitude to distance from centre of Mars
        self.xlat = lat
        self.xlon = lon
        # If the time is a Julian date, convert it to solar longitude and day of the year
        if time_is_JD:
            Ls, Ds = TC.JD_to_Ls(time, JD_Tudat=JD_Tudat)
        else:
            Ls, Ds = time
        # Convert the day of the year to hour of the day
        self.localtime = Ds % 1 * 24
        self.xdate = Ls
        # Call the MCD
        self.call()
        # Return the wind vector
        return [self.merwind, self.zonwind, self.vertwind]