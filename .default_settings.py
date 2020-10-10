parameters = {


    # ================================================================================================
    # SECTION 1: PROGRAM SETTINGS
    # ================================================================================================

    #Pickle allows to save and load classes easily and quicly.
    #Set to False only if not installed or supported, but some functionalities might be compromised.
    'use_pickle'    :   True,

    #Ask mode (n-tot, n-g) when looking for isotopes in the finder.
    'ask_mode'      :   False,

    #Number of mesh point to the left and right of peak edges when plotting peaks.
    'prplot'        :   20,

    #In bar plots, bar scaling. This factor is the ratio (highest bar):(highest sample spectrum peak).
    'bar_scaling'   :   0.8,

    #In bar plots, bar width (us).
    'bar_width'     :   0.7,

    #

    # ================================================================================================
    # SECTION 2: EXPERIMENTAL VALUES
    # ================================================================================================

    #Either 'n-g' or 'n-tot'. It is the energy-ToF conversion mode used if unspecified.
    'default_mode'  :   'n-g',

    #Either 'n-g' or 'n-tot'. It is the default mode for sample importing.
    'default_smode' :   'n-g',

    #L0 lengths (m). L0_g and L0_t provide different 
    'L0_g'          :   22.804,
    'L0_t'          :   23.404,#22.884,

    # ================================================================================================
    # SECTION 3: PEAK DETECTION VALUES
    # ================================================================================================

    #True if horizontal threshhold in peak detection is given in ToF, False if it is given in energy.
    'thr_in_tof'    :   True,

    #If thr_in_tof is True: minimum and maximum value of ToF (us) in peak detection
    #                 False: minimum and maximum value of energy (eV) in peak detection
    'thr_min'       :   5,
    'thr_max'       :   2000,

    #Minimum XS (b) in peak detection (in general).
    'crs_min'       :   80,

    #Exceptions in Minimum XS (b):
    #Make sure to enter the element symbol between single quotes, colon and its own crs_min,
    #followed by a comma.
    #If different values for XS in the n-tot and n-g modes are to be provided, it is possible to do so
    #by means of a tuple (XS(n-tot), XS(n-g).
   
   ##e.g.
   #'crs_exc'       :   {
   #'Fe'    : (20, 5),
   #'Sn'    : 70.4,
   #'Cu'    : 55,
   #},
   
    'crs_exc'       :   {
    'Cu'    :   (20, 5),
    'Fe'    :   (20, 5),

    },

    # ================================================================================================
    # SECTION 4: COMPUTATIONAL PARAMETERS
    # ================================================================================================

    #Maximum value prange can get.
    'prangemax'     :   500,

    #Maximum allowed slope (abolute value) at outermost left side of spectra.
    #i.e. starting from left, everything will be set to 0 until the (unsigned) slope reaches this value.
    'maxleftslope'  :   3000,

    #Maximum allowed slope outside the peaks, as in, far away from them.
    'maxouterslope' :   10,

    #Peak edges is set when its slope has fallen down to this fraction of the one nearby the peak summit.
    'slopedrop'     :   .1,

    #Density of boxes (box/b) for slope computation.
    'dboxes'        :   100,

    #Smoothing iterations on slope derivative for computations.
    'itersmooth'    :   1,

    #Smoothing iterations on sample peak detection.
    'itersmoothsamp':   0,

    #Strip-peaks iterations for background fitting in sample imports.
    'iterspeaks'    :   4,

    #Number of coefficients in smaple background fitting, i.e., polynomial order + 1
    'fitting_coeff' :   8,

    #Default tolerance value (us) for finding nearby peaks in pmatch function.
    'max_match'     :   3.5,

    # ================================================================================================
    # SECTION 0: FUNDAMENTAL CONSTANTS
    # ================================================================================================
    # Hopefully have stayed the same since this code was written, so you probably don't want to touch
    # these much. Precision can be modified. Please handle with care, though.

    #Neutron mass (kg)
    'mn'            :   1.68e-27,

    #Fundamental charge (C)
    'e'             :   1.60e-19,

}
