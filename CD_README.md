This readme describes the changes I have made to the version of CM1 attached to this readme.

FLEXVARS:
* Initial forcing
    * `var1`: Location of initial forcing in x direction. Implemented for warm bubble (`iinit=1`) and cold pool (`iinit=2`) forcings. For the warm bubble, sets the center of the bubble in the x direction. For the cold pool, sets the eastern edge of the cold pool.
    * `var2`: Center of bubble in y direction. Only implemented for the warm bubble forcing.
    * `var3`: Vertical location of initial forcing
        * Warm bubble: Height of the center of the bubble above the ground (m)
        * Cold pool: Depth of cold pool (m)
    * `var4`: Horizontal radius of bubble (m). Warm bubble only
    * `var5`: Vertical radius of bubble (m). Warm bubble only
    * `var6`: Max potential temperature perturbation
        * Warm bubble: Max potential temperature perturbation at bubble center (K)
        * Cold pool: Max potential temperature perturbation at surface (K)
    * `var7`: Flag controlling whether $q_v$ or the relative humidity should be kept constant in the warm bubble forcing. `var7=0` keeps $q_v$ constant, `var7=1` keeps RH constant. Warm bubble only
* Initial profile
    * `var8`: Mixed-layer $q_{v0}$ value for the Weisman-Klemp sounding. I think this was 11e-3 (ie 11 g/kg) in the default CM1 setup. The values in the original Weisman-Klemp 1982 paper were 11-14 g/kg.