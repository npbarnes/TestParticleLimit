<img src="/Experimental Example Image.JPG" width=350/> 

This repo had tools for simulating an ionospheric barium injection experiment similar to CRRES[^1] or KiNET-X[^2][^3].
It includes simple neutral cloud models and barium photochemical models.

The simulation assumes test particles and small gyroradius limits, so it may not be accurate at small spatial or temporal scales or whenever electrodynamics is important.
However, it is well suited for studying how neutral cloud and photochemical models affect overall ion cloud morphology.

See `plotting_script.jl` for a simple example of usage. The result is below, feel free to compare it with the experimental image of the two KiNET-X releases above.

<img src="/Simulated Ion Cloud.png" width=350/>

My PhD dissertation[^3] has a detailed comparison between these test particle simulations, a more sophisticated (but slower) fully electromagnetic simulation, and KiNET-X data.

[^3]: Barnes, Nathan. "Obstacles to plasma flow in an ion kinetic regime: Application to a terrestrial ionospheric active plasma experiment and new horizons observations of the Pluto system." PhD diss., University of Alaska Fairbanks, 2024.

[^2]: Delamere, P. A., K. Lynch, M. Lessard, R. Pfaff, M. Larsen, D. L. Hampton, M. Conde et al. "Alfvén wave generation and electron energization in the KiNET-X sounding rocket mission." Physics of Plasmas 31, no. 11 (2024).

[^1]: Johnson, M. H., and John Kierein. "Combined release and radiation effects satellite (CRRES): Spacecraft and mission." Journal of Spacecraft and Rockets 29, no. 4 (1992): 556-563.
