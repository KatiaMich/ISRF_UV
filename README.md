# ISRF_UV
This code enables to compute the interstellar radiation field at the chosen wavelength and location. Wavelength should be between 1000 and 2000 A, and location within 200 pc from the Sun for a good accuracy.

Chose a directory where you save all the python files, and a directory where you save the catalog. The catalog can be found in .......  The path to these directories has to be written down in the 'Commands' file.
The 'Commands' file call the FluxCalculator and the ReadDustMap files to compute the ISRF. Then it writes down the results in a directory called 'result' or 'resultKral', that you have to create beforehand.
For chosing the wavelength or wavelengths at which you want to compute the ISRF, you have to modify the first lines of the 'Commands' file.

I provide several python scripts:

Commands_Sun.py computes the ISRF at the location of the Sun. A pandas file with the detail of individual fluxes from all contributors in the catalog is created in the 'result' directory. The computation is done separately for Hipparcos stars, StarHorse stars with a known spectral type and StarHorse stars without a known spectral type, so 3 pandas files are created. You can chose the wavelength, you can comment the lines that compute the flux from SH stars if you want to have a fast rought estimate of the ISRF.

Commands_Kral_parallel.py computes the ISRF for each target coordinate from Kral's catalog. You can chose the wavelength. If you want to save the detail of all individual fluxes for each target coordinate, you can uncomment the line saving the pandas files. Otherwise the summed fluxes are saved in 190 numpy arrays.

Commands_spectrum.py computes the spectrum of the ISRF at the location of the Sun (you can change these coordinates to compute the spectrum somewhere else). You can chose a list of wavelengths greater than 1000 Angstrom. A pandas file is saved, s well as the figure.

Commands_Anisotropy.py computes the map of the ISRF at the location of the Sun (you can change these coordinates to compute the map somewhere else). You can chose the wavelength. A 2D array is saved, as well as the galactic map plotted on a figure.

