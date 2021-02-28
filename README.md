![image](https://github.com/patternizer/budyko-calculation/blob/master/sea-ice-1969-03-01-ERA5-ma-12-0.1.png)
![image](https://github.com/patternizer/budyko-calculation/blob/master/budyko-calculation-ERA5-ma-12-0.1.png)

# budyko-calculation

Python code to perform the Budyko calculation of sea-ice boundary median latitude as a function of CO₂ for ERA5 (+BE) and JRA-55 interpolated 1x1 degree sea-ice fields as part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org. This is a recalculation using reanalysis of Fig 3 by M.I. Budyko in his 1972 paper 'The Future of Climate'.

## Contents

* `budyko-calculation.py` - python code to perform the Budyko calculation of sea-ice boundary median latitude as a function of CO₂ for ERA5 (+BE) and JRA-55 interpolated 1x1 degree sea-ice fields
* `make_gif.py` - python code to generate an animated GIF from monthly polar plots of minimum and latitudinal median of the sea ice extent boundary

The first step is to clone the latest budyko-calculation code and step into the check out directory: 

    $ git clone https://github.com/patternizer/budyko-calculation.git
    $ cd budyko-calculation

### Using Standard Python

The code should run with the [standard CPython](https://www.python.org/downloads/) installation and was tested in a conda virtual environment running a 64-bit version of Python 3.8+.

budyko-calculation scripts can be run from sources directly, once the dependencies are satisfied. In addition to the python libraries imported in budyko-calculation.py you will need to have a working version of LaTeX and Ghostscript installed. For generation of the animated GIF you will need to install the python library imageio. Change the flags in the settings as desired.

Run with:

    $ python budyko-calculation.py
    $ python make_gif.py

## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

