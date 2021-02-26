![image](https://github.com/patternizer/budyko-calculation/blob/master/budyko-calculation-ERA5-ma-12-0.8.png)

# budyko-calculation

Python code to perform the Budyko calculation of sea-ice boundary mean latitude as a function of CO₂ for ERA5 (+BE) and JRA-55 interpolated 1x1 degree sea-ice fields as part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org. Recalculation using reanalysis data of Fig 3 by M.I. Budyko in his 1972 paper ' The Future of Climate'.

## Contents

* `budyko-calculation.py` - python code to perform the Budyko calculation of sea-ice boundary mean latitude as a function of CO₂ for ERA5 (+BE) and JRA-55 interpolated 1x1 degree sea-ice fields

The first step is to clone the latest budyko-calculation code and step into the check out directory: 

    $ git clone https://github.com/patternizer/budyko-calculation.git
    $ cd budyko-calculation

### Using Standard Python

The code should run with the [standard CPython](https://www.python.org/downloads/) installation and was tested in a conda virtual environment running a 64-bit version of Python 3.8+.

budyko-calculation scripts can be run from sources directly, once the dependencies are satisfied (both needed python libraries and a working R installation).

Run with:

    $ python budyko-calculation.py

## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

