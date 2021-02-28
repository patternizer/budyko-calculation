#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Quick animated GIF maker
"""

#------------------------------------------------------------------------------
# PROGRAM: make_gif.py
#------------------------------------------------------------------------------
# Version 0.1
# 28 February, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#------------------------------------------------------------------------------

import os, glob
import imageio

#----------------------------------------------------------------------------
# SETTINGS
#----------------------------------------------------------------------------

use_era5 = True
if use_era5 == True:
    reanalysisstr = 'ERA5'
else:
    reanalysisstr = 'JRA55'
gifstr = 'budyko-sea-ice-' + reanalysisstr + '.gif'

#----------------------------------------------------------------------------
# MAKE GIF
#----------------------------------------------------------------------------

images = sorted(glob.glob('sea-ice-*' + reanalysisstr + '*.png'))
var = [imageio.imread(file) for file in images]
imageio.mimsave(gifstr, var, fps = 10)

#----------------------------------------------------------------------------
# CLI --> MAKE GIF & MP4
#----------------------------------------------------------------------------

# PNG --> GIF:
# convert -delay 10 -loop 0 sea-ice*ERA5*.png budyko-sea-ice-ERA5.gif
# convert -delay 10 -loop 0 sea-ice*JRA55*.png budyko-sea-ice-JRA55.gif

# GIF --> MP4
# ffmpeg -i budyko-sea-ice-ERA5.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" budyko-sea-ice-ERA5.mp4
# ffmpeg -i budyko-sea-ice-JRA55.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" budyko-sea-ice-JRA55.mp4


# -----------------------------------------------------------------------------
print('** END')
