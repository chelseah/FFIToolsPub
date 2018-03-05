#!/usr/bin/env python 
# 
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT) 
# 
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 


import logging
import os
from patools.util.configurable import ConfigurableObject

# create logger
logger = logging.getLogger(__name__)



class Photometry(object):

    def __init__(self):
        return

    def __call__(self, frame, xylist, outfile=''):
        return


class Fiphot(ConfigurableObject):
    config_keys = ['Fiphot']
    def __init__(self, gain=2.1, magtoflux=21., skyfit_sigma=3, skyfit_niter=4, disjoint_radius=2, apertures='2.5:4.0:3.0'):
        super(Fiphot, self).__init__()
        self.gain = float(gain)
        self.magtoflux = float(magtoflux)
        self.skyfit_sigma = int(skyfit_sigma)
        self.skyfit_niter = int(skyfit_niter)
        self.disjoint_radius = int(disjoint_radius)
        # need to be more optimally decided
        self.apertures = apertures 
        logger.debug("fiphot configuration: gain=%.2f, magtoflux=%.1f, skyfit_sigma=%d, skyfit_niter=%d, disjoint_radius=%d, apertures=%s", self.gain, self.magtoflux, self.skyfit_sigma, self.skyfit_niter, self.disjoint_radius, self.apertures) 

    def getphotfile(self, frame):
        photfile = os.path.basename(os.path.splitext(frame)[0]) + '.fiphot'
        return photfile
        
    def __call__(self, frame, xylist, outfile='', outdir='', dry=False):
        if outfile == '':
            outfile = self.getphotfile(frame)
        if outdir == '':
            outdir = os.path.dirname(xylist.name)
        # come up with some more general way to generate file num
        try:
            filenum = os.path.basename(frame).split('.')[0].split('_')[1].split('-')[0]
        except:
            filenum = 0
        
        cmdline = "fiphot --input %s --input-list %s --col-id %d --col-xy %d,%d " \
              "--gain %.2f --mag-flux %f,10000 --apertures %s " \
              "--sky-fit 'median,sigma=%d,iterations=%d' " \
              "--disjoint-radius %d --format 'ISXY,MmBbXYs' " \
              "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
              "--comment '--comment' --output %s --serial %s" \
              % (frame, xylist.name, xylist.colid, xylist.colx,
                 xylist.coly, self.gain, self.magtoflux, self.apertures,
                 self.skyfit_sigma, self.skyfit_niter,
                 self.disjoint_radius, outdir+'/'+outfile, filenum)
        logger.debug("Excute fiphot command: %s" % cmdline)
        os.system(cmdline)

