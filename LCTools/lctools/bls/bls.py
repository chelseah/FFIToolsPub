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
# but ITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 


import os
import time
import commands
import logging
import sys

from lctools.util.configurable import ConfigurableObject

# create logger
logger = logging.getLogger(__name__)

class VtBls(ConfigurableObject):
    LENGTH = 20
    config_keys = ["BLS"]
    """
    A python class that configure and run the vartools BLS algorithm.
    """
    def __init__(self, f0=0.0625, f1=2.0, fn=100000, qmin=0.008, qmax=0.08, nbin=200, peaknum=3, outdir='', analfile='.blsanal'):
        super(VtBls, self).__init__()
        self.f0 = float(f0)
        self.f1 = float(f1)
        self.fn = int(fn)
        self.qmin = float(qmin)
        self.qmax = float(qmax)
        self.nbin = int(nbin)
        self.peaknum = int(peaknum)
        self.outdir = outdir
        self.analfile = analfile
        logger.debug("readfrom configure file, f0=%f, f1=%f, fn=%d, qmin=%f, qmax=%f, nbin=%d, peaknum=%d, outdir=%s, analfile=%s", self.f0, self.f1, self.fn, self.qmin, self.qmax, self.nbin, self.peaknum, self.outdir, self.analfile)
        return
    
    def __call__(self, lcfile, replace=False):
        if self.outdir == '':
            outdir = os.path.dirname(lcfile.name) + '/'
        else:
            outdir = self.outdir
        blspath = os.path.dirname(outdir)
        
        if blspath == "":
            blspath = "./"
        modelpath = blspath
        phasepath = blspath
        if (not os.path.exists(lcfile.name)):
            logger.warning("Input file %s do not exists, do nothing", lcfile.name)
            return [None, None, None]

        else:
            cmdline = "vartools -i %s -header " \
                      "-readformat 1 %d %d %d " \
                      "-BLS q %f %f %f %f %d %d %d %d 1 %s 1 %s 0 " \
                      "fittrap nobinnedrms ophcurve %s -0.1 1.1 0.001" \
                      % (lcfile.name, lcfile.cols['jd'], lcfile.cols['mag'],lcfile.cols['magerr'],\
                         self.qmin, self.qmax, 1./self.f1, 1./self.f0, self.fn,\
                         self.nbin, 0, self.peaknum, blspath, modelpath, phasepath)
            logger.info("excuting command line %s", cmdline)
            # print cmdline
            status, output = commands.getstatusoutput(cmdline)
            # print output
            header, result = output.split('\n')
            newheader = header.split()[1:VtBls.LENGTH+1]
            newlines = ""
            paramdics = []
            for i in xrange(self.peaknum):
                paramdics.append({})
                newline = result.split()[1 + i * VtBls.LENGTH:(1 + VtBls.LENGTH * (i + 1))]
                for key, val in zip(newheader, newline):
                    paramdics[i][key.lstrip("BLS_").rstrip("_1_0")] = eval(val)
            blsanal = paramdics
            blsfile = os.path.join(blspath, os.path.basename(lcfile.name)+'.bls')
            blsmodelfile = os.path.join(blspath, os.path.basename(lcfile.name)+'.bls.model')
            return [blsanal, blsfile, blsmodelfile]
        

