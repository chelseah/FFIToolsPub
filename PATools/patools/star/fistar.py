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



class SourceFinder(ConfigurableObject):

    def __init__(self, threshold=200, sdkflag=True):
        super(SourceFinder, self).__init__()
        
        self.threshold = threshold
        self.sdkflag = sdkflag

        logger.debug("in baseclass, threshold=%f, sdkflag=%d", self.threshold, self.sdkflag) 
    def __call__(self, infile, outfile=''):
        raise NotImplementedError


class Fistar(ConfigurableObject):
    # FIXME: configurableObject does not work for inherient class
    config_keys = ['Fistar']
    def __init__(self, threshold=200, sdkflag=False):
        super(Fistar, self).__init__()
        self.threshold = int(threshold)
        self.sdkflag = sdkflag
        self.ext = '.fistar'
        
        logger.debug("Configure fistar: threshold=%f, sdkflag=%d", self.threshold, self.sdkflag) 
    def get_outfile(self, infile):
        outfile = os.path.basename(os.path.splitext(infile)[0])+self.ext
        return outfile
    def __call__(self, infile, outfile='', out_dir=''):
        if outfile == '':
            outfile = self.get_outfile(infile)
        if out_dir == '':
            out_dir = os.path.dirname(infile)
        cmdline = 'fistar %s -o %s -s flux ' \
                  '--model elliptic --flux-threshold %d --algorithm uplink ' \
                  '--iterations symmetric=2,general=1 ' \
                  '--format id,x,y,bg,amp,s,d,k,flux,s/n --comment ' \
                  % (infile, out_dir+'/'+outfile, self.threshold)
        if not self.sdkflag:
            cmdline+='--only-candidates '
        logger.debug("Excuting Fistar commands: %s", cmdline) 
        os.system(cmdline) 
        return
