#
# Copyright (C) 2015 - Zach Berta-Thompson <zkbt@mit.edu> (MIT License)
#               2016 - Chelsea Huang (GPLv3 License)
#               2016 - Massachusetts Institute of Technology (MIT)
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
#
"""
Base class for light curves and light curve generation classes.
"""

import ConfigParser
import logging
import numpy as np
import os

logger = logging.getLogger(__name__)

class FilebasedLightCurve(object):
    """Light curve based on user input file"""

    def __init__(self, name=''):
        self.name = name
        self.data = {'jd': [], # Julian Date
                     'mag': [], # Raw light curve data
                     'magerr': [], #Processed light curve data
                     } 
        self.labels = {'jd': 'BJD', 'mag': 'Relatvie Magnitude',
                'magerr': 'relative magnitude error'}
        return
    
    def load_from_file(self, label='all'):
        """
        Load data from a file.  Allow user to specify load all or specific
        columns.
        """
        raise NotImplementedError

    def write_to_file(self, outfile):
        # Output lightcurve to a file as designed.
        raise NotImplementedError

    def check_data(self, label):
        if len(self.data[label]) == 0:
            try:
                self.load_from_file(label=label)
            except (IndexError, KeyError):
                raise IndexError, "Column %s is required, but cannot found in file %s" % (label, self.name)


    def plot_time_serie(self, show=True, label='mag'):
        # only works for the sub
        import matplotlib
        from matplotlib import pyplot as plt
        self.check_data('jd')
        self.check_data(label)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if self.data['jd'][0] > 2454000:
            ax.plot(self.data['jd'] - 2454000, self.data[label], 'k.')
            ax.set_xlabel("BJD-2454000")
        else:
            ax.plot(self.data['jd'], self.data[label], 'k.')
            ax.set_xlabel("BJD")
        
        ax.set_ylabel(self.labels[label])
        yformatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(yformatter)
        ax.grid(True)
        if show:
            plt.show()
        else:
            pngfile = os.path.splitext(self.name)[0] + '.png'
            plt.savefig(pngfile)



class SimpleLightCurve(FilebasedLightCurve):
    def __init__(self, name=''):
        super(SimpleLightCurve, self).__init__(name)
        self.cols={}    
        
    def set_cols(self, **kwargs):
        for key in kwargs:
            self.cols[key.lstrip('col')] = int(kwargs[key])
