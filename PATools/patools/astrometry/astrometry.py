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



import math
import numpy as np
import scipy as sp
from scipy import linalg
import os
import logging
from patools.util.configurable import ConfigurableObject
from patools.util.dataio import readcolumn 
from patools.phot import Fiphot
from starlist import Starlist
import astropy
import astropy.wcs as wcs
import sklearn
from sklearn.neighbors import NearestNeighbors
# create logger
logger = logging.getLogger(__name__)


def ls_polynxy(x,y,f,n):
    xi = sp.ones(len(x))
    A=(xi)[:,np.newaxis]

    for i in range(1,n+1):
      for j in range(0,i+1):
        A=np.hstack((A,(1./math.factorial(i-j)*1./math.factorial(j)*x**(i-j)*y**j)[:,np.newaxis]))
    c,resid,rank,sigma = sp.linalg.lstsq(A,f)
    return c

def poly_ordern_eval(x,y,m):
    order = int(np.sqrt(2*(len(m)-1)))-1
    #print len(m),order
    z=np.zeros(len(x))
    k=1
    for i in range(1,order+1):
      for j in range(0,i+1):
        #z+=m[k]*math.factorial(i-j)*math.factorial(j)*x**(i-j)*y**j
        z+=m[k]*1./math.factorial(i-j)*1./math.factorial(j)*x**(i-j)*y**j
        k+=1
    z+=m[0]
    return z








class Astrometry(object):

    def __init__(self):
        return

    def __call__(self, infile, outfile=''):
        return

# FIXME: configurable object problem, see fistar.py

class Anet(ConfigurableObject):
    config_keys = ['Anet', 'Catalog']
    def __init__(self, ra=None, dec=None, q=0.01, tweak=2, order=4, radius=10, catfile=None, maxdistance=1):
        self.ra = float(ra)
        self.dec = float(dec)
        self.q = float(q)
        self.tweak = int(tweak)
        self.radius = int(radius)
        self.order = int(order)
        self.maxdistance = int(maxdistance)
        self.catfile = catfile
        # FIXME: correct the debug information
        # logger.debug("order=%d, maxdistance=%d, unitarity=%f, catfile=%s", self.order, self.maxdistance, self.unitarity, self.catfile)
         
        self.catalog = Starlist(self.catfile)
        self.cfg = None

    def get_wcsfile(self, infile): 
        wcsfile = os.path.basename(os.path.splitext(infile)[0]) + '.wcs'
        return wcsfile

    def get_transfile(self, infile): 
        transfile = os.path.basename(os.path.splitext(infile)[0]) + '.trans'
        return transfile

    def get_xylsfile(self, infile):
        xylsfile = os.path.basename(os.path.splitext(infile)[0])+'.xyls'

        return xylsfile

    def polyfit(self, match_arr,transfile,xshift=0,yshift=0):
        i_x=match_arr[:,0]
        i_y=match_arr[:,1]
        r_x=match_arr[:,2]-xshift
        r_y=match_arr[:,3]-yshift
        flag=True
        reslim=0.05
        oldres=0
        #print np.median(i_x-r_x)
        #print np.median(i_y-r_y)
        #plt.plot(i_x,i_x-r_x,'b.')
        #plt.plot(i_y,i_y-r_y,'r.')
        #plt.show()
        while flag:
            m1 = ls_polynxy(r_x, r_y, i_x, self.order)
            m2 = ls_polynxy(r_x, r_y, i_y, self.order)

            i_xe = poly_ordern_eval(r_x,r_y,m1)
            i_ye = poly_ordern_eval(r_x,r_y,m2)

            res = ((i_x-i_xe)**2.+(i_y-i_ye)**2.)**0.5
            resstd=(np.sum(res)/len(res))**0.5
            print resstd
            if resstd<reslim or abs(resstd-oldres)/resstd<0.05:
                break
            else:
                keepindex=res<8*np.median(res)
                i_x=i_x[keepindex]
                i_y=i_y[keepindex]
                r_x=r_x[keepindex]
                r_y=r_y[keepindex]
                oldres=resstd
        print np.median(res)
        #plt.plot(i_x,i_x-i_xe,'b.')
        #plt.plot(i_y,i_y-i_ye,'r.')
        #plt.show()
        self.output_trans(transfile, m1, m2, res)
        return
   
    def output_trans(self, transfile,m1,m2,res):
        fout = open(transfile,mode="w")
        fout.write("# Type: polynomial of order=%d (number of coefficients: %d)\n"% (self.order,len(m1)))
        fout.write("type = polynomial\n")
        fout.write("order = %d\n" % (self.order))
        fout.write("# Initial transformation of (x_img,y_img):\n")
        fout.write("offset = 0, 0\n")
        fout.write("scale = 1\n")
        fout.write("basisshift = 0, 0\n")
        fout.write("# Coefficients of the x fit:\n")
        fout.write("dxfit=")
        for i in xrange(len(m1)-1):
            fout.write("%14.7g," % m1[i])
        fout.write("%14.7g\n" % m1[-1])
        fout.write("# Coefficients of the y fit:\n")
        fout.write("dyfit=")
        for i in xrange(len(m2)-1):
            fout.write("%14.7g," % m2[i])
        fout.write("%14.7g\n" % m2[-1])
        if len(res)<500:
            fout.write("# fit with %d stars,failed\n" % len(res))
        else:
            fout.write("# fit with %d stars\n" % len(res))
        fout.write("# Residual: %f %f\n" % ((np.sum(res)/len(res))**0.5,np.median(res)))
        fout.close()
        return
   
    def read_transfile(self, transfile):
        with open(transfile, mode='r') as fin:
            for line in fin.readlines():
                if line.startswith("order"):
                    self.order = int(line.lstrip("order = "))
                if line.startswith("dxfit"):
                    coeffs = line.lstrip("dxfit=").split(",")
                    m1 = [float(c) for c in coeffs]
                if line.startswith("dyfit"):
                    coeffs = line.lstrip("dyfit=").split(",")
                    m2 = [float(c) for c in coeffs]
        
        return [m1, m2]
    def kd_match(self, src, dist):
        
        nbrs=NearestNeighbors(n_neighbors=1,algorithm='auto').fit(dist)
        #distances, indices are in the shape of src
        distances,indices=nbrs.kneighbors(src)
        index=distances[:,0]<self.maxdistance
        match_dist=dist[indices[index]]
        match_src=src[index,:]
        #print match_dist
        #print match_src
        #print src.shape,dist.shape
        #print match_dist.shape,match_src.shape
    
        u,ind,inv=np.unique(indices[index],return_index=True,return_inverse=True)
        #uniq_mask=np.zeros(match_src.shape[0],dtype=bool)
        #uniq_mask[np.unique(indices[index],return_index=True)[1]]=True
        
        duplicate=indices[index][ind[np.bincount(inv)>1]]
        uniq_mask=ind[np.bincount(inv)==1]
        #print indices[index][uniq_mask]
        match_distance=distances[index]
        final_dist=match_dist[uniq_mask,:,:]
        final_src=match_src[uniq_mask,:]
        #print duplicate.shape
        #print "duplicate",duplicate
        #print indices[index]
        #if 0:
        #print final_dist
        for i in xrange(len(duplicate)):
            #print indices[index]
            duplicate_dist=match_distance[indices[index]==duplicate[i]]
            #print i,duplicate[i],indices[index][:,0]
            duplicate_src=match_src[indices[index][:,0]==duplicate[i],:]
            #print duplicate_src
            
            keepsrc=duplicate_src[duplicate_dist==min(duplicate_dist),:]
            #print keepsrc,dist[duplicate[i]]
            keepsrc=keepsrc[0].reshape([1,2])
            #print keepsrc,keepsrc.shape,final_src.shape 
            final_src=np.concatenate((final_src,keepsrc),axis=0)
            keepdist=dist[duplicate[i]].reshape([1,1,2])
            
            #print keepdist,final_dist.shape,keepdist.shape
            final_dist=np.concatenate((final_dist,keepdist),axis=0)
            #break
        #print final_src.shape,final_dist.shape 
        match_arr=np.hstack([final_dist[:,0,:],final_src])    
        return match_arr




    def proj_sky_to_xy(self, xylist, inputcat='', wcsfile='', transfile='', xlim=[0, 2048], ylim=[0, 2048], out_dir=''):
        if wcsfile == '':
            wcsfile = self.get_wcsfile(xylist)
        if transfile == '':
            transfile = self.get_transfile(xylist)
        if out_dir =='':
            out_dir = os.path.dirname(xylist)
        if inputcat == '':
            inputcat = self.catalog
        # FIXME: read input cat

        ids = []; readcolumn(ids, self.catalog.colid, self.catalog.name, datformat='str'); ids= np.array(ids)
        ra = []; readcolumn(ra, self.catalog.colra, self.catalog.name); ra= np.array(ra)
        dec = []; readcolumn(dec, self.catalog.coldec, self.catalog.name); dec= np.array(dec)
        vmag = []; readcolumn(vmag, self.catalog.colmag, self.catalog.name); vmag= np.array(vmag)
        x0, y0 = self._wcs_sky_to_xy(ra, dec, wcsfile)
        
        x, y = self._poly_xy_to_xy(x0, y0, transfile)
        index = (x>xlim[0]) & (x<xlim[1]) & (y>ylim[0]) & (y<ylim[1])

        ids = ids[index]
        ra = ra[index]
        dec = dec[index]
        vmag = vmag[index]
        x = x[index]
        y = y[index]
        fout = open(xylist, mode="w")
        for i in xrange(len(ids)):
            fout.write("%s %f %f %f %f %f\n" % (ids[i], ra[i], dec[i], vmag[i], x[i], y[i]))
        fout.close()
        return
    
    @staticmethod
    def _wcs_sky_to_xy(ra, dec, wcsfile):

        w = wcs.WCS(wcsfile)
        x, y = w.all_world2pix(ra, dec, 1)
        x-=0.5
        y-=0.5
        return [x, y]


    def _poly_xy_to_xy(self, x0, y0, transfile):
        m1, m2 = self.read_transfile(transfile)
        x = poly_ordern_eval(x0, y0, m1)
        y = poly_ordern_eval(x0, y0, m2)
        return [x, y]

    def __call__(self, starlist, wcsfile='', transfile='', out_dir=''):
        if wcsfile == '':
            wcsfile = self.get_wcsfile(starlist.name) 
        if transfile == '':
            transfile = self.get_transfile(starlist.name) 
        if out_dir =='':
            out_dir = os.path.dirname(starlist.name)
        import tempfile
        infile = tempfile.NamedTemporaryFile(delete=False)
        converting = "text2fits -H 'id X Y bg amp S D K flux SN' -f dfffffffff -n '-' %s %s " % (starlist.name, infile.name)
        logging.debug("convert ascii to fits: %s", converting)
        os.system(converting)
        cmdline = "solve-field --overwrite %s -w 2048 -e 2058 -L 20 -H 30 -u app --ra %f --dec %f --radius %d -p -q %.2f --tweak %d --wcs %s -M none -R none -B none -P none " % (infile.name, self.ra, self.dec, self.radius, self.q, self.tweak, out_dir+'/'+ wcsfile)
        logging.debug("Excute astrometry.net: %s", cmdline)

        os.system(cmdline)
        os.unlink(infile.name)
        logging.debug("Project astrometry.net solution")
        
        # FIXME: read input cat
        ra = []; readcolumn(ra, self.catalog.colra, self.catalog.name); ra= np.array(ra)
        dec = []; readcolumn(dec, self.catalog.coldec, self.catalog.name); dec= np.array(dec)

        # FIXME: read input starlist 
        star_x = []; readcolumn(star_x, 2, starlist.name); star_x= np.array(star_x)
        star_y = []; readcolumn(star_y, 3, starlist.name); star_y= np.array(star_y)
        x0, y0 = self._wcs_sky_to_xy(ra, dec, wcsfile)
        src = np.swapaxes(np.vstack([x0, y0]), 0, 1)
        dist = np.swapaxes(np.vstack([star_x, star_y]), 0, 1)
        match_arr = self.kd_match(src,dist)

        self.polyfit(match_arr, transfile)

        

