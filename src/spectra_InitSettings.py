class Settings:
    def __init__(self):
        from ..settings import parameters
        for param in parameters:
            setattr(self,param,parameters[param])
        
        self.e_min   = self.t2E(self.thr_max,self.default_mode) if self.thr_in_tof else self.thr_min
        self.e_max   = self.t2E(self.thr_min,self.default_mode) if self.thr_in_tof else self.thr_max
        self.e_min_g = self.t2E(self.thr_max,'n-g')             if self.thr_in_tof else self.thr_min
        self.e_max_g = self.t2E(self.thr_min,'n-g')             if self.thr_in_tof else self.thr_max
        self.e_min_t = self.t2E(self.thr_max,'n-tot')           if self.thr_in_tof else self.thr_min
        self.e_max_t = self.t2E(self.thr_min,'n-tot')           if self.thr_in_tof else self.thr_max
        self.tof_min = self.thr_min                             if self.thr_in_tof else self.E2t(self.thr_max,self.default_mode)
        self.tof_max = self.thr_max                             if self.thr_in_tof else self.E2t(self.thr_min,self.default_mode)
        self.tof_min_g = self.thr_min                           if self.thr_in_tof else self.E2t(self.thr_max,'n-g')
        self.tof_max_g = self.thr_max                           if self.thr_in_tof else self.E2t(self.thr_min,'n-g')
        self.tof_min_t = self.thr_min                           if self.thr_in_tof else self.E2t(self.thr_max,'n-tot')
        self.tof_max_t = self.thr_max                           if self.thr_in_tof else self.E2t(self.thr_min,'n-tot')


    def _funcs(self):
        return self.L0, self.E2t, self.t2E, self.dt2dE, self.dE2dt
    
    def L0(self,mode=None):
        if mode is None: mode = self.default_mode
        return self.L0_g if mode=='n-g' else (self.L0_t if mode=='n-tot' else None)
    def E2t(self,En,mode=None):
        return self.L0(mode)*1E6*(0.5*self.mn/(self.e*En))**0.5
    def t2E(self,t,mode=None):
        return 0.5*(self.mn/self.e)*((self.L0(mode)*1E6)/t)**2
    def dt2dE(self,dt,t,mode=None):
        return 10**12*self.mn*self.L0(mode)**2/(self.e*t**3)*dt
    def dE2dt(self,dE,E,mode=None):
        return 0.5*self.L0(mode)*1E6*(0.5*self.mn/self.e)**0.5*E**(-1.5)*dE

    def xbounds(self,tof=False):
        return (self.tof_min, self.tof_max) if tof else (self.e_min, self.e_max)

    def ybounds(self,symb,mode='n-tot'):
        exc = self.crs_exc.get(symb, self.crs_min)
        if isinstance(exc, tuple) or isinstance(exc, list):
            return exc[{'n-tot': 0, 'n-g': 1}.get(mode)]
        elif isinstance(exc, float) or isinstance(exc, int):
            return exc
        else:
            return self.crs_min


    def pack(self,**kwargs):
        return {attr: kwargs.get(attr, getattr(self,attr)) for attr in self.__dict__}



class PeakAttributes:
    totup = lambda ty: lambda inp: tuple([ty(el) for el in inp.split(',')])
    tobool= lambda inp: inp in ['True','true','1',1,'yes']
    attr = (
            (   'fullname'      ,   str             ),      #00   #Name of the isotope they belong to
            (   'num'           ,   int             ),      #01   #Peak label, mainly. Should match with integral_
            (   'center_'       ,   int             ),      #02   #Energy rank of peak center. For example: in an isotope with N peaks, 0 is the first peak in the spectrum and N is the last one.
            (   'center'        ,   float           ),      #03   #Energy value of peak center
            (   'icenter'       ,   int             ),      #04   #Index number of peak center in raw data file
            (   'center_tof'    ,   float           ),      #05   #Time of Flight (us) value of peak center
            (   'integral'      ,   float           ),      #06   #Integral value of peak
            (   'integral_'     ,   int             ),      #07   #Integral rank of peak. For example: in an isotope with N[0] peaks[errors], 0 is the most intense peak and N is the least intense one.
            (   'integral_tof'  ,   float           ),      #08   #Integral value of peak (us)
            (   'width'         ,   float           ),      #09   #Width value of peak (should be xlims[1]-xlims[0])
            (   'width_'        ,   int             ),      #10   #Width rank of peak
            (   'height'        ,   float           ),      #11   #Height value of peak. Computed as the XS (or NC) at peak center minus the arithmetic mean of the XSs (or NCs) values on the peak bounds.
            (   'height_'       ,   int             ),      #12   #Height rank of peak
            (   'fwhm'          ,   float           ),      #13   #Full width at half maximum of peak, i.e. energy distance when the XS (or NC) has dropped to 1/2 along the peak respect to the maximum value
            (   'fwhm_'         ,   int             ),      #14   #Full width at half maximum rank of peak
            (   'ahh'           ,   float           ),      #15   #AHH: Integral divided by squared height
            (   'ahh_'          ,   int             ),      #16   #AHH rank of peak
            (   'ahw'           ,   float           ),      #17   #AHW: Integral divided by height times width
            (   'ahw_'          ,   int             ),      #18   #AHW rank of peak
            (   'xlims'         ,   totup(float)    ),      #19   #Energy values of peak bounds
            (   'ilims'         ,   totup(int)      ),      #20   #Data indices of peak bounds
            (   'outerslope'    ,   totup(float)    ),      #21   #Ideally, slope value far away from the peak and uncertainty (tuple). See documentation in definepeak for further information on this.
            (   'peakreason'    ,   totup(int)      ),      #22   #Left and right boundaries reason (tuple). See documentation in definepeak for further information on this.
            (   'yvals'         ,   totup(float)    ),      #23   #CS (or NC) values at peak boundaries (tuple)
            (   'coords'        ,   totup(float)    ),      #24   #Energy and CS (or NC) values at peak summit
            (   'coords_tof'    ,   totup(float)    ),      #25   #Time of Flight (us) and CS (or NC) values at peak summit
            (   'prange'        ,   int             ),      #26   #prange value. See documentation in definepeak for further information on this.
            (   'user_edited'   ,   tobool          ),      #27   #True if user has edited the peak, False if the peak attributes are all from computation.
            (   'user_defined'  ,   tobool          ),      #28   #True if user has defined the peak, False if the peak is defined from computation.
            )
    dattr = dict(attr)
    size = len(attr)


    def has(self,inp):
        return inp in self.dattr

    def gettup(self,ind):
        if not isinstance(ind,int): return (None,str)
        if ind >= self.size: return (None,str)
        return self.attr[ind]

    def get(self,ind):
        return self.gettup(ind)[0]

    def conv(self,ind,inp):
        return self.gettup(ind)[1](inp)

    def maketuple(self,elements):
        return tuple([self.conv(i,elements[i]) for i in range(self.size)])

    def getlist(self):
        return [self.get(i) for i in range(self.size)]

class ErrorReporter:
    def __init__(self):
        self.start()

    def start(self):
        self.errcount = 0
        self.errdict = dict()

    def add(self,excep,loc,label,num,consoleprint=False):
        if consoleprint:
            print('!!! EXCEPTION RAISE\n\tFunction: {}\n\tSubstance: {}\n\tExeption: {}'.format(loc, label, excep.args[0]))
        self.errcount += 1
        self.errdict[self.errcount] = dict(func=loc,excep=excep,subs=label,peakc=num)

    def get(self):
        return self.errdict

    def present(self):
        if self.errcount == 0:
            print('No errors raised.')
        else:
            print('Errors report:')
            print('{:5s}  {:>12s} {:>4s} {}'.format('#Err.','Substance','#Pos','Description'))
            for el in self.errdict:
                entry = self.errdict[el]
                print('{:>5}: {:>12} {:>4} {}'.format(el, entry['subs'], entry['peakc'], entry['excep'].args[0]))

class Path:
    def __init__(self):
        from ..paths import paths_
        for p in paths_:
            setattr(self, p, paths_[p])
    
    def __get(self,d):
        if hasattr(self,d):
            return getattr(self,d) 
        else:
            import os
            return os.path.join(self.cwd, d)
    
    def path(self, d = None):
        if d is None:
            return self.cwd
        else:
            return self.__get(d)
    
    def isd(self,d):
        import os
        #return os.path.isdir(os.path.join(path,d))
        return os.path.isdir(self.path(d))
    
    def isfx(self,f,d=None):
        import os
        
        return os.path.isfile(os.path.join(self.path(d),f))
        #return os.path.isfile(getattr(self,f))
    
    def isf(self,f,d=None):
        assert 0==1, 'Function was deprecated and must not be used!'
        #return self.isfx(f+'.txt',d)
    
    def join(self, d, f):
        import os
        return os.path.join(self.path(d), f)
    
    

cf = Settings()
peakattr = PeakAttributes()
err = ErrorReporter()
paths = Path()