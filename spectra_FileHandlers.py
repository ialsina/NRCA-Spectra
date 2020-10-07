import os
import sys
from tqdm import tqdm
from datetime import datetime
import time
import shutil
import numpy as np

import spectra_Basics as basics
import spectra_Finders as finder
import spectra_Mixer as mixer
from spectra_InitSettings import cf, peakattr, err

#Notice: some of the functions presented in this file are actually methods of classes presented in spectra_Objects.

def psave(self):
    """This is a method.
    Method that saves the pickle representation of the catalog class to the cwd.
    input:
        - self: catalog instance"""
    assert cf.use_pickle, "Pickle isn't activated"
    import pickle
    with open('spcat.pickle', 'wb') as f:
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
    return None

def pload():
    """Function that loades the pickle representation of the catalog class from the cwd.
    output:
        - class"""
    assert cf.use_pickle, "Pickle isn't activated"
    import pickle
    with open('spcat.pickle', 'rb') as f:
        data = pickle.load(f)
    print('Catalog imported.\nDate of creation: {}\nDate of last modification: {}'.\
        format(data.date_created, data.date_modified))
    return data


def ImportFile(file, check=False):
    """Function that loads a single file of raw data and returns an array of its content.
    Usable for isotope, element, compound and sample files.
    Notice:
        - All the lines starting with a hash (#) will be ignored.
        - All the lines starting with less than two columns will be ignored.
        - Columns are bounded by spaces.
        - Double spaces will be converted into single spaces.
        - First column corresponds to x values, second column to y values. The rest of them are ignored.
        - All the lines with non-numeric data will be ignored.
        - All lines starting with a greater than (>) will considered to carry information about abundance.
        - All spaces at the begining or at the end of the line will be removed.
        - There is a feature for Kaeri files that checks that the filename corresponds with the heading.
            This might be inconvenient in case of format change.
        - Duplicate x values are averaged.
    input:
        - file: str
            whole path to file to import.
        - check: bool
            Whether or not perform the filename-heading check.
            Set to False for every file that doesn't fit the Kaeri format, or either if it has changed at some point.
    output:
        - np.ndarray with x and y values (horizontal)
        - abundances dictionary. Only appliable to elements and compounds. Empty dictionary otherwise."""
    xlist,ylist = [],[]
    with open(file,'r') as iFile:
        
        lines = iFile.readlines()

        # Perform the checking, if asked to
        if check and not 'element' in file and not 'compound' in file:
            wanted_list = basics.InterpretName(file.split('/')[-1].split('.')[0])
            wanted_isotope = wanted_list[1] + '-' + wanted_list[2]
            wanted_mode = wanted_list[3]
            actual_isotope = lines[0].split(' ')[0].replace(')','').split('(')[0]
            actual_mode = lines[0].split(' ')[0].replace(')','').split('(')[1].replace(',','-')
            if actual_isotope != wanted_isotope or actual_mode != wanted_mode:
                print(wanted_isotope,actual_isotope,wanted_mode,actual_mode)
                print('Isotope error in file:', file)
                return None,None

        # lines stores every line in the file.
        # processed gets filled with the lines of usable data.
        processed = []
        abund = dict()
        for line in lines:
            if line[0] == "#": continue
            while line.count('\n') + line.count('  ') != 0: line = line.replace('\n','').replace('  ',' ')
            while line[0] == ' ': line = line[1:]
            while line[-1] == ' ': line = line[:-1]

            # Abundance line here.
            if line[0] == ">":
                if line.count(':') != 1: continue
                key,val = line.replace('>','').replace(' ','').split(':')
                abund[key] = val
                continue
            if not line.replace(' ','').replace('-','').replace('.','')[0].isnumeric(): continue
            if line.count(' ') == 0: continue

            # If it made it this fare, this is usable data.
            processed.append(line)
        
        # lines now stores the usable data
        # processed will get filled with the lines of non-duplicate data (unique x values).
        lines = processed.copy()
        processed = []
        toadd = [0,0]
        summing = 1
        
        for i in range(len(lines)):
            line = lines[i]
            linels = lines[i].split(' ')
            linels_next = lines[i+1].split(' ') if i!=len(lines)-1 else [0,0]
            try:
                linels = [float(el) for el in linels]
                linels_next = [float(el) for el in linels_next]
            except:
                print('Error produced when converting data to float.')
            # Accumulates and averages the data we find...
            toadd[0] = linels[0]
            toadd[1] = (toadd[1]*(summing-1) + linels[1])/summing
            summing += 1
            # ...until the x value changes, or we get to the end of the file. Then reset.
            if linels[0] != linels_next[0] or i == len(lines)-1:
                processed.append(toadd)
                summing = 1
                toadd = [0,0]
        
        # Lists are created and returned, along with abundance
        for i in range(len(processed)):
            xlist.append(processed[i][0])
            ylist.append(processed[i][1])
        iFile.close()
        return np.array([xlist,ylist]), abund


def ImportFileB(file):
    """Same function as ImportFile (see its docstrings), but imports three columns instead of two.
    The third column corresponds to the y error. Appliable to sample files.
    Doesn't have filename-heading check.
    Outputs an array with x, y, err_x, err_y"""
    xlist,ylist,exlist,eylist = [],[],[],[]
    global lines,processed
    with open(file,'r') as iFile:
        lines = iFile.readlines()
        processed = []
        for line in lines:
            if line[0] == "#": continue
            while line.count('\n') + line.count('  ') != 0: line = line.replace('\n','').replace('  ',' ')
            while line[0] == ' ': line = line[1:]
            while line[-1] == ' ': line = line[:-1]
            if line.count(' ') != 2: continue
            processed.append(line)
        lines = processed.copy()
        processed = []
        toadd = [0,0,0]
        summing = 1
        
        for i in range(len(lines)):
            line = lines[i]
            linels = lines[i].split(' ')
            linels_next = lines[i+1].split(' ') if i!=len(lines)-1 else [0,0,0]
            try:
                linels = [float(el) for el in linels]
                linels_next = [float(el) for el in linels_next]
            except:
                print('Error produced when converting data to float')
            toadd[0] = linels[0]
            toadd[1] = (toadd[1]*(summing-1) + linels[1])/summing
            toadd[2] = (toadd[2]*(summing-1) + linels[2])/summing
            summing += 1
            if linels[0] != linels_next[0] or i == len(lines)-1:
                processed.append(toadd)
                summing = 1
                toadd = [0,0,0]
        
        for i in range(len(processed)):
            xlist.append(processed[i][0])
            ylist.append(processed[i][1])
            if i==0:
                exlist.append(processed[1][0]-processed[0][0])
            elif i==len(processed)-1:
                exlist.append(processed[-1][0]-processed[-2][0])
            else:
                exlist.append(max(processed[i+1][0]-processed[i][0], processed[i][0]-processed[i-1][0]))
            eylist.append(processed[i][2])
        iFile.close()
        return np.array([xlist,ylist,exlist,eylist])


def ImportData(directory=None):
    """Imports every file in the specified directory. Intended to be used for the Data files.
    It creates the actual Isotope, Element and Compound instances on the fly.
    How this works:
        - First, it looks for the 'peakprops.txt file. If it exists, asks whether to import it.
        - If importing: propsdict dictionary is loaded.
            Every time a file is imported, propsdict is consulted to get information on the peaks.
            If it is found, the appropriate class is called with this information to use it.
            If not, None is passed so that the class knows that there is work to do.
        - If not importing: propsdict dictionary is an empty dictionary.
            Every time a file is imported, propsdict is consulted to get information on the peaks.
            Obviously there is not, so the appropriate classes are called with this None that
            tells them: get to work!
        - All the Isotope, Element and Compound instances are stored in dictionaries and then returned.
    input:
        - directory: str
            directory path. If unspecified (recomended), imports 'data' file.
    output:
        - isotdict: dictionary of Isotopes
        - elemdict: dictionary of Elements
        - compdict: dictionary of Compounds"""
    from spectra_Objects import Isotope, Element, Compound

    isotdict = dict()
    elemdict = dict()
    compdict = dict()
    if basics.isf('peakprops'):
        inp = input('Import peak properties from file? ([y]/n) >')
        if inp in ['q','quit']:
            return dict(),dict(),dict()
        elif not inp in ['n','no']:
            propsdict = LoadPeaks()
        else:
            propsdict = dict()
    else:
        print('Peak properties file not found.')
        propsdict = dict()
    if propsdict == dict():
        print('Computing peak properties.')
    else:
        print('Importing files into database.')
    time.sleep(0.5)
    isotcount, elemcount, allycount, igncount= 0, 0, 0, 0
    if directory is None: directory = os.path.join(os.getcwd(), 'data')
    for file in tqdm(os.listdir(directory), leave=False):
        filename = os.fsdecode(file)
        if not os.path.getsize(directory+'/'+filename) > 0:
            # EMPTY FILE, ignore
            igncount+=1
            continue
        if filename.endswith(".asm") or filename.endswith(".py"):
            # Weird files? No thanks!
            continue
        else:
            name = filename.replace('element_','').replace('compound_','').split('.')[0]
            arrout, abund = ImportFile(os.path.join(directory, filename), False)
            if filename.startswith('element'):
                # ELEMENT DATAFILE
                elemdict[name] = Element(name, arrout, abund, propsdict.get(name))
                elemcount+=1
                continue
            elif filename.startswith('compound'):
                compdict[name] = Compound(name,arrout, abund, propsdict.get(name))
                allycount+=1
            else:
                # ISOTOPE DATAFILE
                if not arrout is None:
                    name = filename.split('.')[0]
                    isotdict[name] = Isotope(name, arrout, propsdict.get(name))
                    isotcount+=1
                continue
    print(isotcount,'isotope files imported.')
    print(elemcount,'element files imported.')
    print(allycount,'compound files imported.')
    if igncount>0: print(igncount,'empty files ignored.')

    return isotdict, elemdict, compdict

def ImportSamp(directory=None):
    """Imports every file in the specified directory. Intended to be used for the Sample files.
    It creates the actual Sample instances on the fly.
    It works in a very similar way than ImportData (see docstrings on that function). Differences are:
        - It asks the user if he wants to rename the Sample files as they are imported, and gives him
            the chance to do it (would be rude to ask and then not do it).
        - Information on the mode ('n-tot', 'n-g') is gathered and handled.
        - There is no propsdict here.
        - There is no abundances handling here.
    input:
        - directory: str
            directory path. If unspecified (recomended), imports 'data' file.
    output:
        - sampdict: dictionary of Samples"""
    from spectra_Objects import Sample

    sampdict = dict()
    if not basics.isd('samples_n-tot'):
        print('samples_n-tot folder does not exist!')
        return dict()

    if not basics.isd('samples_n-g'):
        print('samples_n-g folder does not exist!')
        return dict()

    print('We are about to load samples spectra.')

    naming = False
    inp = input('Name manually? (y/[n]) >')
    if inp in ['y','yes']:
        naming = True
    elif inp in ['q','quit']:
        return dict()

    sampcount, igncount= 0, 0
    for mode in ['n-tot', 'n-g']:
        print('Importing {} sample files into database'.format(mode))
        if naming: print('Type in the sample names. [] will name them as the filename and [-] will omit the import')
        time.sleep(0.5)
        folder = 'samples_' + mode
        directory = os.path.join(os.getcwd(), folder)
        filelist = os.listdir(directory)
        for file in tqdm(filelist, disable=naming, leave=False):
            filename = os.fsdecode(file)
            if not os.path.getsize(os.path.join(directory,filename)) > 0:
                # EMPTY FILE
                igncount+=1
                continue
            if filename.endswith(".asm") or filename.endswith(".py"):
                # Weird files? No thanks!
                continue
            else:
                arrout = ImportFile(os.path.join(directory,filename), False)[0]
                # SAMPLE DATAFILE
                if not arrout is None:
                    if naming:
                        name = input('{}/{}: {} >'.format(sampcount+igncount+1, len(filelist), filename))
                        if name in ['q','quit']:
                            return dict()
                        elif name in ['','=']:
                            name = filename.split('.')[0]
                        elif name in ['-','x']:
                            continue
                    else:
                        name = filename.split('.')[0]
                    sampdict[name] = Sample(name, arrout, mode)
                    sampcount+=1
                continue
    print(sampcount,'sample files imported.')
    if igncount>0: print(igncount,'empty files ignored.')
    return dict(sampdict)


def LoadPeaks():
    """Reads the peakprops.txt file to return a big dictionary with all of peak properties stored in there.
    output:
        - THE dictionary of peaks."""
    from spectra_Objects import Peak

    peaks = dict()
    if not basics.isf('peakprops'):
        print('peakprops.txt must exist')
        return None
    with open('peakprops.txt','r') as iFile:
        curIsot = ''
        lines = iFile.readlines()
        for line in lines:
            line = line.replace('=','').replace('|','').replace('\n','')
            # In case some stupid human touched it.
            if line.replace(' ','') == '': continue
            while line[0] == ' ': line = line[1:]
            while line[-1] == ' ': line = line[:-1]
            if line[0] in ['#','>']:
                   continue
            elif line[0] == ':':
                curIsot = line.replace(':','')
                peaks[curIsot] = dict()
                if line == '::': break      #End of file
            else:
                while '  ' in line: line = line.replace('  ',' ')
                linels = line.split(';')
                if len(linels) == peakattr.size:
                    peaks[curIsot][int(linels[1])] = Peak(peakattr.maketuple(linels))
                else:
                    print('Bad amount of data separators in peakprops.txt file found on import.\nComputing peak properties instead.')
                    return dict()
    return peaks
    print('File of properties imported')


def CreateEmptyFiles():
    """Creates all the empty files to be filled from Kaeri (for example) after the file 'Isotopes.txt'"""
    if basics.isf('Isotopes'):
        try:
            if basics.isd('EmptyFiles'): shutil.rmtree(basics.path+'/EmptyFiles')
            os.mkdir(basics.path+'/EmptyFiles')
        except Exception as e:
            print('Error in output file creation:',e)
            sys.exit()
        
        with open('Isotopes.txt','r') as iFile:
            for line in iFile.readlines():
                line = line.replace('\n','')
                print(line)
                for suff in ['_n-tot','_n-g']:
                    f = open(basics.path+'/EmptyFiles/'+line+suff+'.txt','w')
                    f.close()
    else:
        print('There must be an Isotopes file')


def MixOut(Dict, container_in):
    """Creates the empty Natural_out.txt file and Compound_out.txt file for the user to fill.
    inputs:
        - container_in: tuple
            tuple of isots, ielems, elems, and non_unique used for the listing."""
    instructions = """Instructions:
# Insert desired abundance of isotopes after the colon (:). Spaces will be ignored.
# A colon (:) with a character string on its left must be used to name the element/compound.
# Allowed punctuation signs and symbols for the name are .,-+_!?[]()%&/ but nothing else.
# It is possible to mix and join different isotopes and/or elements together.
# Only compounds with least two non-empty abundance entries will be considered.
# If the abundances don't sum up to 1, they will be rescaled automatically.
# A period (.) indicates decimal place and an asterisk (*) can be used to indicate mutliplication.
# Only 0 to 9 digits and the above symbols should be used
# Some examples are provided. Either leave abundances empty or delete them.
# Change the file name to 'Compound_in.txt' and load from program
# The file must end with a double colon (::)\n"""
    
    isots, ielems, elems, non_unique, comps = container_in

    with open('Natural_out.txt','w') as oFile:
        oFile.write('# NATURAL ABUNDANCES'+'\n'+instructions)
        oFile.write('#'*70+'\n\n')
        for elem in ielems:
            oFile.write('\n:'+elem+'\n')
            for res in sorted(finder.Seek(Dict, (elem,'n-tot'))):
                isot = res.split('_')[0]
                if elem in non_unique:
                    oFile.write(isot+':\n')
                else:
                    oFile.write(isot+':\n')

        oFile.write('\n::')
    print('Natural_out.txt file exported.')


    with open('Compound_out.txt','w') as oFile:
        oFile.write('# COMPONENT PROPORTION IN COMPOUNDS'+'\n'+instructions)

        oFile.write('#'*70+'\n\n')
        oFile.write('\n##Example:\n#:MyCompound1\n#22-Ti:50\n#29-Cu:49.2\n#24-Cr:.8\n')
        oFile.write('\n##Example:\n#:special_steel\n#6-C:.25\n#26-Fe-56:.75*.85\n#26-Fe-54:.75*.10\n#26-Fe-57:.75*.05\n')
        oFile.write('\n#Isotopes:\n')
        for item in sorted(list(np.unique([el.split('_')[0] for el in isots]))):
            oFile.write(item+'\n')
        oFile.write('\n#Elements:\n')
        for item in sorted(list(np.unique([el.split('_')[0] for el in elems]))):
            oFile.write(item+'\n')
        oFile.write('\n#Compounds:\n')
        for item in sorted(list(np.unique([el.split('_')[0] for el in comps]))):
            oFile.write(item+'\n')
        oFile.write('\n#New compounds: (build below)\n')
        oFile.write('\n'*10+'::')

    
    print('Compound_out.txt file exported.')
    
def MixInFile(filename, Dict, permitted, kind=None):
    """Reads a single file, either Natural_in, or Compounds_in."""
    from spectra_Objects import Isotope, Element, Compound
    dictout = dict()
    dictcomp = dict()
    
    if not basics.isf(filename):
        print('File does not exist')
        return dict()
    else:
        if not input('Importing '+filename+'. Continue?\n'\
                     'Note: this is likely to take a while. ([y]/n) >') in ['y','yes','']: return dict()
    
    mixels = []
    mixname = ""
    mix = np.empty((0))
    skip = False

    #Objects assigned depending on who are we running the function for
    assert kind in ['compound', 'element'], "Bad value for the 'kind' argument"
    obj = Compound if kind=='compound' else Element

    with open(filename,'r') as iFile:
        for line in iFile.readlines():
            if line[0] in ["#",">"]: continue
            line = line.replace('=','').replace('\n','')
            line = line.replace(' ','')
            if line == '': continue
            if line.count(':') == 0: continue
            if line != "::" and line.count(':') != 1:
                print('Misplacement of colon (:) in line:',line)
                return dict()
            if line[0] == ":":
                
                #NEW ITEM
                #If we are not in a bad line: normalize every component, store
                if np.size(mix) >=1 and not skip:
                    dictcomp[mixname] = {x:y for x,y in zip(mixels,mix/np.sum(mix))}
                
                #And restart, in any case:
                skip = False
                mixname = line[1:]
                mixels = []
                mix = np.empty((0))
                continue
            perc = 1
            comp,percstr = line.split(':')
            if percstr == '': continue
            if comp not in permitted:
                skip = True
                print('Warning: Bad component:', comp, '- Skipping import')
                continue
            try:
                perclist = [float(el) for el in percstr.split('*')]
            except:
                skip = True
                print('Warning: Non numerical input in component:', comp, '- Skipping import')
                continue
            for el in perclist:
                perc = perc*el
            mixels.append(comp)
            mix = np.append(mix,perc)

    count = 0
    for mix in dictcomp:
        for suff in ['_n-tot','_n-g']:
            #Compute weighted outcome for every mode
            print('Computing:', mix, suff)
            dictout[mix+suff] = obj(mix+suff,mixer.GetWeighted(Dict,mix,suff,dictcomp[mix]),dictcomp[mix],peaksdict=None)
            ExportWeighted(dictout[mix+suff])
            count+=1
    
    os.rename(filename,filename.split('_')[0]+'_through.txt')
    print(filename, 'file imported.')
    print(count, kind+'(s) computed and exported.')
    
    return dictout

def MixIn(Dict,permitted):
    """Makes the function MixInFile work for Elements and Compounds.
    If the respective files are there and ready, it returns a dictionary with them all computed.
    If not, it returns an empty dictionary."""
    dictelements = MixInFile('Natural_in',Dict,permitted, 'element') if basics.isf('Natural_in') else dict()
    dictcompounds = MixInFile('Compound_in',Dict,permitted, 'compound') if basics.isf('Compound_in') else dict()
    return dictelements, dictcompounds

def ExportWeighted(*args,directory=None):
    if directory is None: directory = os.getcwd() + '/data/'
    for el in args:
        with open(directory+el.kind+'_'+el.fullname+'.txt','w') as oFile:
            oFile.write(el.fullname+'\n')
            for comp in el.abundances:
                oFile.write('>{:>10s}:{:<8.6f}\n'.format(comp, float(el.abundances[comp])))
            for xx,yy in zip(el.spectrum[0],el.spectrum[1]):
                oFile.write('{:12.6e} {:14.8e}\n'.format(xx,yy))


def ExportProps2(Dict):
    """Given a Data dictionary, exports a comprehensive non-human-readable file of ugly values."""
    with open('peakprops.txt','w') as oFile:
        oFile.write('# ISOTOPE PEAK PROPERTIES FILE\n')
        oFile.write('>TOF:{},{};E:{},{};b:{};L0:{},{}'.format(cf.tof_min,cf.tof_max,cf.e_min,cf.e_max,cf.crs_min,cf.L0_t,cf.L0_g)+'\n')
        oFile.write('>Date: '+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+'\n')
        oFile.write('#'+';'.join(peakattr.getlist())+'\n')
        for isot in sorted(Dict):
            oFile.write(':'+Dict[isot].fullname+'\n')
            for peak in Dict[isot].peaks:
                oFile.write(';'.join([str(getattr(Dict[isot].peaks[peak],attrstr)) for attrstr in peakattr.getlist()])\
                    .replace('(','').replace(')','').replace(' ','')+'\n')
        oFile.write('::')
        oFile.close()
    print('peakprops.txt exported')

def infoone(isot):
    """Given a Data instance, squishes all the information that goes into human-readable files out of it."""
    linesout = []
    for peak in range(len(isot.peaks)):
        line = isot.peaks[peak]
        linesout.append('{:>4d}{:1s} {:>10.3e} {:^6s} {:<10.3e} {:10.3e} {:10.3e} {:<6s} {:10.3e} {:<6s}'.\
                    format( line.num,\
                            '*' if line.integral == 0 else '',\
                            line.center,\
                            '('+str(line.center_)+')',\
                            basics.E2t(line.center,isot.mode),\
                            line.integral,\
                            line.width,\
                            '('+str(line.width_)+')',\
                            line.height,
                            '('+str(line.height_)+')'))

    return linesout



def ExportProps(Dict):
    """Given a Data dictionary, exports a human-readable file of nice values."""
    with open('PeakProperties.txt','w') as oFile:
        oFile.write('# ISOTOPE PEAK Properties\n')
        oFile.write('\n# Peaks conditioned to: ToF in [{},{}] us, L0_g = {}, L0_t = {}, C.S. > {} b'.format(cf.tof_min,cf.tof_max,cf.L0_g,cf.L0_t,cf.crs_min))
        oFile.write('\n# Which corresponds to: E in [{},{}] eV, C.S. > {} b'.format(cf.e_min,cf.e_max,cf.crs_min))
#        oFile.write('\n\n>Isotopes:{}\n>Non-empty:{}\n>Complete:{}\n>Peaks:{}\n>Successful:{}\n'.\
#            format(countIsotopes,countNonEmpty,countComplete,countPeaks,countSuccess))
        oFile.write('\n# Computation date and time: '+datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        oFile.write('\n'+'='*93+'\n\n')
        oFile.write('#'+'='*92+'\n')
        oFile.write('#{:>3s}  {:>10s}        {:<10s} {:>10s} {:>17s} {:>17s}\n'.\
                    format('Rk.','Energy (eV)', 'TOF (us)', 'Integral','Peak width','Peak height'))
        oFile.write('#'+'='*92+'\n')
        for isot in sorted(Dict):
            oFile.write('\n:{:=<92}'.format(isot)+'\n')
            for line in infoone(Dict[isot]): oFile.write(line+'\n')
        oFile.close()
    print('PeakProperties file exported')


def BackName():
    """Renames everything with _through in its name back to _in"""
    if basics.isf('Natural_through'): os.rename('Natural_through.txt','Natural_in.txt')
    if basics.isf('Compounds_through'): os.rename('Compounds_through.txt','Compounds_in.txt')

            
            
# UNUSED, ALTHOUGH MIGHT BE USEFUL T SOME POINT.
# =============================================================================
# def BadPeaks(request=None):
#     if os.path.isdir(path+'/widepeaks'): shutil.rmtree(path+'/widepeaks')
#     os.mkdir(path+'/widepeaks')
#     with open('widepeaks.txt','w') as oFile:
#         for isot in tqdm(sorted(Properties)):
#             badone = False
#             for peak in Properties[isot]:
#                 if Properties[isot][peak].get('width',0) > 900:
#                     badone = True
# #                    print(isot,peak,Properties[isot][peak].get('width'))
#                     oFile.write('{}\t{}\t{}\n'.format(isot,peak,Properties[isot][peak].get('width')))
#             if badone:
#                 Define(isot)
#                 PropsOne()
#                 Plot0()
#                 plt.savefig(path+'/widepeaks/'+isot+'.png')
#                 plt.close()
#         oFile.close()    
# =============================================================================