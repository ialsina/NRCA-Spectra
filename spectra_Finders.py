from spectra_Basics import *
from spectra_InitSettings import cf, err

def Query(askmode=cf.ask_mode):
    """Basically, asks the user what to look up.
    input:
        - askmode: bool
            Whether or not ask the mode.
    output:
        - tuple with the desired input to look up and the mode."""
    inp = input('Input a substance >') if cf.ask_mode else input(' >')
    if inp in ['q','quit','exit','']: return None
    qm = input('Input a mode: [n,tot]; [n,g]; [] both or N/A >') if askmode else ''
    if qm in ['n-tot','n,tot','ntot','n-t','n,t','nt']:
        qm = 'n-tot'
    elif qm in ['n-g','n,g','ng','n-gamma','n,gamma','ngamma']:
        qm = 'n-g'
    elif qm in ['',' ','-']:
        qm = ''
    else:
        print('Mode not understood. Showing both')
        qm = ''
    return (inp,qm)


def Seek(Dict, tuplein):
    """Searches a query from the function Query into a dictionary.
    input:
        - Dict: dictionary.
            dictionary to look up the query. Notice that this must be a special dictionary with
            the keys: 'isotopes', 'elements', 'compounds', and 'samples', and then their contents, either
            dictionaries or lists. This is because the criteria to look this up is different in each case.
        - tuplein: results from Query().
    output: a list of every item matching in the dictioanry matching the criteria."""
    if not isinstance(tuplein,tuple): return None
    if len(tuplein) != 2: return None
    inp, qm = tuplein
    if inp is None or qm is None: return None
    outp = []
    
    for isot in inp.replace(' ','').split(','):
        isotl = isot.split('-')

        #different possibilities have different treatments, right?
        if len(isotl) == 0:
            qz,qe,qn = '','',''
        elif len(isotl) == 1:
            if isotl[0].isnumeric():
                qz,qe,qn = isotl[0],'',''
            else:
                qz,qe,qn = '',isotl[0],''
        elif len(isotl) == 2:
            if isotl[0].isnumeric() and isotl[1].isnumeric():
                qz,qe,qn = isotl[0],'',isotl[1]
            elif isotl[0].isnumeric() and not isotl[1].isnumeric():
                qz,qe,qn = isotl[0],isotl[1],''
            elif not isotl[0].isnumeric() and isotl[1].isnumeric():
                qz,qe,qn = '',isotl[0],isotl[1]
            else:
                return None
        elif len(isotl) == 3:
            if isotl[0].isnumeric() and not isotl[1].isnumeric() and isotl[2].isnumeric():
                qz,qe,qn = isotl[0],isotl[1],isotl[2]
            else:
                return None

        # Bank joins everything back together to be looped over.
        Bank = dict(Dict.get('isotopes',dict()), **dict(Dict.get('elements',dict()), **dict(Dict.get('compounds',dict()), **Dict.get('samples',dict()))))
        for isot in Bank:

            #Below, by 'bit' we understand either atomic number, symbol, mass number or mode.
            #if the current page is not in compounds or samples, we require every query-bit of information the user gave to be the same.
            if not isot in dict(Dict.get('compounds'), **Dict.get('samples')):
                z,e,n,m = InterpretName(isot)
                #if the user didn't give the information, both the query-bit and the dictionary-bit are set to 0 or to ''.
                if qz in ['','0']:
                    qz='0'
                    z='0'
                if qe=='':e=''
                if qn in ['','0']:
                    qn='0'
                    n='0'
                if qm=='':m=''
                cond = int(qz) == int(z) and qe == e and int(qn) == int(n) and qm == m
            #if the current page is in samples, we require that the user's query is either the name or in one of the parts between underscores (_)
            elif isot in Dict.get('samples'):
                parts = isot.split('_')
                cond = inp in parts or inp == isot
            #if the current page is in compounds, the input must be contained in the name
            else:
                cond = inp in isot
            #if the condition is satisfied, there you go.
            if cond: outp.append(isot)
    return sorted(outp)



def Select(Dict,recursive=False,ask_if_one=True,askmode=cf.ask_mode,restrict=None):
    """Recursively (or not) asks the user enter a series of queries (a single query) and to
     select from the possible results that match them (it).
    inputs:
        - Dict: dictionary
            dictionary to look up the query. Notice that this must be a special dictionary with
            the keys: 'isotopes', 'elements', 'compounds', and 'samples', and then their contents, either
            dictionaries or lists. This is because the criteria to look this up is different in each case.
        - recursive: bool
            False: only once is asked.
            True: it asks recursively and fills a list and until s/he enters an empty string.
        - ask_if_one: bool
            Should we ask the user if his query provides only one result?
        - restrict: str or None
            str: if the query results do not contain this string, they aren't even presented.
                 (particularly useful to select a single mode)
            None: function disabled.
    outputs:
        - list with all the user's choices.
        """
    outp = []
    while True:
        curquery = Seek(Dict,Query(askmode))
        if curquery is None: break
        if restrict:
            for el in curquery:
                if not restrict in el: curquery.remove(el)
        if len(curquery) != 1 or ask_if_one:
            print('\n'.join(['{}: {}'.format(i+1,curquery[i]) for i in range(len(curquery))]))
            curadd = input('Select above results of your query ([0] all; [] None) >')
            if curadd in ['q','quit']: return None
            for el in curadd.replace('-',',').split(','):
                if el == '': break
                if not el.isnumeric():
                    print('Non numeric input, sorry')
                    return None
                if el == '0':
                    outp.extend(curquery)
                    break
                if int(el) > len(curquery):
                    print(el, 'is too big')
                    break
                outp.append(curquery[int(el)-1])
        else:
            outp.extend(curquery)
        if not recursive: return outp[0] if len(outp) != 0 else None
        print('Current elements:',outp)
    return list(np.unique(outp))