import csv, json, pprint, re, os, numpy as np
from scipy import stats
from tqdm import *

pp = pprint.PrettyPrinter(indent=2)

def extract_acs_sum(fname):
    # print 'Extracting: %s'%fname
    typ = fname[11:12]
    st_nm = fname[17:19]
    seq_nm = fname[-11:-7]
    acs = []
    with open(fname, 'rb') as f:
        for row in csv.reader(f, delimiter=','):
            acs.append(row)
    # Get Header info
    with open('./acs_sas/%s%s_%s.sas'%(typ,st_nm,seq_nm)) as f:
        read_flg = False
        header = []
        for ln in f:
            if ln.strip() not in ('', ';', 'RUN;') and read_flg:
                header.append(re.split('\s*', ln.strip())[0])
            read_flg = True if 'INPUT' in ln else read_flg
    # Get geography
    with open('./acs_geog/g20135%s.txt'%st_nm,'rb') as f:
        geog = {}
        nottract = []
        for row in f:
            logrecno = row[13:20]
            geoid = row[178:218].strip()
            geog[logrecno] = geoid
        for i, row in enumerate(acs):
            if geog[row[5]][0:2] == '14':
                row.append(geog[row[5]])
            else:
                # Remove Non tract geography
                nottract.append(i)
        # Reverse sort to delete without messing index
        nottract.sort(reverse=True)
        for row in nottract:
            del acs[row]
    return [header + ['geoid']] + acs

def build_directory(dirnm, force=False):
    # print 'Loading Variable to Summary Table Mapping'
    if force or not os.path.isfile('directory.json'):
        directory = {}
        for fname in tqdm(os.listdir(dirnm)):
            typ = fname[0:1]
            st_nm = fname[1:3]
            seq_nm = fname[4:8]
            with open(dirnm + fname, 'rb') as f:
                read_flg = False
                for ln in f:
                    if ln.strip() not in ('', ';', 'RUN;') and read_flg:
                        varnm = re.split('\s*', ln.strip())[0]
                        try:
                            directory[varnm].append([typ, st_nm, seq_nm])
                        except KeyError:
                            directory[varnm] = [[typ, st_nm, seq_nm]]
                    read_flg = True if 'INPUT' in ln else read_flg
        for key in ['FILEID', 'FILETYPE', 'STUSAB', 'CHARITER', 'SEQUENCE', 'LOGRECNO']:
            directory.pop(key)
        with open('directory.json', 'wb') as f:
            json.dump(directory, f)
    with open('directory.json', 'rb') as f:
        directory = json.load(f)
    return directory

def load_varlist(fname):
    print 'Loading ACS 5 year Summary File Variable List'
    with open(fname, 'r') as f:
        varlist = json.load(f)['variables']
    for varname in tqdm(varlist.keys()):
        if '0' in varname:
            split = re.split('_', varname)
            corrected = split[0] + split[1][-1].lower() + split[1].strip('0')[0:-1]
            varlist[corrected] = varlist.pop(varname)
    return varlist

def build_custom(varlist, score):
    print '*****************************************************************************************************************************'
    print 'Building Custom List:'
    print '*****************************************************************************************************************************'
    directory = build_directory('./acs_sas/')
    custom = {}
    header = []
    for var in tqdm(varlist):
        header.extend([var])
        if var in score:
            header.extend(['score_' + var])
        for table in directory[var]:
            typ = table[0]
            st_nm = table[1]
            seq_nm = table[2]
            acs_summary = extract_acs_sum('./acs_summ/' + typ + '20135' + st_nm + seq_nm + '000.txt')
            position = [ i for i,x in enumerate(acs_summary[0]) if x == var ]
            if var in score:
                scoreme = [ x[position[0]] for x in acs_summary[1:] ]
                scored = [ stats.percentileofscore(scoreme, x) for x in scoreme ]
            for i, row in enumerate(acs_summary[1:]):
                try:
                    custom[row[-1]].extend([row[position[0]]])
                except KeyError:
                    custom[row[-1]] = [row[position[0]]]
                if var in score:
                    custom[row[-1]].extend([scored[i]])
    position = [ i for i,x in enumerate(header) if 'score' in x ]
    for tract in custom.keys():
        avg = sum([ float(custom[tract][x]) for x in position ])/len(position)
        custom[tract].extend([avg])
    header.extend(['final_score'])
    custom['varnames'] = header
    with open('./custom.json', 'wb') as f:
        json.dump(custom, f)
    return custom

def build_custom_json(varlist, score):
    print '*****************************************************************************************************************************'
    print 'Building Custom List:'
    print '*****************************************************************************************************************************'
    directory = build_directory('./acs_sas/')
    custom = {}
    header = []
    for var in tqdm(varlist):
        header.extend([var])
        if var in score:
            header.extend(['score_' + var])
        for table in directory[var]:
            typ = table[0]
            st_nm = table[1]
            seq_nm = table[2]
            acs_summary = extract_acs_sum('./acs_summ/' + typ + '20135' + st_nm + seq_nm + '000.txt')
            position = [ i for i,x in enumerate(acs_summary[0]) if x == var ]
            if var in score:
                scoreme = [ x[position[0]] for x in acs_summary[1:] ]
                scored = [ stats.percentileofscore(scoreme, x) for x in scoreme ]
            for i, row in enumerate(acs_summary[1:]):
                try:
                    custom[row[-1]][var] = row[position[0]]
                except KeyError:
                    custom[row[-1]] = {}
                    custom[row[-1]][var] = row[position[0]]
                if var in score:
                    custom[row[-1]]['score_'+var] = scored[i]
    position = [ x for x in header if 'score' in x ]
    for tract in custom.keys():
        avg = sum([ float(custom[tract][x]) for x in position ])/len(position)
        custom[tract]['final_score'] = avg
    header.extend(['final_score'])
    # custom['varnames'] = header
    with open('./custom.json', 'wb') as f:
        json.dump(custom, f)
    return custom

if __name__=='__main__':
    # Test case
    # variables = ['B17001e2', 'B22002e2', 'B23025e7']
    # scores = ['B17001e2', 'B22002e2', 'B23025e7']
    # Actual pull
    # variables = ['B00001e1', 'B00002e1', 'B01001e2', 'B01001e26', 'B01001e27', 'B01001e28', 'B01001e29', 'B01001e30', 'B02001e1', 'B02001e2', 'B02001e3', 'B02001e4', 'B02001e5', 'B02001e6', 'B02001e7', 'B02001e8', 'B02001e9', 'B02001e10', 'B17001e1', 'B17001e2', 'B18101e1', 'B18105e1', 'B21100e1', 'B21100e2', 'B21100e3', 'B22001e1', 'B22001e2', 'B22001e3', 'B22001e4', 'B22002e1', 'B22002e2', 'B22002e3', 'B22002e4', 'B22002e5', 'B22002e6', 'B22002e7', 'B22002e8', 'B22002e9', 'B22002e10', 'B22002e11', 'B22002e12', 'B22002e13', 'B22002e14', 'B23024e1', 'B23024e2', 'B23024e3', 'B23024e4', 'B23024e5', 'B23024e6', 'B23024e7', 'B23024e8', 'B23024e9', 'B23025e1', 'B23025e5', 'B23025e7']
    # scores = ['B17001e2', 'B22002e2', 'B23025e7']
    # Tutorial pull
    # Girls under 5, Same house 1 year ago, SNAP benefits the past year
    variables = ['B01001e27', 'B07001e17', 'B22002e2']
    scores = variables
    varlist = load_varlist('./variables.json')
    header = { x: varlist[x] for x in variables }
    with open('./header.json', 'wb') as f:
        json.dump(header, f)
    print '*****************************************************************************************************************************'
    print 'Pulling the following:'
    print '*****************************************************************************************************************************'    
    pp.pprint(header)
    scored = { x: varlist[x] for x in scores }
    with open('./scored.json', 'wb') as f:
        json.dump(scored, f)
    print '*****************************************************************************************************************************'
    print 'Scoring the following:'
    print '*****************************************************************************************************************************'
    pp.pprint(scored)
    custom = build_custom_json(variables, scores)
