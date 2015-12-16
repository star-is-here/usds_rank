import csv, json, pprint, re, os, numpy as np, fiona, shapely
from scipy import stats
from tqdm import *

pp = pprint.PrettyPrinter(indent=2)

def loadinput(fname, tp):
    with open(fname) as data:
        if tp == 'csv':
            loadme = []
            csvloadme = csv.reader(data)
            for row in csvloadme:
                loadme.append(row)
        elif tp == 'json':
            loadme = json.load(data)
        elif tp == 'jsons':
            loadme = json.loads(data)
        else:
            print 'Need defined filetype'
    return loadme

def extract_acs_sum(fname):
    # print 'Extracting: %s'%fname
    typ = fname[9:10]
    st_nm = fname[15:17]
    seq_nm = fname[17:21]
    acs = []
    with open(fname, 'rb') as f:
        for row in csv.reader(f, delimiter=','):
            acs.append(row)
    # Get Header info
    with open('./2013_5yr_SAS/%s%s_%s.sas'%(typ,st_nm,seq_nm)) as f:
        read_flg = False
        header = []
        for ln in f:
            if ln.strip() not in ('', ';', 'RUN;') and read_flg:
                header.append(re.split('\s*', ln.strip())[0])
            read_flg = True if 'INPUT' in ln else read_flg
    # Get geography
    with open('./tab4/sumfile/prod/2009thru2013/geo/g20135%s.txt'%st_nm,'rb') as f:
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

# def build_custom(varlist, score):
#     print '*****************************************************************************************************************************'
#     print 'Building Custom List:'
#     print '*****************************************************************************************************************************'
#     directory = build_directory('./acs_sas/')
#     custom = {}
#     header = []
#     for var in tqdm(varlist):
#         header.extend([var])
#         if var in score:
#             header.extend(['score_' + var])
#         for table in directory[var]:
#             typ = table[0]
#             st_nm = table[1]
#             seq_nm = table[2]
#             acs_summary = extract_acs_sum('./acs_summ/' + typ + '20135' + st_nm + seq_nm + '000.txt')
#             position = [ i for i,x in enumerate(acs_summary[0]) if x == var ]
#             if var in score:
#                 scoreme = [ x[position[0]] for x in acs_summary[1:] ]
#                 scored = [ stats.percentileofscore(scoreme, x) for x in scoreme ]
#             for i, row in enumerate(acs_summary[1:]):
#                 try:
#                     custom[row[-1]].extend([row[position[0]]])
#                 except KeyError:
#                     custom[row[-1]] = [row[position[0]]]
#                 if var in score:
#                     custom[row[-1]].extend([scored[i]])
#     position = [ i for i,x in enumerate(header) if 'score' in x ]
#     for tract in custom.keys():
#         avg = sum([ float(custom[tract][x]) for x in position ])/len(position)
#         custom[tract].extend([avg])
#     header.extend(['final_score'])
#     custom['varnames'] = header
#     with open('./custom.json', 'wb') as f:
#         json.dump(custom, f)
#     return custom

def build_custom_json(varlist, score, force=False):
    with open('header.json') as f:
        existing_score = json.load(f).keys()
    with open('header.json') as f:
        existing_variables = json.load(f).keys()
    if score != existing_score or varlist != existing_variables:
        force = True
    if force or not os.path.isfile('custom.json'):
        print '*****************************************************************************************************************************'
        print 'Building Custom JSON:'
        print '*****************************************************************************************************************************'
        directory = build_directory('./2013_5yr_SAS/')
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
                acs_summary = extract_acs_sum('./group2/' + typ + '20135' + st_nm + seq_nm + '000.txt')
                position = [ i for i,x in enumerate(acs_summary[0]) if x == var ]
                if var in score:
                    scoreme = [ x[position[0]] for x in acs_summary[1:] ]
                    scored = [ stats.percentileofscore(scoreme, x) for x in scoreme ]
                for i, row in enumerate(acs_summary[1:]):
                    try:
                        custom[row[-1][7:]][var] = row[position[0]]
                    except KeyError:
                        custom[row[-1][7:]] = {}
                        custom[row[-1][7:]][var] = row[position[0]]
                    if var in score:
                        custom[row[-1][7:]]['score_'+var] = scored[i]
        position = [ x for x in header if 'score' in x ]
        for tract in custom.keys():
            avg = sum([ float(custom[tract][x]) for x in position ])/len(position)
            custom[tract]['final_score'] = avg
        header.extend(['final_score'])
        # custom['varnames'] = header
        with open('custom.json', 'wb') as f:
            json.dump(custom, f)
    with open('custom.json', 'rb') as f:
        custom = json.load(f)
    return custom

if __name__=='__main__':
    ########################################################################################################################################
    # CFPB Complaints Database
    # ZCTA to Tract Relationship file here: http://www2.census.gov/geo/docs/maps-data/data/rel/zcta_tract_rel_10.txt
    ########################################################################################################################################
    print 'Loading Complaints:'
    complaints = loadinput('./Consumer_Complaints.csv','csv')
    compl_list = []
    for complaint in tqdm(complaints[1:]):
        compl_list.append(complaint[9])
    print 'Complete'

    print 'Counting Complaints:'
    compl_cnt = {}
    # This Craps out for some reason
    # compl_cnt = {x:compl_list.count(x) for x in compl_list}
    for zips in tqdm(compl_list):
        if 'X' not in zips:
            if compl_cnt.has_key(zips):
                compl_cnt[zips] += 1
            else:
                compl_cnt[zips] = 1
    print 'Complete'

    print 'Migrating to Census Tract:'
    ziptract = [ [(x[0],x[4]),float(x[18])/100] for x in loadinput('./zcta_tract_rel_10.txt', 'csv')[1:] ]
    for i, row in tqdm(enumerate(ziptract)):
        zipcode = row[0][0]
        if compl_cnt.has_key(zipcode):
            ziptract[i].append(row[1]*compl_cnt[zipcode])
        else:
            ziptract[i].append(0)
    print 'Complete'

    print 'Counting Census Tract:'
    compl_append = {}
    for row in tqdm(ziptract):
        tract = row[0][1]
        if compl_append.has_key(tract) == False:
            compl_append[tract] = [row[2], 1]
        else:
            compl_append[tract][0] += row[2]
            compl_append[tract][1] += 1
    print 'Complete'
    ########################################################################################################################################
    # ACS 5-year summary files
    ########################################################################################################################################
    # Variable JSON
    # http://api.census.gov/data/2014/acs5/variables.json
    ########################################################################################################################################
    # SAS Table generation code
    # http://www2.census.gov/programs-surveys/acs/summary_file/2013/documentation/user_tools/2013_5yr_SAS.zip
    ########################################################################################################################################
    # Geography Codes
    # http://www2.census.gov/programs-surveys/acs/summary_file/2013/data/5_year_entire_sf/2013_ACS_Geography_Files.zip
    ########################################################################################################################################
    # Tract level ACS 5-year estimates
    # http://www2.census.gov/programs-surveys/acs/summary_file/2013/data/5_year_entire_sf/Tracts_Block_Groups_Only.tar.gz
    ########################################################################################################################################
    # Test pull
    ########################################################################################################################################
    # variables = ['B17001e2', 'B22002e2', 'B23025e7']
    # scores = ['B17001e2', 'B22002e2', 'B23025e7']
    ########################################################################################################################################
    # Big pull
    ########################################################################################################################################
    # variables = ['B00001e1', 'B00002e1', 'B01001e2', 'B01001e26', 'B01001e27', 'B01001e28', 'B01001e29', 'B01001e30', 'B02001e1', 'B02001e2', 'B02001e3', 'B02001e4', 'B02001e5', 'B02001e6', 'B02001e7', 'B02001e8', 'B02001e9', 'B02001e10', 'B17001e1', 'B17001e2', 'B18101e1', 'B18105e1', 'B21100e1', 'B21100e2', 'B21100e3', 'B22001e1', 'B22001e2', 'B22001e3', 'B22001e4', 'B22002e1', 'B22002e2', 'B22002e3', 'B22002e4', 'B22002e5', 'B22002e6', 'B22002e7', 'B22002e8', 'B22002e9', 'B22002e10', 'B22002e11', 'B22002e12', 'B22002e13', 'B22002e14', 'B23024e1', 'B23024e2', 'B23024e3', 'B23024e4', 'B23024e5', 'B23024e6', 'B23024e7', 'B23024e8', 'B23024e9', 'B23025e1', 'B23025e5', 'B23025e7']
    # scores = ['B17001e2', 'B22002e2', 'B23025e7']
    ########################################################################################################################################
    # Tutorial pull
    ########################################################################################################################################
    # Girls under 5, Same house 1 year ago, SNAP benefits the past year
    # Merge against Financial Complaints @ http://www.consumerfinance.gov/complaintdatabase/
    ########################################################################################################################################
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
    ########################################################################################################################################
    # Building a national tract level geojson
    # Different vintages from https://www.census.gov/geo/maps-data/data/cbf/cbf_tracts.html
    # We use 2010 decennial vintage
    ########################################################################################################################################
    rootdir = './tracts/'
    us_features = []
    for subdir, dirs, files in os.walk(rootdir):
        for f in tqdm(files):
            ################################################################################################################################
            # Process only compressed files
            ################################################################################################################################
            if f[-3:] == 'zip':
                features = []
                ############################################################################################################################
                # Census provides single layer tract shapefile so we take the first
                ############################################################################################################################
                with fiona.open('/', vfs='zip://'+os.path.join(subdir, f), layer=0) as src:
                    for feat in src:
                        ####################################################################################################################
                        # Toss out extra data to save space
                        ####################################################################################################################
                        tract = feat['properties']['GEO_ID'][9:]
                        if tract in custom.keys():
                            for popme in [u'NAME', u'LSAD', u'STATE', u'COUNTY', u'TRACT', u'CENSUSAREA']:
                                feat['properties'].pop(popme)
                            for rank in variables:
                                feat['properties']['s'+rank] = custom[tract]['score_'+rank]
                            feat['properties']['r'] = custom[tract]['final_score']
                            try:
                                feat['properties']['f'] = compl_append[tract][0]
                            except KeyError:
                                feat['properties']['f'] = 0
                            features.append(feat)
                            us_features.append(feat)
                rank_map = {
                    'type':'FeatureCollection',
                    'features':features,
                    'crs':{'init': u'epsg:4269'}}
                with open('./tracts/%s.geojson'%f, 'wb') as f:
                    json.dump(rank_map, f)
    rank_map = {
        'type':'FeatureCollection',
        'features':us_features,
        'crs':{'init': u'epsg:4269'}}
    with open('./tracts/us_all_tracts.geojson', 'wb') as f:
        json.dump(rank_map, f)
    os.system("topojson -o ./tracts/us_all_tracts.topojson ./tracts/us_all_tracts.geojson")