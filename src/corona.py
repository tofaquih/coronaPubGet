#!/usr/bin/env python
# coding: utf-8


__author__ = "Tariq Faquih"
__copyright__ = "Copyright 2020, Clinical Epidemiology Department, LUMC"
__credits__ = ["Tariq Faquih", "Linda Nab", "Ype Jong"]
__license__ = "MIT License"
__maintainer__ = "Tariq Faquiih"
__email__ = "t.o.faquih@lumc.nl"
__status__ = "Development"


# # Import modules

# In[1]:


import json, csv , os , sys , datetime
from Bio import Entrez
from Bio import Medline
from datetime import datetime
from datetime import timedelta
import warnings



# # COVID class

# This class objects send pubmed queries through the API and stores the output in json files

# In[56]:



class COVID:
 
    def __init__(self , startD , endD , json_file , dict_file , Eemail , apikey):

        #list of query terms to be used in the search_function

        self.querydict = json.load(open("./jsonfiles/pubmed_search_terms.json"))[0]
        
        self.Eemail = Eemail
        self.apikey = apikey
        
        #json_file stores the proper json format to be used in the googlesheet
        #dict_file stores the output in a dictionary to be loaded in later uses


        #read the stored dictionary file (dict_file) or create a new blank dictionary
        if os.path.isfile(dict_file):
            self.mainDict = json.load(open(dict_file))
        else:
            self.mainDict = {}
        
        #set counter for how many articles are added and create a list to store the log messages
        self.NumNew = 0
        self.Log = []
        
        #for each query term in the querylist, run the search function using the provideed start
        #and end dates
        for K in self.querydict.keys():
            if K != 'MainTerm':
                for DB in ('pubmed' , 'pmc'):
                    self.search_function(K , DB , startD , endD )
            
        #for X in querylist:
        #    self.search_function(X , startD , endD )
        
        #add the total number of added articles to the log list
        self.Log.append(self.NumNew)
        
        #Write the main json file
        with open(json_file , 'w') as fp:   
            Output=[]
            for Key,Item in self.mainDict.items():
                Output.append(Item)
            json.dump(Output, fp)
            
        #Write the exact dict as json file (easily read by the script)
        with open(dict_file , 'w') as fp:   
            json.dump(self.mainDict, fp)
        
        with open('./logs/log_{}_{}.txt'.format(startD.replace('/' , '') , endD.replace('/' , '')) ,'w'  , newline='') as fp:
            W = csv.writer(fp)
            W.writerow(['Search results for range {} to {}'.format(startD , endD)])
            W.writerow(['Number of Records Added: {}'.format(self.Log[-1])])
            for Line in self.Log[:-1]:
                W.writerow([Line])
        

    def search_function (self , MyTerms, DB , startD , endD):

        
        if self.Eemail != None:
            Entrez.email = self.Eemail
            
        if self.apikey != None:
            Entrez.api_key = self.apikey

        #MainTerm = '"COVID-19"[All Fields]'
        DateRange = '"{}"[PDAT] : "{}"[PDAT]'.format(startD , endD)
        if MyTerms == 'Default':
            Query = self.querydict['MainTerm'] + ' AND ' + DateRange
        else:
            Query = self.querydict['MainTerm'] + ' AND ' + self.querydict[MyTerms] + ' AND ' + DateRange
            
        print(Query)
        self.Log.append(Query)
        search_results = Entrez.read(
            Entrez.esearch(
                db=DB, term=Query,  datetype="pdat", usehistory="y" , sort = 'relevance' 
            )
        )
        count = int(search_results["Count"])
        self.Log.append("Found %i results" % count)
        
        print("Found %i results" % count)

        batch_size = 10
        out_handle = open("./pubmed_results/corona_{}_papers.txt".format(MyTerms), "w" , encoding="utf-8")
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            print("Going to download record %i to %i" % (start + 1, end))
            self.Log.append("Going to download record %i to %i" % (start + 1, end))
            fetch_handle = Entrez.efetch(
                db=DB,
                rettype="medline",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=search_results["WebEnv"],
                query_key=search_results["QueryKey"],
            )
            data = fetch_handle.read()

            dataresults = data.split('\nPMID')[1:]
            self.add2dict(dataresults , MyTerms)
            #print(data)
            fetch_handle.close()
            out_handle.write(data)
        out_handle.close()
        
    def FormatAbstract (self, AB):
        Abstract = ''
        if 'AB' in AB.keys():
            Abstract = AB['AB']
            tempAbs = Abstract.split(' ')
            lastN=0
            newabs = []
            for N in range(30, len(tempAbs)+30, 30) :
                newabs.append(' '.join(tempAbs[lastN:N]))
                lastN= N
                
            Abstract = '\n'.join(newabs)
                
        else: 
            Abstract = 'NA'
            
        return(Abstract)
            
    def add2dict(self , dataresults , Q):
        for hit in dataresults:
            m1 = 'PMID' + hit
            parse_res = Medline.read(m1.split('\n'))
            if 'PMID' in parse_res.keys():
                PMID = parse_res['PMID']
                
            elif 'PMCID' in parse_res.keys():
                PMID = parse_res['PMCID']
                
            else: PMID = ''
                
            Tag = Q
            if PMID in self.mainDict.keys():
                print('PMID [{}] exists'.format(PMID))
                self.Log.append('PMID exists')
                if Q not in self.mainDict[PMID]['Tag']:
                    self.mainDict[PMID]['Tag'] = self.mainDict[PMID]['Tag']+' '+Tag
                    
                if self.mainDict[PMID]['Abstract'] == 'NA':
                    NewABS = self.FormatAbstract(parse_res)
                    if NewABS != 'NA':
                        self.mainDict[PMID]['Abstract'] = self.FormatAbstract(parse_res)
                        print('Updated Abstract')
                        self.Log.append('Updated Abstract')
                        

                    
                continue
            else:
                if 'TI' in parse_res.keys(): 
                    Title = parse_res['TI']
                else: 
                    print('PMID {} has no title (yet). Skipping'.format(PMID))
                    continue
                
                if 'JT' in parse_res.keys():
                    JournalName  = parse_res['JT']
                else:
                    JournalName = ''
                    
                #Date of Publication [dp] - Date searching includes both print and electronic dates of publication. 
                #Searching for a single date does not include items when the electronic date of publication is after the print date.
                if 'DP' in parse_res.keys():
                    dateP  = parse_res['DP']
                else: dateP = ''
                
                #The date the citation first entered PubMed. In PMC DEP is used instead
                if 'EDAT' in parse_res.keys():
                    dateC  = datetime.strptime(parse_res['EDAT'] , '%Y/%m/%d %H:%S')
                    dateC = dateC.strftime("%Y%m%d")
                elif 'DEP' in parse_res.keys():
                    dateC = parse_res['DEP']
                else: dateC = ''

                Abstract = self.FormatAbstract(parse_res)
                
                    
                Link= 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(PMID)   


                self.mainDict[PMID] = {'PMID': PMID, 'Title':Title ,
                                'JournalName':JournalName ,
                                  'Date Added':dateC ,
                                'Publication Date':dateP , 
                                'Abstract':Abstract , 
                                'Link':Link,
                                'Tag':Tag   }
                self.NumNew +=1
                print('Added new PMID: {}'.format(PMID))
                self.Log.append('Added new PMID: {}'.format(PMID))


# In[60]:

def MakeTemplate():
    headerslist = []
    for H in ('PMID',
              'Title',
              'Link',
                'JournalName',
                'Date Added',
                'Publication Date',
                'Abstract',
                'Tag'):
        headers = '=ImportJSON("https://raw.githubusercontent.com/tofaquih/coronaPubGet/master/jsonfiles/results4.json", "/{}", "noInherit,noTruncate",$A$1)'.format(H)
        headerslist.append(headers)

    headerslist
    with open('./templates/template.csv' ,'w'  , newline='' ) as fp:
        W = csv.writer(fp, delimiter=';')
        W.writerow(headerslist)


        
if __name__ == '__main__':
    Today = datetime.now()
    StartDate = Today - timedelta(days=3)

    Today = Today.strftime("%Y/%m/%d")
    StartDate = StartDate.strftime("%Y/%m/%d")
    
    json_file ='./jsonfiles/results4.json'
    dict_file ='./jsonfiles/results_dictionary4.json'
    
    if os.path.isfile('credentials.json'):
        entrezinfo = json.load(open('credentials.json'))[0]
        Eemail = entrezinfo['email']
        apikey = entrezinfo['APIKEY']
    else:
        Eemail = None
        apikey = None
        warnings.warn("""File 'credentials.json' is missing. Thereofore no NCBI email will be used.
Please create the file 'credentials.json' at the root directory of the repo.
Example of the file format:
[{'email':"username@email.com" , 
 'APIKEY':'000000000000000000000000000000000000'}]""")
        
    COVID(StartDate , Today ,json_file , dict_file ,Eemail ,apikey )
    MakeTemplate()


    #COVID('2020/01/01' , '2020/03/30',json_file , dict_file )


