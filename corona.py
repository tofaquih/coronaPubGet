#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python


# In[2]:


__author__ = "Tariq Faquih"
__copyright__ = "Copyright 2020, Clinical Epidemiology Department, LUMC"
__credits__ = ["Tariq Faquih", "Linda Nab", "Ype Jong"]
__license__ = "MIT License"
__maintainer__ = "Tariq Faquiih"
__email__ = "t.o.faquih@lumc.nl"
__status__ = "Development"


# # Import modules

# In[3]:


import json, csv , os , sys , datetime
from Bio import Entrez
from Bio import Medline


# # COVID class

# This class objects send pubmed queries through the API and stores the output in json files

# In[4]:



class COVID:
 
    def __init__(self , startD , endD):

        #list of query terms to be used in the search_function

        self.querydict = {'COVID':None , 
                 'Big Five':'(NEJM[journal] OR BMJ[journal] OR lancet[journal] OR nature[journal] OR JAMA[journal])', 
                 'Elderly':'elderly[TITLE]',
                 'Clinical Trial':'clinical trial[Title/Abstract]' , 
                 'Italy':'italy[Title/Abstract]' , 
                 'Netherlands':'netherlands[Title/Abstract]' , 
                 'Case Control':'case control study' , 
                 'Epidemiology':'epidemiology' , 
                 'Mortality':'mortality', 
                 'Pregnant':'pregnant[TITLE]',
                 'Treatment':'(treatment[All Fields] OR drug[All Fields] OR intervention[All Fields] OR recovery[All Fields])' }
        
        #json_file stores the proper json format to be used in the googlesheet
        #dict_file stores the output in a dictionary to be loaded in later uses
        json_file ='results3.json'
        dict_file ='results_dictionary3.json'

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
            self.search_function(K , startD , endD )
            
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
        
        with open('log_{}_{}.txt'.format(startD.replace('/' , '') , endD.replace('/' , '')) ,'w'  , newline='') as fp:
            W = csv.writer(fp)
            W.writerow(['Search results for range {} to {}'.format(startD , endD)])
            W.writerow(['Number of Records Added: {}'.format(self.Log[-1])])
            for Line in self.Log[:-1]:
                W.writerow([Line])
        

    def search_function (self , MyTerms, startD , endD):

        Entrez.email = "tariqf549@gmail.com"
        MainTerm = """(((("coronavirus"[MeSH Terms] OR "coronavirus"[All Fields]) AND ("COVID-19"[All Fields] OR "severe acute respiratory syndrome coronavirus 2"[Supplementary Concept] OR "severe acute respiratory syndrome coronavirus 2"[All Fields] OR "2019-nCoV"[All Fields] OR "SARS-CoV-2"[All Fields] OR "2019nCoV"[All Fields] ))))"""

        #MainTerm = '"COVID-19"[All Fields]'
        DateRange = '"{}"[MHDA] : "{}"[MHDA]'.format(startD , endD)
        if MyTerms == 'COVID':
            Query = MainTerm + ' AND ' + DateRange
        else:
            Query = MainTerm + ' AND ' + self.querydict[MyTerms] + ' AND ' + DateRange
            
        print(Query)
        self.Log.append(Query)
        search_results = Entrez.read(
            Entrez.esearch(
                db="pubmed", term=Query,  datetype="pdat", usehistory="y" , sort = 'relevance' 
            )
        )
        count = int(search_results["Count"])
        self.Log.append("Found %i results" % count)
        
        print("Found %i results" % count)

        batch_size = 10
        out_handle = open("pubmed_results/corona_{}_papers.txt".format(MyTerms), "w")
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            print("Going to download record %i to %i" % (start + 1, end))
            self.Log.append("Going to download record %i to %i" % (start + 1, end))
            fetch_handle = Entrez.efetch(
                db="pubmed",
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
            PMID = parse_res['PMID']
            Tag = Q
            if PMID in self.mainDict.keys():
                print('PMID [{}] exists'.format(PMID))
                self.Log.append('PMID exists')
                if Q not in self.mainDict[PMID]['Tag']:
                    self.mainDict[PMID]['Tag'] = self.mainDict[PMID]['Tag']+';'+Tag
                    
                if self.mainDict[PMID]['Abstract'] == 'NA':
                    NewABS = self.FormatAbstract(parse_res)
                    if NewABS != 'NA':
                        self.mainDict[PMID]['Abstract'] = self.FormatAbstract(parse_res)
                        print('Updated Abstract')
                        self.Log.append('Updated Abstract')
                        

                    
                continue
            else:
                Title = parse_res['TI']
                dateP = parse_res['DP']
                dateC = parse_res['MHDA']

                if 'JT' in parse_res.keys():
                    JournalName  = parse_res['JT']
                else: JournalName = ''

                #if 'LR' in parse_res.keys():
                #    dateMod  = parse_res['LR']
                #else: dateMod = ''

                Abstract = self.FormatAbstract(parse_res)
                
                    
                Link= 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(PMID)   


                self.mainDict[PMID] = {'PMID': PMID, 'Title':Title ,
                                'JournalName':JournalName ,
                                  'Creation Date':dateC ,
                                'Publication Date':dateP , 
                                'Abstract':Abstract , 
                                'Link':Link,
                                'Tag':Tag   }
                self.NumNew +=1
                print('Added new PMID: {}'.format(PMID))
                self.Log.append('Added new PMID: {}'.format(PMID))


# In[5]:


if __name__ == '__main__':
    Today =datetime.datetime.now()
    StartDate = Today - datetime.timedelta(days=3)

    Today = Today.strftime("%Y/%m/%d")
    StartDate = StartDate.strftime("%Y/%m/%d")
    
    COVID(StartDate , Today )



# In[1]:


def MakeTemplate():
    headerslist = []
    for H in ('PMID',
              'Title',
              'Link',
                'JournalName',
              'Creation Date',
                'Publication Date',
                'Abstract',
                'Tag'):
        print(H)
        headers = '=ImportJSON("https://raw.githubusercontent.com/tofaquih/coronaPubGet/master/results3.json", "/{}", "noInherit,noTruncate",$A$1)'.format(H)
        headerslist.append(headers)

    headerslist
    with open('template.csv' ,'w'  , newline='' ) as fp:
        W = csv.writer(fp, delimiter=';')
        W.writerow(headerslist)

