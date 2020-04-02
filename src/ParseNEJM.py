#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests
import json
import csv,os,sys
from datetime import datetime
from bs4 import BeautifulSoup


#content = requests.get(url)
#page = BeautifulSoup(open('Coronavirus (COVID19) _ JAMA Network.html'), "html.parser")

nejmpage = "https://www.nejm.org/coronavirus?query=main_nav_lg"
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:66.0) Gecko/20100101 Firefox/66.0",
    "Accept-Encoding": "*",
    "Connection": "keep-alive"
}
htmlpage = requests.get(nejmpage, headers=headers)

page = BeautifulSoup(htmlpage.text, "html.parser")


# In[2]:


def FormatAbstract (AB):

    Abstract = ''
    tempAbs = AB.split(' ')
    lastN=0
    newabs = []
    for N in range(30, len(tempAbs)+30, 30) :
        newabs.append(' '.join(tempAbs[lastN:N]))
        lastN= N

    Abstract = '\n'.join(newabs)

    return(Abstract)


# In[3]:


txt = []
link = []
ID = []
AB = []
Dates = []
Type = []

mydiv = page.find("div", attrs={"class": "col-md-2-3 o-col o-col--primary o-col--keyline-top"})

for entry in mydiv.find_all("b", attrs={"class": "m-article__title"}):
    txt.append(entry.text)
for entry in mydiv.find_all("a", attrs={"class": "m-article__type"}):
    Type.append(entry.text)
for entry in mydiv.find_all("a", attrs={"class": "m-article__link"}):
    link.append('https://www.nejm.org' + entry.get('href')) 
    articalid = entry.get('href').split('/')[-1].split('?')[0]
    ID.append(articalid)
for entry in mydiv.find_all("span", attrs={"class": "m-article__blurb"}):
    if entry.text != None:
        abstract = FormatAbstract(entry.text.replace('\n' , ''))
        AB.append(abstract)
    else:
        AB.append('NA')

for entry in mydiv.find_all("em", attrs={"class": "m-article__date"}):
    D = datetime.strptime(entry.text + ' 2020', '%b %d %Y')
    Dates.append(D.strftime("%Y%m%d"))

    


# In[4]:


json_file = '../jsonfiles/nejm.json'
dict_file = '../jsonfiles/nejm_dict_file.json'


#read the stored dictionary file (dict_file) or create a new blank dictionary
if os.path.isfile(dict_file):
    jsondict = json.load(open(dict_file))
else:
    jsondict = {}    
    
    
for N in range(0, len(ID)) :
    if ID[N] not in jsondict.keys():
        jsondict[ID[N]] = {'ID' : ID[N] , 
                          'txt' : txt[N],
                          'AB' : AB[N],
                          'link' : link[N],
                           'Dates' : Dates[N],
                           'Type' : Type[N]
                          }
jsondict


#Write the main json file
with open(json_file , 'w') as fp:   
    Output=[]
    for Key,Item in jsondict.items():
        Output.append(Item)
    json.dump(Output, fp)

#Write the exact dict as json file (easily read by the script)
with open(dict_file , 'w') as fp:   
    json.dump(jsondict, fp)


# In[5]:


def MakeTemplate():
    headerslist = []
    for H in ('ID' ,'txt' ,'AB' ,'link' ,'Dates' ,'Type' ):

        print(H)
        headers = '=ImportJSON("https://raw.githubusercontent.com/tofaquih/coronaPubGet/master/jsonfiles/nejm.json", "/{}", "noInherit,noTruncate",$A$1)'.format(H)
        headerslist.append(headers)

    headerslist
    with open('../templates/nejm_template.csv' ,'w'  , newline='' ) as fp:
        W = csv.writer(fp, delimiter=';')
        W.writerow(headerslist)
        


# In[6]:


MakeTemplate()

