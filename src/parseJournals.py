#!/usr/bin/env python
# coding: utf-8
__author__ = "Tariq Faquih"
__copyright__ = "Copyright 2020, Clinical Epidemiology Department, LUMC"
__credits__ = ["Tariq Faquih", "Linda Nab", "Ype Jong"]
__license__ = "MIT License"
__maintainer__ = "Tariq Faquiih"
__email__ = "t.o.faquih@lumc.nl"
__status__ = "Development"
# In[21]:


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


# In[22]:



def BeautifulSoup_parse (URL):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:66.0) Gecko/20100101 Firefox/66.0",
        "Accept-Encoding": "*",
        "Connection": "keep-alive"
    }
    htmlpage = requests.get(URL, headers=headers)
    
    print('Reading URL: {}'.format(URL))
    page = BeautifulSoup(htmlpage.text, "html.parser")
    return(page)


# In[23]:


def WriteJsonOutput(json_file , dict_file , jsondict):
#Write the main json file
    with open(json_file , 'w') as fp:   
        Output=[]
        for Key,Item in jsondict.items():
            Output.append(Item)
        json.dump(Output, fp)

    print('JSON {} write done'.format(json_file))
    
    #Write the exact dict as json file (easily read by the script)
    with open(dict_file , 'w') as fp:   
        json.dump(jsondict, fp)
        
    print('JSON {} write done'.format(dict_file))
    
    return(json_file)


# In[24]:


def nejm_parse(page, json_out_path ):
    txt = []
    link = []
    ID = []
    AB = []
    Dates = []
    Type = []

    print('Scraping NEJMA')
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

        
    json_file = '{}nejm.json'.format(json_out_path)
    dict_file = '{}nejm_dict_file.json'.format(json_out_path)

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
    WriteJsonOutput(json_file , dict_file , jsondict)

    
    return(json_file)


# In[25]:


def jama_parse(page, json_out_path ):

    txt = []
    link = []
    ID = []
    AB = []
    Dates = []
    Type = []
    print('Scraping JAMA')
    for entry in page.find_all("a", attrs={"class": "article--title d-b mb05"}):
        txt.append(entry.text)
        link.append(entry.get('href'))
        ID.append(entry.get('href').split('/')[-1])

    for abstract in page.find_all("div", attrs={"class": "article--excerpt"}): 
        if abstract.p != None:
            AB.append(FormatAbstract(abstract.p.text))
        else:
            AB.append('NA')

    for Date in page.find_all("div", attrs={"class": "article--date meta-item no-wrap"}): 
        D = datetime.strptime(Date.text, '%B %d, %Y')
        Dates.append(D.strftime("%Y%m%d"))

    for typemeta in page.find_all("div", attrs={"class": "article--type meta-item"}): 
        Type.append(typemeta.text)

    
    
    json_file = '{}jama.json'.format(json_out_path)
    dict_file = '{}jama_dict_file.json'.format(json_out_path)


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


    WriteJsonOutput(json_file , dict_file , jsondict)

    
    return(json_file)


# In[26]:


def lancent_parse(page, json_out_path ):

    txt = []
    link = []
    ID = []
    AB = []
    Dates = []
    Type = []
    print('Scraping LANCENT')
    mydiv = page.select('div.widget-body.body.body-regular')[0]

    for entry in mydiv.find_all("h4", attrs={"class": "title"}):
        txt.append(entry.text)

    for entry in mydiv.find_all("div", attrs={"class": "articleType"}):
        Type.append(entry.text)

    for entry in mydiv.find_all("div", attrs={"class": "published-online"}):
        Pdate = entry.text.split('Published: ')[-1]
        Pdate = datetime.strptime(Pdate , '%B %d, %Y')
        Dates.append(Pdate.strftime("%Y%m%d"))

    for entry in mydiv.find_all("div", attrs={"class": "doi"}):

        ID.append(entry.text.split('/')[-1])
        link.append(entry.text.split('DOI: ')[-1])

    json_file = '{}lancent.json'.format(json_out_path)
    dict_file = '{}lancent_dict_file.json'.format(json_out_path)


    #read the stored dictionary file (dict_file) or create a new blank dictionary
    if os.path.isfile(dict_file):
        jsondict = json.load(open(dict_file))
    else:
        jsondict = {}    


    for N in range(0, len(ID)) :
        if ID[N] not in jsondict.keys():
            jsondict[ID[N]] = {'ID' : ID[N] , 
                              'txt' : txt[N],
                              'link' : link[N],
                               'Dates' : Dates[N],
                               'Type' : Type[N]
                              }



    WriteJsonOutput(json_file , dict_file , jsondict)

    
    return(json_file)


# In[27]:


def MakeTemplate(json_file):
    headerslist = []
    for H in ('ID' ,'txt' ,'AB' ,'link' ,'Dates' ,'Type' ):
        headers = '=ImportJSON("https://raw.githubusercontent.com/tofaquih/coronaPubGet/master/jsonfiles/{}", "/{}", "noInherit,noTruncate",$A$1)'.format(json_file , H)
        headerslist.append(headers)
    
    TemplateFileName = './templates/{}_template.csv'.format(json_file.split('/')[-1])
    with open(TemplateFileName ,'w'  , newline='' ) as fp:
        W = csv.writer(fp, delimiter=';')
        W.writerow(headerslist)
    print('Template file {} written'.format(TemplateFileName))    


# In[28]:


import requests
import json
import csv,os,sys
from datetime import datetime
from bs4 import BeautifulSoup

if __name__ == '__main__':
    #content = requests.get(url)
    #page = BeautifulSoup(open('Coronavirus (COVID19) _ JAMA Network.html'), "html.parser")

    nejmpage = "https://www.nejm.org/coronavirus?query=main_nav_lg"
    jamapage = "https://jamanetwork.com/collections/46099/coronavirus-covid19?appId=scweb&fl_ContentType=Article&fl_Categories=Coronavirus+(COVID19)"
    lancetpage = "https://www.thelancet.com/coronavirus"

    json_out_path = './jsonfiles/'


    # In[18]:


    MakeTemplate(nejm_parse(BeautifulSoup_parse(nejmpage) , json_out_path))
    MakeTemplate(jama_parse(BeautifulSoup_parse(jamapage) , json_out_path))
    MakeTemplate(lancent_parse(BeautifulSoup_parse(lancetpage) , json_out_path))

    #nejm_parse(pages)

