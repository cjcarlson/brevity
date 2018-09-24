#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 13:36:31 2018

@author: Casey

This is where I do network things.



"""
import pandas as pd
import networkx as nx
import numpy as np
import itertools

wild_host_association_net = nx.Graph()
dom_host_association_net = nx.Graph()
all_host_association_net = nx.Graph()


olival_associations = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/olival nature 2017/associations.csv')
olival_hosts = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/olival nature 2017/hosts.csv')
olival_viruses = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/olival nature 2017/viruses.csv')
host_data = olival_hosts[['hHostNameFinal','hOrder']]
host_data.columns = ['host_name', 'order']

unique_hosts = np.unique(host_data['host_name'])

unique_orders = np.unique(host_data['order'])
for each_order in unique_orders:
    wild_host_association_net.add_node(each_order)
    dom_host_association_net.add_node(each_order)
    all_host_association_net.add_node(each_order)
    
associations = olival_associations[['vVirusNameCorrected','hHostNameFinal', 'WildDomInReference']]
associations.columns = ['virus', 'host', 'wild_or_dom']

unique_viruses = np.unique(associations['virus'])
host_order_dict = {}
for each_host in unique_hosts[0:10]:
    each_host_data = host_data.loc[host_data['host_name'] == each_host]
    host_order = each_host_data['order']
    host_order.to_string()
    host_order_dict[each_host] = host_order
print host_order_dict
    

#for each_virus in unique_viruses[0:10]:
#    #print each_virus
#    wild_virus_data = associations.loc[(associations['virus'] ==each_virus) & (associations['wild_or_dom']=='wild')]
#    #print wild_virus_data
#    #print wild_virus_data['host']
#    host_combinations = list(itertools.combinations(list(wild_virus_data['host']), 2))
#    for each_combination in host_combinations:
#        for each_host in each_combination:
            
    #dom_virus_data =  associations.loc[(associations['virus'] ==each_virus) & (associations['wild_or_dom']=='domestic')]
    #all_virus_data =  associations.loc[(associations['virus'] ==each_virus)]
    #print wild_virus_data['host']

#print list(itertools.combinations([1, 2, 3], 2))
















