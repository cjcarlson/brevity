#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:20:18 2018

"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

@author: Casey

Take Olival data and create unipartite host association networks for all hosts.
The nodes in the networks are the orders of the hosts i.e.'CARNIVORA' 'CETARTIODACTYLA' 'CHIROPTERA' 'CINGULATA' 'DIDELPHIMORPHIA' 'DIPROTODONTIA' 'EULIPOTYPHLA' 'LAGOMORPHA' 'PERAMELEMORPHIA' 'PERISSODACTYLA' 'PILOSA' 'PRIMATES' 'PROBOSCIDEA' 'RODENTIA' 'SCANDENTIA'
The association data has virus, host, and wild or domestic
I go through each association in the Olival association data. For each virus, I get the orders of all of the hosts that the virus is associated with. I create an edge between all orders that have the same virus. 
If an edge already exists between two orders, I increase the weight of the edge by one for each repeated association between the same two orders
Result:
    Network where nodes are orders, edges are viral associations (a virus is found in both of the orders that are the nodes), the weight of the edge indicates the number of viruses that are shared between two orders
Calculate:
    Degree
    Eigenvector centrality
    Node strength

Note: This file uses Python 2.7 and NetworkX version 1.11

"""
import pandas as pd
import networkx as nx
import numpy as np
import itertools

def read_in_files_and_create_nets():
    """
    create empty network for all hosts
    read in olival data
    input: none
    output: Empty network and olival data (hosts, associations, viruses)
    """
    #create empty nets for 3 types of networks
    all_host_association_net = nx.Graph()    
    #read in olival data
    olival_associations = pd.read_csv('~/Github/brevity/Olival Nature 2017 Raw Data/associations.csv')
    olival_hosts = pd.read_csv('~/Github/brevity/Olival Nature 2017 Raw Data/hosts.csv')
    olival_viruses = pd.read_csv('~/Github/brevity/Olival Nature 2017 Raw Data/viruses.csv')
    
    return all_host_association_net, olival_associations, olival_hosts, olival_viruses

def host_and_order_data(all_host_data, all_association_data, all_host_association_net):
    """
    deal with the host and order data and add nodes to the network for each order
    input: host data, association data, empty network
    output: dictionary of hosts and their orders, list of unique viruses, network with nodes but no edges, association data
    """
    
    host_data = all_host_data[['hHostNameFinal','hOrder']] # just get the host name and the order from the olival host data
    host_data.columns = ['host_name', 'order'] # change column names to be easier
    unique_hosts = np.unique(host_data['host_name']) # get a list of unique hosts in the host data
    unique_orders = np.unique(host_data['order']) # get a list of the unique orders in the host data
    
    for each_order in unique_orders: # loop through each order and create a node in each network for every order
        all_host_association_net.add_node(each_order)
    
    associations = all_association_data[['vVirusNameCorrected','hHostNameFinal', 'WildDomInReference']] # get the virus name, the host name, and whether the host sampled was wild or domestic
    associations.columns = ['virus', 'host', 'wild_or_dom']  # change the column names to be easier
    
    unique_viruses = np.unique(associations['virus']) # get a list of unique viruses in the dataset
    
    host_order_dict = {} # create a dictionary. Will host host:order -  so I can reference what order each host is in
    unique_hosts = list(unique_hosts) # make unique hosts into a list from a dataframe so I can use the index function below
    human_index= unique_hosts.index('Homo_sapiens') # remove humans- ignoring humans in our sharing network
    del(unique_hosts[human_index])
    associations= associations[associations.host!= 'Homo_sapiens'] # remove human data from associations data
    
   # create a dictionary. Will host host:order -  so I can reference what order each host is in
    for each_host in unique_hosts: # go through each host
        each_host_data = host_data.loc[host_data['host_name'] == each_host] # find the data on each host in the host data 
        host_order = each_host_data['order'] # get the order that the host is in 
        host_order = host_order.item() # just the order (get rid of the row index)
        host_order_dict[each_host] = host_order # add the host and order to the dictionary as host: order
    return host_order_dict, unique_viruses, all_host_association_net, associations

def add_edges_to_host_association_net(host_order_dict, unique_viruses, all_host_association_net, associations):
    """
    Add edges to the network and weights to the edges that are already formed
    input: host-order dictionary, list of unique viruses, network with nodes but no edges, association data
    output: dictionary of hosts and their orders, list of unique viruses, network with nodes but no edges, association data
    """
    
   
    for each_virus in unique_viruses:# go through each virus in the data
            
        all_virus_data = associations.loc[associations['virus'] ==each_virus]# get the association data for each virus for all hosts
        
        
        host_order_list = [] # create a list for the order that each host is in
        for each_host in all_virus_data['host']:
            host_order_list.append(host_order_dict[each_host])
        #create all possible combinations of the hosts for each virus- these will be the edges because they all share a virus
        host_combinations = list(itertools.combinations(list(all_virus_data['host']), 2))     # create all possible 2 host combinations of hosts selected above     
        
        for each_combination in host_combinations:# go through each combination
            order_combination = []# an empty list to stick the orders of the hosts in
            for each_host in each_combination:# go through each host in the list of 2 hosts that are associated                   
                each_host_order = host_order_dict[each_host]# get the order of each host from the host- order dictionary  
                order_combination.append(each_host_order)   # put the 2 orders in the order combination list so I can make an edge between them   
               
            if all_host_association_net.has_edge(order_combination[0], order_combination[1]) == False:# if there is no edge betweent the orders
                all_host_association_net.add_edge(order_combination[0], order_combination[1], weight = 1)# create an edge with a weight of 1 between the orders
            if all_host_association_net.has_edge(order_combination[0], order_combination[1]) == True:# if there is already an edge between the 2 orders
                edge_weight = all_host_association_net[order_combination[0]][order_combination[1]]['weight'] # get the current weight of that edge
                updated_edge_weight = edge_weight+1# add 1 to the weight of the edge and save as updated weight
                all_host_association_net.add_edge(order_combination[0], order_combination[1], weight = updated_edge_weight)# redraw the edge with the new weight
    #remove nodes with no edges
    for each_node in nx.nodes(all_host_association_net):
        each_node_neighbors = all_host_association_net.neighbors(each_node)
        if len(each_node_neighbors) == 0:
            all_host_association_net.remove_node(each_node)
    #remove self-loops        
    all_host_association_net.remove_edges_from(all_host_association_net.selfloop_edges())
    
    return all_host_association_net
      
            

if __name__ == '__main__':
    all_host_association_net, olival_associations, olival_hosts, olival_viruses = read_in_files_and_create_nets()
    host_order_dict, unique_viruses, all_host_association_net, associations = host_and_order_data(olival_hosts, olival_associations, all_host_association_net)
    all_host_association_net = add_edges_to_host_association_net(host_order_dict, unique_viruses, all_host_association_net, associations)
    
    #calculating network statistics: degree, eigenvector centrality, average node strength

    for each_node in nx.nodes(all_host_association_net):
        node_degree = all_host_association_net.degree(each_node)
        neighbor_weights = []
        node_neighbors = all_host_association_net.neighbors(each_node)
        for each_node_neighbor in node_neighbors:
            edge_weight = all_host_association_net[each_node][each_node_neighbor]['weight']
            neighbor_weights.append(edge_weight)
        node_strength = sum(neighbor_weights)
        if node_degree >0:
            average_node_strength = float(node_strength)/float(node_degree)
    nx.eigenvector_centrality(all_host_association_net, 100)
            






