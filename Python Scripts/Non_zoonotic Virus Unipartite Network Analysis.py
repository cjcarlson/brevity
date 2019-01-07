#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:47:30 2018

@author: Casey

Take Olival data and create unipartite host association networks for non-zoonotic viruses only
The nodes in the networks are the orders of the hosts i.e.'CARNIVORA' 'CETARTIODACTYLA' 'CHIROPTERA' 'CINGULATA' 'DIDELPHIMORPHIA' 'DIPROTODONTIA' 'EULIPOTYPHLA' 'LAGOMORPHA' 'PERAMELEMORPHIA' 'PERISSODACTYLA' 'PILOSA' 'PRIMATES' 'PROBOSCIDEA' 'RODENTIA' 'SCANDENTIA'
I go through each association in the Olival association data. For each ZOONOTIC virus, I get the orders of all of the hosts that the virus is associated with. I create an edge between all orders that have the same virus. 
If an edge already exists between two orders, I increase the weight of the edge by one for each repeated association between the same two orders
This is only for zoonotic viruses
Result:
    Network where nodes are orders, edges are viral associations (a virus is found in both of the orders that are the nodes), the weight of the edge indicates the number of viruses that are shared between two orders
    For zoonotic viruses only
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
    #create empty nets
    all_host_association_net = nx.Graph()    
    #read in olival data
    #olival_associations = pd.read_csv('~/GitHub/brevity/Olival Nature 2017 Raw Data/associations.csv')
    #olival_hosts = pd.read_csv('~/GitHub/brevity/Olival Nature 2017 Raw Data/hosts.csv')
    #olival_viruses = pd.read_csv('~/GitHub/brevity/Olival Nature 2017 Raw Data/viruses.csv')
    
    olival_associations = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/Olival Nature 2017 Raw Data/associations.csv')
    olival_hosts = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/Olival Nature 2017 Raw Data/hosts.csv')
    olival_viruses = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/Olival Nature 2017 Raw Data/viruses.csv')
    
    
    return all_host_association_net, olival_associations, olival_hosts, olival_viruses

def host_and_order_data(all_host_data, all_association_data, all_virus_data, all_host_association_net):
    """
    deal with the host and order data and add nodes to the network for each order
    input: host data, association data, empty network
    output: dictionary of hosts and their orders, list of unique viruses, network with nodes but no edges, association data
    """
    
    #create dict of virusus and whether they are zoonotic or not
    virus_data = all_virus_data[['vVirusNameCorrected', 'IsZoonotic']]
    virus_data.columns = ['virus', 'is_zoonotic']
    is_virus_zoonotic_dict = {}
    for each_entry in range(0, len(virus_data['virus'])):
        virus_name = virus_data['virus'][each_entry]
        is_virus_zoonotic= virus_data['is_zoonotic'][each_entry]
        is_virus_zoonotic_dict[virus_name] = is_virus_zoonotic
    
    host_data = all_host_data[['hHostNameFinal','hOrder']] # just get the host name and the order from the olival host data
    host_data.columns = ['host_name', 'order'] # change column names to be easier
    unique_hosts = np.unique(host_data['host_name']) # get a list of unique hosts in the host data
    unique_orders = np.unique(host_data['order']) # get a list of the unique orders in the host data
    
    for each_order in unique_orders: # loop thorugh each order and create a node in each network for every order 
        all_host_association_net.add_node(each_order)
    
    associations = all_association_data[['vVirusNameCorrected','hHostNameFinal', 'WildDomInReference']] # get the virus name, the host name, and whether the host sampled was wild or domestic
    associations.columns = ['virus', 'host', 'wild_or_dom']  # change the column names to be easier 
    unique_viruses = np.unique(associations['virus']) # get a list of unique viruses in the dataset
    host_order_dict = {} # create a dictionary. Will host host:order -  so I can reference what order each host is in
    unique_hosts = list(unique_hosts)
    human_index= unique_hosts.index('Homo_sapiens') # remove humans
    del(unique_hosts[human_index])
    associations= associations[associations.host!= 'Homo_sapiens']
   # create a dictionary. Will host host:order -  so I can reference what order each host is in
    for each_host in unique_hosts: # go through each host
        each_host_data = host_data.loc[host_data['host_name'] == each_host] # find the data on each host in the host data 
        host_order = each_host_data['order'] # get the order that the host is in 
        host_order = host_order.item() # just the order (get rid of the row index)
        host_order_dict[each_host] = host_order # add the host and order to the dictionary as host: order
    
    return host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict 


def add_edges_to_host_association_net(host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict ):

    for each_virus in unique_viruses:# go through each virus in the data
        # NON- zoonotic viruses only!!!!
        if is_virus_zoonotic_dict[each_virus] == 0: 
            
            all_virus_data = associations.loc[associations['virus'] ==each_virus]# get the association data for each virus for all hosts
            
            host_order_list = []
            for each_host in all_virus_data['host']:
                
                host_order_list.append(host_order_dict[each_host])
            
            
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
                    
#    #remove nodes with no edges
    for each_node in nx.nodes(all_host_association_net):
        
        each_node_neighbors = all_host_association_net.neighbors(each_node)
        if len(each_node_neighbors) == 0:
            
            all_host_association_net.remove_node(each_node)
#    #remove self-loops        
    all_host_association_net.remove_edges_from(all_host_association_net.selfloop_edges())
    
    return all_host_association_net
          

if __name__ == '__main__':

    all_host_association_net, olival_associations, olival_hosts, olival_viruses = read_in_files_and_create_nets()
    host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict  = host_and_order_data(olival_hosts, olival_associations, olival_viruses, all_host_association_net)
    non_zoonotic_host_association_net = add_edges_to_host_association_net(host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict )
    
    
    
#    id_for_file = range(0, len(nx.nodes(all_host_association_net)))
#    print id_for_file
#    for each_id in id_for_file:
#        node = nx.nodes(non_zoonotic_host_association_net)[each_id]
#        print node
#    
#        filename = '/Users/admin/Dropbox (Bansal Lab)/brevity_project/non_zoonotic_nodelist.csv'
#        f = open(filename, "a+")
#        f.write(str(each_id)+","+str(node)+"\n")
#        f.close()
#        
#    nx.write_edgelist(non_zoonotic_host_association_net, '/Users/admin/Dropbox (Bansal Lab)/brevity_project/non_zoonotic_edgelist.csv', data = ['weight'], delimiter = ',')
    
#    
#    #calculating network statistics: degree, eigenvector centrality, average node strength
#
    for each_node in nx.nodes(non_zoonotic_host_association_net):
        node_degree = non_zoonotic_host_association_net.degree(each_node)
        neighbor_weights = []
        node_neighbors = non_zoonotic_host_association_net.neighbors(each_node)
        
        for each_node_neighbor in node_neighbors:
            edge_weight = non_zoonotic_host_association_net[each_node][each_node_neighbor]['weight']
            #print each_node, each_node_neighbor, edge_weight
            neighbor_weights.append(edge_weight)
        node_strength = sum(neighbor_weights)
        
        if node_degree >0:
            average_node_strength = float(node_strength)/float(node_degree)
            print each_node
            print average_node_strength
            print node_degree
    print nx.eigenvector_centrality(non_zoonotic_host_association_net, 100)


