#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 13:32:15 2019

@author: admin
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 11:10:17 2019



Here: network. code allows for adding of human node, but there won't be a human node for the non-zoonotic network because that doesn't make sense
calculate number of samples in each order for node size
calculate Jaccard's index and use for edge weights



@author: admin
"""

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
    olival_associations = pd.read_csv('~/GitHub/brevity/Olival Nature 2017 Raw Data/associations.csv')
    olival_hosts = pd.read_csv('~/GitHub/brevity/Olival Nature 2017 Raw Data/hosts.csv')
    olival_viruses = pd.read_csv('~/GitHub/brevity/Olival Nature 2017 Raw Data/viruses.csv')

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
    all_host_association_net.add_node('HUMANS')
    
    associations = all_association_data[['vVirusNameCorrected','hHostNameFinal', 'WildDomInReference']] # get the virus name, the host name, and whether the host sampled was wild or domestic
    associations.columns = ['virus', 'host', 'wild_or_dom']  # change the column names to be easier 
    unique_viruses = np.unique(associations['virus']) # get a list of unique viruses in the dataset
    host_order_dict = {} # create a dictionary. Will host host:order -  so I can reference what order each host is in
    unique_hosts = list(unique_hosts)
   # create a dictionary. Will host host:order -  so I can reference what order each host is in
    for each_host in unique_hosts: # go through each host
        each_host_data = host_data.loc[host_data['host_name'] == each_host] # find the data on each host in the host data 
        host_order = each_host_data['order'] # get the order that the host is in 
        host_order = host_order.item() # just the order (get rid of the row index)
        host_order_dict[each_host] = host_order # add the host and order to the dictionary as host: order
    
    return host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict 



def create_dict_of_orders_and_num_appropriate_viruses(host_order_dict, associations, is_virus_zoonotic_dict):
    orders= np.unique(host_order_dict.values())
    order_virus_dict = {}
    for each_order in orders:
        order_virus_dict[each_order] = []
    order_virus_dict['HUMANS'] = []

    num_associations = len(associations['host'])
    for each in range(0, num_associations):
        virus = associations['virus'][each]
        host = associations['host'][each]
       
        is_virus_zoonotic = is_virus_zoonotic_dict[virus]
        if is_virus_zoonotic == 0:
            
            if host == 'Homo_sapiens':
                host_order = 'HUMANS'
                
            else:
                host_order = host_order_dict[host]
            current_order_virus_list = order_virus_dict[host_order]
            
            new_order_virus_list = current_order_virus_list
            new_order_virus_list.append(virus)
        
            order_virus_dict[host_order] = new_order_virus_list
   
    return order_virus_dict



def add_edges_to_host_association_net(host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict, order_virus_dict):
    unique_orders = list(np.unique(host_order_dict.values()))
    unique_orders.append('HUMANS')
   
    
    order_combinations = list(itertools.combinations(unique_orders, 2))
    for each_combination in order_combinations:
        order1 = each_combination[0]
        order2 = each_combination[1]
        all_order1_viruses = order_virus_dict[order1]
        all_order2_viruses = order_virus_dict[order2]
        unique_order1_viruses = list(np.unique(all_order1_viruses))
        unique_order2_viruses = list(np.unique(all_order2_viruses))
        
        num_order1_viruses = len(unique_order1_viruses)
        num_order2_viruses = len(unique_order2_viruses)

        shared_viruses = [each_order1_virus for each_order1_virus in unique_order1_viruses if each_order1_virus in unique_order2_viruses]
        num_shared_viruses = len(shared_viruses)
        
        #jaccards_index
        if num_shared_viruses > 0:
            # C/(N1+N2 - C)
            jaccards_index = float(num_shared_viruses)/float((num_order1_viruses+num_order2_viruses-num_shared_viruses))
            all_host_association_net.add_edge(order1, order2, weight = jaccards_index)
     
    #remove self-loops        
    all_host_association_net.remove_edges_from(all_host_association_net.selfloop_edges())
    
    #remove nodes with no edges
    for each_node in nx.nodes(all_host_association_net):
        each_node_neighbors = all_host_association_net.neighbors(each_node)
        if len(each_node_neighbors) == 0:
            all_host_association_net.remove_node(each_node)  
    
    return all_host_association_net



def get_num_samples_for_node_size(associations, host_order_dict, all_host_association_net, is_virus_zoonotic_dict):
    #print is_virus_zoonotic_dict
    num_associations = len(associations['host'])
    
    order_num_samples_dict = {}
    for each_order in np.unique(host_order_dict.values()):
        order_num_samples_dict[each_order] = 0
    order_num_samples_dict['HUMANS'] = 0
  
    for each in range(0, num_associations):
        virus = associations['virus'][each]
        is_virus_zoonotic = is_virus_zoonotic_dict[virus]
        if is_virus_zoonotic == 0:
            host = associations['host'][each]
            if host == 'Homo_sapiens':
                previous_num_human_samples = order_num_samples_dict['HUMANS']
                current_num_human_samples = previous_num_human_samples+1
                order_num_samples_dict['HUMANS'] = current_num_human_samples
                
            else:
                get_host_order = host_order_dict[host]
                previous_num_order_samples = order_num_samples_dict[get_host_order]
                current_num_order_samples = previous_num_order_samples+1
                order_num_samples_dict[get_host_order] = current_num_order_samples
    
    
    for each_order in order_num_samples_dict.keys():
        if not each_order in nx.nodes(all_host_association_net):
            del order_num_samples_dict[each_order]
   
    
    nx.set_node_attributes(all_host_association_net, 'num_samples', order_num_samples_dict)
    return all_host_association_net
          

if __name__ == '__main__':

    all_host_association_net, olival_associations, olival_hosts, olival_viruses = read_in_files_and_create_nets()
    host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict  = host_and_order_data(olival_hosts, olival_associations, olival_viruses, all_host_association_net)
    order_virus_dict = create_dict_of_orders_and_num_appropriate_viruses(host_order_dict, associations, is_virus_zoonotic_dict)
    non_zoonotic_host_association_net = add_edges_to_host_association_net(host_order_dict, unique_viruses, all_host_association_net, associations, is_virus_zoonotic_dict, order_virus_dict)
    non_zoonotic_host_association_net = get_num_samples_for_node_size(associations, host_order_dict, non_zoonotic_host_association_net, is_virus_zoonotic_dict)

    node_strength_dict = {}
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
            node_strength_dict[each_node] = average_node_strength
            
    
    num_samples_dict = nx.get_node_attributes(non_zoonotic_host_association_net, 'num_samples')
    node_dict = {}
    ids = []
    labels=  []
    num_samples = []
    for each_node in nx.nodes(non_zoonotic_host_association_net):
       ids.append(each_node)
       labels.append(each_node)
       each_num_samples = num_samples_dict[each_node]
       num_samples.append(each_num_samples)
    node_dict = {'Id':ids, 'Label': labels, 'Num_samples':num_samples}
    node_dict_df = pd.DataFrame(node_dict)

