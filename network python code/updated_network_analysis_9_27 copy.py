#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:20:18 2018

@author: admin
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 13:36:31 2018

@author: Casey

Take Olival data and create unipartite host association networks.
Can create three types of networks- just wild hosts, just domestic hosts, or all hosts
The nodes in the networks are the orders of the hosts i.e.'CARNIVORA' 'CETARTIODACTYLA' 'CHIROPTERA' 'CINGULATA' 'DIDELPHIMORPHIA' 'DIPROTODONTIA' 'EULIPOTYPHLA' 'LAGOMORPHA' 'PERAMELEMORPHIA' 'PERISSODACTYLA' 'PILOSA' 'PRIMATES' 'PROBOSCIDEA' 'RODENTIA' 'SCANDENTIA'
The association data has virus, host, and wild or domestic
I go through each association in the Olival association data. For each virus, I get the orders of all of the hosts that the virus is associated with. I create an edge between all orders that have the same virus. 
If an edge already exists between two orders, I increase the weight of the edge by one for each repeated association between the same two orders
I do this jsut for wild hosts, just for domestic hosts, and for all hosts.
Result:
    Network where nodes are orders, edges are viral associations (a virus is found in both of the orders that are the nodes), the weight of the edge indicates the number of viruses that are shared between two orders



"""
import pandas as pd
import networkx as nx
import numpy as np
import itertools





def read_in_files_and_create_nets():
    """
    create empty networks for wild, dom, or all hosts
    read in olival data
    input: none
    output: 3 networks and olival data (hosts, associations, viruses)
    """
    #create empty nets for 3 types of networks
    wild_host_association_net = nx.Graph()
    dom_host_association_net = nx.Graph()
    all_host_association_net = nx.Graph()    
    #read in olival data
    olival_associations = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/olival nature 2017/associations.csv')
    olival_hosts = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/olival nature 2017/hosts.csv')
    olival_viruses = pd.read_csv('/Users/admin/Desktop/GitHub/brevity/olival nature 2017/viruses.csv')
    return wild_host_association_net, dom_host_association_net, all_host_association_net, olival_associations, olival_hosts, olival_viruses





def host_and_order_data(all_host_data, all_association_data, wild_host_association_net, dom_host_association_net, all_host_association_net):
    """
    deal with the host and order data and add nodes to the network for each order
    input: host data, association data, 3 empty networks
    output: dictionary of hosts and their orders, list of unique viruses, 3 networks with nodes but no edges, association data
    """
    
    host_data = all_host_data[['hHostNameFinal','hOrder']] # just get the host name and the order from the olival host data
    host_data.columns = ['host_name', 'order'] # change column names to be easier
    
    unique_hosts = np.unique(host_data['host_name']) # get a list of unique hosts in the host data
    
    unique_orders = np.unique(host_data['order']) # get a list of the unique orders in the host data
    #print unique_orders
    for each_order in unique_orders: # loop thorugh each order and create a node in each network for every order (I am assuming all orders will be present in all 3 networks currently- may need to change)
        wild_host_association_net.add_node(each_order)
        dom_host_association_net.add_node(each_order)
        all_host_association_net.add_node(each_order)
    
    associations = all_association_data[['vVirusNameCorrected','hHostNameFinal', 'WildDomInReference']] # get the virus name, the host name, and whether the host sampled was wild or domestic
    associations.columns = ['virus', 'host', 'wild_or_dom']  # change the column names to be easier 
    #print associations[associations['host']=='Homo_sapiens']
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
   
    return host_order_dict, unique_viruses, wild_host_association_net, dom_host_association_net, all_host_association_net, associations






def add_edges_to_host_association_net(wild_dom_or_all, host_order_dict, unique_viruses, wild_host_association_net, dom_host_association_net, all_host_association_net, associations):
#    if wild_dom_or_all == 'wild': # create wild host network
#        for each_virus in unique_viruses: # go through each virus in the data
#            wild_virus_data = associations.loc[(associations['virus'] ==each_virus) & (associations['wild_or_dom']=='wild')] # get the association data for each virus where the host is wild
#            host_combinations = list(itertools.combinations(list(wild_virus_data['host']), 2)) # create all possible 2 host combinations of hosts selected above   
#            for each_combination in host_combinations: # go through each combination
#                order_combination = [] # an empty list to stick the orders of the hosts in 
#                for each_host in each_combination: # go through each host in the list of 2 hosts that are associated
#                    each_host_order = host_order_dict[each_host] # get the order of each host from the host- order dictionary
#                    order_combination.append(each_host_order)  # put the 2 orders in the order combination list so I can make an edge between them     
#                if wild_host_association_net.has_edge(order_combination[0], order_combination[1]) == False: # if there is no edge betweent the orders
#                    wild_host_association_net.add_edge(order_combination[0], order_combination[1], weight = 1) # create an edge with a weight of 1 between the orders
#                if wild_host_association_net.has_edge(order_combination[0], order_combination[1]) == True: # if there is already an edge between the 2 orders
#                    edge_weight = wild_host_association_net[order_combination[0]][order_combination[1]]['weight'] # get the current weight of that edge
#                    updated_edge_weight = edge_weight+1 # add 1 to the weight of the edge and save as updated weight
#                    wild_host_association_net.add_edge(order_combination[0], order_combination[1], weight = updated_edge_weight) # redraw the edge with the new weight
#    
#    if wild_dom_or_all == 'dom': # create domestic host network
#        for each_virus in unique_viruses:# go through each virus in the data
#            dom_virus_data = associations.loc[(associations['virus'] ==each_virus) & (associations['wild_or_dom']=='domestic')]# get the association data for each virus where the host is domestic
#            host_combinations = list(itertools.combinations(list(dom_virus_data['host']), 2))  # create all possible 2 host combinations of hosts selected above   
#            for each_combination in host_combinations:# go through each combination
#                order_combination = []# an empty list to stick the orders of the hosts in
#                for each_host in each_combination:# go through each host in the list of 2 hosts that are associated
#                    each_host_order = host_order_dict[each_host]# get the order of each host from the host- order dictionary
#                    order_combination.append(each_host_order)       # put the 2 orders in the order combination list so I can make an edge between them     
#                if dom_host_association_net.has_edge(order_combination[0], order_combination[1]) == False:# if there is no edge betweent the orders
#                    dom_host_association_net.add_edge(order_combination[0], order_combination[1], weight = 1)# create an edge with a weight of 1 between the orders
#                if dom_host_association_net.has_edge(order_combination[0], order_combination[1]) == True:# if there is already an edge between the 2 orders
#                    edge_weight = dom_host_association_net[order_combination[0]][order_combination[1]]['weight'] # get the current weight of that edge
#                    updated_edge_weight = edge_weight+1# add 1 to the weight of the edge and save as updated weight
#                    dom_host_association_net.add_edge(order_combination[0], order_combination[1], weight = updated_edge_weight)# redraw the edge with the new weight
    num_viruses_with_more_than_one_host = 0
    num_viruses_with_more_than_one_order = 0
    if wild_dom_or_all == 'all': # create all host network
        for each_virus in unique_viruses:# go through each virus in the data
            #num_viruses = num_viruses+1
            
            all_virus_data = associations.loc[associations['virus'] ==each_virus]# get the association data for each virus for all hosts
            #print all_virus_data['host']
            #print len(all_virus_data['host'])
            if len(all_virus_data['host'])>1:
                    num_viruses_with_more_than_one_host = num_viruses_with_more_than_one_host+1
            host_order_list = []
            for each_host in all_virus_data['host']:
                host_order_list.append(host_order_dict[each_host])
            #print host_order_list
            if len(np.unique(host_order_list))>1:
                num_viruses_with_more_than_one_order = num_viruses_with_more_than_one_order+1
            
            
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
    #print num_viruses_with_more_than_one_host
    #print num_viruses_with_more_than_one_order
    # just want to return the network that I have specified 
    if wild_dom_or_all == 'wild':
        desired_network = wild_host_association_net
    if wild_dom_or_all == 'dom':
        desired_network = dom_host_association_net
    if wild_dom_or_all == 'all':
        desired_network = all_host_association_net
    for each_node in nx.nodes(desired_network):
        if len(desired_network.neighbors(each_node)) == 0:
            desired_network.remove_node(each_node)
    
    return desired_network# just return the desired network
            
            
        
        
            

if __name__ == '__main__':

    wild_host_association_net, dom_host_association_net, all_host_association_net, olival_associations, olival_hosts, olival_viruses = read_in_files_and_create_nets()
    host_order_dict, unique_viruses, wild_host_association_net, dom_host_association_net, all_host_association_net, associations = host_and_order_data(olival_hosts, olival_associations, wild_host_association_net, dom_host_association_net, all_host_association_net)
    what_type_of_net = 'all' # specify 'wild', 'dom', or 'all'
    desired_net = add_edges_to_host_association_net(what_type_of_net, host_order_dict, unique_viruses, wild_host_association_net, dom_host_association_net, all_host_association_net, associations)
    
#    for each_node in nx.nodes(desired_net):
#        each_node_neighbors = desired_net.neighbors(each_node)
#        print each_node
#        print each_node_neighbors
    print nx.number_of_edges(desired_net)
    desired_net.remove_edges_from( desired_net.selfloop_edges())
    print nx.number_of_edges(desired_net)
    
    for each_node in nx.nodes(desired_net):
        node_degree = desired_net.degree(each_node)
        print each_node
        #print node_degree
        neighbor_weights = []
        node_neighbors = desired_net.neighbors(each_node)
        for each_node_neighbor in node_neighbors:
            edge_weight = desired_net[each_node][each_node_neighbor]['weight']
            
            neighbor_weights.append(edge_weight)
        node_strength = sum(neighbor_weights)
        if node_degree >0:
            average_node_strength = float(node_strength)/float(node_degree)
            print average_node_strength

#    print nx.number_of_edges(desired_net)
#    edge_weight_matrix = nx.adjacency_matrix(desired_net, weight = 'weight')
#    filename = "/Users/admin/Dropbox (Bansal Lab)/brevity_project/reanalysis/edge_weight_matrix_all_net.csv"
#    #for each_edge in nx.edges(desired_net):
#        #print each_edge['weight']
#    print nx.get_edge_attributes(desired_net, 'weight')
#        #print each
#        
#    f= open(filename,"a+")
#    #f.write(str(nx.get_edge_attributes(desired_net, 'weight')))
#    #f.close()
#    for each_edge, each_weight in nx.get_edge_attributes(desired_net, 'weight').iteritems():
#        f.write(str(each_edge)+"\n"+str(each_weight)+"\n")
#        #print each_weight
#    f.close()
#    for each_line in edge_weight_matrix:
#        f.write(str(each_line)+"\n")
#    f.close()
        
    #edge_weight_matrix.to_csv('/Users/admin/Dropbox (Bansal Lab)/brevity_project/reanalysis/edge_weight_matrix_all_net.csv', index = False)
    #nx.write_edgelist(desired_net, '/Users/admin/Dropbox (Bansal Lab)/brevity_project/reanalysis/olival_net_edgelist_no_humans_9_27.csv', data = ['weight'], delimiter = ',')


#    for each_node in nx.nodes(desired_net):
#        #node_neighbors = desired_net.neighbors(each_node)
#        #num_node_neighbors = len(node_neighbors)
##        print each_node
##        print node_neighbors
##        print num_node_neighbors
#        print each_node
    print nx.eigenvector_centrality(desired_net, 100)



