import csv
import json
import numpy as np
import networkx as nx
from numba import jit

# JIT cannot handle yet with string/char arrays, that's why
#    I have to use int arrays. The states are:
S = 0
E = 1
I = 2
R = 3
G = 10


class Agent_Country:
    def __init__(self, args, graph):
        self.iterations = 0

        self.args = args
        self.graph = graph
        self.N = len(graph.nodes)
        self.pos = nx.get_node_attributes(self.graph, "pos")

        if("distribution" in args):
            self.distribution = args["distribution"]
            self.num_of_groups = len(self.distribution)
            self.distribution_cumulative = [0] + [sum(self.distribution[:i+1]) for i in range(self.num_of_groups)]
            self.type_dict = {agent: i for i in range(self.num_of_groups) for agent in range(self.distribution_cumulative[i], self.distribution_cumulative[i + 1])}
        else:
            print("Please give the distribution")

        if("probas" in args):
            self.probas = args["probas"]
        else:
            print("Please give the distribution")

        # === Init states ===
        self.init_states()
        
        # === Infections ===
        if("infected_agents" in args):
            self.init_seeds = args["infected_agents"]
            self.seed[self.init_seeds] = self.init_seeds

            self.states[args["infected_agents"]] = I
            
            nodes=list(self.graph.nodes())
            nx.set_node_attributes(self.graph, {nodes[ind]: -1 for ind in args["infected_agents"]}, 'infector')          

            nx.set_node_attributes(self.graph, {nodes[ind]: 0 for ind in args["infected_agents"]}, 'first_inf_time')          

            self.timers[self.init_seeds] = np.random.geometric(self.args["gamma"])+1
            
            indexes = nx.get_node_attributes(graph, "index")
        else:
            print("Please give the infected_agents")

    def barabasi_inf(self, graph):
        ind = np.random.randint(0,self.N)
        inf_nodes = np.array([node for node in graph.neighbors(ind)])[:10]
        self.states[inf_nodes] = np.random.choice(np.array([I,G]), inf_nodes.size)
        self.timers[inf_nodes] = np.random.geometric(self.args["gamma"])+1



    def grid_inf(self):
        # === Infect some agents in the center ===
        # IGI
        # GIG
        # IGI
        n = int(np.sqrt(self.N))
        c = n//2
        inf_nodes = np.array([n*(c-1)+c-1, n*(c-1)+c, n*(c-1)+c+1, n*c+c-1, n*c+c, n*c+c+1, n*(c+1)+c-1, n*(c+1)+c, n*(c+1)+c+1])
        self.states[inf_nodes] = np.array([G,I,G,I,G,I,G,I,G])
        #self.states[inf_nodes] = np.array([1,10,1,10,1,10,1,10,1])
        self.timers[inf_nodes] = np.random.geometric(self.args["gamma"])+1


    def run(self):
        self.log_json()
        for i in range(self.args["max_iteration"]):
            if self.check_stop():
                print("Infection extinction at iteration: {}".format(i))
                break

            self.step()
            self.log_json()
        
        print("")
    
    def check_stop(self):
        return (np.sum(self.states == I)+ np.sum(self.states == G) == 0)

    def check_stop2(self):
        return not ( (self.seed[self.states == I]==self.init_seeds[0]).any() and (self.seed[self.states == I]==self.init_seeds[1]).any() )
    
    def step(self):
        #print("\rStep {}".format(self.iterations), end = '')

        beta = self.args["beta"]
        if "beta2" in self.args:
            beta2 = self.args["beta2"]
        else:
            beta2 = 0.1

        gamma = self.args["gamma"]
        xi = self.args["xi"]
        # === Infections ===
        all_infected = Agent_Country.update_neighs(self.states, self.indexes, self.timers, self.Rtimers,
                              self.neighs_groups, self.infector, self.seed, beta, gamma, xi, self.num_of_groups,
                              self.type_dict, self.distribution_cumulative,  self.probas)
        
        nodes=list(self.graph.nodes())
        nx.set_node_attributes(self.graph, {nodes[ind]: int(self.infector[ind]) for ind in all_infected}, 'infector')          
        nx.set_node_attributes(self.graph, {nodes[ind]: self.iterations+1 for ind in all_infected}, 'first_inf_time')          

        #for ind in all_infected:
        #    print(self.iterations,ind,self.graph.nodes[nodes[ind]])

        self.iterations += 1

    #@jit(nopython = True)
    def update_neighs(states, indexes, timers, Rtimers, neighs, infector, seed,
                      beta, gamma, xi, num_of_groups, type_dict, distribution_cumulative, probas):
        # 0 : Suscepted
        # 1 : Infected
        neigh_arr = neighs[0]
        slices = neighs[1]


        if xi<=0:
            states[timers==1]=R
        elif xi>=1:
            states[timers==1]=S
        else:
            states[timers==1]=R
            Rtimers[timers==1] =np.random.geometric(xi, size=np.sum(timers==1))
            states[Rtimers==1]=S
            Rtimers[Rtimers>0]-=1
            
        seed[timers==1] = -1
        timers[timers>0]-=1

        
        #if (maxInf>0) and (maxInf<int(np.sum(states == I)) + int(np.sum(states == G)):
        #    beta=beta/10
        #    print("Max reached")

        all_infected=[np.int64(x) for x in range(0)]
        infected_dict = {agent: [] for agent in indexes}
        for ind in indexes[states == I]:
            infected_indexes = []
            for infected_type in range(num_of_groups):
                adj_list = neigh_arr[slices[ind][infected_type][0]:slices[ind][infected_type][1]]
                S_adj_list = adj_list[states[adj_list]==0]
                num_infected = np.random.binomial(S_adj_list.size, probas[type_dict[ind], infected_type])
                infected_indexes = np.random.choice(S_adj_list, num_infected, replace = False)
                for infected_index in infected_indexes:
                    infected_dict[infected_index].append(ind)

        infected = []
        for infected_index in infected_dict:
            if len(infected_dict[infected_index]) >= 1:
                infected += infected_index,
                if len(infected_dict[infected_index]) > 1:
                    weights = np.array([probas[type_dict[ind], type_dict[infected_index]] for ind in infected_dict[infected_index]])
                    ind = np.random.choice(infected_dict[infected_index], p = weights/sum(weights))
                    infected_dict[infected_index] = [ind]

                states[infected_index] = -I
                timers[infected_index] = np.random.geometric(gamma, size=1)+1
                seed[infected_index] = seed[infected_dict[infected_index][0]]
                infector[infected_index] = infected_dict[infected_index][0]

        
        all_infected += infected


        for i,s in enumerate(states):
            if s<0:
                states[i]=abs(s)
                
        return all_infected

    # === Helper Functions ===
    def get_neigh_flattened(self, neighs):
        arr = []
        slices = np.zeros(shape = (len(neighs),2))
        ind = 0
        for i,adj_list in enumerate(neighs):
            n = len(adj_list)
            slices[i] = (ind, ind+n)
            arr += adj_list

            ind += n

        return (np.array(arr, dtype = np.int32), np.array(slices, dtype = np.int32))

    def get_neigh_flattened_groups(self, neighs, distribution, type_dict):
        arr = []
        slices = np.zeros(shape = (len(neighs),len(distribution),2), dtype=np.int32)
        ind = 0
        for i,adj_list in enumerate(neighs):
            type_i = type_dict[i]
            for j, group_size in enumerate(distribution):
                n = group_size
                if j == type_i:
                    n = group_size - 1
                slices[i][j] = (ind, ind + n)
                ind += n

            arr += adj_list

        return (np.array(arr, dtype = np.int32), slices)

    def init_states(self):
        #self.states = np.array(["S"]*self.N)
        self.states = np.zeros(self.N, dtype = np.int8)
        self.indexes = np.array(np.arange(self.N, dtype= np.int32))
        self.timers = np.zeros(self.N, dtype = np.int32)
        self.Rtimers = np.zeros(self.N, dtype = np.int32)
        self.seed = -np.ones(self.N, dtype=np.int32)
        
        all_neighs = []
        indexes = nx.get_node_attributes(self.graph, "index")
        
        med_deg = np.median([val for (node, val) in self.graph.degree()])
        for name in self.graph.nodes():
            neighs = [indexes[neigh] for neigh in self.graph[name].keys() if "second_type" not in self.graph[name][neigh]]
            all_neighs.append(neighs)

        self.neighs = self.get_neigh_flattened(all_neighs)
        if("distribution" in self.args):
            self.neighs_groups = self.get_neigh_flattened_groups(all_neighs, self.distribution, self.type_dict)
        self.infector = -np.ones(self.N, dtype=np.int32)

        
    def log_json(self):
        # === Delete previous logs ===
        if self.iterations == 0:
            with open(self.args["logfile"], 'w') as f:
                pass
        
        with open(self.args["logfile"], 'a') as outfile:
            outfile.write(json.dumps(self.get_json())+"\n")

    def get_seed_info(self, index):
        if(index not in self.init_seeds):
            return False
        else:     
            data = {}
            data["I"] = int(np.sum((self.seed == index) & (self.states == I)))
            return data


    def encode_dists(dists, normalize=False):
        if isinstance(dists, float):
            return (0,{0:0})
        if len(dists)==0:
            return (0,{0:0})
        hist, edges = np.histogram(dists, range= (0.0, dists.max()), bins=int(np.ceil(dists.max()))+1)
        if normalize:
            return (np.mean(dists),{str(edges[i]):float(hist[i])/(normalize*max(1,4*np.ceil(edges[i]))) for i in range(len(hist))})
        else:               
            return (np.mean(dists),{str(edges[i]):int(hist[i]) for i in range(len(hist))})


        
    def get_agent_info(self, index, name):
        data = {}
        data["index"] = index
        data["name"] = name
        data["population"] = 1 # This is the original
        data["d_SIR"] = {
            "S": int(self.states[index] == S),
            "I": int(self.states[index] == I),
            "R": int(self.states[index] == R),
        }
        return data

    def get_country_info(self):
        data = {"d_SIR": {} }
        data["d_SIR"]["S"] = int(np.sum(self.states == S))
        data["d_SIR"]["I"] = int(np.sum(self.states == I))
        data["d_SIR"]["R"] = int(np.sum(self.states == R))
        

        return data

