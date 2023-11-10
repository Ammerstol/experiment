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
                              self.neighs, self.neighs2, self.infector, self.seed, beta, beta2, gamma, xi)
        
        nodes=list(self.graph.nodes())
        nx.set_node_attributes(self.graph, {nodes[ind]: int(self.infector[ind]) for ind in all_infected}, 'infector')          
        nx.set_node_attributes(self.graph, {nodes[ind]: self.iterations+1 for ind in all_infected}, 'first_inf_time')          

        #for ind in all_infected:
        #    print(self.iterations,ind,self.graph.nodes[nodes[ind]])

        self.iterations += 1

    #@jit(nopython = True)
    def update_neighs(states, indexes, timers, Rtimers, neighs, neighs2, infector, seed,
                      beta, beta2, gamma, xi):
        # 0 : Suscepted
        # 1 : Infected
        neigh_arr = neighs[0]
        slices = neighs[1]
        neigh_arr2 = neighs2[0]
        slices2 = neighs2[1]


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

        for ind in indexes[states == I]:   
            adj_list = neigh_arr[slices[ind][0]:slices[ind][1]]
            S_adj_list = adj_list[states[adj_list]==0]

            adj_list2 = neigh_arr2[slices2[ind][0]:slices2[ind][1]]
            S_adj_list2 = adj_list2[states[adj_list2]==0]

            # === Choose infected: ===
            new_inf_num = np.random.binomial(S_adj_list.size,min(1,beta))
            infected = np.random.choice(S_adj_list, replace = False, size=new_inf_num)
            all_infected+=list(infected)


            states[infected] = -I
            timers[infected] = np.random.geometric(gamma, size=len(infected))+1
            seed[infected] = seed[ind]
            infector[infected] = ind
        
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

    def init_states(self):
        #self.states = np.array(["S"]*self.N)
        self.states = np.zeros(self.N, dtype = np.int8)
        self.indexes = np.array(np.arange(self.N, dtype= np.int32))
        self.timers = np.zeros(self.N, dtype = np.int32)
        self.Rtimers = np.zeros(self.N, dtype = np.int32)
        self.seed = -np.ones(self.N, dtype=np.int32)
        
        all_neighs = []
        all_neighs2 = []
        indexes = nx.get_node_attributes(self.graph, "index")
        
        med_deg = np.median([val for (node, val) in self.graph.degree()])
        for name in self.graph.nodes():
            neighs = [indexes[neigh] for neigh in self.graph[name].keys() if "second type" not in self.graph[name][neigh]]
            all_neighs.append(neighs)
            neighs2 = [indexes[neigh] for neigh in self.graph[name].keys() if "second type" in self.graph[name][neigh]]
            all_neighs2.append(neighs2)

        self.neighs = self.get_neigh_flattened(all_neighs)
        self.neighs2 = self.get_neigh_flattened(all_neighs2)
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

