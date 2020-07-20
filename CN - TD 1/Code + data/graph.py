""" CR15 graph library """
class Graph(object):

    """def __init__(self, graph_dict={}):
        " initializes a graph object "
        self.__graph_dict = graph_dict"""

    def __init__(self, graph_dict=None):
        """ initializes a graph object """
        if graph_dict:
            self.__graph_dict = graph_dict
        else:
            self.__graph_dict = dict()


    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()





    def __generate_edges(self):
        edges = set() # I use a set structure for easy removing of doublons.
        
        for v1 in self.__graph_dict :
            for v2 in self.__graph_dict[v1] :
                
                ## I check for errors in graph definition
                if v1 == v2: # First error case : loop
                    print("error : a loop is present in the graph")
                    return()
                elif not (v1 in self.__graph_dict[v2]) : # Second error case : asymetry
                    print("error : an edge is asymetrical")
                    return()
                
                ## if no error exist I add the edge to the output
                else:
                    edges.add(frozenset([v1,v2])) # edges are converted to frozenset because we need to arrange them in a set and set() do not allow dynamic structures such as sets

        ### Then I convert output to a list(set) for readability
        edges_list = []
        for e in edges:
            edges_list.append(set(e))
        return edges_list

    def add_vertex(self, vertex):
        if vertex in self.__graph_dict:
            print("The vertex " + vertex + " already was in the graph !")
        else :
            self.__graph_dict[vertex]=[]

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list. No loops or 
        multiple edges. To complete."""

        ## Early conversion to a set assures that the function works the same on different data type
        edge = set(edge) 

        ## first I check edge size
        if (len(edge) > 2) :
            print("Error : edge has more than two elements !")
            return()

        elif (len(edge) == 1):
            print("Error : edge is a loop !")
            return()

        elif (len(edge) < 1):
            print("Error : edge is empty !")
            return()
        
        ## Then I extract vertices
        v1 = edge.pop()
        v2 = edge.pop()

        
        ## Then I check that the edge is properly defined
        if set([v1,v2]) in self.edges():
            print("The edge " + str([v1,v2]) + " already existed")
            return()

        elif not ((v1 in self.vertices()) and (v2 in self.vertices())):
            print("error : undefined vertex !")
            return()

        ## if no error can be detected, I add the edge as two symmetrical entries in the __graph_dict structure"
        self.__graph_dict[v1].append(v2)
        self.__graph_dict[v2].append(v1)

        return()

    
    def vertex_degree(self):
        degree = dict()
        for v in self.__graph_dict.keys():
            degree[v] = len(self.__graph_dict[v])

        return degree


    def find_isolated_vertices(self):
        degree = self.vertex_degree()
        isolated_vertices=set([v for v in degree if degree[v] == 0])

        return isolated_vertices


    def density(self):
        n = len(self.vertices())
        m = len(self.edges())

        # to avoid limit problem, I do not allow density calculation when the number of vertices do not allow for any edge
        if n < 2 :
            print("No density can be computed for graph with n<2")
            return()

        return 2*m/(n*(n-1))
    

    def degree_sequence(self):
        degree = self.vertex_degree()
        deg_seq = list(degree.values())
        deg_seq.sort(reverse=True)
        deg_seq = tuple(deg_seq)

        return(deg_seq)


    def erdos_gallai(self,seq):
        n = len(seq)

        ## check non-increase condition
        for i in range(1,n):
            if seq[i]>seq[i-1]:
                print("Sequance isn't non increasing")
                return(False)

        ## check sum parity
        if sum(seq)%2 == 1:
            print("Sequence has odd sum")
            return(False)

        ## check inequality for all k
        for k in range(1,n):
            l_term = sum(seq[0:k])

            r_term = k*(k-1)
            r_term += sum(min(list(seq[k+1:]),[k for i in range(k+1,n)]))

            if l_term > r_term:
                print("Inequality invalid for k = " + str(k))
                return(False)
            
        ## and if no condition is invalid
        return(True)


    def global_clustering_coefficient(self):
        ## Let's first build the set of all triplets in the graph
        global_triplets = set()

        for v in self.vertices():
            local_triplets = set()
            for i in range(len(self.__graph_dict[v])):
                for j in range(i+1,len(self.__graph_dict[v])):
                    s = frozenset([v,self.__graph_dict[v][i],self.__graph_dict[v][j]])
                    local_triplets.add(s)
                    
            for t in local_triplets :
                    global_triplets.add(t)

        ## Then let's count the number of triplets that are triangles
        t_counter = 0
        for s in global_triplets:
            [v1,v2,v3] = list(s)
            if (v1 in self.__graph_dict[v2]) and (v2 in self.__graph_dict[v3]) and (v3 in self.__graph_dict[v1]): #here we (ab)use the symmetry property to gain a few operations
                t_counter+=1

        ## ...and return the result. By convention, we will take 0 for all graph where no triplet exist.
        if len(global_triplets) > 0 :
            return(t_counter/len(global_triplets))
        else:
            return(0)

            
    def connected_component(self):
        
        ## Initialisation
        visited = dict() # this store both the set of visited vertices and the order in which they were visited
        count = 0 # this is the variable counting visit order
        connected_components_list = []

        for v in self.vertices():
            visited[v] = count # 0 is understood as "not visited"

            
        ## Calculation
        while (0 in visited.values()):
            
            ## Initialisation of a single connected component's exploration
            
            adj_vertices = set() # contains the list of vertices that are adjacent to the visited ensemble
            connected_component = set() # contains the set of elements of the connected component currently active

            v = [key  for (key, value) in visited.items() if value == 0][0] # pick an unvisited vertex
            adj_vertices.add(v)

            ## explore until end of connected component
            while len(adj_vertices) > 0:
                v = adj_vertices.pop()
                connected_component.add(v)
                l = [v_adj for v_adj in self.__graph_dict[v] if visited[v_adj] == 0]
                adj_vertices |= set(l)
                                
                count += 1
                visited[v] = count

            ## add connected component to connected component list
            connected_components_list.append(frozenset(connected_component)) # frozenset is used to allow connected components list to be converted as a set, which is useful to check equality
             
        return(connected_components_list)

    def shortest_path(self):
        distances = dict()

        ## Distances initialisation
        for v1 in self.vertices():
            distances[v1] = dict()
            for v2 in self.vertices():
                distances[v1][v2] = float("inf")

        ## definition of origin node for the BFS
        for v_ini in self.vertices(): 
            d = 0
            adj_vertices = set() ## vertices that are adjacent to the currently explored subset
            adj_vertices.add(v_ini)

            ## explore with BFS
            while len(adj_vertices) > 0: 
                d_class = [v for v in adj_vertices] ## distance class
                
                ## vertices closest to origine are explored first
                for v in d_class :
                    distances[v_ini][v] = d
                    adj_vertices.remove(v)
                    
                ## then we add their child to the adjacent set
                for v1 in d_class :
                    adj_vertices |= set([v2 for v2 in self.__graph_dict[v1] if distances[v_ini][v2] > d])
                d+=1

        return(distances)


    def diameter(self):
        d_max = 0
        distances = self.shortest_path()
        for v1 in self.vertices():
            for v2 in self.vertices():
                d_max = max(distances[v1][v2],d_max)        
        return(d_max)
                
    def biggest_component_diameter(self):
        d_max = 0
        distances = self.shortest_path()
        for v1 in self.vertices():
            for v2 in self.vertices():
                if float("inf") > distances[v1][v2]:
                    d_max = max(distances[v1][v2],d_max)        
        return(d_max)


    def spanning_forest(self):
        ## Initialisation of the spanning forest as a new graph
        spanning_forest = Graph(dict())
        
        
        while [v for v in self.vertices() if v not in spanning_forest.vertices()]:
            
            ## Initialisation of a new connected component
            v_ini = [v for v in self.vertices() if v not in spanning_forest.vertices()][0] # Current spanning tree's central vertex
            spanning_forest.add_vertex(v_ini)
            v_frontier = [v_ini] # list of current spanning trees elements whose neighbours have not been explored yet

            ## Each frontier node is used to explore its neighbours. Then it is removed from frontier nodes.
            while v_frontier:
                v1 = v_frontier[0]
                v1_childs = [v for v in self.__graph_dict[v1] if v not in spanning_forest.vertices()]
                v_frontier.remove(v1)

                ## Each explored node is added to the spanning tree, as well as its parent edge
                for v2 in v1_childs :
                    spanning_forest.add_vertex(v2)
                    spanning_forest.add_edge([v1,v2])
                    v_frontier.append(v2)

        return(spanning_forest)
        
