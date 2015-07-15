class BDFAttributes(object):

    def __init__(self):
        pass

    @property
    def nnodes(self):
        return len(self.nodes)

    @property
    def node_ids(self):
        return self.nodes.keys()

    #def get_nodes(self):
        #nodes = []
        #for nid, node in sorted(iteritems(self.nodes)):
            #nodes.append(node)
        #return nodes

    #--------------------
    # Elements CARDS

    @property
    def nelements(self):
        return len(self.elements)

    @property
    def element_ids(self):
        return self.elements.keys()

    #--------------------
    # Property CARDS

    @property
    def nproperties(self):
        return len(self.properties)

    @property
    def property_ids(self):
        return self.properties.keys()

    #--------------------
    # Material CARDS

    @property
    def material_ids(self):
        return self.materials.keys()

    @property
    def nmaterials(self):
        return len(self.materials)

    #--------------------
    # Coords CARDS

    @property
    def coord_ids(self):
        return self.coords.keys()

    @property
    def ncoords(self):
        return len(self.coords)

    #--------------------

    @property
    def ncaeros(self):
        return len(self.caeros)

    @property
    def caero_ids(self):
        return self.caeros.keys()
