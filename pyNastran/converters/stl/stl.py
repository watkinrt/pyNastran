#pylint:  disable=C0111
from __future__ import print_function
from six import iteritems
from six.moves import zip, range
from copy import deepcopy
#from struct import pack

from numpy import array, zeros, ndarray, cross, where, vstack, unique
from numpy import arctan2, cos, sin, transpose, pi
from numpy.linalg import norm
import scipy

from struct import unpack, Struct, pack

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils import is_binary_file
from pyNastran.utils.log import get_logger


class STL(object):
    #model_type = 'stl'
    #isStructured = False
    #isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')

        self.nodes = None
        self.elements = None
        self.header = None
        self.infilename = None

    def write_stl(self, out_filename, is_binary=False, float_fmt='%6.12f'):
        self.log.info("---writing STL file...%r---" % out_filename)
        assert len(self.nodes) > 0
        solid_name = 'dummy_name'
        if is_binary:
            self.write_binary_stl(out_filename)
        else:
            self.write_stl_ascii(out_filename, solid_name, float_fmt=float_fmt)

    def read_stl(self, stl_filename):
        self.infilename = stl_filename
        self.log.info("---starting reading STL file...%r---" % self.infilename)

        if is_binary_file(stl_filename):
            self.read_binary_stl(stl_filename)
        else:
            self.read_ascii_stl(stl_filename)

        #self.log.info("nPoints=%s  nElements=%s" % (self.nPoints, self.nElements))
        self.log.info("---finished reading STL file...%r---" % self.infilename)
        #assert self.nPoints > 0, 'nPoints=%s' % self.nPoints
        #assert self.nElements > 0, 'nElements=%s' % self.nElements


    def write_binary_stl(self, stl_filename):
        """Write an STL binary file."""
        infile = open(stl_filename, "wb")

        if hasattr(self, 'header'):
            self.header.ljust(80, '\0')
            infile.write(self.header)
        else:
            header = '%-80s' % stl_filename
            infile.write(pack('80s', header))
        #avector = [0., 0., 0.]
        #bvector = [0., 0., 0.]
        #cvector = [0., 0., 0.]
        nelements = self.elements.shape[0]
        infile.write(pack('i', nelements))
        elements = self.elements

        p1 = self.nodes[elements[:, 0], :]
        p2 = self.nodes[elements[:, 1], :]
        p3 = self.nodes[elements[:, 2], :]
        avector = p2 - p1
        bvector = p3 - p1
        n = cross(avector, bvector)
        del avector, bvector
        #n /= norm(n, axis=1)

        s = Struct('12fH')
        for eid, element in enumerate(elements):
            data = s.pack(n[eid, 0], n[eid, 1], n[eid, 2],
                          p1[eid, 0], p1[eid, 1], p1[eid, 2],
                          p2[eid, 0], p2[eid, 1], p2[eid, 2],
                          p3[eid, 0], p3[eid, 1], p3[eid, 2], 0)
            infile.write(data)
        infile.close()

    def read_binary_stl(self, stl_filename):
        """Read an STL binary file."""
        with open(stl_filename, 'rb') as infile:
            data = infile.read()

        self.header = data[:80]
        print('header = %r' % self.header.rstrip())
        nelements, = unpack('i', data[80:84])
        j = 84

        inode = 0
        nodes_dict = {}
        elements = zeros((nelements, 3), 'int32')

        s = Struct('12fH')
        for ielement in range(nelements):
            (nx, ny, nz, ax, ay, az, bx, by, bz,
             cx, cy, cz, i) = s.unpack(data[j:j+50])

            t1 = (ax, ay, az)
            t2 = (bx, by, bz)
            t3 = (cx, cy, cz)
            if t1 in nodes_dict:
                i1 = nodes_dict[t1]
            else:
                i1 = inode
                nodes_dict[t1] = inode
                inode += 1

            if t2 in nodes_dict:
                i2 = nodes_dict[t2]
            else:
                i2 = inode
                nodes_dict[t2] = inode
                inode += 1

            if t3 in nodes_dict:
                i3 = nodes_dict[t3]
            else:
                i3 = inode
                nodes_dict[t3] = inode
                inode += 1
            elements[ielement] = [i1, i2, i3]
            j += 50
        assert inode > 0, inode
        nnodes = inode + 1 # accounting for indexing
        self.elements = elements
        nodes = zeros((nnodes, 3), 'float64')

        for node, inode in iteritems(nodes_dict):
            nodes[inode] = node
        self.nodes = nodes


    def _get_normals_data(self, elements, nodes=None):
        """
        This is intended as a submethod to help handle the problem of bad normals
        """
        if nodes is None:
            nodes = self.nodes
        self.log.debug("get_normals...elements.shape %s" % str(elements.shape))
        p1 = nodes[elements[:, 0]]
        p2 = nodes[elements[:, 1]]
        p3 = nodes[elements[:, 2]]
        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)
        normals_norm = norm(v123, axis=1)
        inan = where(normals_norm == 0)[0]
        return v123, normals_norm, inan

    def remove_elements_with_bad_normals(self, elements, nodes=None):
        v123, normals_norm, inan = self._get_normals_data(elements, nodes=nodes)
        if len(inan):
            inotnan = where(normals_norm != 0)[0]
            self.elements = elements[inotnan, :]
            normals_norm = normals_norm[inotnan]
            v123 = v123[inotnan]
            self.log.info('removing %i elements with coincident nodes' % len(inan))

        normals = v123
        normals[:, 0] /= normals_norm
        normals[:, 1] /= normals_norm
        normals[:, 2] /= normals_norm
        return normals

    def get_area(self, elements, stop_on_failure=True):
        v123, normals_norm, inan = self._get_normals_data(elements, nodes=self.nodes)

        if stop_on_failure:
            msg = 'Failed Elements: %s\n' % inan
            if len(inan) > 0:
                for inani in inan:
                    msg += '  eid=%s nodes=%s\n' % (inani, elements[inani, :])
                    for ni in elements[inani]:
                        msg += '    nid=%s node=%s\n' % (ni, self.nodes[ni, :])
                raise RuntimeError(msg)
        return 0.5 * normals_norm

    def get_normals(self, elements, nodes=None, stop_on_failure=True):
        """
        Parameters
        ----------
        elements : (n, 3) ndarray ints
            the elements to get the normals for
        nodes : (n, ) ndarray; default=None -> all
            a subset of the nodes
        stop_on_failure : bool (default=True)
            True:  crash if there are coincident points
            False: delete elements
        """
        if nodes is None:
            nodes = self.nodes
        v123, normals_norm, inan = self._get_normals_data(elements, nodes=nodes)

        if stop_on_failure:
            msg = 'Failed Elements: %s\n' % inan
            if len(inan) > 0:

                for ifail, inani in enumerate(inan):
                    msg += '  eid=%s nodes=%s\n' % (inani, elements[inani, :])
                    for ni in elements[inani]:
                        msg += '    nid=%s node=%s\n' % (ni, nodes[ni, :])
                    if ifail > 10:
                        break
                msg += 'Failed Elements: %s; n=%s\n' % (inan, len(inan))
                raise RuntimeError(msg)
        else:
            inotnan = where(normals_norm != 0)[0]
            normals_norm[inan] = array([1., 0., 0.])
            elements = elements[inotnan, :]
            normals_norm = normals_norm[inotnan]
            v123 = v123[inotnan]

        # we need to divide our (n,3) array in 3 steps
        normals = v123 # / normals_norm
        normals[:, 0] /= normals_norm
        normals[:, 1] /= normals_norm
        normals[:, 2] /= normals_norm
        return normals

    def flip_normals(self, i=None):
        """
        Flips the normals of the specified elements.

        Parameters
        ----------
        i : (n, ) ndarray ints; default=None -> all
            the indicies to flip
        """
        self.log.info("---flip_normals---")
        if i is None:
            elements = self.elements
        else:
            elements = self.elements[i, :]

        n0, n1, n2 = elements[:, 0], elements[:, 1], elements[:, 2]
        elements2 = elements.copy()
        elements2[:, 0] = n0
        elements2[:, 1] = n2
        elements2[:, 2] = n1
        if i is None:
            self.elements = elements2
        else:
            self.elements[i, :] = elements2 #[i, :]

    def get_normals_at_nodes(self, elements, normals=None, nid_to_eid=None):
        """
        Calculates the normal vector of the nodes based on the average
        element normal.

        Parameters
        ----------
        elements : ????
            The elements...should be removed
        normals : (n, 3) ndarray floats
            The elemental normals
        nid_to_eid : Dict[int] = [int, int, ... ]
            key = node_id
            value = list of element_ids

        Returns
        -------
        normals_at_nodes : (nnodes, 3) ndarray ints
            the normals
        """
        if normals is None:
            nodes = self.nodes
            normals = self.get_normals(elements, nodes=self.nodes)

        if nid_to_eid is None:
            from collections import defaultdict
            nid_to_eid = defaultdict(list)
            eid = 0
            for (n1, n2, n3) in elements:
                nid_to_eid[n1].append(eid)
                nid_to_eid[n2].append(eid)
                nid_to_eid[n3].append(eid)
                eid += 1
            del eid, n1, n2, n3

        normals_at_nodes = zeros(nodes.shape, 'float64')
        eid = 0
        for nid, elementsi in iteritems(nid_to_eid):
            pe = normals[elementsi]
            m = pe.mean(axis=0)
            normals_at_nodes[nid] = m/norm(m)
            eid += 1
        return normals_at_nodes

    def equivalence_nodes(self, tol=1e-5):
        nnodes = self.nodes.shape[0]

        # build the kdtree
        kdt = scipy.spatial.KDTree(self.nodes)

        # find the node ids of interest
        nids_new = unique(self.elements.ravel())
        nids_new.sort()

        # check the closest 10 nodes for equality
        deq, ieq = kdt.query(self.nodes[nids_new, :], k=10, distance_upper_bound=tol)

        # get the ids of the duplicate nodes
        slots = where(ieq[:, 1:] < nnodes)
        replacer = unique(ieq[slots])

        # update the duplcated node id with it's partner
        # we'll pick the minimum ID
        for r in replacer:
            ip = where(ieq[r, :] < nnodes)[0]
            possible = ieq[r, ip]

            # node 11 can become node 10, but node 10 cannot become node 11
            ip2 = where(r > possible)[0]

            if len(ip2):
                # replace the node ids
                possible2 = possible[ip2]
                r_new_nid = possible2.min()
                ireplace = where(self.elements == r)
                self.elements[ireplace] = r_new_nid

    def project_boundary_layer(self, nodes, elements, volume_bdfname):
        """
        Create a boundary layer mesh.

        Parameters
        ----------
        nodes : (n, 3) ndarray floats
            The nodes on the surface.
        elements : (n, 3) ndarray ints
            The elements on the surface.
        volume_bdfname : str
            The CPENTA bdf file to write.

        Returns
        -------
        nodes2 : (n, 3) ndarray floats
            The boundary layer nodes
        elements2 : (n, 6) ndarray ints
            The boundary layer elements
        """
        self.log.info("project_mesh...")

        normals_at_nodes = self.get_normals_at_nodes(elements, normals=None, nid_to_eid=None)

        #print "normals_at_nodes[4]", normals_at_nodes[4]
        #----------- make boundary layer---------------
        # deltaN = a^N * delta
        delta = 0.1
        b = 1.0
        a = 1.1
        N = 13
        r = array(range(N))
        #r = 1100

        deltaNs = b * a**r * delta
        self.log.info('deltaNs = %s' % deltaNs)
        if not isinstance(deltaNs, ndarray):
            deltaNs = array([deltaNs])
        N = len(deltaNs)
        #print "N = ", N

        nid = 0
        #print "deltaNs =", deltaNs
        nnodes = nodes.shape[0]
        nodes2 = zeros((nnodes * (N+1), 3), 'float64')
        nodes2[:nnodes, :] = nodes

        nelements = elements.shape[0]
        elements2 = zeros((nelements * (N+1), 3), 'int32')
        elements2[:nelements, :] = elements
        #print "nodes.shape =", nodes.shape
        #print "nodes2.shape =", nodes2.shape
        #print "elements2.shape =", elements2.shape

        ni = 0
        #print "nelements =", nelements
        cid = None
        nid = 1
        eid2 = 1
        pid = 100
        mid = 100
        with open(volume_bdfname, 'wb') as bdf:
            bdf.write('CEND\nBEGIN BULK\n')
            bdf.write('$NODES in Layer=0\n')
            for (x, y, z) in nodes:
                card = ['GRID', nid, cid, x, y, z]
                bdf.write(print_card_8(card))
                nid += 1

            for deltaN in deltaNs:
                outer_points = nodes + normals_at_nodes * deltaN
                nodes2[ni*nnodes : (ni+1)*nnodes, :] = outer_points

                nnbase = ni * nnodes
                nnshift = (ni+1) * nnodes

                nebase = (ni) * nelements
                neshift = (ni + 1) * nelements
                elements2[neshift : neshift + nelements, :] = elements + nnodes * (ni + 1)

                self.log.info('nodes = %s' % str(nodes))
                self.log.info('deltaN = %s' % str(deltaN))
                self.log.info('normals_at_nodes = %s' % str(normals_at_nodes))
                self.log.info('outer_points = %s' % str(outer_points))
                bdf.write('$NODES in Layer=%i\n' % (ni + 1))
                for x, y, z in outer_points:
                    card = ['GRID', nid, cid, x, y, z]
                    bdf.write(print_card_8(card))
                    nid += 1

                bdf.write('$SOLID ELEMENTS in Layer=%i\n' % (ni + 1))
                for eid in range(nelements):
                    (n1, n2, n3) = elements2[nebase  + eid] + 1
                    (n4, n5, n6) = elements2[neshift + eid] + 1
                    card = ['CPENTA', eid2, pid, n1, n2, n3, n4, n5, n6]
                    bdf.write(print_card_8(card))
                    eid2 += 1

                card = ['PSOLID', pid, mid]
                bdf.write(print_card_8(card))

                E = 1e7
                G = None
                nu = 0.3
                card = ['MAT1', mid, E, G, nu]
                bdf.write(print_card_8(card))

                pid += 1
                mid += 1
                ni += 1

        #print(elements2)
        #for node in nodes:
            #normal = normals_elements[nid]
            #nid += 1
        #print(deltaN)

        #----------- make far field---------------
        nodes3 = nodes2[nnbase:, :]
        nodes3 = nodes2[nnbase:, :]

        elements3 = elements2[nebase:, :]
        elements3 = elements2[nebase:, :]
        self.log.debug("done projecting...")
        return nodes2, elements2


    def write_stl_ascii(self, out_filename, solid_name, float_fmt='%.6f'):
        """
        solid solid_name
         facet normal -6.601157e-001 6.730213e-001 3.336009e-001
          outer loop
            vertex 8.232952e-002 2.722531e-001 1.190414e+001
            vertex 8.279775e-002 2.717848e-001 1.190598e+001
            vertex 8.557653e-002 2.745033e-001 1.190598e+001
          endloop
         endfacet
        end solid
        """
        nodes = self.nodes
        elements = self.elements
        self.log.info("---write_stl_ascii...%r---" % out_filename)
        msg = ''
        node_format = ' facet normal %s %s %s\n' % (float_fmt, float_fmt, float_fmt)
        vertex_format = '     vertex %s %s %s\n' % (float_fmt, float_fmt, float_fmt)
        msg += 'solid %s\n' % solid_name
        out = open(out_filename, 'wb')
        out.write(msg)

        nelements = elements.shape[0]
        normals = self.get_normals(elements)
        for element, normal in zip(elements, normals):
            try:
                p1, p2, p3 = nodes[element]
            except IndexError:
                print(element)
                raise

            #msg  += '  facet normal -6.601157e-001 6.730213e-001 3.336009e-001\n'
            msg = node_format % tuple(normal)
            msg += '   outer loop\n'
            msg += vertex_format % tuple(p1)
            msg += vertex_format % tuple(p2)
            msg += vertex_format % tuple(p3)
            msg += '   endloop\n'
            msg += ' endfacet\n'
            #print(msg)
            out.write(msg)
        msg = 'endsolid\n'
        out.write(msg)
        out.close()


    def read_ascii_stl(self, stl_filename):
        infile = open(stl_filename, 'r')
        line = infile.readline()
        inode = 0
        ielement = 0
        nodes_dict = {}
        elements = []
        while line:
            if 'solid' in line[:6].lower():
                line = infile.readline()  # facet
                while 'facet' in line.strip()[:5].lower():
                    #facet normal -6.665299e-001 6.795624e-001 3.064844e-001
                    #   outer loop
                    #      vertex 8.142845e-002 2.731541e-001 1.190024e+001
                    #      vertex 8.186898e-002 2.727136e-001 1.190215e+001
                    #      vertex 8.467505e-002 2.754588e-001 1.190215e+001
                    #   endloop
                    #endfacet

                    infile.readline() # outer loop
                    L1 = infile.readline().strip()
                    L2 = infile.readline().strip()
                    L3 = infile.readline().strip()

                    v1 = L1.split()[1:]
                    v2 = L2.split()[1:]
                    v3 = L3.split()[1:]
                    infile.readline() # endloop
                    infile.readline() # endfacet
                    t1 = tuple(v1)
                    t2 = tuple(v2)
                    t3 = tuple(v3)

                    assert len(v1) == 3, '%r' % L1
                    assert len(v2) == 3, '%r' % L2
                    assert len(v3) == 3, '%r' % L3

                    if t1 in nodes_dict:
                        i1 = nodes_dict[t1]
                    else:
                        i1 = inode
                        nodes_dict[t1] = inode
                        inode += 1

                    if t2 in nodes_dict:
                        i2 = nodes_dict[t2]
                    else:
                        i2 = inode
                        nodes_dict[t2] = inode
                        inode += 1

                    if t3 in nodes_dict:
                        i3 = nodes_dict[t3]
                    else:
                        i3 = inode
                        nodes_dict[t3] = inode
                        inode += 1
                    element = [i1, i2, i3]
                    elements.append(element)
                    ielement += 1
                    line = infile.readline()  # facet
                #print "end of solid..."
            elif 'endsolid' in line.lower():
                line = infile.readline()
            elif line.strip() == '':
                line = infile.readline()
            else:
                self.log.error(line)
                #line = f.readline()
                raise NotImplementedError('multiple solids are not supported; line=%r' % line)
                #break
        infile.close()

        assert inode > 0, inode
        nnodes = inode + 1 # accounting for indexing
        self.elements = array(elements, 'int32')
        nodes = zeros((nnodes, 3), 'float64')

        for node, inode in iteritems(nodes_dict):
            nodes[inode] = node
        self.nodes = nodes

    def scale_nodes(self, xscale, yscale, zscale):
        self.nodes[:, 0] *= xscale
        self.nodes[:, 1] *= yscale
        self.nodes[:, 2] *= zscale

    def shift_nodes(self, xshift, yshift, zshift):
        self.nodes[:, 0] += xshift
        self.nodes[:, 1] += yshift
        self.nodes[:, 2] += zshift

    def flip_axes(self, axes, scale):
        if axes == 'xy':
            x = deepcopy(self.nodes[:, 0])
            y = deepcopy(self.nodes[:, 1])
            self.nodes[:, 0] = y * scale
            self.nodes[:, 1] = x * scale
        elif axes == 'yz':
            y = deepcopy(self.nodes[:, 1])
            z = deepcopy(self.nodes[:, 2])
            self.nodes[:, 1] = z * scale
            self.nodes[:, 2] = y * scale
        elif axes == 'xz':
            x = deepcopy(self.nodes[:, 0])
            z = deepcopy(self.nodes[:, 2])
            self.nodes[:, 0] = z * scale
            self.nodes[:, 2] = x * scale

    def create_mirror_model(self, xyz, tol):
        """
        Creates a mirror model.

        Parameters
        ----------
        xyz : str {x, y, z}
            the direction of symmetry
        tol: float
            the tolerance for symmetry plane nodes

        .. note:: All elements on the symmetry plane will be removed
        """
        assert xyz in ['x', 'y', 'z'], 'xyz=%r' % xyz
        assert tol >= 0.0, 'tol=%s' % tol

        nnodes = self.nodes.shape[0]
        if xyz == 'x':
            xyzi = 0
        elif xyz == 'y':
            xyzi = 1
        elif xyz == 'z':
            xyzi = 2
        else:
            raise RuntimeError(xyz)

        # the nodes on the symmetry plane
        i = where(self.nodes[:, xyzi] < tol)[0]

        # smash the symmetry nodes to 0.0
        self.nodes[i, xyzi] = 0.
        nodes_sym = deepcopy(self.nodes)
        nodes_sym[:, xyzi] *= -1.

        # we're lazy and duplicating all the nodes
        # but will only write out a subset of them
        nodes = vstack([self.nodes, nodes_sym])

        # create the symmetrical elements
        elements2 = []
        elements3 = []
        for element in self.elements:
            epoints = nodes[element, xyzi][0]
            je = where(epoints <= tol)[0]
            if len(je) < 3:  # not a symmetry element, so we save it
                elements2.append(element)

                # duplicate the node if it's not on the symmetry plane
                element3 = [elementi if elementi in i else (elementi + nnodes)
                            for elementi in element]

                # the normal is now backwards, so we flip it
                element3.reverse()
                elements3.append(element3)

        self.nodes = nodes
        self.elements = array(elements2 + elements3, dtype='int32')


def _rotate_model(stl):  # pragma: no cover
    nodes = stl.nodes
    elements = stl.elements
    if 0:
        # rotate the model
        x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
        #i = where(y > 0.0)[0]
        R = x**2 + y**2
        theta = arctan2(y, x)
        iRz = where(R == 0)[0]
        theta[iRz] = 0.0

        min_theta = min(theta)
        dtheta = max(theta) - pi/4
        theta2 = theta + min_theta

        x2 = R * cos(theta2)
        y2 = R * sin(theta2)
        #print "x.shape", x.shape, y2.shape
        nodes_rotated = transpose(vstack([x2, y2, z]))
        #print "nodes.shape", nodes_rotated.shape
        #print nodes_rotated

    if 0:
        # project the volume
        (nodes2, elements2) = stl.project_mesh(nodes_rotated, elements)

    # write the model
    stl_geom_out = 'rotated.stl'
    stl.write_stl_ascii(stl_geom_out, 'sphere')


if __name__ == '__main__':  # pragma: no cover
    from pyNastran.converters.stl.stl_reshape import main
    main()

