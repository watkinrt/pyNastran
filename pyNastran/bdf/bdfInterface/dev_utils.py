from __future__ import print_function
from six import iteritems
from six.moves import zip
import scipy
from pyNastran.bdf.bdf import BDF
from numpy import array, zeros, unique, where, arange, hstack, vstack, searchsorted

def bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                          neq_max=10):
    """
    Equivalences nodes; keeps the lower node id; creates two nodes with the same

    .. warning:: only handles CQUAD4, CTRIA3
    .. warning:: assumes cid=0
    .. warning:: renumbers nodes
    """
    model = BDF()
    model.read_bdf(bdf_filename, xref=True)

    # quads / tris
    quads = []
    quadmap = []
    tris = []
    trimap = []

    inode = 1
    nid_map = {}
    renumber_nodes = True
    if renumber_nodes:
        for nid, node in sorted(iteritems(model.nodes)):
            node.nid = inode
            nid_map[inode - 1] = nid
            inode += 1
    else:
	    raise NotImplementedError()
    #model.write_bdf('A_' + bdf_filename_out)

    nids = array([node.nid for nid, node in sorted(iteritems(model.nodes))], dtype='int32')
    nnodes = len(nids)
    i = arange(nnodes, dtype='int32')
    #nids2 = vstack([i, nids]).T

    nodes_xyz = array([node.xyz for nid, node in sorted(iteritems(model.nodes))], dtype='float32')

    nids_new = set([])
    for eid, element in sorted(iteritems(model.elements)):
        emap = []

        if element.type == 'CQUAD4':
            nids_quads.append(element.node_ids)
            eids_quads.append(element.eid)
        elif element.type == 'CTRIA3':
            nids_tris.append(element.node_ids)
            eids_tris.append(element.eid)
        else:
            raise NotImplementedError(element.type)

    nids_quads = array(quads, dtype='int32') - 1
    #eids_quads = array(quadmap, dtype='int32')
    nids_tris = array(tris, dtype='int32') - 1
    #eids_tris = array(trimap, dtype='int32')

    # build the kdtree
    try:
        kdt = scipy.spatial.cKDTree(nodes_xyz)
    except RuntimeError:
        print(nodes_xyz)
        raise RuntimeError(nodes_xyz)

    # find the node ids of interest
    nids_new = unique(hstack([
        nids_quads.flatten(), nids_tris.flatten()
    ]))
    nids_new.sort()
    inew = searchsorted(nids, nids_new, side='left')

    # check the closest 10 nodes for equality
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)

    # get the ids of the duplicate nodes
    slots = where(ieq[:, 1:] < nnodes)
    irows, icols = slots
    replacer = unique(ieq[slots])

    skip_nodes = []
    for (irow, icol) in zip(irows, icols):
        nid1 = nid_map[irow]
        nid2 = nid_map[icol]

        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        node2.nid = node1.nid
        node2.xyz = node1.xyz
        assert node2.cp == node1.cp
        assert node2.cd == node1.cd
        assert node2.ps == node1.ps
        assert node2.seid == node1.seid
        skip_nodes.append(nid2)

    #model.remove_nodes = skip_nodes
    #model._write_nodes = _write_nodes
    model.write_bdf(bdf_filename_out)


def _write_nodes(self, outfile, size, is_double):
    """
    Writes the NODE-type cards

    :param self: the BDF object
    """
    if self.spoints:
        msg = []
        msg.append('$SPOINTS\n')
        msg.append(self.spoints.write_card(size, is_double))
        outfile.write(''.join(msg))

    if self.nodes:
        msg = []
        msg.append('$NODES\n')
        if self.gridSet:
            msg.append(self.gridSet.print_card(size))
        for (nid, node) in sorted(iteritems(self.nodes)):
            if nid not in self.remove_nodes:
                msg.append(node.write_card(size, is_double))
        outfile.write(''.join(msg))


def eq2():
    lines = [
        '$pyNastran: version=msc',
        '$pyNastran: punch=True',
        '$pyNastran: encoding=ascii',
        '$NODES',
        '$ Nodes to merge:',
        '$ 5987 10478',
        '$   GRID        5987           35.46     -6.      0.',
        '$   GRID       10478           35.46     -6.      0.',
        '$ 5971 10479',
        '$   GRID        5971           34.92     -6.      0.',
        '$   GRID       10479           34.92     -6.      0.',
        '$ 6003 10477',
        '$   GRID        6003             36.     -6.      0.',
        '$   GRID       10477             36.     -6.      0.',
        'GRID        5971           34.92     -6.      0.',
        'GRID        5972           34.92-5.73333      0.',
        'GRID        5973           34.92-5.46667      0.',
        'GRID        5987           35.46     -6.      0.',
        'GRID        5988           35.46-5.73333      0.',
        'GRID        5989           35.46-5.46667      0.',
        'GRID        6003             36.     -6.      0.',
        'GRID        6004             36.-5.73333      0.',
        'GRID        6005             36.-5.46667      0.',
        'GRID       10476             36.     -6.    -1.5',
        'GRID       10477             36.     -6.      0.',
        'GRID       10478           35.46     -6.      0.',
        'GRID       10479           34.92     -6.      0.',
        'GRID       10561           34.92     -6.    -.54',
        '$ELEMENTS_WITH_PROPERTIES',
        'PSHELL         1       1      .1',
        'CQUAD4      5471       1    5971    5987    5988    5972',
        'CQUAD4      5472       1    5972    5988    5989    5973',
        'CQUAD4      5486       1    5987    6003    6004    5988',
        'CQUAD4      5487       1    5988    6004    6005    5989',
        'PSHELL        11       1      .1',
        'CTRIA3      9429      11   10561   10476   10478',
        'CTRIA3      9439      11   10478   10479   10561',
        'CTRIA3      9466      11   10476   10477   10478',
        '$MATERIALS',
        'MAT1           1      3.              .3',
    ]
    bdf_filename = 'nonunique2.bdf'
    f = open(bdf_filename, 'wb')
    f.write('\n'.join(lines))
    f.close()
    bdf_equivalence_nodes('nonunique2.bdf', 'unique2.bdf', 0.01)

def eq1():
    msg = 'CEND\n'
    msg += 'BEGIN BULK\n'
    msg += 'GRID,1,,0.,0.,0.\n'
    msg += 'GRID,2,,0.,0.,0.5\n'
    msg += 'GRID,3,,0.,0.,0.51\n'
    msg += 'GRID,10,,0.,0.,1.\n'
    msg += 'GRID,11,,0.,0.,1.\n'
    msg += 'CTRIA3,1,1,1,2,11\n'
    msg += 'CTRIA3,2,1,1,2,11\n'
    msg += 'CTRIA3,3,1,2,3,11\n'
    msg += 'CTRIA3,4,1,1,2,10\n'
    msg += 'PSHELL,1,1,0.1\n'
    msg += 'MAT1,1,3.0,, 0.3\n'
    msg += 'ENDDATA'

    bdf_filename = 'nonunique.bdf'

    f = open(bdf_filename, 'wb')
    f.write(msg)
    f.close()

    bdf_filename_out = 'unique.bdf'
    tol = 0.2
    bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol)

if __name__ == '__main__':
    eq1()
    eq2()
