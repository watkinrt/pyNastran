from __future__ import print_function

from six import iteritems
from six.moves import range
import os

from numpy import vstack, amax, amin, arange, ones, zeros, where, unique

from pyNastran.converters.ugrid.surf_reader import TagReader
from pyNastran.converters.ugrid.ugrid_reader import UGRID
from pyNastran.converters.ugrid.ugrid2d_reader import UGRID2D_Reader
from pyNastran.utils import print_bad_path, is_binary_file

#VTK_TRIANGLE = 5

import vtk
from vtk import vtkTriangle, vtkQuad


class UGRID_IO(object):
    def __init__(self):
        pass

    def get_ugrid_wildcard_geometry_results_functions(self):
        data = (
            'AFLR3 Ugrid',
            'AFLR3 Ugrid (*.ugrid)', self.load_ugrid_geometry,
            None, None)
        return data

    def load_ugrid_geometry(self, ugrid_filename, dirname, name='main', plot=True):
        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skip_reading:
        #    return
        if is_binary_file(ugrid_filename):
            model = UGRID(log=self.log, debug=True)
            base, fmt, ext = os.path.basename(ugrid_filename).split('.')
            is_2d = False
        else:
            base, ext = os.path.basename(ugrid_filename).split('.')
            model = UGRID2D_Reader(log=self.log, debug=True)
            is_2d = True

        self.model_type = 'ugrid'
        print('ugrid_filename = %s' % ugrid_filename)


        assert ext == 'ugrid', ugrid_filename
        model.read_ugrid(ugrid_filename)

        if is_2d:
            tris = model.tris
            quads = model.quads
        else:
            tris = model.tris - 1
            quads = model.quads - 1

        #self.nodes = nodes
        #self.tris  = tris
        #self.quads = quads
        #self.pids = pids

        #self.tets = tets
        #self.penta5s = penta5s
        #self.penta6s = penta6s
        #self.hexas = hexas

        nnodes = model.nodes.shape[0]
        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        nelements = ntris + nquads

        nodes = model.nodes
        self.nElements = nelements
        self.nNodes = nnodes

        print("nNodes = %s" % self.nNodes)
        print("nElements = %s" % self.nElements)
        assert nelements > 0, nelements

        self.grid.Allocate(self.nElements, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)
        self.log.info('max = %s' % mmax)
        self.log.info('min = %s' % mmin)

        diff_node_ids = model.check_hanging_nodes(stop_on_diff=False)
        if len(diff_node_ids):
            red = (1., 0., 0.)
            self.create_alternate_vtk_grid('hanging_nodes', color=red, line_width=5, opacity=1., point_size=10, representation='point')
            self._add_ugrid_nodes_to_grid('hanging_nodes', diff_node_ids, nodes)
            self._add_alt_actors(self.alt_grids)


        for inode, node in enumerate(nodes):
            points.InsertPoint(inode, node)

        if ntris:
            for eid, element in enumerate(tris):
                elem = vtkTriangle()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
        if nquads:
            for eid, element in enumerate(quads):
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                elem.GetPointIds().SetId(3, element[3])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())

        self.nElements = nelements
        self.grid.SetPoints(points)
        self.grid.Modified()
        print('update...')
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        #print("updated grid")

        # loadCart3dResults - regions/loads
        self. turn_text_on()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['AFLR UGRID Surface', '']}
        cases = {}
        ID = 1

        if hasattr(model, 'pids'):
            form, cases = self._fill_ugrid3d_case(ugrid_filename, cases, ID, nnodes, nelements, model)
        else:
            form, cases = self._fill_ugrid2d_case(ugrid_filename, cases, ID, nnodes, nelements, model)

        if plot:
            self._finish_results_io2(form, cases)

    def _add_ugrid_nodes_to_grid(self, name, diff_node_ids, nodes):
        """
        based on:
          _add_nastran_nodes_to_grid
        """
        nnodes = nodes.shape[0]
        assert nnodes > 0, nnodes
        # if nnodes == 0:
            # return
        nnodes = len(diff_node_ids)
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        for nid in diff_node_ids:
            node = nodes[nid, :]
            print('nid=%s node=%s' % (nid, node))
            points.InsertPoint(nid, *node)

            if 1:
                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, nid)
            else:
                elem = vtk.vtkSphere()
                sphere_size = self._get_sphere_size(dim_max)
                elem.SetRadius(sphere_size)
                elem.SetCenter(points.GetPoint(nid))

            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        self.alt_grids[name].SetPoints(points)

    def clear_surf(self):
        pass

    # def _load_ugrid_results(self, openfoam_filename, dirname):
        # pass

    def _fill_ugrid2d_case(self, base, cases, ID, nnodes, nelements, model):
        cases_new = []
        results_form = []

        geometry_form = [
            ('ElementID', 0, []),
            ('NodeID', 1, []),
        ]

        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        cases[(ID, 0, 'ElementID', 1, 'centroid', '%i', '')] = eids
        cases[(ID, 1, 'NodeID', 1, 'node', '%i', '')] = nids

        form = [
            ('Geometry', None, geometry_form),
        ]
        #if has_tag_data:
            #form.append(('Tag Data', None, tag_form),)

        #results_form = []
        #if len(results_form):
            #form.append(('Results', None, results_form))
        return form, cases

    def _fill_ugrid3d_case(self, base, cases, ID, nnodes, nelements, model):
        tag_filename = base + '.tags'
        mapbc_filename = base.split('.')[0] + '.mapbc'
        print('mapbc_filename = %r' % mapbc_filename)

        cases_new = []
        has_tag_data = False
        has_mapbc_data = False
        results_form = []
        mapbc_form = []

        geometry_form = [
            #('Region', 0, []),
            ('ElementID', 0, []),
            ('NodeID', 1, []),
            ('SurfaceID', 2, []),
            #('normSpacing', 3, []),
            #('BL_thick', 4, []),
            #('ReconFlag', 5, []),
            #('GridBC', 6, []),
        ]

        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        #nelements = ntris + nquads
        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        #grid_bcs = element_props[:, 2]

        #npids = len(model.pids)
        pids = model.pids
        cases[(ID, 0, 'ElementID', 1, 'centroid', '%i', '')] = eids
        cases[(ID, 1, 'NodeID', 1, 'node', '%i', '')] = nids
        cases[(ID, 2, 'SurfaceID', 1, 'centroid', '%i', '')] = pids

        n = 3
        if os.path.exists(tag_filename):
            #surf_ids = element_props[:, 0]
            #recon_flags = element_props[:, 1]
            #cases[(ID, 2, 'ReconFlag', 1, 'centroid', '%i')] = recon_flags
            #cases[(ID, 3, 'GridBC',    1, 'centroid', '%i')] = grid_bcs

            tagger = TagReader()
            data = tagger.read_tag_filename(tag_filename)

            int_data = ones((nelements, 8), dtype='int32') * -10.
            float_data = zeros((nelements, 2), dtype='float64')
            for key, datai in sorted(iteritems(data)):
                #self.log.info(datai)
                [name, is_visc, is_recon, is_rebuild, is_fixed, is_source, is_trans, is_delete, bl_spacing, bl_thickness, nlayers] = datai
                i = where(pids == key)[0]
                int_data[i, :] = [is_visc, is_recon, is_rebuild, is_fixed, is_source, is_trans, is_delete, nlayers]
                float_data[i, :] = [bl_spacing, bl_thickness]
                self.log.info('data[%i] = %s' % (key, name))

            has_tag_data = True
            tag_form = []
            tag_form.append( ('is_visc',      n, []) )
            tag_form.append( ('is_recon',     n+1, []) )
            tag_form.append( ('is_rebuild',   n+2, []) )
            tag_form.append( ('is_fixed',     n+3, []) )
            tag_form.append( ('is_source',    n+4, []) )
            tag_form.append( ('is_trans',     n+5, []) )
            tag_form.append( ('is_delete',    n+6, []) )
            tag_form.append( ('nlayers',      n+7, []) )
            tag_form.append( ('bl_spacing',   n+8, []) )
            tag_form.append( ('bl_thickness', n+9, []) )

            cases[(ID, n, 'is_visc',      1, 'centroid', '%i', '')] = int_data[:, 0]
            cases[(ID, n + 1, 'is_recon',   1, 'centroid', '%i', '')] = int_data[:, 1]
            cases[(ID, n + 2, 'is_rebuild', 1, 'centroid', '%i', '')] = int_data[:, 2]
            cases[(ID, n + 3, 'is_fixed',   1, 'centroid', '%i', '')] = int_data[:, 3]
            cases[(ID, n + 4, 'is_source',  1, 'centroid', '%i', '')] = int_data[:, 4]
            cases[(ID, n + 5, 'is_trans',   1, 'centroid', '%i', '')] = int_data[:, 5]
            cases[(ID, n + 6, 'is_delete',  1, 'centroid', '%i', '')] = int_data[:, 6]
            cases[(ID, n + 7, 'nlayers',    1, 'centroid', '%i', '')] = int_data[:, 7]

            cases[(ID, n + 8, 'bl_spacing',   1, 'centroid', '%.3e', '')] = float_data[:, 0]
            cases[(ID, n + 9, 'bl_thickness', 1, 'centroid', '%.3e', '')] = float_data[:, 1]
            n += 10
        else:
            self.log_info('tag_filename=%r could not be found' % tag_filename)

        if os.path.exists(mapbc_filename):
            has_mapbc_data = True
            mapbc = open(mapbc_filename, 'r')
            lines = mapbc.readlines()
            lines = [line.strip() for line in lines
                     if not line.strip().startswith('#') and line.strip()]
            npatches = int(lines[0])
            mapbcs = zeros(pids.shape, dtype='int32')
            for ipatch in range(npatches):
                line = lines[ipatch + 1]
                iline, bc_num, name = line.split()
                iline = int(iline)
                bc_num = int(bc_num)
                assert ipatch + 1 == iline, 'line=%r; ipatch=%s iline=%s' % (line, ipatch + 1, iline)
                islot = where(pids == ipatch + 1)[0]
                if len(islot) == 0:
                    upids = unique(pids)
                    msg = 'ipatch=%s not found in pids=%s' % (ipatch + 1, upids)
                    raise RuntimeError(msg)
                mapbcs[islot] = bc_num
                print(line)
            mapbc_form.append(('Map BC', n, []))
            cases[(ID, n, 'Map BC', 1, 'centroid', '%i', '')] = mapbcs
        else:
            self.log_info('mapbc_filename=%r could not be found' % mapbc_filename)


        #norm_spacing = model.node_props[:, 0]
        #bl_thickness = model.node_props[:, 1]
        #cases[(ID, 1, 'normSpacing', 1, 'node', '%.3e', '')] = norm_spacing
        #cases[(ID, 2, 'BL_thick',    1, 'node', '%.3e', '')] = bl_thickness

        form = [
            ('Geometry', None, geometry_form),
        ]
        if has_tag_data:
            form.append(('Tag Data', None, tag_form),)
        if has_mapbc_data:
            form.append(('Map BC Data', None, mapbc_form),)

        results_form = []
        if len(results_form):
            form.append(('Results', None, results_form))
        print(form)
        return form, cases
