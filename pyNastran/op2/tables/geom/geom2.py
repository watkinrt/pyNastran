# pylint: disable=W0612,C0103,C0301,C0302,C0303,W0613,C0111,R0914,C0326,R0201
from struct import unpack, Struct
from six import b
from six.moves import range

from pyNastran.bdf.cards.elements.elements import CGAP
from pyNastran.bdf.cards.elements.damper import (CDAMP1, CDAMP2, CDAMP3,
                                                 CDAMP4, CDAMP5, CVISC)
from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.elements.shell import (CTRIA3, CQUAD4, CTRIA6,
                                                CQUADR, CQUAD8, CQUAD, CQUADX,
                                                CSHEAR)
from pyNastran.bdf.cards.elements.rods import CROD, CTUBE, CONROD
from pyNastran.bdf.cards.elements.bars import CBAR
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2,
                                               CMASS3, CMASS4)
from pyNastran.bdf.cards.elements.solid import (CTETRA4, CTETRA10, CPENTA6,
                                                CPENTA15, CHEXA8, CHEXA20)
from pyNastran.bdf.cards.thermal.thermal import CHBDYG, CONV  # , CONVM, CHBDYP
from pyNastran.bdf.cards.nodes import SPOINTs
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class GEOM2(GeomCommon):

    def _read_geom2_4(self, data, ndata):
        return self._read_geom_4(self._geom2_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom2_map = {
            (2408,   24,  180): ['CBAR', self._read_cbar],       # record 8
            (4001,   40,  275): ['CBARAO', self._read_cbarao],   # record 9  - not done
            (5408,   54,  261): ['CBEAM', self._read_cbeam],     # record 10
            (11401, 114, 9016): ['CBEAMP', self._read_cbeamp],   # record 11 - not done
            (4601,   46,  298): ['CBEND', self._read_cbend],     # record 12 - not done
            (2608,  26,   60): ['CBUSH', self._read_cbush],      # record 13 - not done
            (5608,   56,  218): ['CBUSH1D', self._read_cbush1d], # record 14 - not done
            (2315,   23,  146): ['CCONE', self._read_ccone],     # record 15 - not done
            (201,     2,   69): ['CDAMP1', self._read_cdamp1],      # record 16
            (301,     3,   70): ['CDAMP2', self._read_cdamp2],      # record 17
            (401,     4,   71): ['CDAMP3', self._read_cdamp3],      # record 18
            (501,     5,   72): ['CDAMP4', self._read_cdamp4],      # record 19
            (10608, 106, 404): ['CDAMPS', self._read_cdamp5], # record 20

            (601, 6, 73): ['CELAS1', self._read_celas1],      # record 29
            (701, 7, 74): ['CELAS2', self._read_celas2],      # record 30
            (801, 8, 75): ['CELAS3', self._read_celas3],      # record 31
            (901, 9, 76): ['CELAS4', self._read_celas4],      # record 32
            # record 33
            # record 34
            # record 35
            (8515, 85, 209): ['CFLUID2', self._readCFLUID2],   # record 35 - not done
            (8615, 86, 210): ['CFLUID3', self._readCFLUID3],   # record 36 - not done
            (8715, 87, 211): ['CFLUID4', self._readCFLUID4],   # record 37 - not done
            (1908, 19, 104): ['CGAP', self._read_cgap],       # record 39 - buggy
            # record 40
            # record 41
            # record 42
            (10808, 108, 406): ['CHBDYG', self._readCHBDYG], # record 43
            (10908, 109, 407): ['CHBDYP', self._readCHBDYP], # record 44 - not done
            (7308,  73,  253): ['CHEXA', self._readCHEXA],     # record 45
            (1001,  10,   65): ['CMASS1', self._read_cmass1],    # record 51
            (1101,  11,   66): ['CMASS2', self._read_cmass2],    # record 52
            (1201,  12,   67): ['CMASS3', self._read_cmass3],    # record 53
            (1301,  13,   68): ['CMASS4', self._read_cmass4],    # record 54
            (2508,  25,    0): ['CMFREE', self._readCMFREE],     # record 55 - not done
            (1401,  14,   63): ['CONM1', self._read_conm1],      # record 56 - not done
            (1501,  15,   64): ['CONM2', self._read_conm2],      # record 57
            (1601,  16,   47): ['CONROD', self._read_conrod],    # record 58
            (12701, 127, 408): ['CONV', self._readCONV],     # record 59 - not tested
            (8908,  89, 422): ['CONVM', self._readCONVM],     # record 60 - not tested
            # record 61
            (4108,  41, 280): ['CPENTA', self._readCPENTA],   # record 62
            # record 63
            # record 64
            # record 65
            # record 66
            # record 67
            (9108,   91,  507): ['CQUAD', self._readCQUAD],     # record 68 - not tested
            (2958,   51,  177): ['CQUAD4', self._readCQUAD4],   # record 69 - maybe buggy on theta/Mcsid field
            (13900, 139, 9989): ['CQUAD4', self._readCQUAD4],# record 70 - maybe buggy on theta/Mcsid field
            (4701,   47,  326): ['CQUAD8', self._readCQUAD8],   # record 71 - maybe buggy on theta/Mcsid field
            # record 72
            # record 73
            (8009,   80,  367): ['CQUADR', self._readCQUADR],   # record 74 - not tested
            (9008,   90,  508): ['CQUADX', self._readCQUADX],   # record 75 - not tested
            # record 76
            # record 77
            # record 78
            # record 79
            (3001, 30, 48): ['CROD', self._read_crod],        # record 80
            # record 81
            # record 82
            # record 83
            # record 84
            # record 85
            (12201,122,9013): ['CTETP', self._readCTETP],    # record 86 - not done
            (5508,  55, 217): ['CTETRA', self._readCTETRA],   # record 87
            # record 88
            # record 89
            # record 90
            # record 91
            # record 92
            (5959, 59, 282): ['CTRIA3', self._read_ctria3],   # record 93 - maybe buggy on theta/Mcsid field
            # record 94
            (4801, 48, 327): ['CTRIA6', self._read_ctria6],   # record 95 - buggy
            # record 96
            # record 97
            (9200, 92, 385): ['CTRIAR', self._read_ctriar],   # record 98  - not done
            # record 99
            (6108, 61, 107): ['CTRIAX6', self._read_ctriax6], # record 100 - not done
            # record 101
            # record 102
            (3701, 37, 49): ['CTUBE', self._read_ctube],      # record 103
            (3901, 39, 50): ['CVISC', self._read_cvisc],      # record 104 - not done
            # record 105
            # record 106
            # record 107
            # record 108
            # record 109
            # record 110
            # record 111
            # record 112
            # record 113
            (5201, 52, 11):   ['PLOTEL', self._read_plotel],    # record 114 - not done
            # record 115
            # record 116
            # record 117
            (5551,  49,  105): ['SPOINT', self._read_spoint],   # record 118
            (11601,116, 9942): ['VUBEAM', self._read_vubeam],  # record 119 - not done
            (2108, 21, 224): ['', self._read_fake],
            (3101, 31, 61): ['', self._read_fake],
            (4301, 43, 28): ['', self._read_fake],
            (5601, 56, 296): ['', self._read_fake],
            (6908, 69, 115): ['', self._read_fake],
            (6808, 68, 114): ['', self._read_fake],
            (7409, 74, 9991): ['', self._read_fake],
            (7509, 75, 9992): ['', self._read_fake],
            (7609, 76, 9993): ['', self._read_fake],
            (8100, 81, 381): ['', self._read_fake],
            (8200, 82, 383): ['', self._read_fake],
            (8308, 83, 405): ['', self._read_fake],
            (11201, 112, 9940): ['', self._read_fake],
            (12801, 128, 417): ['', self._read_fake],
            (13900, 139, 9984): ['', self._read_fake],
            (14000, 140, 9990): ['', self._read_fake],
            (16000, 160, 9988): ['', self._read_fake],
            (16100, 161, 9986): ['', self._read_fake],
            (16300, 163, 9989): ['', self._read_fake],
            (16700, 167, 9981): ['', self._read_fake],
            (16800, 168, 9978): ['', self._read_fake],
            (16500, 165, 9987): ['', self._read_fake],
            (2708, 27, 59): ['', self._read_fake],
            (5008, 50, 258): ['', self._read_fake],
            (16400, 164, 9983) : ['', self._read_fake],
            (3201, 32, 478): ['', self._read_fake],
            (11000, 110, 6667): ['', self._read_fake],
            (12301, 123, 9921): ['', self._read_fake],
            (12401, 124, 9922): ['', self._read_fake],
            (12600, 126, 6661): ['', self._read_fake],
            (14700, 147, 6662): ['', self._read_fake],
            (7309, 73, 0): ['', self._read_fake],
            (17200, 172, 6663): ['', self._read_fake],
            (17300, 173, 6664): ['', self._read_fake],
            (11501, 115, 9941): ['', self._read_fake],    # record
            (12501, 125, 9923): ['', self._read_fake],    # record
            (3401, 34, 9600): ['', self._read_fake],    # record
            (2208, 22, 225): ['', self._read_fake],  # record
            (17000, 170, 9980): ['', self._read_fake],  # record
            (7701, 77, 8881): ['', self._read_fake],  # record
            (12901, 129, 482): ['', self._read_fake],  # record
            (7801, 78, 8883): ['', self._read_fake],  # record
            (4408, 44, 227): ['', self._read_fake],  # record
            (17100, 171, 9979): ['', self._read_fake],  # record
            (2901, 29, 9601): ['', self._read_fake],  # record
            (4508, 45, 228): ['', self._read_fake],  # record
            (16600, 166, 9985) : ['', self._read_fake],  # record
            (16200, 162, 9982) : ['', self._read_fake],  # record
            (16900, 169, 9977) : ['', self._read_fake],  # record

            (1701, 17, 980) : ['', self._read_fake],  # record
            (1801, 18, 986) : ['', self._read_fake],  # record
            (8801, 88, 984) : ['', self._read_fake],  # record
            (8401, 84, 985) : ['', self._read_fake],  # record
            (17200, 172, 1000) : ['', self._read_fake],  # record
            (23500, 235, 6662) : ['', self._read_fake],  # record
            (23800, 238, 6665) : ['', self._read_fake],  # record
            (23900, 239, 6666) : ['', self._read_fake],  # record
            (1976, 1, 1996) : ['', self._read_fake],  # record
            (6120, 1, 60434) : ['', self._read_fake],  # record
            (2024, 1001, 2024) : ['', self._read_fake],  # record
            (801, 1, 572) : ['', self._read_fake],  # record

            (5701, 57, 981) : ['', self._read_fake],  # record
            (5801, 58, 982) : ['', self._read_fake],  # record
            (6111, 61, 996) : ['', self._read_fake],  # record
            (6112, 61, 997) : ['', self._read_fake],  # record
            (6113, 61, 998) : ['', self._read_fake],  # record
            (6114, 61, 999) : ['', self._read_fake],  # record
            (3501, 35, 1) : ['', self._read_fake],  # record
            (1001, 100, 10000) : ['', self._read_fake],  # record
            (1118, 1, 1874) : ['', self._read_fake],  # record
            (1801, 18, 986) : ['', self._read_fake],  # record
            (7909, 79, 9946) : ['', self._read_fake],  # record
        }

    def add_element(self, elem, allow_overwrites=True):
        raise RuntimeError('this should be overwritten')

    def addOp2Element(self, elem):
        if elem.eid <= 0:
            self.log.debug(elem)
            return
        self.add_element(elem, allow_overwrites=True)
        #print(str(elem)[:-1])

# 1-AEROQ4 (???)
# AEROT3   (???)
# 1-BEAMAERO (1701,17,0)
# 2-CAABSF (2708,27,59)
# 3-CAXIF2 (2108,21,224)
# 4-CAXIF3 (2208,22,225)
# 5-CAXIF4 (2308,23,226)

    def _read_cbar(self, data, n):
        """
        CBAR(2408,24,180) - the marker for Record 8
        """
        nelements = (len(data) - n) // 64
        for i in range(nelements):
            edata = data[n:n + 64]  # 16*4
            f, = self.struct_i.unpack(edata[28:32])
            if f == 0:
                out = unpack(b(self._endian + '4i3f3i6f'), edata)
                (eid, pid, ga, gb, x1, x2, x3, f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:
                out = unpack(b(self._endian + '4i3f3i6f'), edata)
                (eid, pid, ga, gb, x1, x2, x3, f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:
                out = unpack(b(self._endian + '7if2i6f'), edata)
                (eid, pid, ga, gb, g0, junk, junk, f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a,
                            w2a, w3a, w1b, w2b, w3b], [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBAR.add_op2_data(data_in)
            self.addOp2Element(elem)
            n += 64
        self.card_count['CBAR'] = nelements
        return n

    def _read_cbarao(self, data, n):
        """
        CBARAO(4001,40,275) - the marker for Record 9
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CBARAO in GEOM2\n')
        return n

    def _read_cbeam(self, data, n):
        """
        CBEAM(5408,54,261) - the marker for Record 10
        """
        nelements = (len(data) - n) // 72
        for i in range(nelements):
            edata = data[n:n + 72]  # 18*4
            f, = self.struct_i.unpack(edata[40:44])
            if   f == 0:  # basic cid
                out = unpack(b(self._endian + '6i3f3i6f'), edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:  # global cid
                out = unpack(b(self._endian + '6i3f3i6f'), edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:  # grid option
                out = unpack(b(self._endian + '12i6f'), edata)
                (eid, pid, ga, gb, sa, sb, g0, xx, xx, f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            elem = CBEAM.add_op2_data(data_in)
            self.addOp2Element(elem)
            n += 72
        self.card_count['CBEAM'] = nelements
        return n

    def _read_cbeamp(self, data, n):
        """
        CBEAMP(11401,114,9016) - the marker for Record 11
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CBEAMP in GEOM2\n')
        return n

    def _read_cbend(self, data, n):
        """
        CBEND(4601,46,298) - the marker for Record 12
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CBEND in GEOM2\n')
        return n

    def _read_cbush(self, data, n):
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CBUSH in GEOM2\n')
        return n

    def _read_cbush1d(self, data, n):
        """
        CBUSH1D(5608,56,218) - the marker for Record 14
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CBUSH1D in GEOM2\n')
        return n

    def _read_ccone(self, data, n):
        """
        CCONE(2315,23,0) - the marker for Record 15
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CCONE in GEOM2\n')
        return n

    def _read_cdamp1(self, data, n):
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = unpack(b(self._endian + '6i'), edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP1=%s\n' % str(out))
            (eid, pid, g1, g2, c1, c2) = out
            elem = CDAMP1.add_op2_data(out)
            self.addOp2Element(elem)
            n += 24
        self.card_count['CDAMP1'] = nelements
        return n

    def _read_cdamp2(self, data, n):
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = unpack(b(self._endian + 'if4i'), edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP2=%s\n' % str(out))
            (eid, bdamp, g1, g2, c1, c2) = out
            elem = CDAMP2.add_op2_data(out)
            self.addOp2Element(elem)
            n += 24
        self.card_count['CDAMP2'] = nelements
        return n

    def _read_cdamp3(self, data, n):
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP3=%s\n' % str(out))
            (eid, pid, s1, s2) = out
            elem = CDAMP3.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CDAMP3'] = nelements
        return n

    def _read_cdamp4(self, data, n):
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        s = Struct(b(self._endian + 'ifii'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP4=%s\n' % str(out))
            (eid, bdamp, s1, s2) = out
            elem = CDAMP4(None, out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CDAMP4'] = nelements
        return n

    def _read_cdamp5(self, data, n):
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP5=%s\n' % str(out))
            (eid, pid, s1, s2) = out
            elem = CDAMP5.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CDAMP5'] = nelements
        return n

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7
# CDUM8
# CDUM9

    def _read_celas1(self, data, n):
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        ntotal = 24  # 6*4
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS1=%s\n' % str(out))
            (eid, pid, g1, g2, c1, c2) = out
            elem = CELAS1.add_op2_data(out)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CELAS1'] = nelements
        return n

    def _read_celas2(self, data, n):
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        s1 = Struct(b(self._endian + 'if4iff'))
        ntotal = 32
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+32]
            out = s1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS2=%s\n' % str(out))
            (eid, k, g1, g2, c1, c2, ge, s) = out
            elem = CELAS2.add_op2_data(out)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CELAS2'] = nelements
        return n

    def _read_celas3(self, data, n):
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        ntotal = 16  # 4*4
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+16]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS3=%s\n' % str(out))
            (eid, pid, s1, s2) = out
            elem = CELAS3.add_op2_data(out)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CELAS3'] = nelements
        return n

    def _read_celas4(self, data, n):
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        s = Struct(b(self._endian + 'ifii'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS4=%s\n' % str(out))
            (eid, k, s1, s2) = out
            elem = CELAS4.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CELAS4'] = nelements
        return n

# CFAST
# CFASTP

    def _readCFLUID2(self, data, n):
        """
        CFLUID2(8515,85,209) - the marker for Record 35
        """
        return n

    def _readCFLUID3(self, data, n):
        """
        CFLUID3(8615,86,210) - the marker for Record 36
        """
        return n

    def _readCFLUID4(self, data, n):
        """
        CFLUID4(8715,87,211) - the marker for Record 37
        """
        return n

# CINT

    def _read_cgap(self, data, n):
        """
        CGAP(1908,19,104) - the marker for Record 39
        """
        s1 = Struct(b(self._endian + '4i3fii'))
        nelements = (len(data) - n) // 36
        for i in range(nelements):
            edata = data[n:n + 36]  # 9*4
            out = s1.unpack(edata)
            (eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
            g0 = None
            f2, = self.struct_i.unpack(edata[28:32])
            assert f == f2, 'f=%s f2=%s' % (f, f2)
            if f == 2:
                g0 = self.struct_i.unpack(edata[16:20])
                x1 = None
                x2 = None
                x3 = None

            data_in = [eid, pid, ga, gb, g0, x1, x2, x3, cid]
            elem = CGAP.add_op2_data(data_in)
            self.addOp2Element(elem)
            n += 36
        self.card_count['CGAP'] = nelements
        return n

# CHACAB
# CHACBR
# CHBDYE
# CHBDYG

    def _readCHBDYG(self, data, n):
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        ntotal = 64  # 16*4
        s = Struct(b(self._endian + '16i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (eid, blank, Type, iviewf, iviewb, radmidf, radmidb, blank2,
             g1, g2, g3, g4, g5, g6, g7, g8) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYG=%s\n' % str(out))
            data_in = [eid, Type, iviewf, iviewb, radmidf, radmidb,
                       g1, g2, g3, g4, g5, g6, g7, g8]
            elem = CHBDYG()
            elem.add_op2_data(data_in)
            self.addOp2Element(elem)
            n += ntotal
        return n

    def _readCHBDYP(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping CHBDYP in GEOM2\n')
        return n

    def _readCHEXA(self, data, n):
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        s = Struct(b(self._endian + '22i'))
        ntotal = 88  # 22*4
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+88]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CHEXA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, ]
            big_nodes = [g9, g10, g11, g12, g13, g14, g15, g16,
                         g17, g18, g19, g20]
            if sum(big_nodes) > 0:
                elem = CHEXA20(None, data_in + big_nodes)
            else:
                elem = CHEXA8(None, data_in)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CHEXA'] = nelements
        return n

# CHEXA20F
# CHEXAFD
# CHEXAL
# CHEXP
# CHEXPR

    def _read_cmass1(self, data, n):
        """
        CMASS1(1001,10,65) - the marker for Record 51
        """
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CMASS1.add_op2_data(out)
            self.addOp2Element(elem)
            n += 24
        self.card_count['CMASS1'] = nelements
        return n

    def _read_cmass2(self, data, n):
        """
        CMASS2(1101,11,66) - the marker for Record 52
        """
        s = Struct(b(self._endian + 'if4i'))
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS2=%s\n' % str(out))
            #(eid, m, g1, g2, c1, c2) = out
            elem = CMASS2.add_op2_data(out)
            self.addOp2Element(elem)
            n += 24
        self.card_count['CMASS2'] = nelements
        return n

    def _read_cmass3(self, data, n):
        """
        CMASS3(1201,12,67) - the marker for Record 53
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CMASS3.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CMASS3'] = nelements
        return n

    def _read_cmass4(self, data, n):
        """
        CMASS4(1301,13,68) - the marker for Record 54
        """
        nelements = (len(data) - n) // 16
        s = Struct(b(self._endian + 'ifii'))
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            #(eid, m,s 1, s2) = out
            elem = CMASS4.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CMASS4'] = nelements
        return n

    def _readCMFREE(self, data, n):
        """
        CMFREE(2508,25,0) - the marker for Record 55
        """
        if self.is_debug_file:
            self.binary_debug.write('skipping CMFREE in GEOM2\n')
        return n

    def _read_conm1(self, data, n):
        """
        CONM1(1401,14,63) - the marker for Record 56
        """
        s = Struct(b(self._endian + '3i21f'))
        nelements = (len(data) - n) // 96
        for i in range(nelements):
            edata = data[n:n + 96]  # 24*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM1=%s\n' % str(out))
            (eid, g, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = out
            elem = CONM1()
            elem.add_op2_data(out)
            self.addOp2Element(elem)
            n += 96
        self.card_count['CONM1'] = nelements
        return n

    def _read_conm2(self, data, n):
        """
        CONM2(1501,15,64) - the marker for Record 57
        """
        ntotal = 52  # 13*4
        s = Struct(b(self._endian + '3i10f'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM2=%s\n' % str(out))
            (eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            elem = CONM2.add_op2_data(out)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CONM2'] = nelements
        return n

    def _read_conrod(self, data, n):
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        ntotal = 32  # 8*4
        s = Struct(b(self._endian + '4i4f'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONROD=%s\n' % str(out))
            (eid, n1, n2, mid, a, j, c, nsm) = out
            elem = CONROD.add_op2_data(out)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CONROD'] = nelements
        return n

    def _readCONV(self, data, n):
        """
        CONV(12701,127,408) - the marker for Record 59
        """
        #return
        ntotal = 80  # 20*4
        s = Struct(b(self._endian + '12i8f'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+80]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
            (eid, pconID, flmnd, cntrlnd,
             ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8,
             wt1, wt2, wt3, wt4, wt5, wt6, wt7, wt8) = out
            data_in = [eid, pconID, flmnd, cntrlnd,
                       [ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                       [wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
            elem = CONV(None, data_in)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CONV'] = nelements
        return n

    def _readCONVM(self, data, n):
        """
        CONVM(8908,89,422) - the marker for Record 60
        """
        return n
        ntotal = 28  # 7*4
        s = Struct(b(self._endian + '7i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+28]
            out = unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pconID, flmnd, cntrlnd,
             [ta1, ta2, ta3]) = out
            data_in = [
                eid, pconID, flmnd, cntrlnd,
                [ta1, ta2, ta3]]
            elem = CONVM(None, data_in)  # undefined
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CONVM'] = nelements
        return n

# CPENP

    def _readCPENTA(self, data, n):
        """
        CPENTA(4108,41,280) - the marker for Record 62
        """
        s = Struct(b(self._endian + '17i'))
        nelements = (len(data) - n) // 68
        for i in range(nelements):
            edata = data[n:n + 68]  # 17*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CPENTA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6]
            bigNodes = [g7, g8, g9, g10, g11, g12, g13, g14, g15]
            if sum(bigNodes) > 0:
                elem = CPENTA15(None, data_in + bigNodes)
            else:
                elem = CPENTA6(None, data_in)
            self.addOp2Element(elem)
            n += 68
        self.card_count['CPENTA'] = nelements
        return n

# CPENPR
# CPENT15F
# CPENT6FD
# CQDX4FD
# CQDX9FD

    def _readCQUAD(self, data, n):
        """
        CQUAD(9108,91,507)  - the marker for Record 68
        """
        return self.run_cquad(data, n, CQUAD)

    def run_cquad(self, data, n, Element):
        """common method for CQUAD, CQUADX"""
        s = Struct(b(self._endian + '11i'))
        nelements = (len(data) - n) // 44  # 11*4
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))
        for i in range(nelements):
            edata = data[n:n + 44]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (Element.type, str(out)))
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            data = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9]
            elem = Element(None, data)
            self.addOp2Element(elem)
            n += 44
        self.card_count[Element.type] = nelements
        return n

    def _readCQUAD4(self, data, n):
        """
        CQUAD4(2958,51,177)    - the marker for Record 69
        CQUAD4(13900,139,9989) - the marker for Record 70
        """
        return self.run_cquad4(data, n, CQUAD4)

    def run_cquad4(self, data, n, Element):
        """
        common method for CQUAD4, CQUADR
        """
        nelements = (len(data) - n) // 56
        s = Struct(b(self._endian + '6iffii4f'))
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))
        for i in range(nelements):
            edata = data[n:n + 56]  # 14*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
             t1, t2, t3, t4) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (Element.type, str(out)))
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4))
            data_init = [
                eid, pid, n1, n2, n3, n4, theta, zoffs,
                tflag, t1, t2, t3, t4]
            elem = Element.add_op2_data(data_init)
            self.addOp2Element(elem)
            n += 56
        self.card_count[Element.type] = nelements
        return n

# CQUAD4FD

    def _readCQUAD8(self, data, n):
        """
        CQUAD8(4701,47,326)  - the marker for Record 71
        .. warning:: inconsistent with dmap manual
        """
        #return n
        nelements = (len(data) - n) // 64  # 17*4
        s = Struct(b(self._endian + '10i5fi'))
        for i in range(nelements):
            edata = data[n:n + 64]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CQUAD8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             t3, t4, theta, tflag) = out
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            self.addOp2Element(elem)
            n += 64
        self.card_count['CQUAD8'] = nelements
        return n

# CQUAD9FD
# CQUADP
    def _readCQUADR(self, data, n):
        """
        CQUADR(8009,80,367)  - the marker for Record 74
        """
        return self.run_cquad4(data, n, CQUADR)

    def _readCQUADX(self, data, n):
        """
        CQUADX(9008,90,508)  - the marker for Record 75
        """
        #print("reading CQUADX")
        return self.run_cquad4(data, n, CQUADX)

# CRBAR
# CRBE1
# CRBE3
# CRJOINT

    def _read_crod(self, data, n):
        """
        CROD(3001,30,48)    - the marker for Record 80
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16  # 4*4
        for i in range(nelements):
            edata = data[n:n + 16]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CROD=%s\n' % str(out))
            (eid, pid, n1, n2) = out
            elem = CROD.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CROD'] = nelements
        return n

# CRROD
# CSEAM

    def _readCSHEAR(self, data, n):
        """
        CSHEAR(3101,31,61)    - the marker for Record 83
        """
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // 24  # 6*4
        for i in range(nelements):
            edata = data[n:n + 24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CSHEAR=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4) = out
            elem = CSHEAR(None, out)
            self.addOp2Element(elem)
            n += 24
        self.card_count['CSHEAR'] = nelements
        return n

# CSLOT3
# CSLOT4

    def _readCTETP(self, data, n):
        """
        CTETP(12201,122,9013)    - the marker for Record 86
        .. todo:: create object
        """
        #raise NotImplementedError('needs work...')
        nelements = (len(data) - n) // 108  # 27*4
        s = Struct(b(self._endian + '27i'))
        for i in range(nelements):
            edata = data[n:n+108]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETP=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
             f1, f2, f3, f4, b1, ee1, ee2, ee3, ee4) = out
            #print("out = ",out)
            e = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12]
            f = [f1, f2, f3, f4]
            ee = [ee1, ee2, ee3, ee4]

            #print("e  = ",e)
            #print("f  = ",f)
            #print("b1  = ",b1)
            #print("ee = ",ee)
            data_in = [eid, pid, n1, n2, n2, n3, n4]
            elem = CTETRA4(None, data_in)
            self.addOp2Element(elem)

    def _readCTETRA(self, data, n):
        """
        CTETRA(5508,55,217)    - the marker for Record 87
        """
        s = Struct(b(self._endian + '12i'))
        nelements = (len(data) - n)// 48  # 12*4
        for i in range(nelements):
            edata = data[n:n + 48]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETRA=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = out
            #print("out = ",out)

            data_in = [eid, pid, n1, n2, n3, n4]
            bigNodes = [n5, n6, n7, n8, n9, n10]
            if sum(bigNodes) > 0:
                elem = CTETRA10(None, data_in + bigNodes)
            else:
                elem = CTETRA4(None, data_in)
            self.addOp2Element(elem)
            n += 48
        self.card_count['CTETRA'] = nelements
        return n

# CTETPR
# CTETR10F
# CTETR4FD
# CTQUAD
# CTTRIA

    def _read_ctria3(self, data, n):
        """
        CTRIA3(5959,59,282)    - the marker for Record 93
        """
        ntotal = 52  # 13*4
        s = Struct(b(self._endian + '5iff3i3f'))
        nelements = (len(data) - n)// 52  # 13*4
        for i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3))
            (eid, pid, n1, n2, n3, theta, zoffs, blank1,
             blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA3=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            self.addOp2Element(elem)
            n += ntotal
        self.card_count['CTRIA3'] = nelements
        return n


# CTRIAFD

    def _read_ctria6(self, data, n):
        """
        CTRIA6(4801,48,327)    - the marker for Record 95
        .. warning:: inconsistent with dmap manual
        """
        s = Struct(b(self._endian + '8i4fi'))
        nelements = (len(data) - n) // 52  # 13*4
        for i in range(nelements):
            edata = data[n:n + 52]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, t1, t2, t3, tflag) = out
            elem = CTRIA6.add_op2_data(out)
            self.addOp2Element(elem)
            n += 52
        self.card_count['CTRIA6'] = nelements
        return n

# CTRIA6FD
# CTRIAP

    def _read_ctriar(self, data, n):  # 98
        if self.is_debug_file:
            self.binary_debug.write('skipping CTRIAR in GEOM2\n')
        return n

# CTRIAX

    def _read_ctriax6(self, data, n):  # 100
        if self.is_debug_file:
            self.binary_debug.write('skipping CTRIAX6 in GEOM2\n')
        return n

# CTRIX3FD
# CTRIX6FD

    def _read_ctube(self, data, n):
        """
        CTUBE(3701,37,49)    - the marker for Record 103
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTUBE=%s\n' % str(out))
            (eid, pid, n1, n2) = out
            elem = CTUBE.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CTUBE'] = nelements
        return n

    def _read_cvisc(self, data, n):
        """CVISC(3901,39, 50) - the marker for Record 104"""
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CVISC=%s\n' % str(out))
            #(eid,pid,n1,n2) = out
            elem = CVISC.add_op2_data(out)
            self.addOp2Element(elem)
            n += 16
        self.card_count['CVISC'] = nelements
        return n

    def _read_cweld(self, data, n):  # 105
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELD in GEOM2\n')
        return n

    def _read_cweldc(self, data, n):  # 106
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDC in GEOM2\n')
        return n

    def _read_cweldg(self, data, n):  # 107
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDG in GEOM2\n')
        return n

# CWSEAM
# GENEL
# GMDNDC
# GMBNDS
# GMINTC
# GMINTS
    def _read_plotel(self, data, n):  # 114
        if self.is_debug_file:
            self.binary_debug.write('skipping PLOTEL in GEOM2\n')
        return n
# RADBC
# RADINT
# SINT

    def add_spoint(self, spooint):
        raise RuntimeError('this should be overwritten by the BDF')

    def _read_spoint(self, data, n):
        """
        (5551,49,105)    - the marker for Record 118
        """
        npoints = (len(data) - n) // 4
        fmt = b'%ii' % npoints
        nids = unpack(fmt, data[n:])
        if self.is_debug_file:
            self.binary_debug.write('SPOINT=%s\n' % str(nids))
        spoint = SPOINTs.add_op2_data(list(nids))
        self.add_spoint(spoint)
        self.card_count['SPOINT'] = npoints
        return n

    def _read_vubeam(self, data, n):  # 119
        if self.is_debug_file:
            self.binary_debug.write('skipping VUBEAM in GEOM2\n')
        return n

# VUHEXA
# VUQUAD4
# VUPENTA
# VUTETRA
# VUTRIA
# VUBEAM
# VUHEXA
# VUQUAD4
# WELDP
