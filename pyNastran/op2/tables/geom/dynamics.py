from struct import unpack

from pyNastran.bdf.cards.loads.loads import DAREA
from pyNastran.op2.tables.geom.geom_common import GeomCommon

class DYNAMICS(GeomCommon):
    def _read_dynamics_4(self, data, ndata):
        return self._read_geom_4(self._dynamics_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._dynamics_map = {
            (5307, 53, 379) : ['ACSRCE', self._read_fake], # 1
            #(27, 17, 182): ['DAREA', self._read_darea],  # 2
            (27, 17, 182): ['DAREA', self._read_fake],  # 2

            (37, 18, 183): ['DELAY', self._read_delay],  # 3
            (57, 5, 123): ['DLOAD', self._read_dload],  # 4
            (77, 19, 184): ['DPHASE', self._read_dphase],  # 5
            (107, 1, 86): ['EIGB', self._read_eigb],   # 7
            (207, 2, 87): ['EIGC', self._read_eigc],   # 8
            (257, 4, 158): ['EIGP', self._read_eigp],   # 9
            (307, 3, 85): ['EIGR', self._read_eigr],   # 10

            (308, 8, 348): ['EIGRL', self._read_eigrl],  # 11
            (707, 7, 124): ['EPOINT', self._read_epoint],  # 12
            (1307, 13, 126): ['FREQ', self._read_freq],   # 13
            (1007, 10, 125): ['FREQ1', self._read_freq1],  # 14
            (1107, 11, 166): ['FREQ2', self._read_freq2],  # 15
            (1407, 14, 39): ['FREQ3', self._read_freq3],  # 16
            (1507, 15, 40): ['FREQ4', self._read_freq4],  # 17
            (1607, 16, 41): ['FREQ5', self._read_freq5],  # 18

            (3707, 37, 556) : ['NLRGAP', self._read_fake], # 19

            (3107, 31, 127): ['NONLIN1', self._read_fake], # 20
            (3207, 32, 128): ['NONLIN2', self._read_fake], # 21
            (3207, 33, 129): ['NONLIN3', self._read_fake], # 22
            (3207, 34, 130): ['NONLIN4', self._read_fake], # 23
            (2107, 21, 195): ['RANDPS', self._read_fake], # 24
            (2207, 22, 196): ['RANDT1', self._read_fake], # 25
            (5107, 51, 131): ['RLOAD1', self._read_rload1],  # 26
            (5207, 52, 132): ['RLOAD2', self._read_rload2],  # 27
            (8910, 89, 606): ['ROTORB', self._read_fake],  # 28
            (8210, 82, 599): ['ROTORD', self._read_fake],  # 29
            (8410, 84, 600): ['ROTORG', self._read_fake],  # 30
            (5707, 57, 135): ['SEQEP', self._read_fake],  # 31
            (6207, 62, 136): ['TF', self._read_fake],  # 32
            (6607, 66, 137): ['TIC', self._read_fake],  # 33
            (7107, 71, 138): ['TLOAD1', self._read_tload1],  # 37
            (7207, 72, 139): ['TLOAD2', self._read_tload2],  # 38
            (8307, 83, 142): ['TSTEP', self._read_tstep],  # 39

            (10701, 107, 117) : ['', self._read_fake],
            (10801, 108, 242) : ['', self._read_fake],
            (3807, 38, 505) : ['', self._read_fake],
            (4807, 48, 306) : ['', self._read_fake],
            (11001, 110, 310) : ['', self._read_fake],
            (10901, 109, 260) : ['', self._read_fake],
            (3307, 33, 129) : ['', self._read_fake],
            (11101, 111, 368) : ['', self._read_fake],
        }

#ACSRCE (5307,53,379)

    def _read_area(self, data, n):
        """DAREA(27,17,182) - the marker for Record 2"""
        #print("reading DAREA")
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('DAREA', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iiff', edata)
            #(sid,p,c,a) = out
            darea = DAREA(data=out)
            self.add_DAREA(darea)
            n += ntotal
        return n

    def _read_delay(self, data, n):
        """DELAY(37,18,183) - Record 3"""
        if self.is_debug_file:
            self.binary_debug.write('skipping DELAY in DYNAMICS\n')
        return len(data)

    def _read_dload(self, data, n):
        """DLOAD(57,5,123) - Record 4"""
        if self.is_debug_file:
            self.binary_debug.write('skipping DLOAD in DYNAMICS\n')
        return len(data)

    def _read_dphase(self, data, n):
        """DPHASE(77,19,184) - Record 5"""
        if self.is_debug_file:
            self.binary_debug.write('skipping DPHASE in DYNAMICS\n')
        return len(data)

#DYNRED(4807,48,306)

    def _read_eigb(self, data, n):
        """EIGB(107,1,86) - Record 7"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGB in DYNAMICS\n')
        return len(data)

    def _read_eigc(self, data, n):
        """EIGC(207,2,87) - Record 8"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGC in DYNAMICS\n')
        return len(data)

    def _read_eigp(self, data, n):
        """EIGP(257,4,158) - Record 9"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGP in DYNAMICS\n')
        return len(data)

    def _read_eigr(self, data, n):
        """EIGR(307,3,85) - Record 10"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGR in DYNAMICS\n')
        return len(data)

    def _read_eigrl(self, data, n):
        """EIGRL(308,8,348) - Record 11"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGRL in DYNAMICS\n')
        return len(data)

    def _read_epoint(self, data, n):
        """EPOINT(707,7,124) - Record 12"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EPOINT in DYNAMICS\n')
        return len(data)

    def _read_freq(self, data, n):
        """FREQ(1307,13,126) - Record 13"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ in DYNAMICS\n')
        return len(data)

    def _read_freq1(self, data, n):
        """FREQ1(1007,10,125) - Record 14"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ1 in DYNAMICS\n')
        return len(data)

    def _read_freq2(self, data, n):
        """FREQ2(1107,11,166) - Record 15"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ2 in DYNAMICS\n')
        return len(data)

    def _read_freq3(self, data, n):
        """FREQ3(1407,14,39) - Record 16"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ3 in DYNAMICS\n')
        return len(data)

    def _read_freq4(self, data, n):
        """FREQ4(1507,15,40) - Record 17"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ4 in DYNAMICS\n')
        return len(data)

    def _read_freq5(self, data, n):
        """FREQ5(1607,16,41) - Record 18"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ5 in DYNAMICS\n')
        return len(data)

#NLRSFD
#NOLIN1
#NOLIN2
#NOLIN3
#NOLIN4
#RANDPS
#RANDT1

    def _read_rload1(self, data, n):
        """RLOAD1(5107,51,131) - Record 26"""
        if self.is_debug_file:
            self.binary_debug.write('skipping RLOAD1 in DYNAMICS\n')
        return len(data)

    def _read_rload2(self, data, n):
        """RLOAD2(5107,51,131) - Record 27"""
        if self.is_debug_file:
            self.binary_debug.write('skipping RLOAD2 in DYNAMICS\n')
        return len(data)

#
#RLOAD2(5207,52,132)
#RGYRO
#ROTORG
#RSPINR
#RSPINT
#SEQEP(5707,57,135)
#TF
#TIC
#TIC
#TIC3

    def _read_tload1(self, data, n):
        """TLOAD1(7107,71,138) - Record 37"""
        if self.is_debug_file:
            self.binary_debug.write('skipping TLOAD1 in DYNAMICS\n')
        return len(data)

    def _read_tload2(self, data, n):
        """TLOAD2(7207,72,139) - Record 37"""
        if self.is_debug_file:
            self.binary_debug.write('skipping TLOAD2 in DYNAMICS\n')
        return len(data)

    def _read_tstep(self, data, n):
        """TSTEP(8307,83,142) - Record 38"""
        if self.is_debug_file:
            self.binary_debug.write('skipping TSTEP in DYNAMICS\n')
        return len(data)

#UNBALNC
