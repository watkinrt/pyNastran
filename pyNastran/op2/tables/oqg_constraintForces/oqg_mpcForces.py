from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray


class RealMPCForcesArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class ComplexMPCForcesArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                         C O M P L E X   F O R C E S   O F   M U L T I   P O I N T   C O N S T R A I N T\n', ]
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase, is_sort1)
