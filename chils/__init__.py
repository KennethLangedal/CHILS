import ctypes
import os
import platform

_lib_name = "libCHILS.so"
if platform.system() == "Windows":
    _lib_name = "libCHILS.dll"

_lib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), _lib_name))

_lib.chils_initialize.restype = ctypes.c_void_p
_lib.chils_release.argtypes = [ctypes.c_void_p]
_lib.chils_add_vertex.argtypes = [ctypes.c_void_p, ctypes.c_longlong]
_lib.chils_add_vertex.restype = ctypes.c_int
_lib.chils_add_edge.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
_lib.chils_run_full.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_uint]
_lib.chils_run_local_search_only.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_uint]
_lib.chils_solution_get_size.argtypes = [ctypes.c_void_p]
_lib.chils_solution_get_size.restype = ctypes.c_int
_lib.chils_solution_get_weight.argtypes = [ctypes.c_void_p]
_lib.chils_solution_get_weight.restype = ctypes.c_longlong
_lib.chils_solution_get_time.argtypes = [ctypes.c_void_p]
_lib.chils_solution_get_time.restype = ctypes.c_double
_lib.chils_solution_get_vertex_configuration.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.chils_solution_get_vertex_configuration.restype = ctypes.c_int

class CHILS:
    def __init__(self):
        self.solver = _lib.chils_initialize()

    def __del__(self):
        _lib.chils_release(self.solver)

    def add_vertex(self, weight):
        return _lib.chils_add_vertex(self.solver, weight)

    def add_edge(self, u, v):
        _lib.chils_add_edge(self.solver, u, v)

    def run_full(self, time_limit, n_solutions, seed):
        _lib.chils_run_full(self.solver, time_limit, n_solutions, seed)

    def run_local_search_only(self, time_limit, seed):
        _lib.chils_run_local_search_only(self.solver, time_limit, seed)

    def get_solution_size(self):
        return _lib.chils_solution_get_size(self.solver)

    def get_solution_weight(self):
        return _lib.chils_solution_get_weight(self.solver)

    def get_solution_time(self):
        return _lib.chils_solution_get_time(self.solver)

    def get_solution_vertex_configuration(self, u):
        return _lib.chils_solution_get_vertex_configuration(self.solver, u)
