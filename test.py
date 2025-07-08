import chils

solver = chils.CHILS()

solver.add_vertex(10)
solver.add_vertex(10)
solver.add_vertex(20)

solver.add_edge(0, 1)
solver.add_edge(0, 2)

solver.run_full(20.0, 4, 0)

print(solver.get_solution_weight())