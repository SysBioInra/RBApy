function model = load_model(path_)

model = py.rba.RbaModel()
model = model.from_xml(path_)
matrices = model.get_matrices()
