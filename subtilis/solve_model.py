import os.path
import sys
sys.path.append(os.path.join(sys.path[0],'..'))
import rba

model = rba.RbaModel.from_xml('subtilis')
solver = model.solve('medium_2')
