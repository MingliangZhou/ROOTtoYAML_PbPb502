import pandas as pd

from hepdata_lib import Submission
from hepdata_lib import Table
from hepdata_lib import RootFileReader
from hepdata_lib import Variable
from hepdata_lib import Uncertainty

submission = Submission()

df = pd.read_csv("input.csv")
df = df.fillna("")
for index, fig in df.iterrows():
	print(fig["figure_name"])
	# create a table
	table = Table(fig["figure_name"])
	table.description = fig["description"]
	table.location = fig["paper_location"]

	# read figures
	reader = RootFileReader(fig["file_location"])
	if fig["type_stat"].lower() in ["tgraph", "tgrapherrors", "tgraphasymmerrors"]:
		stat = reader.read_graph(fig["name_stat"])
	elif fig["type_stat"].lower() == "th1":
		stat = reader.read_hist_1d(fig["name_stat"])
	elif fig["type_stat"].lower() == "th2":
		stat = reader.read_hist_2d(fig["name_stat"])
	else:
		print("ERROR: {}, type not recognized!".format(fig["figure_name"]))
	
	if fig["type_syst"].lower() in ["tgraph", "tgrapherrors", "tgraphasymmerrors"]:
		syst = reader.read_graph(fig["name_syst"])
	elif fig["type_syst"].lower() == "th1":
		syst = reader.read_hist_1d(fig["name_syst"])
	elif fig["type_syst"].lower() == "th2":
		syst = reader.read_hist_2d(fig["name_syst"])
	else:
		print("WARNING: {}, systematic errors not found!".format(fig["figure_name"]))

	# read points
	if fig["type_stat"].lower() == "th2":
		x1 = Variable(fig["x1_name"], is_independent=True, is_binned=False, units=fig["x1_units"])
		x1.values = stat["x"]
		x2 = Variable(fig["x2_name"], is_independent=True, is_binned=False, units=fig["x2_units"])
		x2.values = stat["y"]
		y = Variable(fig["y_name"], is_independent=False, is_binned=False, units=fig["y_units"])
		y.values = stat["z"]
	else:
		x1 = Variable(fig["x1_name"], is_independent=True, is_binned=False, units=fig["x1_units"])
		x1.values = stat["x"]
		y = Variable(fig["y_name"], is_independent=False, is_binned=False, units=fig["y_units"])
		y.values = stat["y"]

	if fig["type_stat"].lower() == "tgraphasymmerrors":
		y_stat = Uncertainty("stat. uncertainty", is_symmetric=False)
		y_stat.values = stat["dy"]
		y.add_uncertainty(y_stat)
	elif fig["type_stat"].lower() in ["tgrapherrors", "th1"]:
		y_stat = Uncertainty("stat. uncertainty", is_symmetric=True)
		y_stat.values = stat["dy"]
		y.add_uncertainty(y_stat)

	if fig["type_syst"].lower() == "tgraphasymmerrors":
		y_syst = Uncertainty("syst. uncertainty", is_symmetric=False)
		y_syst.values = syst["dy"]
		y.add_uncertainty(y_syst)
	elif fig["type_syst"].lower() in ["tgrapherrors", "th1"]:
		y_syst = Uncertainty("syst. uncertainty", is_symmetric=True)
		y_syst.values = syst["dy"]
		y.add_uncertainty(y_syst)

	# write table
	if fig["type_stat"].lower() == "th2":
		table.add_variable(x1)
		table.add_variable(x2)
		table.add_variable(y)
	else:
		table.add_variable(x1)
		table.add_variable(y)
	submission.add_table(table)

submission.create_files("output")
