example_scripts_show = {
    "Show mapping": """
#pick your mapping file:
mapping_file = experiment.files[0]
#show it:
figure, axis = mapping_file.show_image()
mapping_file.mapping.show(axis=axis, show_source=True)
figure.show()
""",

    "Show Average Image": """
#show the average image of the selected file
file.show_average_image()
self.update_plots()
""",

    "Show 1D Histogram": """
#Create a selection [WILL BE ADDED TO CURRENT SELECTION RULES]:
#file.create_selection(variable='intensity', channel='red', aggregator='max', operator='>', threshold=3000, name='selection_red')
#file.apply_selections()
#Plot the histogram:
fig, ax = file.show_histogram('intensity', selected=True, bins=100)
fig.show()
""",

    "Show 2D Histogram": """
#show a 2-D histogram of the selected file, with a user-set range
import matplotlib.pyplot as plt
marginal_hist2d_kwargs = {
    'range': ((-10E3, 30E3), (-5E3, 15E3))
}
file.histogram_2D_intensity_per_channel(    bins=100,
    show_marginal=False,
    **marginal_hist2d_kwargs,
)
plt.show()

""",
}

example_scripts_select = {
    "Show selected files": """
for file in experiment.selectedFiles:
    print(file)
""",

    "MyScript": """
None
""",
}

example_scripts_other= {
    "MyScript": """
None
""",
}

example_scripts = {
    "Show": example_scripts_show,
    "Select": example_scripts_select,
    "Analyze": example_scripts_other,
}