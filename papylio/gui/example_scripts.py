example_scripts = {
    "Histogram current file": """
#Create a selection:
file.create_selection(variable='intensity', channel='red', aggregator='max', operator='>', threshold=3000, name='selection_red')
file.apply_selections()
#Plot the histogram:
fig, ax = file.show_histogram('intensity', selected=True, bins=100)
fig.show()
""",

    "Show mapping": """
#pick your mapping file:
mapping_file = experiment.files[0]
#show it:
figure, axis = mapping_file.show_image()
mapping_file.mapping.show(axis=axis, show_source=True)
figure.show()
""",

    "Show selected files": """
for file in experiment.selectedFiles:
    print(file)
""",

    "Show average image": """
file.show_average_image()
""",

    "Preset_B": """
# Any code
print("Hello world")
""",

    "Preset_C": """
# Any code
print("Hello world")
""",


}