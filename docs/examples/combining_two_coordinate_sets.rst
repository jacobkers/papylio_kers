Combining two coordinate sets
========================================

Say, you have two file sets that is measured under different lasers: files_g and files_r
1. Find coordinates from files_g with a config setting
2. Find coordinates from files_r with another config setting
3. Merge this these two coordinates
4. extract traces

.. code-block:: python

    # --- find coordinates with config 1
    config = exp.configuration['find_coordinates']
    config['channels'] = ['donor']
    config['peak_finding']['minimum_intensity_difference'] = 400

    files_g.find_coordinates(configuration=config)

    # --- find coordinates with config 2
    config = exp.configuration['find_coordinates']
    config['channels'] = ['acceptor']
    config['peak_finding']['minimum_intensity_difference'] = 800

    files_r.find_coordinates(configuration=config)


    # merge coordinates
    for file_g, file_r in tqdm.tqdm(zip(files_g, files_r)):   # the first movie was done manually already
        # merge the two coordinates
        new_coord = xr.concat([file_g.coordinates, file_r.coordinates], dim='molecule')

        # [optional] make custom labels of your choice (This may be useful in your downstream analysis)
        laser_per_molecule = xr.DataArray([0]*file_g.number_of_molecules + [1]*file_r.number_of_molecules, dims=('molecule'), name='laser_per_molecule')

        # assign the new combined coordinates to all files
        file_g.coordinates = file_r.coordinates = new_coord
        laser_per_molecule.to_netcdf(file_g.absoluteFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
        laser_per_molecule.to_netcdf(file_r.absoluteFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')


    # extract traces
    files_g.extract_traces()
    files_r.extract_traces()