Link single-molecule and sequencing data
===========================================

To link single-molecule and sequencing data, go through the following steps:

`1. Find coordinates and extract traces from the single-molecule movies`_

`2. Generate tile mappings`_

`3. Generate and fine-tune sequencing matches`_

`4. Link single-molecules and sequences`_

`5. Combine all nc files into a single dataset`_


1. Find coordinates and extract traces from the single-molecule movies
----------------------------------------------------------------------
There are no changes here compared to a regular single-molecule analysis (See: :doc:`/Getting started`).

2. Generate tile mappings
----------------------------------------------------
Import the sam file with all sequences aligned to your reference sequence(s) (See: :doc:`/examples/sequence_alignment` on how to generate the sam file) as shown below. With ``sequence_subset``
you can indicate the indices of the nucleotides of interest.

.. code-block:: python

    aligned_sam_filepath = r'path_to_sam_file'
    sequence_subset = [12, 13, 14, 15, 16, 17, 18]
    exp.import_sequencing_data(Path(aligned_sam_filepath), remove_duplicates=True,
                           add_aligned_sequence=True, extract_sequence_subset=sequence_subset)

Now you can generate the tile mappings. This means that the transformation for each sequencing tile with respect to the
single-molecule fluorescence data is determined and applied. As ``files_of_interest`` you give the files used in
`step 1 <#find-coordinates-and-extract-traces-from-the-single-molecule-movies>`_ and as ``mapping_sequence_name`` you
give the name of the reference sequence that corresponds to the molecules detected in `step 1 <#find-coordinates-and-extract-traces-from-the-single-molecule-movies>`_
Set ``surface`` to 0 for objective-type TIRF and to 1 for prism-type TIRF. As ``scale`` and ``rotation``
for the ``AffineTransform`` you enter the previously determined values for this specific fluorescence set-up and
sequencer combination. Finally, to determine the translation for each sequencing tile, cross-correlation is performed
between the single-molecule fluorescence data and sequencing data.

.. code-block:: python

    mapping_sequence_name = 'my_sequence'
    exp.generate_tile_mappings(files_of_interest, mapping_sequence_name=mapping_sequence_name, surface=0)
    exp.tile_mappings.transformation = AffineTransform(scale=[29.53, -29.53], rotation=0.0039328797312210935)
    exp.tile_mappings.serial.cross_correlation(divider=20, kernel_size=7, gaussian_sigma=1, plot=False)

To evaluate whether tile mapping was successful you can plot the translation in x and y for each tile:

.. code-block:: python

    exp.tile_mappings.scatter_parameters('translation', 'translation', 'x', 'y', save=False)

In case of successful tile mapping, the tiles should form a line separated ~50-100 in x and ~25000 in y (MiSeq coordinates).

In case the tile mapping did not work well, try again with different settings for the cross-correlation step.
When you are satisfied with the tile mappings, you can save them:

.. code-block:: python

    exp.tile_mappings.save()


3. Generate and fine-tune sequencing matches
------------------------------
In this step, each single-molecule image is matched to the sequencing tiles. For each single-molecule image there might
be slight deviations in the mapping due to image aberrations or inaccuracies in the stage position. Therefore, the
translation, rotation and scaling for each single-molecule image are fine-tuned. Make sure to exclude the molecules that
you do not expect to see in the fluorescence images (e.g. the calibration sequence).

.. code-block:: python

    # Generate sequencing matches
    files_of_interest.parallel_processing_kwargs['require'] = 'sharedmem'
    files_of_interest.get_sequencing_data(margin=5)
    files_of_interest.parallel_processing_kwargs.pop('require')
    files_of_interest.generate_sequencing_match(overlapping_points_threshold=25,
                                        excluded_sequence_names=['*', 'CalSeq'])
    sequencing_matches = exp.sequencing_matches(files_of_interest)

    # Fine-tuning translation using cross-correlation
    sequencing_matches.parallel.cross_correlation(divider=1/5, gaussian_sigma=1.3, crop=True, plot=False)

    # Further fine-tuning using kernel-correlation
    bounds = ((0.99, 1.01), (-0.01, 0.01), (-1, 1), (-1, 1))
    sequencing_matches.kernel_correlation(bounds, sigma=0.06, crop=True,
                                         strategy='best1bin', maxiter=1000, popsize=50, tol=0.001,
                                         mutation=0.25, recombination=0.7, seed=None, callback=None,
                                         disp=False, polish=True, init='sobol', atol=0,
                                         updating='immediate', workers=1, constraints=())

4. Link single-molecules and sequences
--------------------------------------
Finally, it is time to link the single-molecules and sequences. To this end, the ``destination_distance_threshold`` has
to be set. This threshold indicates the maximum distance (in micrometers) between a single-molecule and a sequencing cluster
for them to be linked.

.. code-block:: python

    sequencing_matches.destination_distance_threshold = 0.2
    sequencing_matches.determine_matched_pairs()

To evaluate how well the matching process worked, you can plot the result. Here, green represents the single-molecules,
red represents the sequencing clusters and blue represents sequence-linked single-molecules. In case the matching process
was successful, all tiles should appear mainly blue.

When satisfied with the matches, you can save them and insert the sequencing data into the single-molecule files datasets:

.. code-block:: python

    sequencing_matches.save()
    files_of_interest.insert_sequencing_data_into_file_dataset()

5. Combine all nc files into a single dataset
---------------------------------------------
For further analysis of the sequence-linked traces, it is convenient to combine all them all into a single dataset:

.. code-block:: python

    import xarray as xr

    files_of_interest = exp.files[exp.files.relativeFilePath.str.regex('Scan')]
    save_path = r'save_path'
    dataset_name = 'complete_dataset.nc'
    ds = xr.open_mfdataset([file.relativeFilePath.with_suffix('.nc') for file in files_of_interest
                            if 'sequence_tile' in file.dataset.data_vars], combine='nested',
                            concat_dim='molecule', data_vars='minimal', coords='minimal',
                            compat='override', engine='h5netcdf', parallel=False)
    ds.to_netcdf(save_path, engine='h5netcdf', mode='w')

To open the dataset:

.. code-block:: python

    ds = xr.open_dataset(, engine='h5netcdf')

