# Quantification of NAS-Score Components
### General Information:
This project was designed to analyze whole slide images (WSI) for liver biopsies in order to gain information about histological features. Data that is gathered can be used to analyze: 

* Macrosteatosis (fatty vesicle count and area)
* Nucleus count and distribution
* Collagen proportionate area

 The analysis uses tissue sections stained in H&E and Masson's Trichrome. Supported datatypes include mrxs and czi files. For more information, please see the corresponding [publication](https://link.springer.com/10.1007/s00292-024-01298-6}). 


### Fat and Nuclei Detection
This code processes a WSI in patches in order to detect fatty vesicles. Optionally, if analyzing a H&E WSI, nuclei can be detected using [HoVer-Net](https://github.com/vqdang/hover_net). Nucleus detection is optional, and can only be performed on H&E stained tissue sections. Fat detection can be performed on both Masson's Trichrome and H&E stained WSIs. 

#### Environment
* Numpy v 1.20.2
* openslide-python v 3.3
* scikit-image v 0.18.1
* pandas v 1.2.4
* findmaxima2d v 0.0.25
* scipy v 1.6.2
* opencv v 4.5.2
* PIL v 8.2.0
* Czifile v 2019.7.2
* json5 v 0.9.6
* Matplotlib v 3.3.4
* Pytorch v 2.2.1
* docopt v 0.6.2
* tqdm v 4.66.2
* imgaug v 0.4.0
* termcolor v 2.4.0
* HoVer-Net

#### Usage:
To start: call run_workflow() in run_fat_and_nuclei_detection.py.
#####  Parameters:
* `output-folder-name` (str)
    * name of folder to store outputs under.
* `output-dir` (str) 
    * where you want the outputs to be stored.
* `input-dir` (str) 
    * where the files you want to process are located.
* `hovernet-dir` (str)
    * path to HoVer-Net code.
* `hovernet-weights-file` (str)
    * path to HoVer-Net weights.
* `WSI-file-type` (str) 
    * can be “.mrxs” or “.czi”.
* `WSI-queue` (list of str) 
    * A list of comma separated file names as a string, ie 'file_a,file_b' [no spaces!] of the WSIs to be processed in the input directory. When not specified, all files with the given file_type in the directory are processed. Default = "".
* `PATCH_SIZE` (int) 
    * patch size in micrometers. Patches are always processed in square. Default is 1500x1500 µm<sup>2</sup>.
* `PATCH_LEVEL` (int) 
    * level at which patches are processed. Default is 0.
* `MIN_AREA` (int) 
    * in µm<sup>2</sup>. Smallest size of a macrovesicular fat vesicle. Default is 40 µm<sup>2</sup>.
* `MAX_AREA` (int)
    * in µm<sup>2</sup>. Maximum area of fat vesicle. Default is 10000 µm<sup>2</sup>. Keep in mind that fatty objects can be connected into one big object.
* `MIN_EXTENT` (int)
    * minimum extent of a fat vesicle. Default is 0.5.
* `FILL_FG_SIZE` (int)
    * fills in objects in the foreground. Default: is calculated from MIN_AREA.
* `FILL_BG_SIZE` (int) 
    * fills in objects in the background. Default: is calculated from MAX_AREA.
* `MAX_AXIS_RATIO` (int) 
    * maximum axis ratio of a fat vesivle. Default: 2.
* `SUBPATCH_SIZES_FACTOR` (int)
    * int representing how many subpatches a patch will be divided into. Default=20.
* `DISK`(float, in µm)
    * Kernel used image preprocessing. Default: 8.3478513356562
* Nucleus Detection Parameters (see [HoVer-Net Repo](https://github.com/vqdang/hover_net) for more information).
    * `run-hovernet` (bool)
    * `hovernet_gpu` (str of an int)
        * which gpus should be used.
    * `hovernet_num_nuc_types` (str of an int)
        * how many nuclei types are to be detected. Depends on weights being used.
    * `hovernet_model_mode` (str)
        * can be “original” or “fast”.
    * `hovernet-run-function` (str) 
        * Default: “./run_infer.py”.
    * `hovernet-batch-size` (str of an int)
        * what batch size should be used.
    * `data_saver_mode` (bool) 
        * if data_saver_mode = True, temporary outputs are automatically deleted. Default = False.

#### Example:
    python .\src\run_fat_and_nuclei_detection.py\
    --input-dir "F:\test\test_img\"\
    --output-dir "F:\test\test_cli"\
    --output-folder-name "fat_and_nuclei_detection"\
    --hovernet-dir "F:\test\hover_net-master"\
    --hovernet-weights-file "F:\test\hover_net-master\hovernet_original_kumar_notype_tf2pytorch.tar"`\
    --WSI-file-type ".mrxs"\
    --WSI-queue "test_img_HE.mrxs"\
    --patch-size "1500"\
    --patch-level "0"\
    --min-fat-vesicle-area "40"\
    --max-fat-vesicle-area "10000"\
    --min-extent "0.5"\
    --max-axis-ratio "2"\
    --subpatch-size-factor "20"\
    --disk "8.35"\
    --run-hovernet\
    --hovernet-gpu "0"\
    --hovernet-num-nuc-types "0"\
    --hovernet-model-mode "original"\
    --hovernet-run-function "./run_infer.py"\
    --hovernet-batch-size "128"\
    --data-saver-mode


#### Output:
Each file is processed, and output for each file is stored in various dataframes. These dataframes contain information about tissue area, fat area, fat objects, mpp, and, if wished, nuclei information. Additionally, binary masks of detected tissue and fatty vesicles are stored on a patchwise basis. For more information about the contents of these dataframes, see fat_and_nucleus_data_analysis_example.ipynb. 
* location: [`output-dir`]/[`output-folder-name`]/[`file_name`]/dataframes/
    * [file_name]_data.csv (dataframe)
        * contains data about each potential fat object that was detected, as well as information about the processed patch.
    * [file_name]_subpatch_df.pkl (dataframe)
        * splits the processed patch into smaller patches, to allow for more granular data analysis.
    * If nucleus information is gathered:
        * [file_name]_raw_global.pkl (dataframe)
            * contains all detected nuclei, even those in unscanned areas.
        *[file_name] _relevant_nuclei.pkl (dataframe)
            * contains only nuclei that lie in regions with relevant detected tissue.
        * [file_name]_fat_and_nucs.pkl (dataframe)
            * contains combined information about fat, tissue area and detected nuclei. This dataframe is required in order to run nucleus cluster analysis. To generate _fat_and_nucs_df.pkl, the WSI must also be processed as subpatches (that means `SUBPATCH_SIZES` must be a list with a length of at least one).
* location: [`output-dir`]/[`output-folder-name`]/[`file_name`]/saved_patches
    * /fat_patches/
        * contains a binary mask of objects determined to be fatty vesicles.
    * /binary_patches/
        * contains a binary mask of detected tissue.
* location: [`output-dir`]/[`output-folder-name`]/[`file_name`]/temporary_data
    * /[`int_of_batch_number`]/
        * /generated_output/
            * location of the generated nucleus information.
        * /original_patches/
            * contains the .tif files of the processed patch. To save memory, these images are automatically deleted if data-saver-mode is activated.


### Nucleus Cluster Analysis

This code takes information about nuclei location (information gained by enabling the HoVer-Net option in the Fat and Nucleus Detection script), and analyzes the nuclei distributions to find nearest neighbors. The script returns information about nuclei clusters, and conducts an area analysis of clusters of nuclei. 


#### Environment
* Numpy v 1.20.2
* openslide-python v 3.3
* scikit-image v 0.18.1
* pandas v 1.2.4
* scipy v 1.6.2
* PIL v 8.2.0
* json5 v 0.9.6
* Matplotlib v 3.3.4
* networkX v 2.5.1


#### Usage:
To start: = run_nucleus_cluster_analysis.py. Fat and Nucleus Detection Script must be executed first!

#####  Parameters:
* `data-dir` (str) 
    * path to output of Fat and Nuclei Detection Script, i.e [`output-dir`]/[`output-folder-name`]. This folder should contain a folder for each processed WSI, as described in the output description of the Fat and Nucleus Detection Readme.
* `max-distance` (int) 
    *  maximum distance (in µm) between nuclei to still be considered neighbors. Default: 15 µm.
* `output-folder-name` (str) 
    * name of experiment, used to generate a folder to store results.

#### Example:
    python run_nucleus_cluster_analysis.py\
    --data-dir "F:\test\test_cli\fat_and_nuclei_detection_results\"\
    --max-distance 15\
    --output-folder-name "test_cluster"



#### Output
For each WSI that is processed, nuclei are clustered, and information about clusters is gathered ("graph stats"). Additionally, information about nucleus location and cluster area (after removing fatty objects that lay within cluster bounds) is stored. For more information about how to use output, see nucleus_cluster_analysis_data_analysis_example.ipynb
* location: [`data-dir`]/[`file_name`]/nucleus_cluster_analysis/[`output-folder-name`]/
    * [file_name]_[`output-folder-name`]_clusters.pkl (lst of lsts).
        * a list of components. Each component contains a list of each nucleus coordinate in that cluster.
    * [file_name]_[`output-folder-name`]_graph_stats.pkl (lst of lsts)
        * structured as a list of lists of nucleus [information](https://networkx.org/documentation/stable/reference/functions.html). Note that graph_stats[1-8] lists each have the length = number of relevant nuclei detected:
            * graph_stats[0] = average clustering coefficient (list of a float with the length of one)<br>
            * graph_stats[1] = clustering coefficients (list)<br>
            * graph_stats[2] = number of neighbors (list)<br>
            * graph_stats[3] = number of non_neighbors (list)<br>
            * graph_stats[4] = average distance (in pixels) to neighbors (list)<br> 
            * graph_stats[5] = average distance (in pixels) to non-neighbor (list)<br>
            * graph_stats[6] = number of common neighbors (list)<br>
            * graph_stats[7] = degree (list)<br>
            * graph_stats[8] = is nucleus in a cluster (list)<br>
    
    * [file_name]_[`output-folder-name`]_global_CCA.pkl (dataframe)
        * this dataframe puts nuclei information on a global (ie WSI) scale. This dataframe contains information about fatty area, nuclei, and area of nuclei clusters.  



### Additional Nucleus Analysis
This script allows an analysis of the number of nuclei clusters in a given (square-shaped) field of view (FOV). To use this, make sure to first run the Nucleus Cluster Analysis script.

#### Environment
* see Nucleus Cluster Analysis Script environment

#### Usage: 
* To start: run_inflammation_by_fov.py
    * By adjusting function “find_20x_fov” you can change the size of FOV. In the script, a 20xFOV is ~1120 x 1120 µm<sup>2</sup>.

    * ##### Parameters
        * `data-dir` (str) 
            * path to output of Fat and Nuclei Detection Script, i.e [`output-dir`]/[`output-folder-name`]. 
        * `CCA-folder-name` (str)
            * name of folder under which the Nucleus Cluster Analysis output is stored.
        * `FOV_output` (str)
            * Folder name to store the FOV results.

#### Example:
    python run_inflammation_by_fov.py\
    --data-dir "F:\test\test_cli\fat_and_nuclei_detection_results\"\
    --CCA-folder-name "test_cluster"\
    --FOV_output "20xFOV"\


#### Output
For each WSI that is processed, a dataframe containing the coordiantes of the FOV and the number of centroids within the FOV is generated.
*  location: [`data-dir`]/[`FOV_output`]/
    * [file_name]_foki_per_20xfov.pkl 


### CPA Analysis
This script analyzes the collagen proportionate area (CPA) in Masson's Trichrome stained WSI. 

#### Environment
* see Fat and Nuclei Detection environment
* scikit-learn v 0.24.2

#### Usage:
* To start: call detect_CPA() in run_CPA_anlaysis.py

#####  Parameters:
* `output-folder-name` (str)
    * name of folder to store outputs under.
* `output-dir` (str) 
    * where you want the outputs to be stored.
* `input-dir` (str) 
    * where the files you want to process are located.
* `WSI-file-type` (str) 
    * can be “.mrxs” or “.czi”.
* `WSI-queue` (list of str) 
    * a list of all the images in the file_loc that you wish to process. If imgs_to_process = [], all files with FILE_TYPE in file_loc are processed.
* `PATCH_SIZE` (int) 
    * patch size in micrometers. Patches are always processed in square. Default is 1500x1500 µm<sup>2</sup>.
* `PATCH_LEVEL` (int) 
    * level at which patches are processed. Default is 0.
* `MIN_AREA` (int) 
    * in µm<sup>2</sup>. Smallest size of a macrovesicular fat vesicle. Default is 40 µm<sup>2</sup>.
* `MAX_AREA` (int)
    * in µm<sup>2</sup>. Maximum area of fat vesicle. Default is 10000 µm<sup>2</sup>. Keep in mind that fatty objects can be connected into one big object.
* `MIN_EXTENT` (int)
    * minimum extent of a fat vesicle. Default is 0.5.
* `FILL_FG_SIZE` (int)
    * fills in objects in the foreground. Default: is calculated from MIN_AREA.
* `FILL_BG_SIZE` (int) 
    * fills in objects in the background. Default: is calculated from MAX_AREA.
* `MAX_AXIS_RATIO` (int) 
    * maximum axis ratio of a fat vesicle. Default: 2.
* `DISK`(float, in µm)
    * Kernel used image preprocessing. Default: 8.3478513356562.

#### Example:

    python .\src\run_CPA_analysis.py \
    --input-dir "F:\test\test_img\" \
    --output-dir "F:\test\test_cli"` \
    --output-folder-name "CPA" \
    --WSI-file-type ".mrxs" \
    --WSI-queue "test_img_HE.mrxs" \
    --patch-size "1500" \
    --patch-level "0" \
    --min-fat-vesicle-area "40" \
    --max-fat-vesicle-area "10000"\
    --min-extent "0.5" \
    --max-axis-ratio "2" \
    --disk "8.35" \

#### Output
For each WSI that is processed, a dataframe with the analysis data is saved. Additionally, binary masks of tissue and fat objects are saved, which can be useful for plotting and further analysis. To see examples how to use generated data, see CPA_data_analysis_example.ipynb.
* [`output-dir`]/[`output-folder-name`]/dataframes
    * [file_name]_CPA.pkl (dataframe)
        * contains the following information on a patch-wise basis:
            * WSI information:
                * `name` name of file processed.
                * `mpps` micrometers per pixels of WSI.
            * Patch information: 
                * `original_patch_key`: patch key.
                * `global_x_coords(pxs)`: x coordinate of the patch on a global (WSI) scale, in pixels.
                * `global_y_coords(pxs)`: y coordinate of hte patch on a global (WSI) scale, in pixels.
                * `num_blue_pixels`: number of fibrotic pixels.
                * `total_tissue_pixels(with_fat)`: total number of tissue pixels detected with fatty pixels.
                * `total number of tissues without fatty pixels`: total number of tissue pixels detected with fatty pixels.
                *  `num_fat_pixels`: number of fatty pixels.



We hope you find our work helpful! If you use our work, please cite the following paper:

```
@article{darling_herausforderungen_2024,
	title = {Herausforderungen der {Automation} bei der quantitativen {Auswertung} von {Leberbiopsien}: {Automatische} {Quantifizierung} von {Leberverfettung}},
	volume = {45},
	issn = {2731-7188, 2731-7196},
	shorttitle = {Herausforderungen der {Automation} bei der quantitativen {Auswertung} von {Leberbiopsien}},
	url = {https://link.springer.com/10.1007/s00292-024-01298-6},
	doi = {10.1007/s00292-024-01298-6},
	language = {de},
	number = {2},
	urldate = {2024-03-27},
	journal = {Die Pathologie},
	author = {Darling, Jessica and Abedin, Nada and Ziegler, Paul K. and Gretser, Steffen and Walczak, Barbara and Barreiros, Ana Paula and Schulze, Falko and Reis, Henning and Wild, Peter J. and Flinner, Nadine},
	month = March,
	year = {2024},
	pages = {115--123},
}
```



### Thank you for the support!

```
@software{reback2020pandas,
    author       = {The pandas development team},
    title        = {pandas-dev/pandas: Pandas},
    month        = feb,
    year         = 2020,
    publisher    = {Zenodo},
    version      = {1.2.4},
    doi          = {10.5281/zenodo.3509134},
    url          = {https://doi.org/10.5281/zenodo.3509134}
    note         = [license](https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html)
}

@InProceedings{ mckinney-proc-scipy-2010,
  author    = { {W}es {M}c{K}inney },
  title     = { {D}ata {S}tructures for {S}tatistical {C}omputing in {P}ython },
  booktitle = { {P}roceedings of the 9th {P}ython in {S}cience {C}onference },
  pages     = { 56 - 61 },
  year      = { 2010 },
  editor    = { {S}t\'efan van der {W}alt and {J}arrod {M}illman },
  doi       = { 10.25080/Majora-92bf1922-00a }
  note         = [license](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
}

@ARTICLE{2020SciPy-NMeth,
  author  = {Virtanen, Pauli and Gommers, Ralf and Oliphant, Travis E. and
            Haberland, Matt and Reddy, Tyler and Cournapeau, David and
            Burovski, Evgeni and Peterson, Pearu and Weckesser, Warren and
            Bright, Jonathan and {van der Walt}, St{\'e}fan J. and
            Brett, Matthew and Wilson, Joshua and Millman, K. Jarrod and
            Mayorov, Nikolay and Nelson, Andrew R. J. and Jones, Eric and
            Kern, Robert and Larson, Eric and Carey, C J and
            Polat, {\.I}lhan and Feng, Yu and Moore, Eric W. and
            {VanderPlas}, Jake and Laxalde, Denis and Perktold, Josef and
            Cimrman, Robert and Henriksen, Ian and Quintero, E. A. and
            Harris, Charles R. and Archibald, Anne M. and
            Ribeiro, Ant{\^o}nio H. and Pedregosa, Fabian and
            {van Mulbregt}, Paul and {SciPy 1.0 Contributors}},
  title   = {{{SciPy} 1.0: Fundamental Algorithms for Scientific
            Computing in Python}},
  journal = {Nature Methods},
  year    = {2020},
  volume  = {17},
  pages   = {261--272},
  adsurl  = {https://rdcu.be/b08Wh},
  doi     = {10.1038/s41592-019-0686-2},
note         = [license](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
}

@InProceedings{SciPyProceedings_11,
  author =       {Aric A. Hagberg and Daniel A. Schult and Pieter J. Swart},
  title =        {Exploring Network Structure, Dynamics, and Function using NetworkX},
  booktitle =   {Proceedings of the 7th Python in Science Conference},
  pages =     {11 - 15},
  address = {Pasadena, CA USA},
  year =      {2008},
  editor =    {Ga\"el Varoquaux and Travis Vaught and Jarrod Millman},
  note = [license](https://networkx.org/documentation/stable/)
}

@Article{Hunter:2007,
  Author    = {Hunter, J. D.},
  Title     = {Matplotlib: A 2D graphics environment},
  Journal   = {Computing in Science \& Engineering},
  Volume    = {9},
  Number    = {3},
  Pages     = {90--95},
  abstract  = {Matplotlib is a 2D graphics package used for Python for
  application development, interactive scripting, and publication-quality
  image generation across user interfaces and operating systems.},
  publisher = {IEEE COMPUTER SOC},
  doi       = {10.1109/MCSE.2007.55},
  year      = 2007
  note         = [license](https://matplotlib.org/stable/users/project/license.html)
}

@misc{paszke2019pytorch,
      title={PyTorch: An Imperative Style, High-Performance Deep Learning Library}, 
      author={Adam Paszke and Sam Gross and Francisco Massa and Adam Lerer and James Bradbury and Gregory Chanan and Trevor Killeen and Zeming Lin and Natalia Gimelshein and Luca Antiga and Alban Desmaison and Andreas Köpf and Edward Yang and Zach DeVito and Martin Raison and Alykhan Tejani and Sasank Chilamkurthy and Benoit Steiner and Lu Fang and Junjie Bai and Soumith Chintala},
      year={2019},
      eprint={1912.01703},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
      note=[license](https://github.com/pytorch/pytorch?tab=License-1-ov-file#readme)
}

@article{da Costa-Luis2019, 
doi = {10.21105/joss.01277}, 
url = {https://doi.org/10.21105/joss.01277}, 
year = {2019}, 
publisher = {The Open Journal}, 
volume = {4}, number = {37},
pages = {1277},
author = {Casper O. da Costa-Luis},
title = {`tqdm`: A Fast, Extensible Progress Meter for Python and CLI},
journal = {Journal of Open Source Software}
note = [license](https://github.com/tqdm/tqdm?tab=License-1-ov-file#readme)
} 

@misc{imgaug,
  author = {Alexander Jung},
  title = {imgaug},
  howpublished = "{https://imgaug.readthedocs.io/en/latest/}",
  year = {2020}, 
  note = [license](https://github.com/aleju/imgaug-doc/blob/master/LICENSE)
}

@misc{termcolor 2.4.0,
  author = {Konstantin Lepa},
  title = {termcolor 2.4.0},
  howpublished = "{https://pypi.org/project/termcolor/}",
  year = {2023}, 
  note = [license](https://pypi.org/project/termcolor/)
}

@Article{         harris2020array,
 title         = {Array programming with {NumPy}},
 author        = {Charles R. Harris and K. Jarrod Millman and St{\'{e}}fan J.
                 van der Walt and Ralf Gommers and Pauli Virtanen and David
                 Cournapeau and Eric Wieser and Julian Taylor and Sebastian
                 Berg and Nathaniel J. Smith and Robert Kern and Matti Picus
                 and Stephan Hoyer and Marten H. van Kerkwijk and Matthew
                 Brett and Allan Haldane and Jaime Fern{\'{a}}ndez del
                 R{\'{i}}o and Mark Wiebe and Pearu Peterson and Pierre
                 G{\'{e}}rard-Marchant and Kevin Sheppard and Tyler Reddy and
                 Warren Weckesser and Hameer Abbasi and Christoph Gohlke and
                 Travis E. Oliphant},
 year          = {2020},
 month         = sep,
 journal       = {Nature},
 volume        = {585},
 number        = {7825},
 pages         = {357--362},
 doi           = {10.1038/s41586-020-2649-2},
 publisher     = {Springer Science and Business Media {LLC}},
 url           = {https://doi.org/10.1038/s41586-020-2649-2}
 notes         = [license](https://numpy.org/doc/stable/license.html)
}

@ARTICLE{Goode2013-cz,
  title     = "{OpenSlide}: A vendor-neutral software foundation for digital
               pathology",
  author    = "Goode, Adam and Gilbert, Benjamin and Harkes, Jan and Jukic,
               Drazen and Satyanarayanan, Mahadev",
  abstract  = "Although widely touted as a replacement for glass slides and
               microscopes in pathology, digital slides present major
               challenges in data storage, transmission, processing and
               interoperability. Since no universal data format is in
               widespread use for these images today, each vendor defines its
               own proprietary data formats, analysis tools, viewers and
               software libraries. This creates issues not only for
               pathologists, but also for interoperability. In this paper, we
               present the design and implementation of OpenSlide, a
               vendor-neutral C library for reading and manipulating digital
               slides of diverse vendor formats. The library is extensible and
               easily interfaced to various programming languages. An
               application written to the OpenSlide interface can transparently
               handle multiple vendor formats. OpenSlide is in use today by
               many academic and industrial organizations world-wide, including
               many research sites in the United States that are funded by the
               National Institutes of Health.",
  journal   = "J. Pathol. Inform.",
  publisher = "Elsevier BV",
  volume    =  4,
  number    =  1,
  pages     = "27",
  month     =  sep,
  year      =  2013,
  keywords  = "Diamond; OpenDiamond; PathFind; digital imaging and
               communications in medicine; digital slide; scanner; whole-slide
               image",
  copyright = "https://creativecommons.org/licenses/by-nc-sa/3.0/",
  language  = "en"  
  note      = [license](https://openslide.org/license/)
}

@article{scikit-image,
 title = {scikit-image: image processing in {P}ython},
 author = {van der Walt, {S}t\'efan and {S}ch\"onberger, {J}ohannes {L}. and
           {Nunez-Iglesias}, {J}uan and {B}oulogne, {F}ran\c{c}ois and {W}arner,
           {J}oshua {D}. and {Y}ager, {N}eil and {G}ouillart, {E}mmanuelle and
           {Y}u, {T}ony and the scikit-image contributors},
 year = {2014},
 month = {6},
 keywords = {Image processing, Reproducible research, Education,
             Visualization, Open source, Python, Scientific programming},
 volume = {2},
 pages = {e453},
 journal = {PeerJ},
 issn = {2167-8359},
 url = {https://doi.org/10.7717/peerj.453},
 doi = {10.7717/peerj.453}
 note = [license](https://scikit-image.org/docs/stable/license.html)
}

@article{opencv_library,
    author = {Bradski, G.},
    citeulike-article-id = {2236121},
    journal = {Dr. Dobb's Journal of Software Tools},
    keywords = {bibtex-import},
    posted-at = {2008-01-15 19:21:54},
    priority = {4},
    title = {{The OpenCV Library}},
    year = {2000}
    note = [license](https://github.com/opencv/opencv/blob/master/LICENSE)
}

@misc{clark2015pillow,
  title={Pillow (PIL Fork) Documentation},
  author={Clark, Alex},
  year={2015},
  publisher={readthedocs},
 url={https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf}
  note = [license](https://github.com/python-pillow/Pillow/blob/main/LICENSE)
}

@misc{graham2019hovernet,
      title={HoVer-Net: Simultaneous Segmentation and Classification of Nuclei in Multi-Tissue Histology Images}, 
      author={Simon Graham and Quoc Dang Vu and Shan E Ahmed Raza and Ayesha Azam and Yee Wah Tsang and Jin Tae Kwak and Nasir Rajpoot},
      year={2019},
      eprint={1812.06499},
      archivePrefix={arXiv},
      primaryClass={cs.CV}
      note = [license](https://github.com/vqdang/hover_net/blob/master/LICENSE)
      
}

@article{scikit-learn,
  title={Scikit-learn: Machine Learning in {P}ython},
  author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
          and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
          and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
          Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
  journal={Journal of Machine Learning Research},
  volume={12},
  pages={2825--2830},
  year={2011}
  note = [license](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING)
}
```

#### Github Repos:
* MaximaFiner, dwaithe, May 26, 2021. https://github.com/dwaithe/MaximaFinder?tab=readme-ov-file. https://github.com/dwaithe/MaximaFinder/blob/master/LICENSE

* czifile, Christoph Gohlke, Laboratory for Fluorescence Dynamics. University of California, Irvine, July 3, 2019. https://github.com/cgohlke/czifile .https://github.com/cgohlke/czifile/blob/master/LICENSE

* docopt, contributors, https://github.com/docopt/docopt. https://github.com/docopt/docopt?tab=MIT-1-ov-file#readme

* json5, Aseem Kishore and contributors, https://github.com/json5/json5, https://github.com/json5/json5/blob/main/LICENSE.md


#### Stackoverflow Answers:
* https://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask, Isaac Sutherland, Sep 17, 2010.

* https://stackoverflow.com/questions/31562534/scipy-centroid-of-convex-hull,  Answered by Alec Day, Sept 11 2019. Edited Jul 18, 2021 by klwire.

* https://stackoverflow.com/questions/35402609/point-on-circle-base-on-given-angle/35402676,  Quentin Pradet, Feb 15, 2016.


